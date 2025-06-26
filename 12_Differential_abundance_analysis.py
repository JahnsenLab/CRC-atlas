"""
Script for performing differential abundance analysis using Milo on single-cell RNA sequencing data.

This script loads the per-cell-type AnnData object generated in 08_Subclustering.py, and performs KNN graph construction, neighborhood identification,
and differential abundance testing using the Milo package. It visualizes results via UMAP and diagnostic plots, and generates 
summary statistics and plots to assess the differential abundance between specified conditions.

"""

######### Load the required packages ######### 
import os
import numpy as np # version 1.26.6
import pandas as pd # version 2.2.0
import seaborn as sns # version 0.13.2
import scanpy as sc # version 1.9.8
import anndata as an # version 0.10.5.post1
import milopy # version 0.1.1
import milopy.utils # version 0.1.1
import milopy.plot as milopl # version 0.1.1
import matplotlib.pyplot as plt # version 0.1.1
import matplotlib # version 3.8.3


######### Specify the paths needed #########

# Path to load the AnnData object
objectPath = "Path/To/Your/Anndata/Object/" 

# Path to store generated results and plots
resultsPath = "Path/To/Store/Results/"  

######### Define cell type of interest #########
celltype = 'Cell_type_of_interest'  # The cell type analyzed, used for creating cell type specific folders and file name

######### Spesify what to analyse (based on metadata) #########

# Column in the metadata for annotation
annotationColumn = "Your_Annotation_Column"  

# Define which conditions to analyse 
analysed = 'Conditions_analysed'  # Example: 'tumor_vs_normal' or 'stage_comparison'

# Set up the design and order of conditions for the analysis
conditionColumn = 'Name_of_metadata_column_containing_conditions_to_analyse'  # Replace with actual column name (e.g., 'Treatment', 'Stage')
condition_1 = 'First_condition_to_analyse'  # Replace with the first condition (e.g., 'Tumor') 
condition_2 = 'Second_condition_to_analyse'  # Replace with the second condition (e.g., 'Normal')

# Build the design formula
design_used = '~' + conditionColumn  # Specify the experimental design, usually the condition column

# Build the contrast for the analysis (compare the two conditions)
contrasts_used = conditionColumn + condition_1 + '-' + conditionColumn + condition_2  # e.g., 'Tumor - Normal'

######### Code #########

# Rescale the figures for better visualization
plt.rcParams['figure.figsize'] = (8, 8) 
sc.settings.verbosity = 3  # Set verbosity level for Scanpy output

# Define K and D values for analysis
k = 10  # Number of neighbors for KNN graph
d = 30  # Number of principal components to use

### LOAD THE MAIN OBJECT ###
adata = an.read_h5ad(objectPath + 'adata_' + celltype + '.h5ad')  # Read the AnnData object

# Visualize the UMAP of the data, colored by annotation and SampleType
sc.pl.umap(adata, color=[annotationColumn], legend_loc="on data")
sc.pl.umap(adata, color=["SampleType"])

# Build the KNN graph with specified neighbors (k) and PCs (d)
sc.pp.neighbors(adata, n_neighbors=k, n_pcs=d)

# Create a folder for the specific cell type and K value to store plots
miloPlotPath = resultsPath + 'Milo/' + 'k_' + str(k) + '/' + celltype + '/'

# Create the folder if it doesn't exist
if not os.path.exists(miloPlotPath):
    os.makedirs(miloPlotPath)

# Set the newly created folder as working directory
os.chdir(miloPlotPath)

# Construct neighborhoods using Milo (assigns cells to representative neighborhoods in the KNN graph)
milo.make_nhoods(adata, prop=0.1)

# Check the median number of cells in each neighborhood
nhood_size = adata.obsm["nhoods"].toarray().sum(0)
plt.hist(nhood_size, bins=20)
plt.xlabel("# cells in neighbourhood")
plt.ylabel("# neighborhoods")
plt.savefig('Neighborhood_distribution_' + celltype + '_k_' + str(k) + '_d_' + str(d) + '_' + analysed + '.pdf', format='pdf', bbox_inches='tight')
np.median(nhood_size)

# Count cells in neighborhoods by sample type
milo.count_nhoods(adata, sample_col="SampleID")

# Differential abundance testing using GLM, with specified experimental design
milo.DA_nhoods(adata, design=design_used, model_contrasts=contrasts_used)

# Check p-values from the differential abundance test
adata.uns["nhood_adata"].obs['PValue']

# Define a function for plotting diagnostic results from the differential abundance test
def plot_milo_diagnostics(adata):
    alpha = 0.1  # significance threshold
    
    with matplotlib.rc_context({"figure.figsize": [12, 12]}):
        # P-value histogram
        plt.subplot(2, 2, 1)
        plt.hist(adata.uns["nhood_adata"].obs["PValue"], bins=20)
        plt.xlabel("Uncorrected P-value")

        # Scatter plot of uncorrected P-value vs SpatialFDR
        plt.subplot(2, 2, 2)
        plt.scatter(adata.uns["nhood_adata"].obs["PValue"], adata.uns["nhood_adata"].obs["SpatialFDR"], s=3)
        plt.xlabel("Uncorrected P-value")
        plt.ylabel("SpatialFDR")

        # Volcano plot
        plt.subplot(2, 2, 3)
        plt.scatter(adata.uns["nhood_adata"].obs["logFC"], -np.log10(adata.uns["nhood_adata"].obs["SpatialFDR"]), s=3)
        plt.axhline(y=-np.log10(alpha), color="red", linewidth=1, label=f"{int(alpha*100)} % SpatialFDR")
        plt.legend()
        plt.xlabel("log-Fold Change")
        plt.ylabel("- log10(SpatialFDR)")

        # MA plot
        plt.subplot(2, 2, 4)
        df = adata.uns["nhood_adata"].obs
        emp_null = df[df["SpatialFDR"] >= alpha]["logFC"].mean()
        df["Sig"] = df["SpatialFDR"] < alpha
        plt.scatter(df["logCPM"], df["logFC"], c=df["Sig"], cmap="coolwarm")
        plt.axhline(y=0, color="grey", linewidth=1)
        plt.axhline(y=emp_null, color="purple", linewidth=1)
        plt.legend(title=f"< {int(alpha*100)} % SpatialFDR")
        plt.xlabel("Mean log-counts")
        plt.ylabel("log-Fold Change")
        plt.tight_layout()

    plt.show()

# Plot diagnostics
plot_milo_diagnostics(adata)

# Save diagnostic plot as a PDF
plt.savefig('Significance_plots_' + celltype + '_k_' + str(k) + '_d_' + str(d) + '_' + analysed + '_new.pdf', format='pdf', bbox_inches='tight')

# Visualize results on UMAP embedding by plotting differential abundance results
milopy.utils.build_nhood_graph(adata)
milopl.plot_nhood_graph(adata, alpha=0.01, min_size=2)
plt.savefig('DA_log_fold_change_umap_' + celltype + '_k_' + str(k) + '_d_' + str(d) + '_' + analysed + '_new.pdf', format='pdf', bbox_inches='tight')

# Annotate neighborhoods with cell type labels
milopy.utils.annotate_nhoods(adata, anno_col=annotationColumn)

# Visualize the fraction of cell types in neighborhoods
plt.hist(adata.uns['nhood_adata'].obs["nhood_annotation_frac"])
plt.xlabel("Celltype fraction")
plt.savefig('Celltype_fraction_' + celltype + '_k_' + str(k) + '_d_' + str(d) + '_' + analysed + '_new.pdf', format='pdf', bbox_inches='tight')

# Rename mixed neighborhoods where less than 60% of the cells have the top label
adata.uns['nhood_adata'].obs["nhood_annotation"] = adata.uns['nhood_adata'].obs["nhood_annotation"].astype(str)
adata.uns['nhood_adata'].obs.loc[adata.uns['nhood_adata'].obs["nhood_annotation_frac"] < 0.6, "nhood_annotation"] = "Mixed"

# Sort neighborhoods by annotation
adata.uns['nhood_adata'].obs = adata.uns['nhood_adata'].obs.sort_values("nhood_annotation")

# Plot the log fold change across different neighborhood annotations
sc.pl.violin(adata.uns['nhood_adata'], "logFC", groupby="nhood_annotation", rotation=90, show=False)
plt.axhline(y=0, color='black', linestyle='--')
plt.ylim(-9, 9)  # Adjust y-axis range
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel('logFC', fontsize=16)
plt.savefig('Violin_log_fold_change_' + celltype + '_k_' + str(k) + '_d_' + str(d) + '_' + analysed + '.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.show()
