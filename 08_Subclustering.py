"""
Script for analyzing single-cell RNA sequencing data by subsetting a specific cell type, 
performing clustering, and identifying differentially expressed genes (DEGs).

This script loads the manually annotated AnnData object generated in 07_Manual_annotation.py,
subsets the data based on a specified cell type, scales the gene expression, performs PCA, clustering,
and generates UMAP visualizations. It also computes average gene expression per cluster and performs
DEG analysis using the Wilcoxon test, saving results in Excel files and updating the AnnData object.

"""

######### Load the required packages ######### 

import os
import scanpy as sc # version 1.9.8
import anndata as an # version 0.10.5.post1
import pandas as pd # version 2.2.0

######### Specify the paths needed #########

# Path to the directory where figures and results will be saved
resultsPath = "/Path/To/Save/Plots/And/Figures/"

# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"

# Path to the directory where the per cell type subdirectories will be stored
technicalPath = "/Path/To/Store/Per/Cell/Type/Subdirectories/"

######### Code #########

# Load the main AnnData object containing cell annotations
filename = 'integrated_Manual_annotation'
adata_main = an.read_h5ad(integrationPath + 'adata_' + filename + '.h5ad')

# Define the cell type of interest for further analysis
celltype = "Cell_type_of_interest"
celltypePath = technicalPath + celltype + "/"

# Create a folder for this specific cell type if it doesn't exist
if not os.path.exists(celltypePath):
   os.makedirs(celltypePath)

# Subset the main object based on the selected cell type
adata = adata_main[adata_main.obs['Manual_annotation'].isin([celltype.replace('_', ' ')])] 
# `.replace('_', ' ')` ensures consistency in cell type names

# Check the unique values in the annotation to confirm successful subsetting
adata.obs['Manual_annotation'].unique()

# Scale gene expression to unit variance (clip extreme values at max 10)
sc.pp.scale(adata, max_value=10)

# Perform Principal Component Analysis (PCA) for dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Set the clustering resolution
res = 0.5

# Create a subfolder to store results for this resolution
resFolder = celltypePath + 'res_' + str(res) + '/'
if not os.path.exists(resFolder):
   os.makedirs(resFolder)

# Change the working directory to the resolution-specific folder
os.chdir(resFolder)

# Create a subfolder within resFolder to store Excel files
resExcelPath = resFolder + 'Excel_files/'
if not os.path.exists(resExcelPath):
   os.makedirs(resExcelPath)

# Compute nearest neighbors using the scVI representation
sc.pp.neighbors(adata, use_rep="X_scVI")

# Compute UMAP embedding for visualization
sc.tl.umap(adata)

# Perform clustering using Leiden algorithm at the specified resolution
sc.tl.leiden(adata, resolution=res)

# Visualize the UMAP embeddings colored by Leiden clusters
sc.pl.umap(adata, color="leiden", frameon=False, 
    save='_' + celltype + '_leiden' + '_res_' + str(res) + '.pdf')

# Visualize UMAP with a custom metadata column called "Paper"
sc.pl.umap(adata, color="Paper", frameon=False, 
    save='_' + celltype + '_paper' + '_res_' + str(res) + '.pdf')


### Compute average and fraction gene expression per cluster ###
# Extract gene names and clusters
gene_ids = adata.raw.var.index.values
clusters = adata.obs['leiden'].cat.categories

# Convert raw counts to a Pandas DataFrame for easier manipulation
obs = adata.raw[:, gene_ids].X.toarray()
obs = pd.DataFrame(obs, columns=gene_ids, index=adata.obs['leiden'])

# Compute the average gene expression per cluster
average_obs = obs.groupby(level=0).mean()

# Transpose and save the average expression to Excel files
average_obs.T.to_excel(resExcelPath + celltype + '_res_' + str(res) + '_average_expression.xlsx')

# Save the processed AnnData object
results_file = resFolder + 'adata_' + celltype + '_res_' + str(res) + '.h5ad'
adata.write(results_file)

### Identify differentially expressed genes (DEGs) ###

# Perform differential gene expression analysis using Wilcoxon test
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# Convert DEG results into a DataFrame
DEG_df = sc.get.rank_genes_groups_df(adata, group=None)

# Filter to retain only positive markers (upregulated genes)
positive_DEG_df = DEG_df[DEG_df['logfoldchanges'] > 0]

# Save the positive DEGs to an Excel file
positive_DEG_df.to_excel(resExcelPath + celltype + '_res_' + str(res) + '_positive_DEG.xlsx')

### Compute the average expression of DEGs per cluster ###

# Extract the list of differentially expressed genes
DEG_genes = positive_DEG_df['names']

# Filter the average expression DataFrame to retain only DEGs
average_obs = average_obs.T  # Transpose to match gene index format
filtered_average_obs = average_obs[average_obs.index.isin(DEG_genes)]

# Save the filtered average expression of DEGs per cluster to an Excel file
excel_file = resExcelPath + celltype + '_res_' + str(res) + '_DEGs_average_expression_per_clusters.xlsx'
filtered_average_obs.to_excel(excel_file)

# Save the AnnData object with DEGs
results_file = resFolder + 'adata_' + celltype + '_res_' + str(res) + '_DEG.h5ad'
adata.write(results_file)
