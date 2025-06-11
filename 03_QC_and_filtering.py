"""
Script for performing quality control (QC) on single-cell RNA sequencing data.

This script loads the  filtered AnnData object generated in 02_Doublet_filtration.py, computes QC metrics,
identifies low-quality cells using Median Absolute Deviation (MAD) thresholding, and removes outliers.
It then normalizes the data, selects highly variable genes, and regresses out unwanted sources of variation
before saving the processed dataset.

"""

######### Load the required packages #########   
import os  
import scanpy as sc # version 1.9.8
import seaborn as sns # version 0.13.2
import numpy as np # version 1.26.6     
from scipy.stats import median_abs_deviation # version 1.12.0

######### Specify the paths needed #########

# Path to the directory where figures and results will be saved
resultsPath = "/Path/To/Save/Plots/And/Figures/"

# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"


######### Code #########

# Set figure parameters: size and frame visibility
sc.set_figure_params(figsize=(6, 6), frameon=False)

# Set the number of parallel jobs to use (2 jobs in this case)
sc.settings.n_jobs = 2

# Set verbosity level for Scanpy's logging (verbosity: errors (0), warnings (1), info (2), hints (3))
sc.settings.verbosity = 3
sc.logging.print_header()

# Set figure parameters: DPI and background color
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Change the current working directory to save plots in the correct folder
os.chdir(resultsPath)

# The file that stores the input object
input_filename = 'pre-integration_filtered'

# Read the filtered AnnData object after doublet filtration
adata = sc.read_h5ad(integrationPath + 'adata_' + input_filename + '.h5ad')

# Make sure gene names are unique by appending a suffix to duplicate gene names
adata.var_names_make_unique()

# Annotate mitochondrial genes (genes starting with "MT-")
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Annotate ribosomal genes (genes starting with "RPS" or "RPL")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

# Annotate heatshock protein genes (genes starting with "HSP")
adata.var["hsp"] = adata.var_names.str.startswith(("HSP"))

# Annotate hemoglobin genes (genes starting with "HB" and not followed by "P")
adata.var["hb"] = adata.var_names.str.startswith(("HB[^(P)]"))

# Calculate QC metrics based on mitochondrial (mt), ribosomal (ribo), heatshock protein (hsp), and hemoglobin (hb) genes
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hsp", "hb"], percent_top=[20], log1p=True, inplace=True)

# Plot QC metrics: total counts, mitochondrial percentage, and gene count per cell
p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# Perform QC filtering using MAD (Median Absolute Deviation) thresholding
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# Apply MAD-based outlier detection to QC metrics
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()

# Filter mitochondrial outliers and cells with >10% mitochondrial genes
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 10)
adata.obs.mt_outlier.value_counts()

# Filter the AnnData object by removing outliers (both general and mitochondrial)
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low-quality cells: {adata.n_obs}")

# Plot filtered data (scatter plot of total counts vs. genes per cell)
p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# Filter genes not detected in at least 20 cells
sc.pp.filter_genes(adata, min_cells=20)

# Filter out cells with fewer than 200 genes
sc.pp.filter_cells(adata, min_genes=200)

# Keep the raw counts intact for later normalization
adata.layers["counts"] = adata.X.copy()

# Normalize the data (scale total counts to 10,000 per cell)
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# Identify highly variable genes based on the data (top 4000 genes)
sc.pp.highly_variable_genes(adata, n_top_genes=4000)

# Store the raw counts in the "raw" attribute of the AnnData object
adata.raw = adata

# Filter the dataset to keep only the highly variable genes
adata = adata[:, adata.var.highly_variable]

# Regress out the effects of total counts and mitochondrial genes
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])

# Scale each gene to unit variance
sc.pp.scale(adata, max_value=10)

# Save the final processed AnnData object
result_filename = 'pre-integration_quality_control'
adata.write(integrationPath + 'adata_' + result_filename + '.h5ad')
