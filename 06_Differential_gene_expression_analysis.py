"""
Script for performing differential gene expression analysis on integrated single-cell RNA sequencing data.

This script loads the integrated and clustered AnnData object generated in 05_Visualization_of_integrated_data.py,
performs DGE analysis using the Wilcoxon rank-sum test, extracts differentially expressed genes (DEGs),
and saves the results as an Excel file. The updated AnnData object is also saved with DEG results.

"""


######### Load the required packages ######### 
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
import openpyxl

######### Specify the paths needed #########

# Path to the directory where figures and results will be saved
resultsPath = "/Path/To/Save/Plots/And/Figures/"

# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"

# Path to the directory where tables will be saved
excelPath = "/Path/To/Save/Tables/"

######### Code #########

# Set figure parameters (size and appearance)
sc.set_figure_params(figsize = (6, 6), frameon = False)

# Set the number of parallel jobs for Scanpy operations
sc.settings.n_jobs = 2

# Set verbosity level (0 = errors, 1 = warnings, 2 = info, 3 = hints)
sc.settings.verbosity = 3             
sc.logging.print_header()  

# Set additional figure parameters
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Set the working directory to ensure figures are saved in the correct folder
os.chdir(resultsPath)

# Define the filename for the input AnnData object
filename = 'integrated_umap'

# Load the preprocessed AnnData object
adata = an.read_h5ad(integrationPath + 'adata_' + filename + '.h5ad')


# Perform differential gene expression analysis
# This ranks genes based on their expression differences across cell types
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key='DEG_celltype')

# Extract differentially expressed genes (DEGs) into a dataframe
DEG_df = sc.get.rank_genes_groups_df(adata, group=None)

# Save the DEGs as an Excel file for further analysis
DEG_df.to_excel(excelPath + 'DEG_integrated.xlsx')

# Define the filename for saving the updated AnnData object
filename = 'adata_integrated_DEG'

# Save the updated AnnData object, now containing DEG results
adata.write(integrationPath + filename + '.h5ad')
