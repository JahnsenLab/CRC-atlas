"""
Script for manually annotating cell clusters in integrated single-cell RNA sequencing data.

This script loads the integrated AnnData object, assigns manual cell type annotations based on Leiden clustering 
and differentially expressed genes (DEGs), generates UMAP visualizations, and saves the annotated dataset.

"""

######### Load the required packages ######### 
import os
import openpyxl # version 3.1.2
import numpy as np # version 1.26.6
import pandas as pd # version 2.2.0
import scanpy as sc # version 1.9.8
import anndata as an # version 0.10.5.post1

######### Specify the paths needed #########

# Path to the directory where figures and results will be saved
resultsPath = "/Path/To/Save/Plots/And/Figures/"

# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"


######### Code #########

# Set figure parameters for visualizations
sc.set_figure_params(figsize=(6, 6), frameon=False)
sc.settings.n_jobs = 2  # Set the number of parallel jobs for computations

# Set verbosity level for logging: 0: Errors, 1: Warnings, 2: Info, 3: Hints
sc.settings.verbosity = 3
sc.logging.print_header()

# Adjust figure settings for better visualization
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Change working directory to where the data is stored
os.chdir(integrationPath)

# Define the file name for the input AnnData object
input_filename = 'adata_integrated_DEG'

# Load the AnnData object 
adata = an.read_h5ad(integrationPath + input_filename + '.h5ad')

# Define the name for manual annotations in the metadata
result_filename = 'Manual_annotation'

# Define cluster annotations based on Leiden clustering results and differentially expressed genes (DEGs)
annotation = {
    "0": "T cells",
    "1": "Epithelial cells",
    "2": "Macrophages",
    "3": "B cells",
    "4": "Fibroblasts",
    "5": "Plasma cells",
    "6": "Endothelial cells",
    "7": "Mast cells",
    "8": "Neutrophil granulocytes",
    "9": "Glia cells",
    "10": "pDCs",
    "11": "Tuft cells"
}

# Map Leiden clusters to their respective annotations
adata.obs[result_filename] = adata.obs.leiden.map(annotation)

# Set the working directory to ensure figures are saved in the correct folder
os.chdir(resultsPath)

# Plot UMAP with the manual annotations
sc.pl.umap(adata, color=[result_filename], frameon=False, save='_integrated_' + result_filename + '.pdf')

# Save the annotated dataset 
adata.write(integrationPath + 'adata_' + result_filename + '.h5ad')
