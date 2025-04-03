"""
Script for computing UMAP visualization and clustering on integrated single-cell RNA sequencing data.

This script loads the integrated AnnData object generated in 04_scVI_integration.py,
computes a neighborhood graph using the scVI latent representation, performs UMAP dimensionality reduction,
applies Leiden clustering, and saves the updated dataset along with UMAP visualizations.

"""

######### Load the required packages ######### 
import scanpy as sc     
import os               
import adata as an

######### Specify the paths needed #########

# Path to the directory where figures and results will be saved
resultsPath = "/Path/To/Save/Plots/And/Figures/"

# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"


######### Code #########

# Load the integrated AnnData object 
adata = an.read_h5ad(integrationPath + "adata_integrated.h5ad")

# Set the working directory to ensure figures are saved in the correct folder
os.chdir(resultsPath)  

# Compute the neighborhood graph using the latent representation from scVI
sc.pp.neighbors(adata, use_rep="X_scVI")  

# Compute UMAP for visualization
sc.tl.umap(adata)  

# Perform clustering using the Leiden algorithm with resolution 0.1
sc.tl.leiden(adata, resolution=0.1)  

# Define the filename for saving UMAP plots
filename = 'integrated_umap'

# Generate UMAP plots colored by "Paper" (dataset source)
sc.pl.umap(adata, color="Paper",
    frameon=False, save='_' + filename + '_paper.pdf')

# Generate UMAP plots colored by "leiden" (clustering result)
sc.pl.umap(adata, color="leiden",
    frameon=False, save='_' + filename + '_leiden.pdf')

# Define the filename for saving the updated AnnData object
filename = 'integrated_umap'

# Save the updated AnnData object containing UMAP embeddings and Leiden clustering results
adata.write(integrationPath + 'adata_' + filename + '.h5ad')
