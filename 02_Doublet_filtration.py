"""
Script for detecting and filtering doublets in single-cell RNA sequencing data using scVI and SOLO.

This script loads the concatenated AnnData object generated in 01_Concatination_of_objects.py,
applies the scVI model for latent representation learning, and uses SOLO to predict doublets.
The doublet predictions are saved, and the dataset is filtered to retain only singlets.

"""

######### Load the required packages ######### 
import os
import scanpy as sc
import anndata as ad

######### Specify the paths needed #########
# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Directory/For/Integration/Results/"

# Path to the directory where technical results (e.g., doublet predictions) will be saved
technicalPath = "/Path/To/Directory/For/Technical/Results/"


######### Code #########
# The file that will store the analysis results
filename = 'pre-integration'

# Read the previously concatenated AnnData object 
adata = sc.read_h5ad(integrationPath + filename + "_Concatenated.h5ad")

# Set up the scVI model for single-cell data analysis
scvi.model.SCVI.setup_anndata(adata)

# Initialize the scVI model using the AnnData object and train it
vae = scvi.model.SCVI(adata)  # Create the scVI model object
vae.train()  # Train the model on the AnnData data

# Identify doublets using the SOLO function from scVI
# SOLO is a method that predicts doublets 
solo = scvi.external.SOLO.from_scvi_model(vae)  # Initialize SOLO with the trained scVI model
solo.train()  # Train SOLO to predict doublets

# Get the predictions from SOLO
# The predictions give the likelihood of doublets in each cell
doubletDF = solo.predict()  # Returns soft predictions with numerical values
doubletDF['prediction'] = solo.predict(soft=False)  # Gives hard predictions (i.e., 'doublet' or 'singlet')

# Save the doublet prediction results as a CSV file
doubletDF.to_csv(technicalPath + 'doublet_data.csv')  

# Save the AnnData object with doublet predictions included
adata.write_h5ad(integrationPath + 'adata_pre-integration_doublet_prediction.h5ad')

# Filter the AnnData object to retain only singlets (non-doublets)
adata = adata[adata.obs.prediction == 'singlet']

# Save the filtered AnnData object (only singlets) to a new file
adata.write_h5ad(integrationPath + 'adata_pre-integration_filtered.h5ad')
