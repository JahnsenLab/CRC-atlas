"""
Script for exporting cell annotations for a specific cell type and clustering resolution.

This script loads the AnnData object for a specified cell type and clustering resolution, generated in 08_Subclustering.py,
extracts the Leiden cluster assignments, and creates an annotation DataFrame. The DataFrame includes cell identifiers,
cluster assignments, and the main cell type, and it is then saved as a CSV file for further analysis.

"""

######### Load the required packages ######### 

import os
import scanpy as sc # version 1.9.8        
import anndata as an # version 0.10.5.post1
import pandas as pd # version 2.2.0
import numpy as np # verion 1.26.6

######### Specify the paths needed #########

#  Path to the directory where the per cell type resolution-specific directories are stored
resfolder = "Path/To/Per/Cell/Type/Resolution/Specific/Directory/"

# Path to the directory where the annotation files will be saved 
annotationPath = "Path/To/Store/Annotation/Files/"

######### Code #########

# Define the cell type of interest for further analysis
celltype = "Cell_type_of_interest"

# Set the resolution used in the clustering 
res = 0.5

# Load the cell type and resolution specific AnnData object
adata = an.read_h5ad(resFolder + 'adata_' + celltype + '_res_' + str(res) + '_DEG.h5ad')
        
# Extract Leiden cluster assignments and create a DataFrame
cleanedAnno = adata.obs[['leiden']].copy()  # Create a copy to avoid warnings

# Reset the index to maintain unique cell identifiers
cleanedAnno.reset_index(inplace=True)

# Rename the index column to 'Cell_ids' for clarity
cleanedAnno.rename(columns={'index': 'Cell_ids'}, inplace=True)

# Add a column specifying the main cell type (replace underscores with spaces)
cleanedAnno['Cell_type'] = celltype.replace('_', ' ')

# Define the cluster type
analysed = 'leiden'

# Save the annotation DataFrame as a CSV file
cleanedAnno.to_csv(annotationPath + celltype + '_' + analysed + '_annotation.csv', index=False)
