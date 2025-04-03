"""
Script for integrating multiple .h5ad files into a single AnnData object.

This script loads multiple .h5ad files from a specified directory, concatenates them into a single AnnData object, 
performs basic preprocessing, ensuring unique variable names and filtering cells and genes, and saves 
the concatenated object to a new .h5ad file.

"""

######### Load the required packages ######### 
import os
import scanpy as sc
import anndata as ad

######### Specify the paths needed #########

# Path to the directory containing .h5ad files for integration
paperObjPath = "/Path/To/Directory/Containing/h5ad_Files/"

# Path to the directory where the final integrated object will be saved
integrationPath = "/Path/To/Save/Integrated/h5ad_Object/"


######### Code ######### 

# Set the name of the file that will store the analysis results
filename = 'pre-integration'

# Initialize an empty list to store AnnData objects
adata_list = []

# Iterate over the files in the directory specified by paperObjPath
# Read each .h5ad file into an AnnData object and append the AnnData object to the list
substring = '.h5ad'
for x in os.listdir(paperObjPath):
    if substring in x:
        adata = an.read_h5ad(paperObjPath + x)  # Read the .h5ad file into an AnnData object
        adata_list.append(adata)  # Add the AnnData object to the list

# Concatenate all the AnnData objects in the list into one large object
adata = an.concat(adata_list, join = "outer")

# Make sure the variable names in the AnnData object are unique
# This avoids conflicts when variables have the same name across different datasets
adata.var_names_make_unique()

# Perform basic filtering of the data:
    # Remove cells with fewer than 200 genes detected
    # Remove genes that are detected in fewer than 10 cells
sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 10)

# Save the filtered AnnData object
adata.write_h5ad(integrationPath + 'adata_' + filename + "_Concatenated.h5ad")
