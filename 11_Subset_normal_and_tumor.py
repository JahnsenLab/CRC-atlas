"""
Script for subsetting tumor and normal cells from the integrated AnnData object.

This script loads the integrated AnnData object with metadata generated in 10_Add_metadata_and_cell_type_annotation_to_main_object.py,
subsets the data into tumor and normal cell populations based on the 'SampleType' column,and saves these subsets as separate AnnData objects for further analysis.

"""

######### Load the required packages ######### 
import anndata as an

######### Specify the paths needed #########

# Path to the directory where integrated .h5ad files are stored and objects will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"

######### Code #########

# Load the main AnnData object with metadata
adata_main = an.read_h5ad(integrationPath + 'adata_with_metadata.h5ad')

# Subset only the tumor cells from the original dataset
adata_T = adata_main[adata_main.obs['SampleType'].isin(["Tumor"])]
# Save the tumor-specific data to a new file
adata_T.write(integrationPath + 'adata_Tumor.h5ad')

# Subset only the normal cells from the original dataset
adata_N = adata_main[adata_main.obs['SampleType'].isin(["Normal"])]
# Save the normal-specific data to a new file
adata_N.write(integrationPath + 'adata_Normal.h5ad')
