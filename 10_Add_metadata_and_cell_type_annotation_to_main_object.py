"""
Script for merging per-cell-type annotations and metadata into the main AnnData object.

This script loads the main AnnData object generated in 07_Manual_annotation.py,
reads in per-cell-type annotation files generated in 09_Make_annotation_files_per_cell_type.py, 
merges them into a unified DataFrame, and adds the annotations (including Leiden subclusters and cell types) 
back into the main AnnData object. It also incorporates external metadata based on the PatientID column 
and saves the updated AnnData object with the added annotations and metadata.

"""

######### Load the required packages ######### 
import os
import pandas as pd # version 1.26.6
import scanpy as sc # version 1.9.8
import anndata as an # version 0.10.5.post1

######### Specify the paths needed #########

# Path to the directory where figures and results will be saved
resultsPath = "/Path/To/Save/Plots/And/Figures/"

# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"

# Path to the directory where the annotation files are saved 
annotationPath = "Path/To/Annotation/Files/"

# Path to the directory where the metadata file is saved 
metadataPath = "Path/To/Metadata/File/"  

######### Code #########

# Read in the main object
adata = an.read_h5ad(integrationPath + 'adata_Manual_annotation.h5ad')

### Merge all per-cell-type annotation files ###

# Ensure the annotation directory exists
if not os.path.exists(annotationPath):
    raise FileNotFoundError(f"Annotation path {annotationPath} does not exist.")

# Get a list of annotation files in the directory
file_list = os.listdir(annotationPath)

# Initialize an empty list to store DataFrames
anno = []

# Loop through each file and read in the annotation data
for file in file_list:
    if file.endswith('.csv') and 'annotation' in file:  # Process only relevant CSV annotation files
        file_path = os.path.join(annotationPath, file)
        df = pd.read_csv(file_path)
        print(file + ' n_cells: ' + str(len(df)))  # Print file name and number of cells
        anno.append(df)  # Append DataFrame to list

# Concatenate all annotation DataFrames into one
annotation_df = pd.concat(anno, ignore_index=True)

# Rename 'leiden' column to 'leiden_subclusters' for clarity
annotation_df = annotation_df.rename(columns={"leiden": "leiden_subclusters"})

# Ensure 'Cell_ids' is set as the index for merging
annotation_df = annotation_df.set_index('Cell_ids')

### ADD FINAL ANNOTATION TO MAIN OBJECT ### 

# Convert adata.obs to a DataFrame for merging
adata_df = adata.obs.copy()

# Merge annotations into the main object based on cell IDs
merged_df = adata_df.merge(annotation_df[['leiden_subclusters', 'Cell_type']],
                           how='left', left_index=True, right_index=True)

# Add the new annotation columns back to the main object
adata.obs['leiden_subclusters'] = merged_df['leiden_subclusters']
adata.obs['Cell_type'] = merged_df['Cell_type']


# Read in the cleaned metadata from an Excel file (e.g. Table_4_CRC_paper.xlsx)
metadata = pd.read_excel(metadataPath + 'metadata_file.xlsx')

# Merge the metadata with the existing adata.obs based on the "SampleID" column
new_adata_obs = adata.obs.merge(metadata, left_on="PatientID", right_on="PatientID", how="left")

# Ensure that the cell barcodes remain intact by resetting the index to the original one
new_adata_obs.index = adata.obs.index
adata.obs = new_adata_obs  # Assign the updated metadata to the AnnData object

# Save the updated AnnData object
filename = 'with_metadata'
adata.write(integrationPath + 'adata_' + filename + '.h5ad')
