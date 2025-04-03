"""
Script for stratifying patients into prognostic groups based on CD8 T cell and macrophage densities in single-cell RNA sequencing data.

This script loads the AnnData object containing only tumor data generated in 11_Subset_normal_and_tumor.py, subsets the data by cell type,
calculates cell counts and fractions per patient and cluster, and performs normalization of cell counts. It groups patients based on normalized cell count data, 
filters patients based on predefined thresholds, and exports the results to Excel files.

""" 

######### Load the required packages ######### 
import scanpy as an
import os
import pandas as pd
import numpy as np


######### Specify the paths needed #########

# Path to the directory where integrated .h5ad files are stored 
integrationPath = "/Path/To/Store/Integrated/AnnData/"

# Path to the folder where the .h5ad files are located
objectFolder = '/Path/To/h5ad_files/'  

# Path to the folder where the folders per prognostic group will be created
groupingPath = '/Path/To/Grouping_folders/'

# Path to the folder where the Excel files per prognostic group will be saved
groupExcelPath = groupingPath + '/Excel_files/'  

# Path to the folder where the objects per prognostic group will be saved
groupSubObjPath = groupingPath + '/Objects_subfolder/'  

# Path to the folder where Excel files are stored and generated excel files will be saved 
excelPath = '/path/to/excel_files_to_read/'


######### Code #########

# Specigy the sample type analyzed 
sample = 'Tumor' 

# Read the main AnnData object containing only tumor data 
adata_main = an.read_h5ad(integrationPath + 'adata_' + sample + '.h5ad')

# Get the unique patient IDs from the main dataset
patient_ids = adata_main.obs['PatientID'].unique()

# Initialize an empty list to store the results of cell counts per patient
results = []

# Iterate through each patient ID and calculate the number of cells
for patient_id in patient_ids:
    num_cells = len(adata_main[adata_main.obs['PatientID'] == patient_id])
    results.append({'PatientID': patient_id, 'Total_number_of_cells': num_cells})

# Create a DataFrame containing total cells per patient 
total_cells = pd.DataFrame(results)

# Cell types to analyze
celltypesList = ['Macrophages', 'T_NK_ILCs']

# Loop over each cell type in the list
for celltype in celltypesList:
    
    # Read the AnnData object
    adata = an.read_h5ad(objectFolder + 'adata_' + celltype + '_' + sample + '.h5ad')

    # Create directories for the current cell type
    groupingPath_cell = groupingPath + '/' + celltype + '/'
    groupExcelPath_cell = groupExcelPath + '/' + celltype + '/'
    groupSubObjPath_cell = groupSubObjPath + '/' + celltype + '/'

    # Check and create directories if they do not exist
    for path in [groupingPath_cell, groupExcelPath_cell, groupSubObjPath_cell]:
        if not os.path.exists(path):
            os.makedirs(path)

    ### GET THE NUMBER/FRACTION OF CELLS PER PATIENT PER CLUSTER ###
    # Create a DataFrame to store the counts of cells per patient and cluster
    cluster_counts = pd.DataFrame(index=adata.obs['PatientID'].unique(), columns=adata.obs['leiden'].unique())

    # Iterate through patients and clusters to count cells
    for patient in cluster_counts.index:
        for cluster in cluster_counts.columns:
            count = ((adata.obs['PatientID'] == patient) & (adata.obs['leiden'] == cluster)).sum()
            cluster_counts.at[patient, cluster] = count

    # Calculate the total number of cells per patient
    total_cells_per_patient = cluster_counts.sum(axis=1)

    # Calculate the fraction of cells in each cluster for each patient
    fraction_cluster_counts = cluster_counts.divide(total_cells_per_patient, axis=0)

    # Get unique clusters, sort them, and convert to strings for correct ordering
    sorted_clusters = sorted(cluster_counts.columns, key=lambda x: int(x))
    sorted_clusters = [str(cluster) for cluster in sorted_clusters]

    # Reorder the columns (clusters) for both DataFrames
    cluster_counts = cluster_counts[sorted_clusters]
    cluster_counts.insert(0, 'PatientID', cluster_counts.index)

    fraction_cluster_counts = fraction_cluster_counts[sorted_clusters]
    fraction_cluster_counts.insert(0, 'PatientID', fraction_cluster_counts.index)

    # Save to Excel 
    count_file = excelPath + celltype + '_' + sample + '_number_of_cells_per_patient_per_cluster.xlsx'
    cluster_counts.to_excel(count_file, index=False)

    fraction_file = excelPath + celltype + '_' + sample + '_fraction_of_cells_per_patient_per_cluster.xlsx'
    fraction_cluster_counts.to_excel(fraction_file, index=False)

    
    # Sum up the number of cells per patient per cluster 
    if celltype == 'Macrophages':
        
        # Sum up the number of cells per patient for ALL clusters
        clusters = list(cluster_counts) 
        clusters.remove('PatientID')
        cluster_counts['Sum'] = cluster_counts[clusters].sum(axis=1, skipna=True)
        cluster_counts = cluster_counts[['PatientID', 'Sum']]

    if not celltype == 'T_NK_ILCs':
        
        # Sum up the number of cells per patient ONLY FOR CD8+ clusters only
        clusters = ['3', '4', '7']  # clusters for CD8+
        cluster_counts['Sum'] = cluster_counts[clusters].sum(axis=1, skipna=True)
        cluster_counts = cluster_counts[['PatientID', 'Sum']]
        

    # Merge with the total cells DataFrame
    merged_df = pd.merge(cluster_counts, total_cells, on='PatientID')
    merged_df.to_excel(groupExcelPath_cell + celltype + '_' + sample + '_Number_of_cells_per_patient_per_cluster_and_total_cells.xlsx')

    # Filter out patients with 'Total_number_of_cells' < 100 and 'Sum' < 10
    merged_df_filtered = merged_df[(merged_df['Total_number_of_cells'] >= 100) & (merged_df['Sum'] >= 10)]
    merged_df_filtered.to_excel(groupExcelPath_cell + celltype + '_' + sample + '_Number_of_cells_per_patient_per_cluster_and_total_cells_filtered.xlsx')

    ### NORMALIZE THE COUNT VALUES ###
    # Define a normalization function
    def normalize_column(row, column):
        return row[column] / row['Total_number_of_cells']

    # Apply the normalization function
    merged_df_filtered['Sum_normalized'] = merged_df_filtered.apply(lambda row: normalize_column(row, 'Sum'), axis=1)

    # Save the normalized values to Excel
    merged_df_filtered.to_excel(groupExcelPath_cell + celltype + '_' + sample + '_normalized_counts_filtered.xlsx', index=False)

    # Group patients into 'high'/'low' categories based on the median normalized sum
    median = merged_df_filtered['Sum_normalized'].median()

    # Create a new column 'Grouping' with default value 'low'
    merged_df_filtered['Grouping'] = 'low'  # Initially assign 'low' to all rows in the 'Grouping' column

    # Update the 'Grouping' column to 'high' where 'Sum_normalized' is greater than the median value
    merged_df_filtered.loc[merged_df_filtered['Sum_normalized'] > median, 'Grouping'] = 'high'  

    # Save the data frame with the grouping information
    merged_df_filtered.to_excel(groupExcelPath_cell + celltype + '_' + sample + '_grouping.xlsx', index=False)
