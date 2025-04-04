"""
Script for generating four prognostic groups based on CD8 T cell and macrophage densities.

This script merges the high/low groupings for CD8 T cells and macrophages for each patient,
assigns a combined prognostic group based on their status in both cell types,and exports the results to Excel files.
It also calculates the number of patients in each prognostic group and provides a summary.

"""

######### Load the required packages ######### 
import os
import pandas as pd

######### Specify the paths needed #########

# Path to the folder where the Excel files per prognostic group are stored and will be saved
groupExcelPath = '/Path/To/Grouping_folders/Excel_files/'  


######### Code #########

# Specigy the sample type analyzed 
sample = 'Tumor'

### Load the files containing the grouping for each cell type ###
# Read in the grouping data for T_NK_ILCs and Macrophages cell types from Excel files
CD8Tcells_df = pd.read_excel(groupExcelPath + 'T_NK_ILCs/' + 'T_NK_ILCs_' + sample + '_grouping.xlsx')
macrophages_df = pd.read_excel(groupExcelPath + 'Macrophages/' + 'Macrophages_' + sample + '_grouping.xlsx')

# Merge the two DataFrames on 'PatientID', creating new columns with suffixes '_T' for CD8 T cells and '_M' for Macrophages
merged_df = pd.merge(CD8Tcells_df, macrophages_df, on='PatientID', suffixes=('_T', '_M'))

# A function that creates a new column 'Prognostic_group' based on combinations of high/low grouping for CD8 T cells and Macrophages
def get_group(row):
    # Check different combinations of 'Grouping_T' (CD8 T cells) and 'Grouping_M' (Macrophages) to assign the 'Prognostic_group'
    if row['Grouping_T'] == 'high' and row['Grouping_M'] == 'low':
        return 'CD8hiTAMlow'  # CD8 high, Macrophages low
    elif row['Grouping_T'] == 'high' and row['Grouping_M'] == 'high':
        return 'CD8hiTAMhi' # CD8 high, Macrophages high
    elif row['Grouping_T'] == 'low' and row['Grouping_M'] == 'low':
        return 'CD8lowTAMlow'# CD8 low, Macrophages low
    elif row['Grouping_T'] == 'low' and row['Grouping_M'] == 'high':
        return 'CD8lowTAMhi' # CD8 low, Macrophages high
    else:
        return 'Unknown'  # Any unexpected combination

# Apply the function to the merged DataFrame to create a new column 'Prognostic_group'
merged_df['Prognostic_group'] = merged_df.apply(get_group, axis=1)

# Create a new DataFrame that only contains the 'PatientID' and the new 'Prognostic_group' columns
merged_df_grouping = merged_df[['PatientID', 'Prognostic_group']]

# Save the new DataFrame to an Excel file, which will contain the prognostic groups
merged_df_grouping.to_excel(groupExcelPath + 'T_NK_ILCs_Macrophages_' + sample + '_prognostic_grouping.xlsx', index=False)

# Calculate the number of patients in each prognostic group
group_counts = merged_df['Prognostic_group'].value_counts()

# Print the number of patients in each prognostic group
print("Number of patients in each group:")
for group, count in group_counts.items():
    print(f"{group}: {count}")
