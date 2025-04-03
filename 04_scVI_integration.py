"""
Script for integrating single-cell RNA sequencing data using scVI.

This script loads the quality-controlled AnnData object generated in 03_QC_and_filtering.py,
sets up and trains an scVI model for batch correction and dimensionality reduction,
extracts the learned latent representation, and saves the integrated dataset.

"""

######### Load the required packages ######### 
import scanpy as sc     
import os               
import scvi 
import adata as an

######### Specify the paths needed #########

# Path to the directory where figures and results will be saved
resultsPath = "/Path/To/Save/Plots/And/Figures/"

# Path to the directory where integrated .h5ad files are stored and will be saved
integrationPath = "/Path/To/Store/Integrated/AnnData/"


######### Code #########

# Setting figure parameters and verbosity
sc.set_figure_params(figsize=(6, 6), frameon=False)  # Setting the figure size and removing frame
sc.settings.n_jobs = 2  # Number of parallel jobs for computations

# Verbosity: controlling log messages (0 = errors, 1 = warnings, 2 = info, 3 = hints)
sc.settings.verbosity = 3
sc.logging.print_header()  # Printing the header of the logging output

# Setting DPI and facecolor for the figures
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Set working directory to save the figures/plots in the correct folder
os.chdir(resultsPath)  

# Load the AnnData object after pre-integration quality control steps
adata = an.read_h5ad(filename=integrationPath + "adata_pre-integration_quality_control.h5ad")

# Setup the scVI model with necessary parameters
# Here we specify the layer for counts, the batch key (SampleID), and the covariates
scvi.model.SCVI.setup_anndata(adata, 
                              layer="counts",  # Layer containing the raw counts
                              batch_key="SampleID",  # Batch covariate
                              categorical_covariate_keys=["SampleID", "Paper"],  # Categorical covariates
                              continuous_covariate_keys=["pct_counts_mt", "total_counts"])  # Continuous covariates

# Initialize the scVI model with the specified number of layers, latent dimensions, and gene likelihood
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

# Train the model
vae.train()

# Get the latent representation (the learned embeddings of cells)
latent = vae.get_latent_representation()

# Print the shape of the latent representation (number of cells by number of latent dimensions)
print(latent.shape)

# Store the latent representation in the AnnData object under the "obsm" attribute
adata.obsm["X_scVI"] = latent  

# Get the normalized expression for each gene using scVI and add it as a new layer
adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size=1e4)

# Define the filename for saving the updated AnnData object
filename = 'integrated'  

# Define the full path for saving the integrated AnnData object
results_file = integrationPath + 'adata_' + filename + '.h5ad'

# Save the updated AnnData object with the latent representation and normalized expression
adata.write(results_file)
