library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(scales)

# Load utilities
source("src/1. data preprocessing/training datasets/data_utils.R")
source("src/1. data preprocessing/training datasets/normalization_utils.R")
source("src/1. data preprocessing/training datasets/subset_utils.R")


# Available methods: cp10k, cp10k_ctnorm, absolute, raw
METHODS_TO_RUN <- c(
  "cp10k",
  "cp10k_ctnorm", 
  "absolute",
  "raw"
)

# Configuration constants
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
GENE_SET_TYPES <- c("metabolic", "nonmetabolic", "random")
CLUSTER_DEFINITIONS <- list(
  lcam_hi = c(44, 9, 17, 28, 46, 11, 42),
  lcam_lo = c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16),
  lcam_both = c(c(44, 9, 17, 28, 46, 11, 42), c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16)),
  macrophages = c(5, 8, 10, 11, 25, 32, 33, 35, 38, 42, 47, 54, 55, 57),
  all_clusters = NULL
)

# Supported normalization methods
NORMALIZATION_METHODS <- c("cp10k", "cp10k_ctnorm", "absolute", "raw")

# Main preprocessing pipeline function
run_preprocessing_pipeline <- function(normalization_method) {
  message(sprintf("\nStarting %s preprocessing", normalization_method))
  
  # Set output directory
  base_output_dir <- file.path("output", "1. data preprocessing", "training datasets", normalization_method)
  
  # Load all datasets
  datasets <- load_all_datasets()
  lung_ldm <- datasets$lung_ldm
  table_s1 <- datasets$table_s1
  annots_list <- datasets$annots_list
  hsa01100_genes <- datasets$hsa01100_genes
  
  # Branch based on normalization method
  if (normalization_method == "cp10k") {
    # Prepare cell metadata for CP10K normalization
    cell_metadata_filtered_df <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)
    
    # Filter UMI table based on cell metadata
    umitab_filtered <- filter_umitab(lung_ldm$dataset$umitab, cell_metadata_filtered_df)
    
    # Clear redundant leader datasets
    rm(lung_ldm)
    gc()
    
    # Apply CP10K normalization (includes aggregation)
    counts_normalized <- normalize_cp10k(cell_metadata_filtered_df, umitab_filtered, annots_list, DOUBLETS)
    
    # Clear intermediate objects
    rm(cell_metadata_filtered_df, umitab_filtered)
    gc()
    
  } else if (normalization_method == "cp10k_ctnorm") {
    # Prepare cell metadata for CP10K with cell type normalization
    cell_metadata_filtered_df <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)
    
    # Filter UMI table based on cell metadata
    umitab_filtered <- filter_umitab(lung_ldm$dataset$umitab, cell_metadata_filtered_df)
    
    # Clear redundant leader datasets
    rm(lung_ldm)
    gc()
    
    # Apply CP10K with cell type normalization (includes aggregation)
    counts_normalized <- normalize_cp10k_celltype(cell_metadata_filtered_df, umitab_filtered, annots_list, DOUBLETS)
    
    # Clear intermediate objects
    rm(cell_metadata_filtered_df, umitab_filtered)
    gc()
    
  } else if (normalization_method == "absolute") {
    # Array-based processing for absolute normalization
    counts_filtered <- prepare_counts_array(lung_ldm, table_s1, DOUBLETS)
    
    # Clear lung_ldm after array preparation
    rm(lung_ldm)
    gc()
    
    counts_long <- array_to_long(counts_filtered)
    
    # Clear array after conversion
    rm(counts_filtered)
    gc()
    
    # Apply absolute normalization
    counts_normalized <- normalize_absolute(counts_long)
    
    # Clear long format data
    rm(counts_long)
    gc()
    
  } else if (normalization_method == "raw") {
    # Array-based processing for raw normalization
    counts_filtered <- prepare_counts_array(lung_ldm, table_s1, DOUBLETS)
    
    # Clear lung_ldm after array preparation
    rm(lung_ldm)
    gc()
    
    counts_long <- array_to_long(counts_filtered)
    
    # Clear array after conversion
    rm(counts_filtered)
    gc()
    
    # Apply raw normalization
    counts_normalized <- normalize_raw(counts_long)
    
    # Clear long format data
    rm(counts_long)
    gc()
    
  } else {
    stop(paste("Unsupported normalization method:", normalization_method))
  }
  
  # Common post-processing for all methods
  counts_log_transformed <- apply_log_transform(counts_normalized)
  gc()
  
  # Create feature identifiers
  counts_processed_long <- create_features(counts_log_transformed, annots_list)
  gc()
  
  # This creates all ablation groups
  process_all_subsets(counts_processed_long, 
                      CLUSTER_DEFINITIONS,
                      GENE_SET_TYPES,
                      hsa01100_genes,
                      table_s1,
                      base_output_dir)
  
  message(sprintf("%s preprocessing complete", normalization_method))
}


# Main execution
message("Starting preprocessing pipeline.")
message("Methods to run:", paste(METHODS_TO_RUN, collapse = ", "))

# Check if we support the normalization methods
invalid_methods <- setdiff(METHODS_TO_RUN, NORMALIZATION_METHODS)
if (length(invalid_methods) > 0) {
  stop(paste("Invalid normalization requested:", paste(invalid_methods, collapse = ", "),
             "\nAvailable methods:", paste(NORMALIZATION_METHODS, collapse = ", ")))
}

# Run each requested method
for (method in METHODS_TO_RUN) {
  tryCatch({
    run_preprocessing_pipeline(method)
  }, error = function(e) {
    message(sprintf("Error in %s preprocessing: %s", method, e$message))
    message("Continuing with next method...")
  })
}

message("\nPreprocessing pipeline complete.")
