library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(scales)

# Load utilities
source("src/1. data preprocessing/training datasets/data_utils.R")
source("src/1. data preprocessing/training datasets/normalization_utils.R")
source("src/1. data preprocessing/training datasets/subset_utils.R")

# Available methods: cp10k, cp10k_ctnorm, cp10k_ctnorm_relative, proportion_ctnorm, proportion_ctnorm_relative, sample_depth
METHODS_TO_RUN <- c(
  "cp10k",
  "cp10k_ctnorm",
  "cp10k_ctnorm_relative",
  "proportion_ctnorm", 
  "proportion_ctnorm_relative",
  "sample_depth"
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
NORMALIZATION_METHODS <- c("cp10k", 
                            "cp10k_ctnorm", 
                            "cp10k_ctnorm_relative", 
                            "proportion_ctnorm", 
                            "proportion_ctnorm_relative", 
                            "sample_depth")

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
  
  # Determine processing approach
  is_relative_method <- grepl("_relative$", normalization_method)
  
  if (is_relative_method) {
    run_relative_approach(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir)
  } else {
    run_global_approach(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir)
  }
  
  message(sprintf("%s preprocessing complete", normalization_method))
}

# Normalize all data then create subsets
run_global_approach <- function(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir) {
  message("Using global normalization approach...")
  
  # Apply normalization to all data
  counts_normalized <- apply_normalization_method(normalization_method, lung_ldm, table_s1, DOUBLETS)
  
  # Apply common post processing
  counts_processed_long <- apply_global_postprocessing(counts_normalized, annots_list)
  
  # Create and save all subsets using the subset processor
  process_all_subsets_unified(counts_processed_long, CLUSTER_DEFINITIONS, GENE_SET_TYPES, 
                              hsa01100_genes, table_s1, base_output_dir, is_relative = FALSE)
}

# Create subsets first then normalize each one
run_relative_approach <- function(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir) {
  message("Using relative normalization approach...")
  
  # Prepare base data for filtering
  base_data <- prepare_base_data_for_relative(normalization_method, lung_ldm, table_s1, DOUBLETS)
  
  # Process each cell type and gene combination with its own normalization
  process_all_subsets_unified(base_data, CLUSTER_DEFINITIONS, GENE_SET_TYPES, 
                              hsa01100_genes, table_s1, base_output_dir, is_relative = TRUE, 
                              normalization_method = normalization_method, annots_list = annots_list)
}

# Aply normalization based on method
apply_normalization_method <- function(normalization_method, lung_ldm, table_s1, doublets) {
  if (normalization_method == "cp10k") {
    cell_metadata_filtered_df <- prepare_cell_metadata(lung_ldm, table_s1, doublets)
    umitab_filtered <- filter_umitab(lung_ldm$dataset$umitab, cell_metadata_filtered_df)
    counts_normalized <- normalize_cp10k(cell_metadata_filtered_df, umitab_filtered, NULL, doublets)
    rm(lung_ldm, cell_metadata_filtered_df, umitab_filtered)
    gc()
    
  } else if (normalization_method == "cp10k_ctnorm") {
    cell_metadata_filtered_df <- prepare_cell_metadata(lung_ldm, table_s1, doublets)
    umitab_filtered <- filter_umitab(lung_ldm$dataset$umitab, cell_metadata_filtered_df)
    counts_normalized <- normalize_cp10k_celltype(cell_metadata_filtered_df, umitab_filtered, NULL, doublets)
    rm(lung_ldm, cell_metadata_filtered_df, umitab_filtered)
    gc()
    
  } else if (normalization_method == "proportion_ctnorm") {
    counts_filtered <- prepare_counts_array(lung_ldm, table_s1, doublets)
    counts_long <- array_to_long(counts_filtered)
    counts_normalized <- apply_proportion_celltype_normalization(counts_long)
    rm(lung_ldm, counts_filtered, counts_long)
    gc()
    
  } else if (normalization_method == "sample_depth") {
    counts_filtered <- prepare_counts_array(lung_ldm, table_s1, doublets)
    counts_long <- array_to_long(counts_filtered)
    counts_normalized <- normalize_sample_depth(counts_long)
    rm(lung_ldm, counts_filtered, counts_long); gc()
    
  } else {
    stop(paste("Unsupported normalization method:", normalization_method))
  }
  
  return(counts_normalized)
}

# Prepare base data for relative processing
prepare_base_data_for_relative <- function(normalization_method, lung_ldm, table_s1, doublets) {
  base_normalization_method <- gsub("_relative$", "", normalization_method)
  
  if (grepl("^cp10k", base_normalization_method)) {
    # Use raw umitab data for CP10K methods
    cell_metadata_filtered_df <- prepare_cell_metadata(lung_ldm, table_s1, doublets)
    umitab_filtered <- filter_umitab(lung_ldm$dataset$umitab, cell_metadata_filtered_df)
    result <- list(
      type = "cell_level",
      cell_metadata = cell_metadata_filtered_df,
      umitab = umitab_filtered
    )
    rm(lung_ldm); gc()
    
  } else {
    # Use aggregated counts per cell type data
    counts_filtered <- prepare_counts_array(lung_ldm, table_s1, doublets)
    counts_long <- array_to_long(counts_filtered)
    result <- list(
      type = "aggregated",
      counts_long = counts_long
    )
    rm(lung_ldm, counts_filtered); gc()
  }
  
  return(result)
}

# Apply post processing
apply_global_postprocessing <- function(counts_normalized, annots_list) {
  counts_log_transformed <- apply_log_transform(counts_normalized)
  counts_processed_long <- create_features(counts_log_transformed, annots_list)
  rm(counts_log_transformed) 
  gc()
  return(counts_processed_long)
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
