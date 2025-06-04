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
  "cp10k"
  # "cp10k_ctnorm",
  # "cp10k_ctnorm_relative",
  # "proportion_ctnorm", 
  # "proportion_ctnorm_relative",
  # "sample_depth"
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
NORMALIZATION_METHODS <- c(
                            "cp10k", 
                            "cp10k_ctnorm", 
                            "cp10k_ctnorm_relative", 
                            "proportion_ctnorm", 
                            "proportion_ctnorm_relative", 
                            "sample_depth")

# Load datasets once at the beginning
message("Loading all datasets once for all methods...")
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
table_s1 <- datasets$table_s1
annots_list <- datasets$annots_list
hsa01100_genes <- datasets$hsa01100_genes


# Main preprocessing pipeline function
run_preprocessing_pipeline <- function(normalization_method) {
  message(sprintf("\nStarting %s preprocessing", normalization_method))
  
  # Set output directory
  base_output_dir <- file.path("output", "1. data preprocessing", "training datasets", normalization_method)
  
  # Determine processing approach
  is_relative_method <- grepl("_relative$", normalization_method)
  
  if (is_relative_method) {
    run_relative_processing(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir)
  } else {
    run_global_processing(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir)
  }
  
  message(sprintf("%s preprocessing complete", normalization_method))
}


# Normalize all data then create subsets
run_global_processing <- function(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir) {
  message("Using global normalization approach...")
  
  # Prepare data for normalization
  prepared_data <- prepare_data_for_normalization(normalization_method, lung_ldm, table_s1, DOUBLETS)
  
  # Apply normalization to prepared data
  counts_normalized <- apply_normalization_method(normalization_method, prepared_data, DOUBLETS)
  
  # Apply log1p
  counts_log_transformed <- apply_log_transform(counts_normalized)

  # Clean up intermediate data
  rm(counts_normalized, prepared_data)
  gc()
  
  # Create and save all subsets
  process_all_subsets(counts_log_transformed, CLUSTER_DEFINITIONS, GENE_SET_TYPES, 
                              hsa01100_genes, table_s1, base_output_dir, 
                              is_relative = FALSE, annots_list = annots_list)
}


# Create subsets first then normalize each one
run_relative_processing <- function(normalization_method, lung_ldm, table_s1, annots_list, hsa01100_genes, base_output_dir) {
  message("Using relative normalization approach...")
  
  # Prepare data for filtering
  prepared_data <- prepare_data_for_normalization(normalization_method, lung_ldm, table_s1, DOUBLETS)
  
  # Process each cell type and gene combination with its own normalization
  process_all_subsets(prepared_data, CLUSTER_DEFINITIONS, GENE_SET_TYPES, 
                              hsa01100_genes, table_s1, base_output_dir, is_relative = TRUE, 
                              normalization_method = normalization_method, annots_list = annots_list)
}

# Apply normalization based on method
apply_normalization_method <- function(normalization_method, prepared_data, doublets) {
  if (normalization_method == "cp10k") {
    if (prepared_data$type == "cell_level") {
      counts_normalized <- normalize_cp10k(prepared_data$cell_metadata, prepared_data$umitab, NULL, doublets)
    } else {
      stop("CP10K normalization requires cell_level data")
    }
    
  } else if (normalization_method == "cp10k_ctnorm") {
    if (prepared_data$type == "cell_level") {
      subset_normalized <- normalize_cp10k(prepared_data$cell_metadata, prepared_data$umitab, NULL, doublets)
      counts_normalized <- apply_celltype_normalization(subset_normalized)
    } else {
      stop("CP10K cell type normalization requires cell_level data")
    }
    
  } else if (normalization_method == "proportion_ctnorm") {
    if (prepared_data$type == "cell_type_level") {
      counts_normalized <- apply_proportion_celltype_normalization(prepared_data$counts_long)
    } else {
      stop("Proportion cell type normalization requires cell_type_level data")
    }
    
  } else if (normalization_method == "sample_depth") {
    if (prepared_data$type == "cell_type_level") {
      counts_normalized <- normalize_sample_depth(prepared_data$counts_long)
    } else {
      stop("Sample depth normalization requires cell_type_level data")
    }
    
  } else {
    stop(paste("Unsupported normalization method:", normalization_method))
  }
  
  return(counts_normalized)
}


# Check if we support the normalization methods
invalid_methods <- setdiff(METHODS_TO_RUN, NORMALIZATION_METHODS)
if (length(invalid_methods) > 0) {
  stop(paste("Invalid normalization requested:", paste(invalid_methods, collapse = ", "),
             "\nAvailable methods:", paste(NORMALIZATION_METHODS, collapse = ", ")))
}

# Run each normalization method
for (method in METHODS_TO_RUN) {
  tryCatch({
    run_preprocessing_pipeline(method)
  }, error = function(e) {
    message(sprintf("Error in %s preprocessing: %s", method, e$message))
    message("Continuing with next method...")
  })
}

message("\nPreprocessing pipeline complete.")
