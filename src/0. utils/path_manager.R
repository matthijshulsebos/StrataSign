# Utility functions for retrieving and organizing dataset paths

create_model_output_dir <- function(norm_type, cell_type, gene_subset, model_name) { # Parameters renamed
  # New structure: norm_type/cell_type/gene_subset/model_name
  output_dir <- file.path("output", "2. models", norm_type, cell_type, gene_subset, model_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  return(output_dir)
}


create_model_output_paths <- function(output_dir, cell_type, gene_subset, model_name) { # Parameters renamed
  list(
    model = file.path(output_dir, paste0(model_name, "_model_", cell_type, "_", gene_subset, ".rds")),
    feature_importance = file.path(output_dir, paste0("feature_importance_", cell_type, "_", gene_subset, ".csv")),
    predictions = file.path(output_dir, paste0("predictions_", cell_type, "_", gene_subset, ".csv")),
    feature_plot = file.path(output_dir, paste0("feature_importance_plot_", cell_type, "_", gene_subset, ".png"))
  )
}


get_model_results <- function(norm_type, cell_type, gene_subset) { # Parameters renamed
  # Base path for the specific combination of norm_type, cell_type, and gene_subset
  base_results_dir <- file.path("output", "2. models", norm_type, cell_type, gene_subset)
  
  results <- list()
  
  if (!dir.exists(base_results_dir)) {
    # warning(paste("Base results directory not found:", base_results_dir))
    return(results) # Return empty list if the base path for this combo doesn't exist
  }

  # List all model directories within this specific combination path
  # These directories are named by model_name
  model_name_dirs <- list.dirs(base_results_dir, full.names = TRUE, recursive = FALSE)
  
  for (model_specific_dir in model_name_dirs) {
    model_name <- basename(model_specific_dir)
    
    # Construct path to feature importance file within the model_specific_dir
    # The file name itself still uses cell_type and gene_subset as per create_model_output_paths
    feature_file_path <- file.path(model_specific_dir, 
                                  paste0("feature_importance_", cell_type, "_", gene_subset, ".csv"))
    
    if (file.exists(feature_file_path)) {
      results[[model_name]] <- feature_file_path
    } else {
      # Optional: warning if a model directory exists but the feature file is missing
      # warning(paste("Feature importance file not found for model", model_name, "at", feature_file_path))
    }
  }
  
  return(results)
}
