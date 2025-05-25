# Utility functions for retrieving and organizing dataset paths

get_dataset_paths <- function(dataset_type,
                              include_cell_type_ablation = TRUE,
                              include_gene_ablation = TRUE,
                              base_path = NULL) {
  
  defined_dataset_types <- c("absolute", "relative", "raw")

  if (dataset_type == "all") {
    types_to_process <- defined_dataset_types
  } else if (dataset_type %in% defined_dataset_types) {
    types_to_process <- c(dataset_type)
  } else {
    valid_options_message <- paste(c(defined_dataset_types, "all"), collapse = "', '")
    stop(paste0("Invalid dataset_type. Must be one of: '", valid_options_message, "'. Provided: '", dataset_type, "'."))
  }
  
  # Determine the root path for dataset_type specific folders
  effective_root_path <- if (is.null(base_path)) {
    file.path("output", "1. data preprocessing", "training datasets")
  } else {
    base_path 
  }
  
  # Define all possible dataset names and versions
  all_dataset_names_const <- c("all_clusters", "lcam_hi", "lcam_lo", "lcam_both", "macrophages")
  all_versions_const <- c("metabolic", "nonmetabolic", "random")
  
  # Filter dataset names and versions based on ablation parameters
  dataset_names_to_process <- if (include_cell_type_ablation) {
    all_dataset_names_const
  } else {
    c("all_clusters")
  }
  
  versions_to_process <- if (include_gene_ablation) {
    all_versions_const
  } else {
    c("metabolic")
  }
  
  final_all_paths <- list() # Accumulator for paths from all processed types
  
  for (current_processing_type in types_to_process) {
    
    current_type_specific_base_path <- file.path(effective_root_path, current_processing_type)
    
    # Create all combinations of dataset paths for the current_processing_type
    dataset_combinations <- expand.grid(
      dataset_name = dataset_names_to_process,
      version = versions_to_process,
      stringsAsFactors = FALSE
    )
    
    for (i in 1:nrow(dataset_combinations)) {
      dataset_name_iter <- dataset_combinations$dataset_name[i]
      version_iter <- dataset_combinations$version[i]
      
      path_set <- list(
        X_train = file.path(current_type_specific_base_path, dataset_name_iter, version_iter, paste0("X_train_", dataset_name_iter, "_", version_iter, ".csv")),
        X_test = file.path(current_type_specific_base_path, dataset_name_iter, version_iter, paste0("X_test_", dataset_name_iter, "_", version_iter, ".csv")),
        y_train = file.path(current_type_specific_base_path, dataset_name_iter, version_iter, paste0("y_train_", dataset_name_iter, "_", version_iter, ".csv")),
        y_test = file.path(current_type_specific_base_path, dataset_name_iter, version_iter, paste0("y_test_", dataset_name_iter, "_", version_iter, ".csv")),
        metadata = file.path(current_type_specific_base_path, dataset_name_iter, version_iter, paste0("metadata_", dataset_name_iter, "_", version_iter, ".csv"))
      )
      
      # Check if all essential files exist
      files_exist <- all(sapply(path_set[c("X_train", "X_test", "y_train", "y_test")], file.exists))
      
      if (files_exist) {
        path_set$dataset_name <- dataset_name_iter
        path_set$version <- version_iter
        path_set$dataset_type_processed <- current_processing_type # Store the actual type processed
        path_set$dataset_key <- paste(current_processing_type, dataset_name_iter, version_iter, sep = "_")
        final_all_paths[[path_set$dataset_key]] <- path_set
      }
    }
  } # End loop over types_to_process
  
  return(list(
    paths = final_all_paths,
    dataset_names = dataset_names_to_process,
    versions = versions_to_process
  ))
}


create_model_output_dir <- function(dataset_type, model_name, dataset_name, version) {
  output_dir <- file.path("output", "2. models", dataset_type, model_name, version, dataset_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  return(output_dir)
}


create_model_output_paths <- function(output_dir, dataset_name, version, model_name) {
  list(
    model = file.path(output_dir, paste0(model_name, "_model_", dataset_name, "_", version, ".rds")),
    feature_importance = file.path(output_dir, paste0("feature_importance_", dataset_name, "_", version, ".csv")),
    predictions = file.path(output_dir, paste0("predictions_", dataset_name, "_", version, ".csv")),
    feature_plot = file.path(output_dir, paste0("feature_importance_plot_", dataset_name, "_", version, ".png"))
  )
}


get_model_results <- function(dataset_type, cell_type, gene_type) {
  # Base path for models
  models_dir <- file.path("output", "2. models", dataset_type)
  
  # List all model directories
  model_dirs <- list.dirs(models_dir, full.names = TRUE, recursive = FALSE)
  
  results <- list()
  
  for (model_dir in model_dirs) {
    model_name <- basename(model_dir)
    
    # Construct path to feature importance file
    feature_file_path <- file.path(model_dir, gene_type, cell_type, 
                                  paste0("feature_importance_", cell_type, "_", gene_type, ".csv"))
    
    if (file.exists(feature_file_path)) {
      results[[model_name]] <- feature_file_path
    }
  }
  
  return(results)
}
