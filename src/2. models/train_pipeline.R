library(yaml)


# Base path of the training datasets
BASE_DATA_PATH <- file.path("output", "1. data preprocessing", "training datasets")

# Load helper functions for managing file paths 
source(file.path("src", "0. utils", "path_manager.R"))

# Defines the path to the configuration file
config_path <- file.path("src", "2. models", "config.yaml") 

# Read config into a list
config <- yaml::read_yaml(config_path)


# === HELPER FUNCTIONS ===

# Extract names of true flags from config list
extract_active_items <- function(config_list, item_name) {
  if (is.null(config_list) || !is.list(config_list)) {
    stop(paste(item_name, "is missing or not a list in config.yaml under training_parameters."))
  }
  # Filters for items marked as true in the config
  active_items <- names(config_list)[unlist(config_list)] 
  if (length(active_items) == 0) {
    warning(paste("No active", item_name, "found."))
  }
  return(active_items)
}


# Load datasets from dataset directory
load_dataset_csvs <- function(x_train_path, x_test_path, y_train_path, y_test_path) {
  X_train_df <- data.table::fread(x_train_path, check.names = FALSE, data.table = FALSE)
  X_test_df  <- data.table::fread(x_test_path, check.names = FALSE, data.table = FALSE)
  y_train_df <- data.table::fread(y_train_path, data.table = FALSE)
  y_test_df  <- data.table::fread(y_test_path, data.table = FALSE)
  return(list(X_train_df = X_train_df, X_test_df = X_test_df, y_train_df = y_train_df, y_test_df = y_test_df))
}


# Saves various model outputs to specified files
results_to_output <- function(model_results, output_paths_list) { 
  if (is.null(model_results)) {
    warning("Model results are empty.")
    return()
  }

  # Save model object
  if (!is.null(model_results$model_object)) {
    # Ensure directory exists
    dir.create(dirname(output_paths_list$model), recursive = TRUE, showWarnings = FALSE) 
    tryCatch({
      # Save the trained model object
      saveRDS(model_results$model_object, output_paths_list$model) 
    }, error = function(e) {
      warning(paste("Error saving model object:", e$message))
    })
  } else {
    warning("Model object is empty.")
  }

  if (!is.null(model_results$predictions_df)) {
    message(paste("Saving predictions to", output_paths_list$predictions))

    # Ensure directory exists
    dir.create(dirname(output_paths_list$predictions), recursive = TRUE, showWarnings = FALSE) 
    
    # Ensure the required columns are present
    required_prediction_cols <- c("y_test", "y_pred")
    if (!all(required_prediction_cols %in% colnames(model_results$predictions_df))) {
      warning("Predictions output does not contain the required columns.")
    } else {
      write.csv(model_results$predictions_df, output_paths_list$predictions, row.names = FALSE) 
    }
  } else {
    warning("Predictions output or predictions output path is empty.")
  }

  if (!is.null(model_results$raw_feature_importance)) {
    message(paste("Saving feature importance to", output_paths_list$feature_importance))

    # Ensure directory exists
    dir.create(dirname(output_paths_list$feature_importance), recursive = TRUE, showWarnings = FALSE)
    
    if (!all(c("Feature", "Value") %in% colnames(model_results$raw_feature_importance))) {
      warning("Feature importance output does not have the required columns.")
    }
    
    # Save feature importance to CSV
    write.csv(model_results$raw_feature_importance, output_paths_list$feature_importance, row.names = FALSE)
  } else {
    warning("Feature importance or feature importance output path is empty.")
  }
}


# === CORE FUNCTIONS ===

# Function to source all model training scripts from a specified directory
source_all_model_scripts <- function(modeling_dir_path) {
  message(paste("Sourcing all training scripts from:", modeling_dir_path))

  if (dir.exists(modeling_dir_path)) {
    # List all R files ending with in the modeling directory
    r_scripts <- list.files(path = modeling_dir_path, pattern = "\\.[Rr]$", full.names = TRUE, recursive = FALSE)
    
    if (length(r_scripts) > 0) {
      message(paste("Found", length(r_scripts), "training scripts in", modeling_dir_path, ":"))
      for (script_path in r_scripts) {
        source(script_path)
      }
    } else {
      warning("No training scripts found.")
    }
  } else {
    warning("Modeling directory not found.")
  }
}

# Main function orchestrating the training pipeline
run_training_pipeline <- function(config_obj) { 
  message("Starting training pipeline.")
  
  # Source all model training scripts
  modeling_dir_path <- file.path("src", "2. models", "modeling") 
  source_all_model_scripts(modeling_dir_path)

  # Get active normalization types from config
  norm_types <- extract_active_items(config_obj$training_parameters$norm_type, "norm_type") 
  # Get active cell types from config
  cell_types <- extract_active_items(config_obj$training_parameters$cell_type, "cell_type") 
  # Get active gene types from config
  gene_types <- extract_active_items(config_obj$training_parameters$gene_type, "gene_type") 
  # Get list of models to run from config
  models_to_run_list <- config_obj$training_parameters$models 

  # Validate extracted parameters
  if (length(norm_types) == 0) {
    stop("No norm_type selected for processing in config.yaml")
  }
  if (length(cell_types) == 0) {
    stop("No cell_type selected for processing in config.yaml")
  }
  if (length(gene_types) == 0) {
    stop("No gene_type selected for processing in config.yaml")
  }
  if (is.null(models_to_run_list) || length(models_to_run_list) == 0) {
    stop("No models specified or model list is empty in config.yaml")
  }

  # Create a grid of all dataset parameter combinations to iterate over
  iteration_grid <- expand.grid(
    norm_type = norm_types,
    cell_type = cell_types,
    gene_type = gene_types,
    # Ensures character vectors are not converted to factors
    stringsAsFactors = FALSE 
  )

  # Loop through each combination of dataset parameters
  for (i in 1:nrow(iteration_grid)) {
    current_norm_type <- iteration_grid$norm_type[i]
    current_cell_type <- iteration_grid$cell_type[i]
    current_gene_type <- iteration_grid$gene_type[i]

    # Construct an identifier for the current iteration for logging
    current_iter_id <- paste0("Norm=", current_norm_type, "_Cell=", current_cell_type, "_Gene=", current_gene_type)

    # Cleaner message for processing combination
    message(paste0("\nProcessing combination: Normalization method=", current_norm_type, 
                   ", Cell type=", current_cell_type, 
                   ", Gene type=", current_gene_type))

    # Construct path to the specific dataset folder for the current combination
    current_dataset_folder <- file.path(BASE_DATA_PATH, current_norm_type, current_cell_type, current_gene_type)

    # Skip to the next iteration if dataset folder doesn't exist
    if (!dir.exists(current_dataset_folder)) {
      warning(paste("Dataset folder not found:", current_dataset_folder, "for", current_iter_id))
      next 
    }

    # Define paths to individual data files
    X_train_file <- file.path(current_dataset_folder, paste0("X_train_", current_cell_type, "_", current_gene_type, ".csv"))
    X_test_file  <- file.path(current_dataset_folder, paste0("X_test_", current_cell_type, "_", current_gene_type, ".csv"))
    y_train_file <- file.path(current_dataset_folder, paste0("y_train_", current_cell_type, "_", current_gene_type, ".csv"))
    y_test_file  <- file.path(current_dataset_folder, paste0("y_test_", current_cell_type, "_", current_gene_type, ".csv"))

    # Check if all required data files exist for the current combination
    required_files <- c(X_train_file, X_test_file, y_train_file, y_test_file)

    if (!all(sapply(required_files, file.exists))) {
      missing_files_str <- paste(required_files[!sapply(required_files, file.exists)], collapse = ", ")
      stop(paste("One or more data files missing for", current_iter_id, ". Missing:", missing_files_str))
      next 
    }
    
    # Load datasets
    message("Loading datasets for current combination.")
    datasets <- load_dataset_csvs(X_train_file, X_test_file, y_train_file, y_test_file)

    message("Loaded all datasets for current combination.")

    # Loop through each model specified in the configuration
    for (model_entry in models_to_run_list) {
      # Check if the train flag for the model is true
      if (isTRUE(model_entry$train)) { 
        model_name <- model_entry$name
        message(paste0("\nTraining model: ", model_name))

        # Init output variables
        model_results <- NULL
        model_train_successful <- FALSE

        # Train the model
        tryCatch({
          # Construct the model training function name
          model_function_name <- paste0("train_", model_name, "_model") 
          
          # Check if the dynamically constructed function name exists
          if (!exists(model_function_name, mode = "function")) {
            stop(paste("Training function '", model_function_name, "not found for model", model_name, ".")) 
          }
          
          # Get the actual function object
          model_func <- get(model_function_name) 
          
          # Call the model specific training function
          model_results <- model_func(
            X_train_df = datasets$X_train_df,
            X_test_df = datasets$X_test_df,
            y_train_df = datasets$y_train_df,
            y_test_df = datasets$y_test_df
          )
        }, error = function(e) {
          stop(paste("Error training model", model_name, "for", current_iter_id, ":", e$message))
        })
        
        # Create output directory for the current model and dataset combination
        output_dir_rel_to_project_root <- create_model_output_dir(
          norm_type = current_norm_type,
          cell_type = current_cell_type,
          gene_subset = current_gene_type,
          model_name = model_name
        )
        
        # Construct full paths for various output files directly
        model_artifact_path <- file.path(output_dir_rel_to_project_root, paste0(model_name, "_model_", current_cell_type, "_", current_gene_type, ".rds"))
        feature_importance_path <- file.path(output_dir_rel_to_project_root, paste0("feature_importance_", current_cell_type, "_", current_gene_type, ".csv"))
        predictions_path <- file.path(output_dir_rel_to_project_root, paste0("predictions_", current_cell_type, "_", current_gene_type, ".csv"))

        output_paths_list <- list(
          model = model_artifact_path,
          feature_importance = feature_importance_path,
          predictions = predictions_path
        )

        results_to_output( 
          model_results = model_results,
          output_paths_list = output_paths_list
        )
        
        message(paste("Finished model:", model_name))
      } else {
        # Skip training if the train flag is not true
      }
    }
  }
}


# === EXECUTE PIPELINE ===

# Run the training pipeline
run_training_pipeline(
  config = config 
)

message("\nPipeline execution finished.")
