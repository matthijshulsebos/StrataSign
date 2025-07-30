
# Creates a base directory for model outputs
create_model_output_dir <- function(norm_type, cell_type, gene_subset, model_name) {
  # Define output directory path
  output_dir <- file.path("output", "2. models", norm_type, cell_type, gene_subset, model_name)

  # Create the directory if it does not exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  return(output_dir)
}


# Creates paths for model object, feature importance, predictions, and feature plot files
create_model_output_paths <- function(output_dir, cell_type, gene_subset, model_name) {
  list(
    model = file.path(output_dir, paste0(model_name, "_model_", cell_type, "_", gene_subset, ".rds")),
    feature_importance = file.path(output_dir, paste0("feature_importance_", cell_type, "_", gene_subset, ".csv")),
    predictions = file.path(output_dir, paste0("predictions_", cell_type, "_", gene_subset, ".csv"))
  )
}
