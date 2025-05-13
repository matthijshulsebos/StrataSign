# Load necessary libraries
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(data.table)
library(tibble)
library(ggplot2)
library(yaml) # Added

# Path relative to the project root directory from where the script is executed
config_path <- "src/2. models/config.yaml" 
if (!file.exists(config_path)) {
  stop("Error: config.yaml not found at path: ", config_path)
}
config <- yaml::read_yaml(config_path)

dataset_type <- config$training_parameters$dataset_type 
include_gene_ablation <- config$training_parameters$include_gene_ablation 
include_cell_type_ablation <- config$training_parameters$include_cell_type_ablation

if (!dataset_type %in% c("absolute", "relative")) {
  stop("Invalid dataset_type in config.yaml. Must be 'absolute' or 'relative'.")
}

base_path <- file.path("output", "1. data preprocessing", "training datasets", dataset_type) 
output_base_dir <- file.path("output", "2. models", dataset_type, "elasticnet") 
summary_dir <- output_base_dir 

# Create base output directory for the model type
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Define cell types and genes to based on config
all_dataset_names <- c("all_clusters", "lcam_hi", "lcam_lo", "lcam_both")
all_versions <- c("metabolic", "nonmetabolic", "random")

if (include_cell_type_ablation) {
  dataset_names_to_process <- all_dataset_names
} else {
  dataset_names_to_process <- c("all_clusters") 
}

if (include_gene_ablation) {
  versions_to_process <- all_versions
} else {
  versions_to_process <- c("metabolic") 
}

# Function to create file paths
create_dataset_paths <- function(dataset_name, version) {
  list(
    X_train = file.path(base_path, dataset_name, version, paste0("X_train_", dataset_name, "_", version, ".csv")),
    X_test = file.path(base_path, dataset_name, version, paste0("X_test_", dataset_name, "_", version, ".csv")),
    y_train = file.path(base_path, dataset_name, version, paste0("y_train_", dataset_name, "_", version, ".csv")),
    y_test = file.path(base_path, dataset_name, version, paste0("y_test_", dataset_name, "_", version, ".csv")),
    metadata = file.path(base_path, dataset_name, version, paste0("metadata_", dataset_name, "_", version, ".csv"))
  )
}

# Init performance summary
performance_summary <- tibble(
  Dataset = character(),
  Accuracy = numeric(),
  Precision = numeric(),
  Recall = numeric(),
  F1_score = numeric(),
  ROC_AUC = numeric()
)

# Main loop 
for (dataset_name in dataset_names_to_process) {
  for (version in versions_to_process) {
    
    print(paste("Processing dataset:", dataset_name, "| version:", version, "| type:", dataset_type))
    
    paths <- create_dataset_paths(dataset_name, version)
    
    files_exist <- all(sapply(paths[c("X_train", "X_test", "y_train", "y_test")], file.exists)) 
    if (!files_exist) {
      warning(paste("Skipping:", dataset_name, version, "- Essential input file(s) not found. Searched in:", dirname(paths$X_train)))
      next 
    }

    # Create output directory for this specific version and dataset_name combination
    current_output_dir <- file.path(output_base_dir, version, dataset_name) 
    dir.create(current_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Load datasets
    X_train <- fread(paths$X_train)
    X_test <- fread(paths$X_test)
    y_train <- fread(paths$y_train)
    y_test <- fread(paths$y_test)

    # Label encoding
    y_train$x <- factor(y_train$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    y_test$x <- factor(y_test$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    
    # Convert to numeric for model training
    y_train_numeric <- as.numeric(as.character(y_train$x))
    y_test_numeric <- as.numeric(as.character(y_test$x))
    
    # Remove zero-variance columns
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, !nzv$zeroVar, with = FALSE]
    
    # Define alpha values to explore
    alpha_values <- c(0.01, 0.05, 0.1, 0.3, 0.5)
    
    # Initialize variables to store best model
    best_cv_fit <- NULL
    best_alpha <- NULL
    best_cvm <- Inf
    
    # Test different alpha values
    for (alpha in alpha_values) {
      print(paste("Testing alpha =", alpha))
      set.seed(42)
      cv_fit <- cv.glmnet(
        as.matrix(X_train), 
        y_train_numeric,
        alpha = alpha,
        family = "binomial",
        type.measure = "deviance",
        nfolds = 5
      )
      
      if (min(cv_fit$cvm) < best_cvm) {
        best_cvm <- min(cv_fit$cvm)
        best_cv_fit <- cv_fit
        best_alpha <- alpha
      }
    }
    
    print(paste("Best alpha =", best_alpha))
    
    # Make predictions
    y_prob <- predict(best_cv_fit, newx = as.matrix(X_test), s = "lambda.min", type = "response")[,1]
    y_pred <- ifelse(y_prob > 0.5, 1, 0)
    
    # Ensure factors have the same levels for confusion matrix
    y_pred <- factor(y_pred, levels = c(0, 1))
    y_test_factor <- factor(y_test_numeric, levels = c(0, 1))
    
    # Calculate performance metrics
    confusion <- confusionMatrix(y_pred, y_test_factor)
    accuracy <- confusion$overall["Accuracy"]
    precision <- confusion$byClass["Pos Pred Value"]
    recall <- confusion$byClass["Sensitivity"]
    f1_score <- confusion$byClass["F1"]
    roc_auc <- as.numeric(roc(y_test_numeric, y_prob)$auc)
    
    # Store performance metrics
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, sep = "_"),
              Accuracy = accuracy, 
              Precision = precision,
              Recall = recall,
              F1_score = f1_score,
              ROC_AUC = roc_auc)
    
    # Get feature importance (coefficients)
    coef_matrix <- as.matrix(coef(best_cv_fit, s = "lambda.min"))
    feature_importance <- data.frame(
      Feature = rownames(coef_matrix),
      Value = coef_matrix[,1]
    ) %>%
      filter(Feature != "(Intercept)") %>%
      arrange(desc(abs(Value)))
    
    # Save feature importance
    write_csv(
      feature_importance,
      file.path(current_output_dir, 
                paste0("feature_importance_", dataset_name, "_", version, ".csv"))
    )
    
    # Save predictions
    write_csv(
      data.frame(
        y_test = y_test$x,
        y_pred = y_pred
      ),
      file.path(current_output_dir, 
                paste0("predictions_", dataset_name, "_", version, ".csv"))
    )
    
    # Save ElasticNet model
    saveRDS(best_cv_fit, file.path(current_output_dir, 
            paste0("elasticnet_model_", dataset_name, "_", version, ".rds")))
  }
}

# Save performance summary
summary_filename <- file.path(summary_dir, "elasticnet_performance_summary.csv") 
write_csv(performance_summary, summary_filename)
cat("Saved performance summary:", summary_filename, "\n")

print(paste("ElasticNet analysis completed successfully for dataset:", dataset_type))
