# Load necessary libraries
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(data.table)
library(tibble)
library(ggplot2)
library(yaml)

# Load the path manager utility
source("src/0. utils/path_manager.R")

# Path relative to the project root directory from where the script is executed
config_path <- "src/2. models/config.yaml" 
if (!file.exists(config_path)) {
  stop("Error: config.yaml not found at path: ", config_path)
}
config <- yaml::read_yaml(config_path)

# Extract config parameters
dataset_type <- config$training_parameters$dataset_type 
include_gene_ablation <- config$training_parameters$include_gene_ablation 
include_cell_type_ablation <- config$training_parameters$include_cell_type_ablation

# Validate dataset_type
if (!dataset_type %in% c("absolute", "relative", "raw")) {
  stop("Invalid dataset_type in config.yaml. Must be 'absolute', 'relative', or 'raw'.")
}

# Use path_manager to get dataset paths
paths_data <- get_dataset_paths(
  dataset_type, 
  include_cell_type_ablation, 
  include_gene_ablation
)

# Extract needed data from paths_data
dataset_paths <- paths_data$paths
dataset_names_to_process <- paths_data$dataset_names
versions_to_process <- paths_data$versions

# Define base output directory
output_base_dir <- file.path("output", "2. models", dataset_type, "elasticnet") 
summary_dir <- output_base_dir 

# Create base output directory for the model type
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

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
    
    # Get path key for the current dataset
    dataset_key <- paste(dataset_name, version, sep = "_")
    
    # Skip if this dataset combination doesn't exist
    if (!dataset_key %in% names(dataset_paths)) {
      warning(paste("Skipping:", dataset_name, version, "- Dataset path not found."))
      next
    }
    
    # Get paths for this dataset
    current_paths <- dataset_paths[[dataset_key]]
    
    # Create output directory using path_manager
    current_output_dir <- create_model_output_dir(dataset_type, "elasticnet", dataset_name, version)
    
    # Get output file paths
    output_paths <- create_model_output_paths(current_output_dir, dataset_name, version, "elasticnet")
    
    # Load datasets
    X_train <- fread(current_paths$X_train)
    X_test <- fread(current_paths$X_test)
    y_train <- fread(current_paths$y_train)
    y_test <- fread(current_paths$y_test)

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
    
    # Convert to matrix format for glmnet
    X_train_matrix <- as.matrix(X_train)
    X_test_matrix <- as.matrix(X_test)
    
    # Set up cross-validation for alpha and lambda
    alpha_values <- c(0.01, 0.05, 0.1, 0.3, 0.5)
    
    # Create grid of alpha and lambda values for tuning
    param_grid <- expand.grid(
      alpha = alpha_values,
      lambda = 10^seq(-5, 3, length.out = 100)
    )
    
    # Set up cross-validation
    set.seed(42)
    cv_folds <- createFolds(y_train_numeric, k = 5, returnTrain = TRUE)
    
    # Tune parameters
    print("Tuning elasticnet parameters...")
    best_cvm <- Inf
    best_alpha <- NA
    best_lambda <- NA
    
    for (alpha in alpha_values) {
      # Fit cross-validated elastic net model
      cv_fit <- cv.glmnet(
        x = X_train_matrix, 
        y = y_train_numeric, 
        family = "binomial", 
        alpha = alpha,
        foldid = sample(1:5, size = length(y_train_numeric), replace = TRUE)
      )
      
      # Check if this alpha gives better results
      if (min(cv_fit$cvm) < best_cvm) {
        best_cvm <- min(cv_fit$cvm)
        best_alpha <- alpha
        best_lambda <- cv_fit$lambda.min
      }
    }
    
    print(paste("Best alpha:", best_alpha, "| Best lambda:", best_lambda))
    
    # Train the final model with tuned parameters
    elasticnet_model <- glmnet(
      x = X_train_matrix, 
      y = y_train_numeric, 
      family = "binomial", 
      alpha = best_alpha, 
      lambda = best_lambda
    )
    
    # Save the model
    saveRDS(elasticnet_model, output_paths$model)
      # Make predictions on the test set
    predictions_prob <- predict(
      elasticnet_model, 
      newx = X_test_matrix, 
      type = "response"
    )
    
    # Ensure predictions_prob is a vector (not a matrix)
    predictions_prob <- as.numeric(predictions_prob)
    
    predictions_class <- ifelse(predictions_prob > 0.5, 1, 0)
      # Performance metrics
    conf_matrix <- confusionMatrix(
      factor(predictions_class, levels = c(0, 1)), 
      factor(y_test_numeric, levels = c(0, 1))
    )
    
    accuracy <- conf_matrix$overall["Accuracy"]
    precision <- conf_matrix$byClass["Precision"]
    recall <- conf_matrix$byClass["Sensitivity"]
    f1_score <- 2 * (precision * recall) / (precision + recall)
    roc_obj <- roc(y_test_numeric, as.numeric(predictions_prob))
    roc_auc <- as.numeric(roc_obj$auc)
    
    # Print performance metrics
    print(paste("Accuracy:", round(accuracy, 4)))
    print(paste("Precision:", round(precision, 4)))
    print(paste("Recall:", round(recall, 4)))
    print(paste("F1 Score:", round(f1_score, 4)))
    print(paste("ROC AUC:", round(roc_auc, 4)))
    
    # Calculate feature importance (coefficient magnitudes)
    coef_matrix <- as.matrix(coef(elasticnet_model))
    feature_importance <- data.frame(
      Feature = rownames(coef_matrix),
      Value = coef_matrix[,1]
    ) %>%
      filter(Feature != "(Intercept)") %>%
      arrange(desc(abs(Value)))
    
    # Save feature importance
    write.csv(feature_importance, output_paths$feature_importance, row.names = FALSE)
    
    # Save predictions
    predictions_df <- data.frame(
      Actual = y_test_numeric,
      Predicted_Class = predictions_class,
      Predicted_Probability = predictions_prob
    )
    write.csv(predictions_df, output_paths$predictions, row.names = FALSE)
    
    # Add to performance summary
    performance_summary <- performance_summary %>%
      add_row(
        Dataset = paste(dataset_name, version, sep = "_"),
        Accuracy = accuracy,
        Precision = precision,
        Recall = recall,
        F1_score = f1_score,
        ROC_AUC = roc_auc
      )
  }
}

# Save overall performance summary
performance_summary_path <- file.path(summary_dir, "performance_summary.csv")
write.csv(performance_summary, performance_summary_path, row.names = FALSE)
