# Load necessary libraries
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(ggplot2)

# Define base path and output directories
base_path <- "output/training_data/raw"
models_dir <- "output/models/elasticnet"
predictions_dir <- "output/models/elasticnet"
feature_importance_dir <- "output/models/elasticnet"
summary_dir <- "output/models/elasticnet"

# Create directories if they don't exist
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Simplified dataset path function - no conditionals needed
create_dataset_paths <- function(dataset_name, version) {
  list(
    X_train = file.path(base_path, dataset_name, version, paste0("X_train_", dataset_name, "_", version, ".csv")),
    X_test = file.path(base_path, dataset_name, version, paste0("X_test_", dataset_name, "_", version, ".csv")),
    y_train = file.path(base_path, dataset_name, version, paste0("y_train_", dataset_name, "_", version, ".csv")),
    y_test = file.path(base_path, dataset_name, version, paste0("y_test_", dataset_name, "_", version, ".csv")),
    metadata = file.path(base_path, dataset_name, version, paste0("metadata_", dataset_name, "_", version, ".csv"))
  )
}

# Define dataset groups and versions
datasets <- list(
  all_clusters = list(),
  lcam_hi = list(),
  lcam_lo = list(),
  lcam_both = list()
)

# Update to use new dataset names
versions <- c("metabolic", "nonmetabolic", "random")

# Populate datasets with paths
for (dataset_name in names(datasets)) {
  for (version in versions) {
    datasets[[dataset_name]][[version]] <- create_dataset_paths(dataset_name, version)
  }
}

# Initialize performance summary
performance_summary <- tibble(
  Dataset = character(),
  Accuracy = numeric(),
  Precision = numeric(),
  Recall = numeric(),
  F1_score = numeric(),
  ROC_AUC = numeric()
)

for (dataset_name in names(datasets)) {
  for (version in names(datasets[[dataset_name]])) {
    print(paste("Processing dataset:", dataset_name, "version:", version))
    
    # Create directories for current dataset and version
    predictions_output_dir <- file.path(predictions_dir, dataset_name, version, "predictions")
    feature_importance_output_dir <- file.path(feature_importance_dir, dataset_name, version, "feature_importance")
    model_output_dir <- file.path(models_dir, dataset_name, version)
    
    dir.create(predictions_output_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(feature_importance_output_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Load datasets
    X_train <- fread(datasets[[dataset_name]][[version]]$X_train)
    X_test <- fread(datasets[[dataset_name]][[version]]$X_test)
    y_train <- fread(datasets[[dataset_name]][[version]]$y_train)
    y_test <- fread(datasets[[dataset_name]][[version]]$y_test)
    metadata <- fread(datasets[[dataset_name]][[version]]$metadata)

    # Ensure proper label encoding
    y_train$x <- factor(y_train$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    y_test$x <- factor(y_test$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    
    # Convert to numeric for model fitting
    y_train_numeric <- as.numeric(as.character(y_train$x))
    y_test_numeric <- as.numeric(as.character(y_test$x))
    
    # Remove constant/zero-variance columns
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, !nzv$zeroVar, with = FALSE]
    
    # Define alpha values to try for elastic net
    alpha_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    
    # Initialize variables to store best model
    best_cv_fit <- NULL
    best_alpha <- NULL
    best_cvm <- Inf
    
    # Try different alpha values
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
    f1_score <- 2 * (precision * recall) / (precision + recall)
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
    
    # Save feature importance with the exact same format as RandomForest
    write_csv(
      feature_importance,
      file.path(feature_importance_output_dir, 
                paste0("feature_importance_", dataset_name, "_", version, ".csv"))
    )
    
    # Save predictions with the exact same format as RandomForest
    write_csv(
      data.frame(
        y_test = y_test$x,
        y_pred = y_pred
      ),
      file.path(predictions_output_dir, 
                paste0("predictions_", dataset_name, "_", version, ".csv"))
    )
    
    # Save ElasticNet model with consistent naming convention
    saveRDS(best_cv_fit, file.path(model_output_dir, 
            paste0("elasticnet_model_", dataset_name, "_", version, ".rds")))
  }
}

# Save performance summary
write_csv(performance_summary, file.path(summary_dir, "elasticnet_performance_summary.csv"))
cat("Saved performance summary:", file.path(summary_dir, "elasticnet_performance_summary.csv"), "\n")

print("ElasticNet analysis completed successfully.")