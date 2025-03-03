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
models_dir <- "output/models/lasso"
predictions_dir <- "output/models/lasso"
feature_importance_dir <- "output/models/lasso"
summary_dir <- "output/models/lasso"

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
    
    # Convert to numeric for LASSO fitting
    y_train_numeric <- as.numeric(as.character(y_train$x))
    y_test_numeric <- as.numeric(as.character(y_test$x))
    
    # Remove constant/zero-variance columns
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, !nzv$zeroVar, with = FALSE]
    
    # Fit LASSO model with LOOCV
    print("Training LASSO model with LOOCV...")
    cv_fit <- cv.glmnet(as.matrix(X_train), y_train_numeric, 
                        alpha = 1,                    # LASSO penalty
                        family = "binomial",          # Binary classification
                        nfolds = nrow(X_train),       # LOOCV: number of folds equals number of samples
                        type.measure = "deviance"     # Measure to evaluate model
    )

    # Get predictions on test set using optimal lambda
    y_prob <- predict(cv_fit, as.matrix(X_test), s = "lambda.min", type = "response")[,1]
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
    
    # Get feature importance (coefficients) using optimal lambda
    coef_matrix <- as.matrix(coef(cv_fit, s = "lambda.min"))
    feature_importance <- data.frame(
      Feature = rownames(coef_matrix),
      Coefficient = coef_matrix[,1]
    ) %>%
      filter(Feature != "(Intercept)") %>%
      mutate(
        Abs_Coefficient = abs(Coefficient),
        Direction = ifelse(Coefficient > 0, "Positive", "Negative"),
        Gene = sapply(strsplit(Feature, "@"), `[`, 1),
        Cluster = sapply(strsplit(Feature, "@"), `[`, 2)
      ) %>%
      arrange(desc(abs(Coefficient)))
    
    # Save feature importance 
    write_csv(
      feature_importance %>% 
        select(Feature, Coefficient) %>%
        rename(Value = Coefficient),
      file.path(feature_importance_output_dir, 
                paste0("feature_importance_", dataset_name, "_", version, ".csv"))
    )
    
    # Save predictions with proper format
    write_csv(
      data.frame(
        y_test = y_test$x,  # This will be 0/1 factors
        y_pred = y_pred     # These are our 0/1 predictions
      ), 
      file.path(predictions_output_dir, 
             paste0("predictions_", dataset_name, "_", version, ".csv"))
    )

    # Save LASSO model
    saveRDS(cv_fit, file.path(model_output_dir, 
            paste0("lasso_model_", dataset_name, "_", version, ".rds")))
  }
}

# Save performance summary
write_csv(performance_summary, file.path(summary_dir, "lasso_performance_summary.csv"))
cat("Saved performance summary:", file.path(summary_dir, "lasso_performance_summary.csv"), "\n")

print("LASSO analysis completed successfully.")