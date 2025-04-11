# Load necessary libraries
library(tidyverse)
library(caret)
library(randomForest)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(ggplot2)

# Define base path and output directories
base_path <- "output/training_data/raw"
models_dir <- "output/models/randomforest"
predictions_dir <- "output/models/randomforest"
feature_importance_dir <- "output/models/randomforest"
summary_dir <- "output/models/randomforest"

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

    # Ensure proper label encoding - keep as "Normal" and "Tumor" 
    y_train$x <- factor(y_train$x, levels = c("Normal", "Tumor"))
    y_test$x <- factor(y_test$x, levels = c("Normal", "Tumor"))
    
    # Remove constant/zero-variance columns
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, !nzv$zeroVar, with = FALSE]
    
    # Tune hyperparameters using cross-validation
    print("Tuning RandomForest hyperparameters...")
    set.seed(42)

    # Define the tuning grid for ranger
    tuning_grid <- expand.grid(
      mtry = c(floor(sqrt(ncol(X_train))), floor(ncol(X_train)/3), floor(ncol(X_train)/2)),
      splitrule = "gini",  # Required for ranger
      min.node.size = c(1, 3, 5)
    )

    # Define control parameters for tuning
    control <- trainControl(
      method = "cv",
      number = 5,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      verboseIter = TRUE
    )

    X_tune <- X_train
    y_tune <- y_train$x

    # Tune model with ranger 
    rf_tune <- train(
      x = as.matrix(X_tune),
      y = y_tune,
      method = "ranger",
      num.trees = 200,
      importance = "permutation",
      tuneGrid = tuning_grid,
      trControl = control,
      metric = "ROC",
      replace = TRUE,
      sample.fraction = 0.7  # Fixed sample fraction
    )

    # Get best parameters
    print("Best tuning parameters:")
    print(rf_tune$bestTune)

    # Use best parameters for final model
    best_mtry <- rf_tune$bestTune$mtry
    best_node_size <- rf_tune$bestTune$min.node.size
    best_sample_fraction <- 0.7  # Use the same value defined in train() call

    # Train final model with best parameters
    print("Training final RandomForest with optimized parameters...")
    rf_fit <- randomForest(
      x = as.matrix(X_train),
      y = y_train$x,  # Use original factor with "Normal"/"Tumor" levels
      ntree = 500,
      mtry = best_mtry,
      nodesize = best_node_size,
      sampsize = round(best_sample_fraction * nrow(X_train)),
      importance = TRUE,
      replace = TRUE
    )
    
    # Make predictions
    y_prob <- predict(rf_fit, newdata = as.matrix(X_test), type = "prob")[, "Tumor"]
    y_pred <- ifelse(y_prob > 0.5, "Tumor", "Normal")
    y_pred <- factor(y_pred, levels = c("Normal", "Tumor"))
    
    # Ensure factors have the same levels for confusion matrix
    y_test_factor <- factor(y_test$x, levels = c("Normal", "Tumor"))
    
    # Calculate performance metrics
    confusion <- confusionMatrix(y_pred, y_test$x)
    accuracy <- confusion$overall["Accuracy"]
    precision <- confusion$byClass["Pos Pred Value"]
    recall <- confusion$byClass["Sensitivity"]
    f1_score <- 2 * (precision * recall) / (precision + recall)
    roc_auc <- as.numeric(roc(y_test$x == "Tumor", y_prob)$auc)
    
    # Store performance metrics
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, sep = "_"),
              Accuracy = accuracy, 
              Precision = precision,
              Recall = recall,
              F1_score = f1_score,
              ROC_AUC = roc_auc)
    
    # Get feature importance - using MeanDecreaseGini for consistency
    imp <- importance(rf_fit, type = 2)  # Mean Decrease in Gini index

    # Debug: Check what's in the importance matrix
    print("Importance matrix structure:")
    print(str(imp))
    print("Importance matrix column names:")
    print(colnames(imp))
    print("First few importance values:")
    print(head(imp))

    # Check if there are any non-zero values
    print("Number of non-zero importance values:")
    print(sum(imp > 0))

    # Try using both importance measures available in RandomForest
    feature_importance <- data.frame(
      Feature = rownames(imp),
      Value = imp[, "MeanDecreaseGini"]  # This might be the issue if column names don't match
    ) %>%
      arrange(desc(Value))

    # Print debug info about non-zero features
    print(paste("Features with non-zero importance:", sum(feature_importance$Value > 0), 
                "out of", nrow(feature_importance)))

    # Save feature importance with exactly the same format as LASSO
    write_csv(
      feature_importance,
      file.path(feature_importance_output_dir, 
                paste0("feature_importance_", dataset_name, "_", version, ".csv"))
    )
    
    # Save predictions with identical format to LASSO
    write_csv(
      data.frame(
        y_test = y_test$x,
        y_pred = y_pred
      ), 
      file.path(predictions_output_dir, 
             paste0("predictions_", dataset_name, "_", version, ".csv"))
    )

    # Save RandomForest model
    saveRDS(rf_fit, file.path(model_output_dir, 
            paste0("randomforest_model_", dataset_name, "_", version, ".rds")))
  }
}

# Save performance summary
write_csv(performance_summary, file.path(summary_dir, "randomforest_performance_summary.csv"))
cat("Saved performance summary:", file.path(summary_dir, "randomforest_performance_summary.csv"), "\n")

print("Random Forest analysis completed successfully.")