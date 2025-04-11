# Load necessary libraries
library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(ggplot2)

# Define base path and output directories
base_path <- "output/training_data/raw"
models_dir <- "output/models/svm_linear"
predictions_dir <- "output/models/svm_linear"
feature_importance_dir <- "output/models/svm_linear"
summary_dir <- "output/models/svm_linear"

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
    
    # Remove constant/zero-variance columns
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, !nzv$zeroVar, with = FALSE]
    
    # Apply PCA for dimensionality reduction
    pca_model <- prcomp(as.matrix(X_train), center = TRUE, scale. = TRUE)
    
    # Calculate explained variance
    explained_variance <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)
    
    # Determine number of components for linear kernel
    num_components_linear <- which(explained_variance >= 0.90)[1]
    print(paste("Number of components (Linear):", num_components_linear))
    
    # Transform the data for linear kernel
    X_train_pca_linear <- predict(pca_model, as.matrix(X_train))[, 1:num_components_linear]
    X_test_pca_linear <- predict(pca_model, as.matrix(X_test))[, 1:num_components_linear]
    
    # Save PCA model
    saveRDS(pca_model, file.path(model_output_dir, paste0("pca_model_", dataset_name, "_", version, ".rds")))
    
    # Train SVM model with linear kernel
    print("Training SVM model with linear kernel...")
    svm_model_linear <- svm(X_train_pca_linear, y_train$x, kernel = "linear", probability = TRUE)
    
    # Save linear SVM model
    saveRDS(svm_model_linear, file.path(model_output_dir, paste0("svm_model_linear_", dataset_name, "_", version, ".rds")))

    # Predict on test set with linear kernel
    y_pred_linear <- predict(svm_model_linear, X_test_pca_linear, probability = TRUE)
    y_prob_linear <- attr(y_pred_linear, "probabilities")[,2]

    # Ensure consistent handling of factor levels
    y_pred_linear_factor <- factor(as.numeric(as.character(y_pred_linear)), 
                                  levels = c(0, 1), 
                                  labels = c("Normal", "Tumor"))
    y_test_factor <- factor(as.numeric(as.character(y_test$x)), 
                           levels = c(0, 1), 
                           labels = c("Normal", "Tumor"))

    # Calculate confusion matrix with explicit positive class
    confusion_linear <- confusionMatrix(y_pred_linear_factor, y_test_factor, positive = "Tumor")

    # Calculate metrics with safety checks
    accuracy_linear <- confusion_linear$overall["Accuracy"]
    precision_linear <- ifelse(is.na(confusion_linear$byClass["Pos Pred Value"]), 0, 
                              confusion_linear$byClass["Pos Pred Value"])
    recall_linear <- ifelse(is.na(confusion_linear$byClass["Sensitivity"]), 0, 
                           confusion_linear$byClass["Sensitivity"])

    # Calculate F1 score safely
    if (precision_linear == 0 && recall_linear == 0) {
        f1_score_linear <- 0
    } else {
        f1_score_linear <- 2 * (precision_linear * recall_linear) / (precision_linear + recall_linear)
    }

    # Calculate ROC AUC consistently
    roc_auc_linear <- as.numeric(roc(y_test_factor == "Tumor", y_prob_linear)$auc)

    # Store performance metrics for linear kernel
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, sep = "_"), 
              Accuracy = accuracy_linear, 
              Precision = precision_linear, 
              Recall = recall_linear, 
              F1_score = f1_score_linear, 
              ROC_AUC = roc_auc_linear)
    
    # Save predictions with proper format
    write_csv(
      data.frame(
        y_test = y_test$x,
        y_pred = y_pred_linear
      ), 
      file.path(predictions_output_dir, 
             paste0("predictions_", dataset_name, "_", version, ".csv"))
    )
    
    # Generate SHAP values for linear kernel
    print("Calculating SHAP values for linear kernel...")
    predictor_linear <- Predictor$new(svm_model_linear, data = as.data.frame(X_train_pca_linear), y = y_train$x)
    
    # Initialize matrix to store SHAP values for all test samples
    all_shap_values_linear <- matrix(0, nrow = num_components_linear, ncol = nrow(X_test_pca_linear))
    
    # Calculate SHAP values for each test sample
    for(i in 1:nrow(X_test_pca_linear)) {
      if(i %% 10 == 0) print(paste("Processing sample", i, "of", nrow(X_test_pca_linear)))
      shapley_linear <- Shapley$new(predictor_linear, x.interest = as.data.frame(X_test_pca_linear[i, , drop = FALSE]))
      all_shap_values_linear[, i] <- shapley_linear$results$phi[1:num_components_linear]
    }
    
    # Calculate mean SHAP values across all test samples
    shap_values_linear <- data.frame(phi = rowMeans(all_shap_values_linear))
    
    # Map SHAP values back to original features using linear components
    pca_loadings_linear <- as.matrix(pca_model$rotation[, 1:num_components_linear])
    
    # Calculate feature contributions (standard way)
    original_feature_contributions_linear <- as.data.frame(t(pca_loadings_linear %*% as.matrix(shap_values_linear)))
    colnames(original_feature_contributions_linear) <- colnames(X_train)
    
    # Extract raw values for regular contributions
    feature_values <- as.numeric(original_feature_contributions_linear[1,])
    feature_names <- colnames(original_feature_contributions_linear)

    # Sort features by absolute importance
    sorted_idx <- order(abs(feature_values), decreasing = TRUE)

    # Keep only top 500 features
    keep_top <- 500
    if (length(sorted_idx) > keep_top) {
      keep_mask <- rep(FALSE, length(feature_values))
      keep_mask[sorted_idx[1:keep_top]] <- TRUE
      feature_values[!keep_mask] <- 0
      
      # Report stats
      message(paste0("Keeping top ", keep_top, " features (", 
                    round(keep_top/length(feature_values)*100, 1), "% of total)"))
    }

    # Create dataframe with raw contributions
    feature_importance_raw <- tibble(
      Feature = feature_names,
      Value = feature_values
    ) %>% arrange(desc(abs(Value)))  # Sort by absolute magnitude
    
    # Calculate TRUE absolute contributions by summing absolute PC contributions
    # Initialize vector with the correct length
    abs_contributions <- numeric(length(feature_names))
    names(abs_contributions) <- feature_names
    
    # For each feature in PCA loadings, calculate absolute contributions
    for (i in 1:ncol(X_train)) {
      feature_name <- colnames(X_train)[i]
      if (feature_name %in% rownames(pca_loadings_linear)) {
        # Extract loadings for this feature across all PCs
        feature_loadings <- pca_loadings_linear[feature_name,]
        # Calculate contribution for each PC and take absolute sum
        pc_contributions <- feature_loadings * shap_values_linear$phi
        abs_contributions[feature_name] <- sum(abs(pc_contributions))
      }
    }
    
    # Create absolute importance dataframe
    feature_importance_abs <- tibble(
      Feature = names(abs_contributions),
      Value = abs_contributions
    ) %>% arrange(desc(Value))
    
    # Save both files - ONLY these two files, no backward compatibility file
    write_csv(
      feature_importance_raw,
      file.path(feature_importance_output_dir, 
                paste0("feature_importance_raw_", dataset_name, "_", version, ".csv"))
    )
    
    write_csv(
      feature_importance_abs,
      file.path(feature_importance_output_dir, 
                paste0("feature_importance_abs_", dataset_name, "_", version, ".csv"))
    )
  }
}

# Save performance summary
write_csv(performance_summary, file.path(summary_dir, "svm_linear_performance_summary.csv"))
cat("Saved performance summary:", file.path(summary_dir, "svm_linear_performance_summary.csv"), "\n")

print("SVM Linear analysis completed successfully.")