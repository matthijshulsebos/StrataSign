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
models_dir <- "output/models/svm_rbf"
predictions_dir <- "output/models/svm_rbf"
feature_importance_dir <- "output/models/svm_rbf"
summary_dir <- "output/models/svm_rbf"

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

# Update the versions to match our new names
versions <- c("metabolic", "matched_nonmetabolic", "matched_random_a", "matched_random_b")

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
    
    # Determine number of components for RBF kernel
    num_components_rbf <- which(explained_variance >= 0.70)[1]
    print(paste("Number of components (RBF):", num_components_rbf))
    
    # Transform the data for RBF kernel
    X_train_pca_rbf <- predict(pca_model, as.matrix(X_train))[, 1:num_components_rbf]
    X_test_pca_rbf <- predict(pca_model, as.matrix(X_test))[, 1:num_components_rbf]
    
    # Save PCA model
    saveRDS(pca_model, file.path(model_output_dir, paste0("pca_model_", dataset_name, "_", version, ".rds")))
    
    # Train SVM model with RBF kernel
    print("Training SVM model with RBF kernel...")
    svm_model_rbf <- svm(X_train_pca_rbf, y_train$x, kernel = "radial", probability = TRUE)
    
    # Save RBF SVM model
    saveRDS(svm_model_rbf, file.path(model_output_dir, paste0("svm_model_rbf_", dataset_name, "_", version, ".rds")))

    # Predict on test set with RBF kernel
    y_pred_rbf <- predict(svm_model_rbf, X_test_pca_rbf, probability = TRUE)
    y_prob_rbf <- attr(y_pred_rbf, "probabilities")[,2]

    # First, properly encode the raw predictions
    y_pred_rbf <- as.numeric(as.character(y_pred_rbf))  # Convert to numeric first
    y_test_raw <- as.numeric(as.character(y_test$x))    # Get raw test values
    
    # Then create factors with consistent levels
    y_pred_factor <- factor(y_pred_rbf, levels = c(0, 1), labels = c("Normal", "Tumor"))
    y_test_factor <- factor(y_test_raw, levels = c(0, 1), labels = c("Normal", "Tumor"))
    
    # Calculate confusion matrix with explicit positive class
    confusion_rbf <- confusionMatrix(y_pred_factor, y_test_factor, positive = "Tumor")
    
    # Calculate metrics ensuring no NA values
    accuracy_rbf <- confusion_rbf$overall["Accuracy"]
    precision_rbf <- ifelse(is.na(confusion_rbf$byClass["Pos Pred Value"]), 0, 
                           confusion_rbf$byClass["Pos Pred Value"])
    recall_rbf <- ifelse(is.na(confusion_rbf$byClass["Sensitivity"]), 0, 
                        confusion_rbf$byClass["Sensitivity"])
    
    # Calculate F1 score safely
    if (precision_rbf == 0 && recall_rbf == 0) {
        f1_score_rbf <- 0
    } else {
        f1_score_rbf <- 2 * (precision_rbf * recall_rbf) / (precision_rbf + recall_rbf)
    }
    
    # Calculate ROC AUC
    roc_auc_rbf <- as.numeric(roc(y_test_factor == "Tumor", y_prob_rbf)$auc)

    # Store performance metrics for RBF kernel
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, sep = "_"), 
              Accuracy = accuracy_rbf, 
              Precision = precision_rbf, 
              Recall = recall_rbf, 
              F1_score = f1_score_rbf, 
              ROC_AUC = roc_auc_rbf)
    
    # Save predictions with proper format
    write_csv(
      data.frame(
        y_test = y_test$x,
        y_pred = y_pred_rbf
      ), 
      file.path(predictions_output_dir, 
             paste0("predictions_", dataset_name, "_", version, ".csv"))
    )

    # Generate SHAP values for RBF kernel
    print("Calculating SHAP values for RBF kernel...")
    predictor_rbf <- Predictor$new(svm_model_rbf, data = as.data.frame(X_train_pca_rbf), y = y_train$x)
    
    # Initialize matrix to store SHAP values for all test samples
    all_shap_values_rbf <- matrix(0, nrow = num_components_rbf, ncol = nrow(X_test_pca_rbf))
    
    # Calculate SHAP values for each test sample
    for(i in 1:nrow(X_test_pca_rbf)) {
      if(i %% 10 == 0) print(paste("Processing sample", i, "of", nrow(X_test_pca_rbf)))
      shapley_rbf <- Shapley$new(predictor_rbf, x.interest = as.data.frame(X_test_pca_rbf[i, , drop = FALSE]))
      all_shap_values_rbf[, i] <- shapley_rbf$results$phi[1:num_components_rbf]
    }
    
    # Calculate mean SHAP values across all test samples
    shap_values_rbf <- data.frame(phi = rowMeans(all_shap_values_rbf))
    
    # Map SHAP values back to original features using RBF components
    pca_loadings_rbf <- as.matrix(pca_model$rotation[, 1:num_components_rbf])
    
    # Calculate feature contributions (standard way - net effect)
    original_feature_contributions_rbf <- as.data.frame(t(pca_loadings_rbf %*% as.matrix(shap_values_rbf)))
    colnames(original_feature_contributions_rbf) <- colnames(X_train)

    # Extract raw values for regular contributions
    feature_values <- as.numeric(original_feature_contributions_rbf[1,])
    feature_names <- colnames(original_feature_contributions_rbf)

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
      if (feature_name %in% rownames(pca_loadings_rbf)) {
        # Extract loadings for this feature across all PCs
        feature_loadings <- pca_loadings_rbf[feature_name,]
        # Calculate contribution for each PC and take absolute sum
        pc_contributions <- feature_loadings * shap_values_rbf$phi
        abs_contributions[feature_name] <- sum(abs(pc_contributions))
      }
    }

    # Create absolute importance dataframe
    feature_importance_abs <- tibble(
      Feature = names(abs_contributions),
      Value = abs_contributions
    ) %>% arrange(desc(Value))

    # Save both files
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
write_csv(performance_summary, file.path(summary_dir, "svm_rbf_performance_summary.csv"))
cat("Saved performance summary:", file.path(summary_dir, "svm_rbf_performance_summary.csv"), "\n")

print("SVM RBF analysis completed successfully.")