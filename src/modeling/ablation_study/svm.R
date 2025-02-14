# Load necessary libraries
library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(ggplot2)

# Define base path
base_path <- "src/modeling/ablation_study"

# Define datasets structure
create_dataset_paths <- function(dataset_name, version) {
  list(
    X_train = file.path(base_path, "datasets", dataset_name, version, paste0("X_train_", dataset_name, ifelse(version == "default", "", paste0("_", version)), ".csv")),
    X_test = file.path(base_path, "datasets", dataset_name, version, paste0("X_test_", dataset_name, ifelse(version == "default", "", paste0("_", version)), ".csv")),
    y_train = file.path(base_path, "datasets", dataset_name, version, paste0("y_train_", dataset_name, ifelse(version == "default", "", paste0("_", version)), ".csv")),
    y_test = file.path(base_path, "datasets", dataset_name, version, paste0("y_test_", dataset_name, ifelse(version == "default", "", paste0("_", version)), ".csv")),
    metadata = file.path(base_path, "datasets", dataset_name, version, paste0("metadata_", dataset_name, ifelse(version == "default", "", paste0("_", version)), ".csv"))
  )
}

# Define dataset groups and versions
datasets <- list(
  all_clusters = list(),
  lcam_hi = list(),
  lcam_lo = list(),
  lcam_both = list()
)

versions <- c("default", "random1", "random2", "random3")

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

# Ensure the output directory exists
output_dir <- "src/modeling/ablation_study/model_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (dataset_name in names(datasets)) {
  for (version in names(datasets[[dataset_name]])) {
    print(paste("Processing dataset:", dataset_name, "version:", version))
    
    # Load datasets
    X_train <- fread(datasets[[dataset_name]][[version]]$X_train)
    X_test <- fread(datasets[[dataset_name]][[version]]$X_test)
    
    # Load y_train and y_test
    y_train <- fread(datasets[[dataset_name]][[version]]$y_train)
    y_test <- fread(datasets[[dataset_name]][[version]]$y_test)
    
    # Load metadata
    metadata <- fread(datasets[[dataset_name]][[version]]$metadata)

    # Ensure Proper Label Encoding
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
    
    # Determine number of components for each kernel
    num_components_linear <- which(explained_variance >= 0.90)[1]
    num_components_rbf <- which(explained_variance >= 0.70)[1]
    
    print(paste("Number of components (Linear):", num_components_linear))
    print(paste("Number of components (RBF):", num_components_rbf))
    
    # Transform the data for linear kernel
    X_train_pca_linear <- predict(pca_model, as.matrix(X_train))[, 1:num_components_linear]
    X_test_pca_linear <- predict(pca_model, as.matrix(X_test))[, 1:num_components_linear]
    
    # Transform the data for RBF kernel
    X_train_pca_rbf <- predict(pca_model, as.matrix(X_train))[, 1:num_components_rbf]
    X_test_pca_rbf <- predict(pca_model, as.matrix(X_test))[, 1:num_components_rbf]
    
    # Save PCA scores with metadata (using linear components for visualization)
    pca_scores <- as.data.frame(X_train_pca_linear)
    pca_scores <- cbind(metadata, pca_scores)
    dir.create(paste0(output_dir, "/", dataset_name, "/", version), recursive = TRUE, showWarnings = FALSE)
    write_csv(pca_scores, paste0(output_dir, "/", dataset_name, "/", version, "/pca_scores_", dataset_name, "_", version, ".csv"))
    
    # Train SVM model with linear kernel
    print("Training SVM model with linear kernel...")
    svm_model_linear <- svm(X_train_pca_linear, y_train$x, kernel = "linear", probability = TRUE)

    # Predict on test set with linear kernel
    y_pred_linear <- predict(svm_model_linear, X_test_pca_linear, probability = TRUE)
    y_prob_linear <- attr(y_pred_linear, "probabilities")[,2]

    # Save predictions and true labels for linear kernel
    write_csv(data.frame(y_test = y_test$x, y_pred = y_pred_linear), paste0(output_dir, "/", dataset_name, "/", version, "/predictions_linear_", dataset_name, "_", version, ".csv"))

    # Calculate performance metrics for linear kernel
    confusion_linear <- confusionMatrix(y_pred_linear, y_test$x)
    accuracy_linear <- confusion_linear$overall["Accuracy"]
    precision_linear <- confusion_linear$byClass["Pos Pred Value"]
    recall_linear <- confusion_linear$byClass["Sensitivity"]
    f1_score_linear <- 2 * (precision_linear * recall_linear) / (precision_linear + recall_linear)
    roc_auc_linear <- as.numeric(roc(y_test$x, y_prob_linear)$auc)

    # Store performance metrics for linear kernel
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, "Linear", sep = "_"), Accuracy = accuracy_linear, Precision = precision_linear, Recall = recall_linear, F1_score = f1_score_linear, ROC_AUC = roc_auc_linear)
    
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
    
    # Save average SHAP values for principal components
    write_csv(shap_values_linear, paste0(output_dir, "/", dataset_name, "/", version, "/shap_values_linear_", dataset_name, "_", version, ".csv"))
    
    # Map SHAP values back to original features using linear components
    pca_loadings_linear <- as.matrix(pca_model$rotation[, 1:num_components_linear])
    original_feature_contributions_linear <- as.data.frame(t(pca_loadings_linear %*% as.matrix(shap_values_linear)))
    colnames(original_feature_contributions_linear) <- colnames(X_train)
    
    # Save original feature contributions for linear kernel
    write_csv(original_feature_contributions_linear, paste0(output_dir, "/", dataset_name, "/", version, "/original_feature_contributions_linear_", dataset_name, "_", version, ".csv"))

    # Train SVM model with RBF kernel
    print("Training SVM model with RBF kernel...")
    svm_model_rbf <- svm(X_train_pca_rbf, y_train$x, kernel = "radial", probability = TRUE)

    # Predict on test set with RBF kernel
    y_pred_rbf <- predict(svm_model_rbf, X_test_pca_rbf, probability = TRUE)
    y_prob_rbf <- attr(y_pred_rbf, "probabilities")[,2]

    # Save predictions and true labels for RBF kernel
    write_csv(data.frame(y_test = y_test$x, y_pred = y_pred_rbf), paste0(output_dir, "/", dataset_name, "/", version, "/predictions_rbf_", dataset_name, "_", version, ".csv"))

    # Calculate performance metrics for RBF kernel
    confusion_rbf <- confusionMatrix(y_pred_rbf, y_test$x)
    accuracy_rbf <- confusion_rbf$overall["Accuracy"]
    precision_rbf <- confusion_rbf$byClass["Pos Pred Value"]
    recall_rbf <- confusion_rbf$byClass["Sensitivity"]
    f1_score_rbf <- 2 * (precision_rbf * recall_rbf) / (precision_rbf + recall_rbf)
    roc_auc_rbf <- as.numeric(roc(y_test$x, y_prob_rbf)$auc)

    # Store performance metrics for RBF kernel
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, "RBF", sep = "_"), Accuracy = accuracy_rbf, Precision = precision_rbf, Recall = recall_rbf, F1_score = f1_score_rbf, ROC_AUC = roc_auc_rbf)

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
    
    # Save average SHAP values for principal components
    write_csv(shap_values_rbf, paste0(output_dir, "/", dataset_name, "/", version, "/shap_values_rbf_", dataset_name, "_", version, ".csv"))

    # Map SHAP values back to original features using RBF components
    pca_loadings_rbf <- as.matrix(pca_model$rotation[, 1:num_components_rbf])
    original_feature_contributions_rbf <- as.data.frame(t(pca_loadings_rbf %*% as.matrix(shap_values_rbf)))
    colnames(original_feature_contributions_rbf) <- colnames(X_train)
    
    # Save original feature contributions for RBF kernel
    write_csv(original_feature_contributions_rbf, paste0(output_dir, "/", dataset_name, "/", version, "/original_feature_contributions_rbf_", dataset_name, "_", version, ".csv"))

    # Save SHAP values for plotting in the corresponding subdirectory
    saveRDS(shapley_linear, paste0(output_dir, "/", dataset_name, "/", version, "/shapley_linear_", dataset_name, "_", version, ".rds"))
    saveRDS(shapley_rbf, paste0(output_dir, "/", dataset_name, "/", version, "/shapley_rbf_", dataset_name, "_", version, ".rds"))
  }
}

# Save performance summary
write_csv(performance_summary, paste0(output_dir, "/model_performance_summary.csv"))
cat("Saved performance summary:", paste0(output_dir, "/model_performance_summary.csv"), "\n")

print("Script completed successfully.")
