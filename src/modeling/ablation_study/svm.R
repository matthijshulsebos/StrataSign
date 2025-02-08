# Load necessary libraries
library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(ggplot2)

# Define datasets
datasets <- list(
  all_clusters = list(
    X_train = "src/modeling/ablation_study/datasets/X_train_all_clusters.csv",
    X_test = "src/modeling/ablation_study/datasets/X_test_all_clusters.csv",
    y_train = "src/modeling/ablation_study/datasets/y_train_all_clusters.csv",
    y_test = "src/modeling/ablation_study/datasets/y_test_all_clusters.csv",
    results = "src/modeling/ablation_study/model_output/results_all_clusters.csv"
  ),
  lcam_hi = list(
    X_train = "src/modeling/ablation_study/datasets/X_train_lcam_hi.csv",
    X_test = "src/modeling/ablation_study/datasets/X_test_lcam_hi.csv",
    y_train = "src/modeling/ablation_study/datasets/y_train_lcam_hi.csv",
    y_test = "src/modeling/ablation_study/datasets/y_test_lcam_hi.csv",
    results = "src/modeling/ablation_study/model_output/results_lcam_hi.csv"
  ),
  lcam_lo = list(
    X_train = "src/modeling/ablation_study/datasets/X_train_lcam_lo.csv",
    X_test = "src/modeling/ablation_study/datasets/X_test_lcam_lo.csv",
    y_train = "src/modeling/ablation_study/datasets/y_train_lcam_lo.csv",
    y_test = "src/modeling/ablation_study/datasets/y_test_lcam_lo.csv",
    results = "src/modeling/ablation_study/model_output/results_lcam_lo.csv"
  ),
  lcam_both = list(
    X_train = "src/modeling/ablation_study/datasets/X_train_lcam_both.csv",
    X_test = "src/modeling/ablation_study/datasets/X_test_lcam_both.csv",
    y_train = "src/modeling/ablation_study/datasets/y_train_lcam_both.csv",
    y_test = "src/modeling/ablation_study/datasets/y_test_lcam_both.csv",
    results = "src/modeling/ablation_study/model_output/results_lcam_both.csv"
  )
)

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
  print(paste("Processing dataset:", dataset_name))
  
  # Load datasets
  X_train <- fread(datasets[[dataset_name]]$X_train)
  X_test <- fread(datasets[[dataset_name]]$X_test)
  
  # Load y_train and y_test
  y_train <- fread(datasets[[dataset_name]]$y_train)
  y_test <- fread(datasets[[dataset_name]]$y_test)

  # Ensure Proper Label Encoding
  y_train$x <- factor(y_train$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
  y_test$x <- factor(y_test$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
  
  # Remove constant/zero-variance columns
  nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
  X_train <- X_train[, !nzv$zeroVar, with = FALSE]
  X_test <- X_test[, !nzv$zeroVar, with = FALSE]
  
  # Apply PCA for dimensionality reduction
  pca_model <- prcomp(as.matrix(X_train), center = TRUE, scale. = TRUE)
  
  # Determine the number of components to retain (90% variance)
  explained_variance <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)
  num_components <- which(explained_variance >= 0.90)[1]
  
  # Transform the data using the selected components
  X_train_pca <- predict(pca_model, as.matrix(X_train))[, 1:num_components]
  X_test_pca <- predict(pca_model, as.matrix(X_test))[, 1:num_components]
  
  # Train SVM model with linear kernel
  print("Training SVM model with linear kernel...")
  svm_model_linear <- svm(X_train_pca, y_train$x, kernel = "linear", probability = TRUE)

  # Predict on test set with linear kernel
  y_pred_linear <- predict(svm_model_linear, X_test_pca, probability = TRUE)
  y_prob_linear <- attr(y_pred_linear, "probabilities")[,2]

  # Save predictions and true labels for linear kernel
  write_csv(data.frame(y_test = y_test$x, y_pred = y_pred_linear), paste0(output_dir, "/predictions_linear_", dataset_name, ".csv"))

  # Calculate performance metrics for linear kernel
  confusion_linear <- confusionMatrix(y_pred_linear, y_test$x)
  accuracy_linear <- confusion_linear$overall["Accuracy"]
  precision_linear <- confusion_linear$byClass["Pos Pred Value"]
  recall_linear <- confusion_linear$byClass["Sensitivity"]
  f1_score_linear <- 2 * (precision_linear * recall_linear) / (precision_linear + recall_linear)
  roc_auc_linear <- as.numeric(roc(y_test$x, y_prob_linear)$auc)

  # Store performance metrics for linear kernel
  performance_summary <- performance_summary %>%
    add_row(Dataset = paste(dataset_name, "Linear", sep = "_"), Accuracy = accuracy_linear, Precision = precision_linear, Recall = recall_linear, F1_score = f1_score_linear, ROC_AUC = roc_auc_linear)
  
  # Generate SHAP values for linear kernel
  predictor_linear <- Predictor$new(svm_model_linear, data = as.data.frame(X_train_pca), y = y_train$x)
  shapley_linear <- Shapley$new(predictor_linear, x.interest = as.data.frame(X_test_pca[1, , drop = FALSE]))
  
  # Save SHAP values for principal components
  shap_values_linear <- as.data.frame(shapley_linear$results$phi)
  shap_values_linear <- shap_values_linear[1:num_components, , drop = FALSE]  # Filter SHAP values to include only retained components
  write_csv(shap_values_linear, paste0(output_dir, "/shap_values_linear_", dataset_name, ".csv"))

  # Map SHAP values back to original features
  pca_loadings <- as.matrix(pca_model$rotation[, 1:num_components])
  original_feature_contributions_linear <- as.data.frame(t(pca_loadings %*% as.matrix(shap_values_linear)))
  colnames(original_feature_contributions_linear) <- colnames(X_train)
  
  # Save original feature contributions for linear kernel
  write_csv(original_feature_contributions_linear, paste0(output_dir, "/original_feature_contributions_linear_", dataset_name, ".csv"))

  # Train SVM model with RBF kernel
  print("Training SVM model with RBF kernel...")
  svm_model_rbf <- svm(X_train_pca, y_train$x, kernel = "radial", probability = TRUE)

  # Predict on test set with RBF kernel
  y_pred_rbf <- predict(svm_model_rbf, X_test_pca, probability = TRUE)
  y_prob_rbf <- attr(y_pred_rbf, "probabilities")[,2]

  # Save predictions and true labels for RBF kernel
  write_csv(data.frame(y_test = y_test$x, y_pred = y_pred_rbf), paste0(output_dir, "/predictions_rbf_", dataset_name, ".csv"))

  # Calculate performance metrics for RBF kernel
  confusion_rbf <- confusionMatrix(y_pred_rbf, y_test$x)
  accuracy_rbf <- confusion_rbf$overall["Accuracy"]
  precision_rbf <- confusion_rbf$byClass["Pos Pred Value"]
  recall_rbf <- confusion_rbf$byClass["Sensitivity"]
  f1_score_rbf <- 2 * (precision_rbf * recall_rbf) / (precision_rbf + recall_rbf)
  roc_auc_rbf <- as.numeric(roc(y_test$x, y_prob_rbf)$auc)

  # Store performance metrics for RBF kernel
  performance_summary <- performance_summary %>%
    add_row(Dataset = paste(dataset_name, "RBF", sep = "_"), Accuracy = accuracy_rbf, Precision = precision_rbf, Recall = recall_rbf, F1_score = f1_score_rbf, ROC_AUC = roc_auc_rbf)

  # Generate SHAP values for RBF kernel
  predictor_rbf <- Predictor$new(svm_model_rbf, data = as.data.frame(X_train_pca), y = y_train$x)
  shapley_rbf <- Shapley$new(predictor_rbf, x.interest = as.data.frame(X_test_pca[1, , drop = FALSE]))
  
  # Save SHAP values for principal components
  shap_values_rbf <- as.data.frame(shapley_rbf$results$phi)
  shap_values_rbf <- shap_values_rbf[1:num_components, , drop = FALSE]  # Filter SHAP values to include only retained components
  write_csv(shap_values_rbf, paste0(output_dir, "/shap_values_rbf_", dataset_name, ".csv"))

  # Map SHAP values back to original features
  original_feature_contributions_rbf <- as.data.frame(t(pca_loadings %*% as.matrix(shap_values_rbf)))
  colnames(original_feature_contributions_rbf) <- colnames(X_train)
  
  # Save original feature contributions for RBF kernel
  write_csv(original_feature_contributions_rbf, paste0(output_dir, "/original_feature_contributions_rbf_", dataset_name, ".csv"))

  # Save SHAP values for plotting
  saveRDS(shapley_linear, paste0(output_dir, "/shapley_linear_", dataset_name, ".rds"))
  saveRDS(shapley_rbf, paste0(output_dir, "/shapley_rbf_", dataset_name, ".rds"))
}

# Save performance summary
write_csv(performance_summary, paste0(output_dir, "/model_performance_summary.csv"))
cat("Saved performance summary:", paste0(output_dir, "/model_performance_summary.csv"), "\n")

print("Script completed successfully.")
