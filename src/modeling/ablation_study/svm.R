# Load necessary libraries
library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(pROC)  # Add this line to load the pROC package
library(ggplot2)

# Define datasets
datasets <- list(
  all_clusters = list(
    X_train = "src/modeling/ablation_study/X_train_all_clusters.csv",
    X_test = "src/modeling/ablation_study/X_test_all_clusters.csv",
    y_train = "src/modeling/ablation_study/y_train_all_clusters.csv",
    y_test = "src/modeling/ablation_study/y_test_all_clusters.csv",
    results = "src/modeling/ablation_study/results_all_clusters.csv"
  ),
  lcam_hi = list(
    X_train = "src/modeling/ablation_study/X_train_lcam_hi.csv",
    X_test = "src/modeling/ablation_study/X_test_lcam_hi.csv",
    y_train = "src/modeling/ablation_study/y_train_lcam_hi.csv",
    y_test = "src/modeling/ablation_study/y_test_lcam_hi.csv",
    results = "src/modeling/ablation_study/results_lcam_hi.csv"
  ),
  lcam_lo = list(
    X_train = "src/modeling/ablation_study/X_train_lcam_lo.csv",
    X_test = "src/modeling/ablation_study/X_test_lcam_lo.csv",
    y_train = "src/modeling/ablation_study/y_train_lcam_lo.csv",
    y_test = "src/modeling/ablation_study/y_test_lcam_lo.csv",
    results = "src/modeling/ablation_study/results_lcam_lo.csv"
  ),
  lcam_both = list(
    X_train = "src/modeling/ablation_study/X_train_lcam_both.csv",
    X_test = "src/modeling/ablation_study/X_test_lcam_both.csv",
    y_train = "src/modeling/ablation_study/y_train_lcam_both.csv",
    y_test = "src/modeling/ablation_study/y_test_lcam_both.csv",
    results = "src/modeling/ablation_study/results_lcam_both.csv"
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

# Iterate over datasets
for (dataset_name in names(datasets)) {
  print(paste("Processing dataset:", dataset_name))
  
  # Load datasets
  print("Loading X_train and X_test...")
  X_train <- fread(datasets[[dataset_name]]$X_train)
  X_test <- fread(datasets[[dataset_name]]$X_test)
  
  # Load y_train and y_test
  print("Loading y_train and y_test...")
  y_train <- fread(datasets[[dataset_name]]$y_train)
  y_test <- fread(datasets[[dataset_name]]$y_test)

  # Ensure Proper Label Encoding
  print("Encoding labels...")
  y_train$x <- factor(y_train$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
  y_test$x <- factor(y_test$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
  
  # Remove constant/zero-variance columns
  print("Removing constant/zero-variance columns...")
  nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
  X_train <- X_train[, !nzv$nzv, with = FALSE]
  X_test <- X_test[, !nzv$nzv, with = FALSE]
  
  # Apply PCA for dimensionality reduction
  print("Applying PCA for dimensionality reduction...")
  pca_model <- prcomp(as.matrix(X_train), center = TRUE, scale. = TRUE)
  X_train_pca <- predict(pca_model, as.matrix(X_train))
  X_test_pca <- predict(pca_model, as.matrix(X_test))
  
  # Train SVM model with linear kernel
  print("Training SVM model with linear kernel...")
  svm_model_linear <- svm(X_train_pca, y_train$x, kernel = "linear", probability = TRUE)

  # Predict on test set with linear kernel
  print("Predicting on test set with linear kernel...")
  y_pred_linear <- predict(svm_model_linear, X_test_pca, probability = TRUE)
  y_prob_linear <- attr(y_pred_linear, "probabilities")[,2]

  # Calculate performance metrics for linear kernel
  print("Calculating performance metrics for linear kernel...")
  confusion_linear <- confusionMatrix(y_pred_linear, y_test$x)
  accuracy_linear <- confusion_linear$overall["Accuracy"]
  precision_linear <- confusion_linear$byClass["Pos Pred Value"]
  recall_linear <- confusion_linear$byClass["Sensitivity"]
  f1_score_linear <- 2 * (precision_linear * recall_linear) / (precision_linear + recall_linear)
  roc_auc_linear <- as.numeric(roc(y_test$x, y_prob_linear)$auc)

  # Store performance metrics for linear kernel
  performance_summary <- performance_summary %>%
    add_row(Dataset = paste(dataset_name, "Linear", sep = "_"), Accuracy = accuracy_linear, Precision = precision_linear, Recall = recall_linear, F1_score = f1_score_linear, ROC_AUC = roc_auc_linear)
  
  # Generate confusion matrix plot for linear kernel
  print("Generating confusion matrix plot for linear kernel...")
  confusion_matrix_plot_linear <- as.data.frame(confusion_linear$table)
  p <- ggplot(confusion_matrix_plot_linear, aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(paste("Confusion Matrix for Linear Kernel -", dataset_name))
  ggsave(paste0("src/modeling/ablation_study/confusion_matrix_linear_", dataset_name, ".png"), plot = p)

  # Generate SHAP values for linear kernel
  print("Generating SHAP values for linear kernel...")
  predictor_linear <- Predictor$new(svm_model_linear, data = as.data.frame(X_train_pca), y = y_train$x)
  shapley_linear <- Shapley$new(predictor_linear, x.interest = as.data.frame(X_test_pca[1, , drop = FALSE]))
  shapley_plot_linear <- shapley_linear$plot() +
    ggtitle(paste("SHAP Values for Linear Kernel -", dataset_name))
  ggsave(paste0("src/modeling/ablation_study/shap_linear_", dataset_name, ".png"), plot = shapley_plot_linear)
  
  # Train SVM model with RBF kernel
  print("Training SVM model with RBF kernel...")
  svm_model_rbf <- svm(X_train_pca, y_train$x, kernel = "radial", probability = TRUE)

  # Predict on test set with RBF kernel
  print("Predicting on test set with RBF kernel...")
  y_pred_rbf <- predict(svm_model_rbf, X_test_pca, probability = TRUE)
  y_prob_rbf <- attr(y_pred_rbf, "probabilities")[,2]

  # Calculate performance metrics for RBF kernel
  print("Calculating performance metrics for RBF kernel...")
  confusion_rbf <- confusionMatrix(y_pred_rbf, y_test$x)
  accuracy_rbf <- confusion_rbf$overall["Accuracy"]
  precision_rbf <- confusion_rbf$byClass["Pos Pred Value"]
  recall_rbf <- confusion_rbf$byClass["Sensitivity"]
  f1_score_rbf <- 2 * (precision_rbf * recall_rbf) / (precision_rbf + recall_rbf)
  roc_auc_rbf <- as.numeric(roc(y_test$x, y_prob_rbf)$auc)

  # Store performance metrics for RBF kernel
  performance_summary <- performance_summary %>%
    add_row(Dataset = paste(dataset_name, "RBF", sep = "_"), Accuracy = accuracy_rbf, Precision = precision_rbf, Recall = recall_rbf, F1_score = f1_score_rbf, ROC_AUC = roc_auc_rbf)

  # Generate confusion matrix plot for RBF kernel
  print("Generating confusion matrix plot for RBF kernel...")
  confusion_matrix_plot_rbf <- as.data.frame(confusion_rbf$table)
  p <- ggplot(confusion_matrix_plot_rbf, aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(paste("Confusion Matrix for RBF Kernel -", dataset_name))
  ggsave(paste0("src/modeling/ablation_study/confusion_matrix_rbf_", dataset_name, ".png"), plot = p)

  # Generate SHAP values for RBF kernel
  print("Generating SHAP values for RBF kernel...")
  predictor_rbf <- Predictor$new(svm_model_rbf, data = as.data.frame(X_train_pca), y = y_train$x)
  shapley_rbf <- Shapley$new(predictor_rbf, x.interest = as.data.frame(X_test_pca[1, , drop = FALSE]))
  shapley_plot_rbf <- shapley_rbf$plot() +
    ggtitle(paste("SHAP Values for RBF Kernel -", dataset_name))
  ggsave(paste0("src/modeling/ablation_study/shap_rbf_", dataset_name, ".png"), plot = shapley_plot_rbf)
}

# Save performance summary
write_csv(performance_summary, "src/modeling/ablation_study/model_performance_summary.csv")

print("Script completed successfully.")