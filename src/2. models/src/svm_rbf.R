library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(ggplot2)
library(yaml)

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
if (!dataset_type %in% c("absolute", "relative")) {
  stop("Invalid dataset_type in config.yaml. Must be 'absolute' or 'relative'.")
}

# Define input paths based on config
base_path <- file.path("output", "1. data preprocessing", "training datasets", dataset_type) 
output_base_dir <- file.path("output", "2. models", dataset_type, "svm_rbf") 
summary_dir <- output_base_dir 

# Create base output directory for the model type
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Define datasets based on config
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

# Initialize performance summary
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

    # Create specific output directory for this version and dataset_name combination
    current_output_dir <- file.path(output_base_dir, version, dataset_name) 
    dir.create(current_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Load datasets
    X_train <- fread(paths$X_train)
    X_test <- fread(paths$X_test)
    y_train <- fread(paths$y_train)
    y_test <- fread(paths$y_test)

    # Ensure proper label encoding
    y_train$x <- factor(y_train$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    y_test$x <- factor(y_test$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    
    # Remove zero-variance columns
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, colnames(X_train), with = FALSE]
    
    # Apply PCA for dimensionality reduction
    pca_model <- prcomp(as.matrix(X_train), center = TRUE, scale. = TRUE)
    
    # Calculate explained variance
    explained_variance <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)
    
    # Determine number of components for RBF kernel
    num_components_rbf <- which(explained_variance >= 0.7)[1]
    print(paste("Number of components (RBF):", num_components_rbf))
    
    # Transform the data for RBF kernel
    X_train_pca_rbf <- predict(pca_model, as.matrix(X_train))[, 1:num_components_rbf]
    X_test_pca_rbf <- predict(pca_model, as.matrix(X_test))[, 1:num_components_rbf]
    
    # Save PCA model
    saveRDS(pca_model, file.path(current_output_dir, paste0("pca_model_", dataset_name, "_", version, ".rds")))
    
    # Train SVM model with RBF kernel
    print("Training SVM model with RBF kernel...")
    svm_model_rbf <- svm(X_train_pca_rbf, y_train$x, kernel = "radial", probability = TRUE)
    
    # Save RBF SVM model
    saveRDS(svm_model_rbf, file.path(current_output_dir, paste0("svm_model_rbf_", dataset_name, "_", version, ".rds")))

    # Predict on test set with RBF kernel
    y_pred_rbf <- predict(svm_model_rbf, X_test_pca_rbf, probability = TRUE)
    y_prob_rbf <- attr(y_pred_rbf, "probabilities")[,2]
    
    # Then create factors with consistent levels
    y_pred_factor <- factor(y_pred_rbf, levels = c(0, 1), labels = c("Normal", "Tumor"))
    y_test_factor <- factor(y_test$x, levels = c(0, 1), labels = c("Normal", "Tumor"))
    
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
      file.path(current_output_dir, 
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
    
    # Calculate feature contributions
    original_feature_contributions_rbf <- as.data.frame(t(pca_loadings_rbf %*% as.matrix(shap_values_rbf)))
    colnames(original_feature_contributions_rbf) <- colnames(X_train)

    # Extract raw values for regular contributions
    feature_values <- as.numeric(original_feature_contributions_rbf[1,])
    feature_names <- colnames(original_feature_contributions_rbf)

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

    # Create dataframe with raw contributions (this will be the main feature importance)
    feature_importance <- tibble(
      Feature = feature_names,
      Value = feature_values
    ) %>% arrange(desc(abs(Value)))  # Sort by absolute magnitude

    # Save feature importance file
    write_csv(
      feature_importance,
      file.path(current_output_dir, 
                paste0("feature_importance_", dataset_name, "_", version, ".csv"))
    )
  }
}

# Save performance summary
summary_filename <- file.path(summary_dir, "svm_rbf_performance_summary.csv") 
write_csv(performance_summary, summary_filename)
cat("Saved performance summary:", summary_filename, "\n")

print(paste("SVM RBF analysis completed successfully for dataset type:", dataset_type))
