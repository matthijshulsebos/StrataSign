library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(iml)
library(data.table)
library(tibble)
library(ggplot2)
library(yaml)

# Path relative to the project root directory
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

# Define paths based on config
base_path <- file.path("output", "1. data preprocessing", "training datasets", dataset_type)
output_base_dir <- file.path("output", "2. models", dataset_type, "svm_linear")
summary_dir <- output_base_dir

# Create base output directory for the model type
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Define datasets and versions to process based on config
all_dataset_names <- c("all_clusters", "lcam_hi", "lcam_lo", "lcam_both", "macrophages")
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
    print(paste("Processing SVM Linear for dataset:", dataset_name, "| version:", version, "| type:", dataset_type))

    # Get file paths for current dataset and version
    paths <- create_dataset_paths(dataset_name, version)

    # Check for essential input files
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

    # Label encoding
    y_train$x <- factor(y_train$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    y_test$x <- factor(y_test$x, levels = c("Normal", "Tumor"), labels = c(0, 1))

    # Remove near-zero variance features
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, colnames(X_train), with = FALSE]

    # Apply PCA
    pca_model <- prcomp(as.matrix(X_train), center = TRUE, scale. = TRUE)

    # Calculate explained variance
    explained_variance <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)

    # Determine number of PCA components for linear kernel (>=90% variance)
    num_components_linear <- which(explained_variance >= 0.9)[1]
    print(paste("Number of components (Linear):", num_components_linear))

    # Transform data using PCA for linear kernel
    X_train_pca_linear <- predict(pca_model, as.matrix(X_train))[, 1:num_components_linear]
    X_test_pca_linear <- predict(pca_model, as.matrix(X_test))[, 1:num_components_linear]

    # Save PCA model
    saveRDS(pca_model, file.path(current_output_dir, paste0("pca_model_", dataset_name, "_", version, ".rds")))

    # Train SVM model with linear kernel
    print("Training SVM model with linear kernel...")
    svm_model_linear <- svm(X_train_pca_linear, y_train$x, kernel = "linear", probability = TRUE)

    # Save linear SVM model
    saveRDS(svm_model_linear, file.path(current_output_dir, paste0("svm_model_linear_", dataset_name, "_", version, ".rds")))

    # Make predictions on test set with linear kernel
    y_pred_linear <- predict(svm_model_linear, X_test_pca_linear, probability = TRUE)
    y_prob_linear <- attr(y_pred_linear, "probabilities")[,2]

    # Prepare factors for confusion matrix
    y_pred_linear_factor <- factor(y_pred_linear,
                                  levels = c(0, 1),
                                  labels = c("Normal", "Tumor"))
    y_test_factor <- factor(y_test$x,
                           levels = c(0, 1),
                           labels = c("Normal", "Tumor"))

    # Calculate confusion matrix
    confusion_linear <- confusionMatrix(y_pred_linear_factor, y_test_factor, positive = "Tumor")

    # Calculate performance metrics
    accuracy_linear <- confusion_linear$overall["Accuracy"]
    precision_linear <- ifelse(is.na(confusion_linear$byClass["Pos Pred Value"]), 0,
                              confusion_linear$byClass["Pos Pred Value"])
    recall_linear <- ifelse(is.na(confusion_linear$byClass["Sensitivity"]), 0,
                           confusion_linear$byClass["Sensitivity"])

    if (precision_linear == 0 && recall_linear == 0) {
        f1_score_linear <- 0
    } else {
        f1_score_linear <- 2 * (precision_linear * recall_linear) / (precision_linear + recall_linear)
    }
    roc_auc_linear <- as.numeric(roc(y_test_factor == "Tumor", y_prob_linear)$auc)

    # Store performance metrics
    performance_summary <- performance_summary %>%
      add_row(Dataset = paste(dataset_name, version, sep = "_"),
              Accuracy = accuracy_linear,
              Precision = precision_linear,
              Recall = recall_linear,
              F1_score = f1_score_linear,
              ROC_AUC = roc_auc_linear)

    # Save predictions
    write_csv(
      data.frame(
        y_test = y_test$x,
        y_pred = y_pred_linear
      ),
      file.path(current_output_dir,
             paste0("predictions_", dataset_name, "_", version, ".csv"))
    )

    # Calculate SHAP values for feature importance
    print("Calculating SHAP values for linear kernel...")
    # Create predictor object
    predictor_linear <- Predictor$new(svm_model_linear, data = as.data.frame(X_train_pca_linear), y = y_train$x)

    # Initialize matrix for SHAP values (PCs x Samples)
    all_shap_values_linear <- matrix(0, nrow = num_components_linear, ncol = nrow(X_test_pca_linear))

    # Calculate SHAP values per test sample for each PC
    for(i in 1:nrow(X_test_pca_linear)) {
      if(i %% 10 == 0) print(paste("Processing sample", i, "of", nrow(X_test_pca_linear)))
      shapley_linear <- Shapley$new(predictor_linear, x.interest = as.data.frame(X_test_pca_linear[i, , drop = FALSE]))
      all_shap_values_linear[, i] <- shapley_linear$results$phi[1:num_components_linear]
    }

    # Calculate mean SHAP values across all test samples for each PC
    shap_values_linear <- data.frame(phi = rowMeans(all_shap_values_linear))

    # Map PC SHAP values back to original features
    pca_loadings_linear <- as.matrix(pca_model$rotation[, 1:num_components_linear])

    # Calculate raw feature contributions from SHAP values and PCA loadings
    original_feature_contributions_linear <- as.data.frame(t(pca_loadings_linear %*% as.matrix(shap_values_linear)))
    colnames(original_feature_contributions_linear) <- colnames(X_train)

    feature_values <- as.numeric(original_feature_contributions_linear[1,])
    feature_names <- colnames(original_feature_contributions_linear)

    # Sort features by absolute importance for filtering
    sorted_idx <- order(abs(feature_values), decreasing = TRUE)

    # Filter raw feature importance to top n features
    keep_top <- 500
    if (length(sorted_idx) > keep_top) {
      keep_mask <- rep(FALSE, length(feature_values))
      keep_mask[sorted_idx[1:keep_top]] <- TRUE
      feature_values[!keep_mask] <- 0

      message(paste0("Keeping top ", keep_top, " features (",
                    round(keep_top/length(feature_values)*100, 1), "% of total)"))
    }

    # Create dataframe for raw feature importance
    feature_importance <- tibble(
      Feature = feature_names,
      Value = feature_values
    ) %>% arrange(desc(abs(Value)))

    # Save feature importance file
    write_csv(
      feature_importance,
      file.path(current_output_dir,
                paste0("feature_importance_", dataset_name, "_", version, ".csv")) # Standard filename
    )
  }
}

# Save performance summary
summary_filename <- file.path(summary_dir, "svm_linear_performance_summary.csv")
write_csv(performance_summary, summary_filename)
cat("Saved performance summary:", summary_filename, "\n")

print(paste("SVM Linear analysis completed successfully for dataset:", dataset_type))
