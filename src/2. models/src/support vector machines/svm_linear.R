library(tidyverse)
library(caret)
library(e1071)
library(pROC)
library(data.table)
library(yaml)

# Load config file
config_path <- "src/2. models/config.yaml"
if (!file.exists(config_path)) stop("Error: config.yaml not found at path: ", config_path)
config <- yaml::read_yaml(config_path)

# Set values from config
dataset_type <- config$training_parameters$dataset_type
include_gene_ablation <- config$training_parameters$include_gene_ablation
include_cell_type_ablation <- config$training_parameters$include_cell_type_ablation

if (!dataset_type %in% c("absolute", "relative", "raw")) {
  stop("Invalid dataset_type in config.yaml. Must be 'absolute', 'relative', or 'raw'.")
}

# Set and create base input paths
base_path <- file.path("output", "1. data preprocessing", "training datasets", dataset_type)
output_base_dir <- file.path("output", "2. models", dataset_type, "svm_linear")
summary_dir <- output_base_dir
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Set variables for dataset permutations
all_dataset_names <- c("all_clusters", "lcam_hi", "lcam_lo", "lcam_both", "macrophages")
all_versions <- c("metabolic", "nonmetabolic", "random")

# Set filters based on config
dataset_names_to_process <- if (include_cell_type_ablation) all_dataset_names else c("all_clusters")
versions_to_process <- if (include_gene_ablation) all_versions else c("metabolic")

# Helper function to create dataset paths
create_dataset_paths <- function(dataset_name, version) {
  list(
    X_train = file.path(base_path, dataset_name, version, paste0("X_train_", dataset_name, "_", version, ".csv")),
    X_test = file.path(base_path, dataset_name, version, paste0("X_test_", dataset_name, "_", version, ".csv")),
    y_train = file.path(base_path, dataset_name, version, paste0("y_train_", dataset_name, "_", version, ".csv")),
    y_test = file.path(base_path, dataset_name, version, paste0("y_test_", dataset_name, "_", version, ".csv"))
  )
}

# Initialize performance_summary
performance_summary <- tibble(
  Dataset = character(),
  Accuracy = numeric(),
  Precision = numeric(),
  Recall = numeric(),
  F1_score = numeric(),
  ROC_AUC = numeric()
)

# Helper function to calculate F-value for ANOVA
calculate_f_value <- function(feature_values, class_labels) {
  # Skip features with near zero variance
  if (sd(feature_values, na.rm = TRUE) < 1e-6) {
    return(0)
  }
  
  # Try to perform ANOVA test
  tryCatch(
    expr = {
      # One-way ANOVA between feature values and class labels
      anova_result <- aov(feature_values ~ class_labels)
      # Extract F-value from first row of ANOVA table
      return(summary(anova_result)[[1]]$"F value"[1])
    },
    error = function(e) {
      warning(paste("ANOVA failed for a feature:", conditionMessage(e)))
      return(0)
    }
  )
}

# Loop over datasets
for (dataset_name in dataset_names_to_process) {
  for (version in versions_to_process) {

    # Indication of current training run
    cat(paste("Processing SVM Linear for dataset:", dataset_name, "| version:", version, "| type:", dataset_type, "\n"))
    
    # Set paths for current dataset
    paths <- create_dataset_paths(dataset_name, version)

    # Check if essential files exist
    if (!all(sapply(paths[c("X_train", "X_test", "y_train", "y_test")], file.exists))) {
      warning(paste("Skipping:", dataset_name, version, "- Essential input file(s) not found."))
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"), Accuracy = NA, Precision = NA, Recall = NA, F1_score = NA, ROC_AUC = NA); next
    }

    # Create output directory for current dataset
    current_output_dir <- file.path(output_base_dir, version, dataset_name)
    dir.create(current_output_dir, recursive = TRUE, showWarnings = FALSE)

    # Load data
    X_train_dt <- fread(paths$X_train, colClasses = "numeric")
    X_test_dt <- fread(paths$X_test, colClasses = "numeric")
    y_train_dt <- fread(paths$y_train)
    y_test_dt <- fread(paths$y_test)

    # Save original features (because ANOVA selection will remove most)
    all_original_features <- colnames(X_train_dt)

    # Create factor labels for target variable
    y_train <- factor(y_train_dt$x, levels = c("Normal", "Tumor"), labels = c(0, 1))
    y_test <- factor(y_test_dt$x, levels = c("Normal", "Tumor"), labels = c(0, 1))

    # Near zero variance filter
    nzv_indices <- nearZeroVar(X_train_dt, saveMetrics = FALSE)
    cols_to_keep_names <- if (length(nzv_indices) > 0) colnames(X_train_dt)[-nzv_indices] else colnames(X_train_dt)

    # This does not really happen but just in case
    if (length(cols_to_keep_names) == 0) {
      warning(paste("All columns removed by NZV for", dataset_name, version, ". Skipping."))
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"), Accuracy = NA, Precision = NA, Recall = NA, F1_score = NA, ROC_AUC = NA); next
    }

    # Set train matrix with only non-NZV columns
    X_train_mat <- as.matrix(X_train_dt[, ..cols_to_keep_names])
    
    # Create empty test matrix with aligned columns
    X_test_mat_aligned <- matrix(0.0, nrow = nrow(X_test_dt), ncol = length(cols_to_keep_names), dimnames = list(NULL, cols_to_keep_names))
    
    # Sets common columns in test set to match train set
    common_cols_in_test <- intersect(cols_to_keep_names, colnames(X_test_dt))

    # Populate empty matrix with the test data
    if(length(common_cols_in_test) > 0) {
        X_test_mat_aligned[, common_cols_in_test] <- as.matrix(X_test_dt[, ..common_cols_in_test])
    }
    
    # Fit scaler on train data
    preProc_scaler <- preProcess(X_train_mat, method = c("center", "scale"))

    # Use scaler to transform both train and test data
    X_train_scaled <- predict(preProc_scaler, X_train_mat)
    X_test_scaled <- predict(preProc_scaler, X_test_mat_aligned)
    
    # This should never happen but just in case
    X_train_scaled[is.na(X_train_scaled) | !is.finite(X_train_scaled)] <- 0
    X_test_scaled[is.na(X_test_scaled) | !is.finite(X_test_scaled)] <- 0

    # ANOVA F-statistic feature selection should select 5000 features
    num_features_to_select <- min(5000, ncol(X_train_scaled))

    # This should never happen but just in case
    if (ncol(X_train_scaled) == 0) {
        warning(paste("X_train_scaled has 0 columns for", dataset_name, version, ". Skipping."))
        performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"), Accuracy = NA, Precision = NA, Recall = NA, F1_score = NA, ROC_AUC = NA); next
    }

    # Apply ANOVA to each feature
    f_values <- apply(X_train_scaled, 2, function(col) {
      calculate_f_value(col, y_train)
    })

    # Replace any NA values with zeros
    f_values[is.na(f_values)] <- 0

    # If all F-values are zero catch the error and skip
    if (all(f_values == 0)) {
      warning(paste("All F-values are zero for", dataset_name, version))
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"), Accuracy = NA, Precision = NA, Recall = NA, F1_score = NA, ROC_AUC = NA); next
    }

    # Order features by F-value and select top features
    ordered_indices <- order(f_values, decreasing = TRUE)
    num_to_actually_select <- min(num_features_to_select, sum(f_values > 0))

    # If all F-values are zero throw an error and skip
    if (num_to_actually_select == 0) {
      warning(paste("Could not select any features with non-zero F-value for", dataset_name, version))
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"), Accuracy = NA, Precision = NA, Recall = NA, F1_score = NA, ROC_AUC = NA); next
    }
    
    # Feature names of the selected features
    selected_feature_names <- colnames(X_train_scaled)[ordered_indices[1:num_to_actually_select]]
    
    # Filter the train and test data to keep only the selected features
    X_train_selected_df <- as.data.frame(X_train_scaled[, selected_feature_names, drop = FALSE])
    X_test_selected_df <- as.data.frame(X_test_scaled[, selected_feature_names, drop = FALSE])

    # Train SVM
    svm_model_linear <- svm(X_train_selected_df, y_train, kernel = "linear", probability = TRUE, scale = FALSE)
    
    # Save the model file
    saveRDS(svm_model_linear, file.path(current_output_dir, paste0("svm_model_linear_", dataset_name, "_", version, ".rds")))

    # Predict on the test set and include probabilities
    y_pred_svm_prob <- predict(svm_model_linear, X_test_selected_df, probability = TRUE)

    # Probabilities for sample to be the positive class which is tumour
    y_prob_values <- attr(y_pred_svm_prob, "probabilities")[, "1"]

    # Create factors out of predidctions and test labels
    y_pred_factor <- factor(y_pred_svm_prob, levels = c(0, 1), labels = c("Normal", "Tumor"))
    y_test_labels_factor <- factor(y_test, levels = c(0, 1), labels = c("Normal", "Tumor"))

    # Create confusion matrix with the factors
    confusion_res <- confusionMatrix(y_pred_factor, y_test_labels_factor, positive = "Tumor")

    # Extract performance metrics from confusion matrix
    accuracy <- confusion_res$overall["Accuracy"]
    precision <- ifelse(is.na(confusion_res$byClass["Pos Pred Value"]), 0, confusion_res$byClass["Pos Pred Value"])
    recall <- ifelse(is.na(confusion_res$byClass["Sensitivity"]), 0, confusion_res$byClass["Sensitivity"])
    f1_score <- if (precision == 0 && recall == 0) 0 else 2 * (precision * recall) / (precision + recall)
    roc_auc <- tryCatch(as.numeric(roc(response = y_test, predictor = y_prob_values, quiet = TRUE)$auc),
                        error = function(e) { warning(paste("ROC AUC calculation failed:", conditionMessage(e))); NA })

    # Add row to performance summary
    performance_summary <- performance_summary %>%
      add_row(Dataset = paste(dataset_name, version, sep = "_"), Accuracy = accuracy,
              Precision = precision, Recall = recall, F1_score = f1_score, ROC_AUC = roc_auc)

    # Save predictions to csv
    write_csv(data.frame(y_test = as.integer(as.character(y_test)), y_pred = as.integer(as.character(y_pred_svm_prob))),
              file.path(current_output_dir, paste0("predictions_", dataset_name, "_", version, ".csv")))

    # Extract SVM coefficients directly
    cat("Extracting SVM coefficients for linear model...\n")
    
    # Multiply support vector coefficients with support vectors to get feature importance
    w <- t(svm_model_linear$coefs) %*% svm_model_linear$SV
    
    # Store feature importance in a data frame
    coef_feature_importance <- data.frame(
      Feature = colnames(X_train_selected_df),
      Value = as.numeric(w),
      stringsAsFactors = FALSE
    )
    
    # Sort the feature importances by magnitude
    coef_feature_importance$Abs_Value <- abs(coef_feature_importance$Value)
    coef_feature_importance <- coef_feature_importance[order(coef_feature_importance$Abs_Value, decreasing = TRUE), ]
    
    # List of non selected features
    non_selected_features <- setdiff(all_original_features, selected_feature_names)
    
    # If there are non selected features add them and assign 0 importance
    if (length(non_selected_features) > 0) {
      zero_importance_df <- data.frame(
        Feature = non_selected_features,
        Value = 0,
        Abs_Value = 0,
        stringsAsFactors = FALSE
      )
      
      # Append zero importance features to feature importance data frame
      coef_feature_importance <- rbind(coef_feature_importance, zero_importance_df)
    }
    
    # Save complete coefficient-based feature importance to file
    write_csv(coef_feature_importance[, c("Feature", "Value")], 
              file.path(current_output_dir, paste0("feature_importance_", dataset_name, "_", version, ".csv")))
    
    cat(paste("Finished processing:", dataset_name, version, "\n\n"))
  }
}

# Save performance summary to CSV
summary_filename <- file.path(summary_dir, "svm_linear_performance_summary.csv")
write_csv(performance_summary, summary_filename)
cat("Saved performance summary:", summary_filename, "\n")
cat(paste("SVM Linear analysis completed successfully for dataset type:", dataset_type, "\n"))
