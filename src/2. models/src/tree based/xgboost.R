library(tidyverse)
library(caret)
library(xgboost)
library(pROC)
library(data.table)
library(tibble)
library(yaml)
library(parallel)

# Some internal calls throw deprecation warnings
options(xgboost.check.deprecation = FALSE)

# Read config file
config_path <- "src/2. models/config.yaml" 
if (!file.exists(config_path)) {
  stop("Error: config.yaml not found at path: ", config_path)
}
config <- yaml::read_yaml(config_path)

# Extract config parameters
dataset_type <- config$training_parameters$dataset_type 
include_gene_ablation <- config$training_parameters$include_gene_ablation 
include_cell_type_ablation <- config$training_parameters$include_cell_type_ablation

# Validate dataset type
if (!dataset_type %in% c("absolute", "relative", "raw")) {
  stop("Invalid dataset_type in config.yaml. Must be 'absolute', 'relative', or 'raw'.")
}

# Define paths based on config
base_path <- file.path("output", "1. data preprocessing", "training datasets", dataset_type) 
output_base_dir <- file.path("output", "2. models", dataset_type, "xgboost") 
summary_dir <- output_base_dir 

# Create base output directory for the model type
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Define datasets based on config
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

# Set up parallel processing
n_cores <- 6  # Leave one core free
cat("Using", n_cores, "cores for processing\n")

# Main loop
for (dataset_name in dataset_names_to_process) {
  for (version in versions_to_process) {
    
    cat(paste("Processing XGBoost for dataset:", dataset_name, "| version:", version, "| type:", dataset_type, "\n"))
    
    paths <- create_dataset_paths(dataset_name, version)
    
    files_exist <- all(sapply(paths[c("X_train", "X_test", "y_train", "y_test")], file.exists)) 
    if (!files_exist) {
      warning(paste("Skipping:", dataset_name, version, "- Essential input file(s) not found."))
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"), 
                                                            Accuracy = NA, Precision = NA, Recall = NA,
                                                            F1_score = NA, ROC_AUC = NA)
      next 
    }

    # Create specific output directory for this version and dataset_name combination
    current_output_dir <- file.path(output_base_dir, version, dataset_name) 
    dir.create(current_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Load datasets
    X_train <- fread(paths$X_train, colClasses = "numeric")
    X_test <- fread(paths$X_test, colClasses = "numeric")
    y_train <- fread(paths$y_train)
    y_test <- fread(paths$y_test)
    
    # Save original feature names for later
    all_original_features <- colnames(X_train)
    
    # Label encoding target variable (1 is tumor)
    y_train_numeric <- as.numeric(y_train$x == "Tumor")
    y_test_numeric <- as.numeric(y_test$x == "Tumor")
    
    # Remove near-zero variance columns
    # Note: X_train is a data.table
    original_colnames <- colnames(X_train) # Save original column names for checking
    
    nzv_indices <- nearZeroVar(X_train, saveMetrics = FALSE) # Get indices of NZV columns TO REMOVE
    
    # Determine the names of columns to keep.
    # If nzv_indices is empty (integer(0)), original_colnames[-integer(0)] correctly returns all original_colnames.
    # Otherwise, it removes columns at the nzv_indices positions.
    cols_to_keep_names <- original_colnames[-nzv_indices]
    
    # Check if all columns were removed (and there were columns to begin with)
    if (length(cols_to_keep_names) == 0 && length(original_colnames) > 0) {
      warning(paste("All columns removed by NZV for", dataset_name, version, ". Skipping."))
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"),
                                                           Accuracy = NA, Precision = NA, Recall = NA,
                                                           F1_score = NA, ROC_AUC = NA)
      next
    }
    
    # Update X_train by selecting only the columns to keep.
    # If cols_to_keep_names is character(0) (e.g., X_train was initially empty),
    # X_train will become/remain an empty data.table with 0 columns.
    X_train <- X_train[, ..cols_to_keep_names]
    
    # Align test set columns with train set
    X_test_aligned <- X_test[, intersect(colnames(X_train), colnames(X_test)), with = FALSE]
    
    # Ensure column order matches
    X_test_aligned <- X_test_aligned[, colnames(X_train), with = FALSE]
    
    # Convert to matrix format
    tryCatch({
      dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train_numeric)
      dtest <- xgb.DMatrix(data = as.matrix(X_test_aligned), label = y_test_numeric)
    }, error = function(e) {
      warning(paste("Error creating XGBoost matrices for", dataset_name, version, ":", e$message))
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"),
                                                           Accuracy = NA, Precision = NA, Recall = NA,
                                                           F1_score = NA, ROC_AUC = NA)
      next
    })
    
    # Tune hyperparameters
    cat("Tuning XGBoost hyperparameters...\n")
    set.seed(42)

    # Define the grid for hyperparameter tuning
    tuning_grid <- expand.grid(
      nrounds = c(100),
      max_depth = c(3, 6, 9),
      eta = c(0.1, 0.3),
      gamma = c(0),
      colsample_bytree = c(0.8, 1.0),
      min_child_weight = c(1),
      subsample = c(0.8)
    )

    # Define tuning parameters
    control <- trainControl(
      method = "cv",
      number = 3,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      verboseIter = TRUE,
      allowParallel = TRUE
    )

    # Data type conversion
    X_tune <- as.matrix(X_train)
    y_tune <- factor(y_train$x, levels = c("Normal", "Tumor"))
    
    # Tune model
    xgb_tune <- tryCatch({
      train(
        x = X_tune,
        y = y_tune,
        method = "xgbTree",
        trControl = control,
        tuneGrid = tuning_grid,
        metric = "ROC",
        verbose = TRUE,
        nthread = 6
      )
    }, error = function(e) {
      warning(paste("Error tuning XGBoost for", dataset_name, version, ":", e$message))
      return(NULL)
    })
    
    if (is.null(xgb_tune)) {
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"),
                                                           Accuracy = NA, Precision = NA, Recall = NA,
                                                           F1_score = NA, ROC_AUC = NA)
      next
    }

    # Get best parameters
    cat("Best tuning parameters:\n")
    print(xgb_tune$bestTune)

    # Extract best parameters
    best_params <- list(
      max_depth = xgb_tune$bestTune$max_depth,
      eta = xgb_tune$bestTune$eta,
      gamma = xgb_tune$bestTune$gamma,
      colsample_bytree = xgb_tune$bestTune$colsample_bytree,
      min_child_weight = xgb_tune$bestTune$min_child_weight,
      subsample = xgb_tune$bestTune$subsample,
      nrounds = xgb_tune$bestTune$nrounds,
      objective = "binary:logistic",
      eval_metric = "logloss",
      alpha = 0,
      nthread = n_cores
    )

    # Train final model with best parameters
    cat("Training final XGBoost with optimized parameters...\n")
    xgb_model <- tryCatch({
      xgb.train(
        params = best_params,
        data = dtrain,
        nrounds = best_params$nrounds,
        watchlist = list(train = dtrain, test = dtest),
        early_stopping_rounds = 10,
        verbose = 0
      )
    }, error = function(e) {
      warning(paste("Error training final XGBoost model for", dataset_name, version, ":", e$message))
      return(NULL)
    })
    
    if (is.null(xgb_model)) {
      performance_summary <- performance_summary %>% add_row(Dataset = paste(dataset_name, version, sep = "_"),
                                                           Accuracy = NA, Precision = NA, Recall = NA,
                                                           F1_score = NA, ROC_AUC = NA)
      next
    }
    
    # Make predictions
    suppressWarnings({
      y_prob <- predict(
        xgb_model, 
        dtest
      )
    })
    
    # Convert probabilities to class labels
    y_pred <- ifelse(y_prob > 0.5, "Tumor", "Normal")

    # Create factors for the confusion matrix
    y_pred <- factor(y_pred, levels = c("Normal", "Tumor"))
    y_test_factor <- factor(y_test$x, levels = c("Normal", "Tumor"))
    
    # Create confusion matrix
    confusion <- confusionMatrix(y_pred, y_test_factor, positive = "Tumor")

    # Extract performance metrics
    accuracy <- confusion$overall["Accuracy"]
    precision <- ifelse(is.na(confusion$byClass["Pos Pred Value"]), 0, confusion$byClass["Pos Pred Value"])
    recall <- ifelse(is.na(confusion$byClass["Sensitivity"]), 0, confusion$byClass["Sensitivity"])
    f1_score <- if (precision == 0 && recall == 0) 0 else 2 * (precision * recall) / (precision + recall)
    roc_auc <- tryCatch(
      as.numeric(roc(y_test$x == "Tumor", y_prob, quiet = TRUE)$auc),
      error = function(e) { warning(paste("ROC AUC calculation failed:", e$message)); NA }
    )
    
    # Store performance metrics
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, sep = "_"),
              Accuracy = accuracy, 
              Precision = precision,
              Recall = recall,
              F1_score = f1_score,
              ROC_AUC = roc_auc)
    
    # Get feature importance
    suppressWarnings({
      importance_matrix <- tryCatch({
        xgb.importance(
          feature_names = colnames(X_train),
          model = xgb_model
        )
      }, error = function(e) {
        warning(paste("Error extracting feature importance for", dataset_name, version, ":", e$message))
        return(NULL)
      })
    })
    
    if (!is.null(importance_matrix)) {
      feature_importance <- data.frame(
        Feature = importance_matrix$Feature,
        Value = importance_matrix$Gain
      ) %>%
        arrange(desc(abs(Value))) 
      
      # Add zero importance for features not used by the model
      non_selected_features <- setdiff(all_original_features, feature_importance$Feature)
      if (length(non_selected_features) > 0) {
        zero_importance_df <- data.frame(
          Feature = non_selected_features,
          Value = 0,
          stringsAsFactors = FALSE
        )
        # Append zero importance features
        feature_importance <- rbind(feature_importance, zero_importance_df)
      }
      
      # Save feature importance to file
      write_csv(
        feature_importance[, c("Feature", "Value")],
        file.path(current_output_dir, paste0("feature_importance_", dataset_name, "_", version, ".csv"))
      )
    }
    
    # Save predictions
    write_csv(
      data.frame(
        y_test = y_test$x,
        y_pred = y_pred
      ), 
      file.path(current_output_dir, paste0("predictions_", dataset_name, "_", version, ".csv"))
    )

    # Save model
    tryCatch({
      xgb.save(xgb_model, file.path(current_output_dir, paste0("xgboost_model_", dataset_name, "_", version, ".model")))
    }, error = function(e) {
      warning(paste("Error saving model for", dataset_name, version, ":", e$message))
    })
    
    cat(paste("Finished processing:", dataset_name, version, "\n\n"))
  }
}

# Save performance summary
summary_filename <- file.path(summary_dir, "xgboost_performance_summary.csv")
write_csv(performance_summary, summary_filename)
cat("Saved performance summary:", summary_filename, "\n")
cat(paste("XGBoost analysis completed successfully for dataset type:", dataset_type, "\n"))
