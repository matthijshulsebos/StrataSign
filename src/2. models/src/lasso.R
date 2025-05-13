# Load necessary libraries
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
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

# Define base input path based on config (relative to project root)
base_path <- file.path("output", "1. data preprocessing", "training datasets", dataset_type) 

# Define base output directories based on config (relative to project root)
output_base_dir <- file.path("output", "2. models", dataset_type, "lasso")
summary_dir <- output_base_dir 

# Create base output directory for the model type
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Define cell types and genes to based on config
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

# Init performance summary
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
    
    # Convert to numeric for model training
    y_train_numeric <- as.numeric(as.character(y_train$x))
    y_test_numeric <- as.numeric(as.character(y_test$x))
    
    # Remove zero-variance columns
    nzv <- nearZeroVar(X_train, saveMetrics = TRUE)
    X_train <- X_train[, !nzv$zeroVar, with = FALSE]
    X_test <- X_test[, !nzv$zeroVar, with = FALSE]
    
    print("Training LASSO model with bootstrapping...")
    n_bootstraps <- 100 
    n_samples <- nrow(X_train)
    all_coefs <- list() 
    all_test_probs <- matrix(NA, nrow = nrow(X_test), ncol = n_bootstraps)
    
    # Get potential feature names from the full training set once
    # This ensures we track all features, even if not selected in some bootstraps
    initial_fit <- glmnet(as.matrix(X_train), y_train_numeric, alpha = 1, family = "binomial", nlambda = 1) 
    feature_names <- rownames(coef(initial_fit))
    
    # Progress bar because bootstrapping takes ages
    pb <- txtProgressBar(min = 0, max = n_bootstraps, style = 3)

    for (i in 1:n_bootstraps) {
      # Ensure different samples each time
      set.seed(42 + i)
      
      # Create bootstrap sample indices
      boot_indices <- sample(1:n_samples, n_samples, replace = TRUE)
      
      # Create bootstrap data
      X_train_boot <- X_train[boot_indices, ]
      y_train_boot <- y_train_numeric[boot_indices]
      
      # Fit LASSO model with CV on the bootstrap sample
      cv_fit_boot <- cv.glmnet(as.matrix(X_train_boot), y_train_boot, 
                               alpha = 1,            # LASSO penalty
                               family = "binomial",  # Binary classification
                               nfolds = 5,           # 5-fold CV for lambda selection
                               type.measure = "deviance",
                               standardize = FALSE # Assuming data is already scaled or handled
      )
      
      # Store coefficients for lambda.min as a named vector for easier access
      coef_vector_boot <- as.vector(coef(cv_fit_boot, s = "lambda.min"))
      names(coef_vector_boot) <- feature_names # Assign standard names
      all_coefs[[i]] <- coef_vector_boot
      
      # Predict on the original test set
      y_prob_boot <- predict(cv_fit_boot, as.matrix(X_test), s = "lambda.min", type = "response")[,1]
      all_test_probs[, i] <- y_prob_boot
      
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    # Close progress bar
    close(pb) 
    print("Bootstrap loop finished.")
    
    # Average test probabilities
    print("Averaging test probabilities...")
    y_prob <- rowMeans(all_test_probs, na.rm = TRUE) 
    y_pred <- ifelse(y_prob > 0.5, 1, 0)
    print("Averaging complete.")
    
    # Calculate feature importance as the average of non-zero coefficients
    print("Starting feature importance calculation...")
    feature_importance_list <- list()
    feature_counter <- 0
    
    for (feature in feature_names) {
        feature_counter <- feature_counter + 1
        if (feature_counter %% 1000 == 0) {
            print(paste("  Processing feature importance", feature_counter, "/", length(feature_names))) 
        }
        
        # Skip intercept
        if (feature == "(Intercept)") next 
        
        # Extract coefficients for this feature across all bootstraps
        feature_coeffs <- sapply(all_coefs, function(coef_vec) {
            # Check if the feature is present in the coefficient vector
            val <- coef_vec[feature]
            if (is.na(val) || val == 0) {
                # Use NA to signify zero or non-selection for averaging non-zeros
                return(NA) 
            } else {
                return(val)
            }
        })
        
        # Filter out NAs and calculate mean
        non_zero_coeffs <- feature_coeffs[!is.na(feature_coeffs)]
        
        if (length(non_zero_coeffs) > 0) {
            avg_non_zero_coef <- mean(non_zero_coeffs)
        } else {
            # Assign 0 if never selected
            avg_non_zero_coef <- 0 
        }
        
        feature_importance_list[[feature]] <- data.frame(
            Feature = feature,
            Avg_NonZero_Coefficient = avg_non_zero_coef
        )
    }
    print("Finished feature importance loop.")
    
    feature_importance <- bind_rows(feature_importance_list) %>%
      mutate(
        Abs_Avg_Coefficient = abs(Avg_NonZero_Coefficient)
      ) %>%
      arrange(desc(Abs_Avg_Coefficient)) # Order by magnitude
    print("Finished feature importance calculation and processing.")

    # Ensure factors have the same levels for confusion matrix
    print("Calculating performance metrics...")
    y_pred <- factor(y_pred, levels = c(0, 1))
    y_test_factor <- factor(y_test_numeric, levels = c(0, 1))
    
    # Calculate performance metrics
    confusion <- confusionMatrix(y_pred, y_test_factor)
    accuracy <- confusion$overall["Accuracy"]
    precision <- confusion$byClass["Pos Pred Value"]
    recall <- confusion$byClass["Sensitivity"]
    f1_score <- confusion$byClass["F1"]
    roc_auc <- as.numeric(roc(y_test_numeric, y_prob)$auc)
    print("Performance metrics calculated.")
    
    # Store performance metrics
    performance_summary <- performance_summary %>% 
      add_row(Dataset = paste(dataset_name, version, sep = "_"),
              Accuracy = accuracy, 
              Precision = precision,
              Recall = recall,
              F1_score = f1_score,
              ROC_AUC = roc_auc)
    
    # Save feature importance 
    print("Saving feature importance...")
    write_csv(
      feature_importance %>% 
        select(Feature, Avg_NonZero_Coefficient) %>% 
        rename(Value = Avg_NonZero_Coefficient),     
      file.path(current_output_dir, 
                paste0("feature_importance_", dataset_name, "_", version, ".csv"))
    )
    print("Feature importance saved.")
    
    # Save aggregated predictions 
    print("Saving predictions...")
    write_csv(
      data.frame(
        y_test = y_test$x,  
        y_pred = y_pred     
      ), 
      file.path(current_output_dir, 
             paste0("predictions_", dataset_name, "_", version, ".csv"))
    )
    print("Predictions saved.")
    
    print(paste("  Finished post-processing for:", dataset_name, version))
    
  }
}

# Save performance summary
print("Saving final performance summary...") 
summary_filename <- file.path(summary_dir, "lasso_performance_summary.csv") 
write_csv(performance_summary, summary_filename)
cat("Saved performance summary:", summary_filename, "\n")

print(paste("LASSO analysis completed successfully for dataset:", dataset_type))
