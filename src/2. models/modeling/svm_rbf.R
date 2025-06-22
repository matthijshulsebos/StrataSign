library(tidyverse)
library(caret)
library(e1071)
library(pROC)

# Helper function to calculate F-value for ANOVA feature selection
calculate_f_value_svmrbf <- function(feature_values, class_labels) {
  # Skip features with near zero variance
  if (sd(feature_values, na.rm = TRUE) < 1e-6) {
    return(0)
  }
  tryCatch(
    expr = {
      anova_result <- aov(feature_values ~ class_labels)
      return(summary(anova_result)[[1]]$"F value"[1])
    },
    error = function(e) {
      return(0)
    }
  )
}

# Function to preprocess data for SVM RBF
preprocess_data_svmrbf <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  # Set training data as dataframes
  X_train_df <- as.data.frame(X_train_df)
  X_test_df <- as.data.frame(X_test_df)

  # Get name of target variable
  y_train_col_name <- names(y_train_df)[1]
  y_test_col_name <- names(y_test_df)[1]

  # Store original feature names
  all_original_features_train <- colnames(X_train_df)

  # Define expected levels consistently
  expected_levels <- c("Normal", "Tumor")
  
  # Convert target variables to factors
  y_train_factor <- factor(y_train_df[[y_train_col_name]], levels = expected_levels)
  y_test_factor <- factor(y_test_df[[y_test_col_name]], levels = expected_levels)
  
  # Set the positive class label
  positive_class_label_factor <- "Tumor"
  
  # Convert target to numeric where tumor is 1 and normal is 0
  y_test_numeric_output <- as.numeric(y_test_factor == positive_class_label_factor)

  # Remove NZV filtering: do not filter features by near zero variance, rely on ANOVA feature selection only
  # ANOVA feature selection to reduce dimensionality
  num_features_to_select_anova <- min(1000, ncol(X_train_df))
  selected_feature_names_anova <- character(0)

  if (ncol(X_train_df) > 0 && num_features_to_select_anova > 0) {
      # Calculate F-values for all features
      f_values <- apply(X_train_df, 2, function(col) calculate_f_value_svmrbf(col, y_train_factor))
      f_values[is.na(f_values)] <- 0

      if (all(f_values == 0)) {
        warning("All ANOVA F-values are zero. Using all features.")
        selected_feature_names_anova <- colnames(X_train_df)
      } else if (sum(f_values > 0) == 0) {
         warning("No features with non-zero F-value.")
         selected_feature_names_anova <- colnames(X_train_df)
      } else {
        # Select top features based on F-values
        ordered_indices <- order(f_values, decreasing = TRUE)
        num_to_actually_select <- min(num_features_to_select_anova, sum(f_values > 0))
        selected_feature_names_anova <- colnames(X_train_df)[ordered_indices[1:num_to_actually_select]]
      }
  } else { 
      selected_feature_names_anova <- colnames(X_train_df) 
  }
  
  # Apply feature selection to both training and test sets
  X_train_selected_df <- as.data.frame(X_train_df[, selected_feature_names_anova, drop = FALSE])
  X_test_selected_df <- as.data.frame(X_test_df[, selected_feature_names_anova, drop = FALSE])

  return(list(
    X_train_processed = X_train_selected_df,
    X_test_processed = X_test_selected_df,
    y_train_factor = y_train_factor,
    y_test_numeric = y_test_numeric_output,
    selected_feature_names = selected_feature_names_anova,
    all_original_features = all_original_features_train,
    positive_class_label = positive_class_label_factor
  ))
}

# Function to train an SVM RBF model
train_svmrbf_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  # Preprocess data
  processed_data <- preprocess_data_svmrbf(X_train_df, X_test_df, y_train_df, y_test_df)
  
  # Set preprocessed variables
  X_train_selected_df <- processed_data$X_train_processed
  X_test_selected_df <- processed_data$X_test_processed
  y_train_factor <- processed_data$y_train_factor
  y_test_numeric_output <- processed_data$y_test_numeric
  selected_feature_names <- processed_data$selected_feature_names
  all_original_features <- processed_data$all_original_features
  positive_class_label_factor <- processed_data$positive_class_label
  
  # Convert to matrix format for SVM
  X_train_selected_mat <- as.matrix(X_train_selected_df)
  X_test_selected_mat <- as.matrix(X_test_selected_df)

  # Bootstrapping parameters
  n_bootstraps <- 25
  n_samples_train <- nrow(X_train_selected_df)
  
  # Initialize variables for bootstrapping
  all_importance_list <- list()
  best_auc <- -Inf
  best_model <- NULL
  best_predictions <- NULL

  for (i in 1:n_bootstraps) {
    if (i %% 10 == 0 || i == n_bootstraps) {
      message(paste("  Bootstrap iteration", i, "/", n_bootstraps))
    }

    # Bootstrap sampling
    set.seed(42 + i)
    boot_indices <- sample(1:n_samples_train, n_samples_train, replace = TRUE)
    X_train_boot <- X_train_selected_df[boot_indices, , drop = FALSE]
    y_train_boot <- y_train_factor[boot_indices]

    # Train SVM model on bootstrap sample
    svm_model_boot <- svm(as.matrix(X_train_boot), y_train_boot, kernel = "radial", probability = TRUE, scale = TRUE)

    # Predict on test data
    y_pred_svm_factor <- predict(svm_model_boot, X_test_selected_mat, probability = TRUE)
    y_prob_values_matrix <- attr(y_pred_svm_factor, "probabilities")
    y_pred_prob_vector <- y_prob_values_matrix[, positive_class_label_factor]

    # Calculate AUC for the current bootstrap model
    roc_obj <- roc(response = y_test_numeric_output, predictor = y_pred_prob_vector, quiet = TRUE, levels = c(0, 1), direction = "<")
    current_auc <- as.numeric(roc_obj$auc)

    # Track the best model
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_model <- svm_model_boot
      best_predictions <- data.frame(
        y_test = y_test_numeric_output,
        y_pred = as.numeric(y_pred_svm_factor == positive_class_label_factor),
        y_pred_prob = as.numeric(y_pred_prob_vector)
      )
    }

    # Permutation-based feature importance for the current bootstrap
    perm_importance_values <- numeric(length(selected_feature_names))
    names(perm_importance_values) <- selected_feature_names

    for (feature_name in selected_feature_names) {
      X_test_permuted_df <- X_test_selected_df
      original_col_values <- X_test_permuted_df[[feature_name]]
      X_test_permuted_df[[feature_name]] <- sample(original_col_values)
      X_test_permuted_mat <- as.matrix(X_test_permuted_df)

      perm_pred_svm_factor <- predict(svm_model_boot, X_test_permuted_mat, probability = TRUE)
      perm_prob_matrix <- attr(perm_pred_svm_factor, "probabilities")
      perm_pred_prob_vector <- perm_prob_matrix[, positive_class_label_factor]

      permuted_roc_obj <- roc(response = y_test_numeric_output, predictor = perm_pred_prob_vector, quiet = TRUE, levels = c(0, 1), direction = "<")
      permuted_auc <- as.numeric(permuted_roc_obj$auc)
      perm_importance_values[feature_name] <- current_auc - permuted_auc
    }

    # Save importance values for this bootstrap iteration
    all_importance_list[[i]] <- perm_importance_values
  }

  # Aggregate feature importance across bootstraps
  importance_matrix <- do.call(rbind, all_importance_list)
  avg_importance <- colMeans(importance_matrix, na.rm = TRUE)

  feature_importance_df <- data.frame(
    Feature = names(avg_importance),
    Value = avg_importance
  )

  # Add features not used in modeling with zero importance
  features_not_in_model <- setdiff(all_original_features, feature_importance_df$Feature)
  if (length(features_not_in_model) > 0) {
    zero_importance_df <- data.frame(
      Feature = features_not_in_model,
      Value = 0
    )
    feature_importance_df <- rbind(feature_importance_df, zero_importance_df)
  }

  # Sort feature importance by absolute value
  feature_importance_df <- feature_importance_df %>% arrange(desc(abs(Value)))

  # Retrain the final model on the entire dataset
  final_model <- svm(
    as.matrix(X_train_selected_df), y_train_factor,
    kernel = "radial",
    probability = TRUE,
    scale = TRUE
  )

  # Use the final model to make predictions on the test set
  final_predictions <- predict(final_model, X_test_selected_mat, probability = TRUE)
  final_prob_matrix <- attr(final_predictions, "probabilities")
  final_y_pred_prob_vector <- final_prob_matrix[, positive_class_label_factor]
  final_y_pred_numeric <- as.numeric(final_predictions == positive_class_label_factor)

  # Update predictions dataframe with final model predictions
  predictions_output_df <- data.frame(
    y_test = y_test_numeric_output,
    y_pred = final_y_pred_numeric,
    y_pred_prob = as.numeric(final_y_pred_prob_vector)
  )

  # Return the final model and updated predictions
  return(list(
    model_object = final_model,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
