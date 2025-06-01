library(tidyverse)
library(caret)
library(e1071)

# Helper function to calculate F-value for ANOVA feature selection
calculate_f_value_svm <- function(feature_values, class_labels) {
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

# Function to preprocess data for SVM Linear
preprocess_data_svm <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
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

  # Near zero variance filter
  nzv_features <- nearZeroVar(X_train_df, saveMetrics = FALSE, names = TRUE)
  cols_to_keep_after_nzv <- setdiff(colnames(X_train_df), nzv_features)
  
  # Apply NZV filter to both sets
  X_train_nzv_filtered <- X_train_df[, cols_to_keep_after_nzv, drop = FALSE]
  X_test_nzv_filtered <- X_test_df[, cols_to_keep_after_nzv, drop = FALSE]

  # ANOVA feature selection to reduce dimensionality
  num_features_to_select_anova <- min(5000, ncol(X_train_nzv_filtered))
  selected_feature_names_anova <- character(0)

  if (ncol(X_train_nzv_filtered) > 0 && num_features_to_select_anova > 0) {
      # Calculate F-values for all features
      f_values <- apply(X_train_nzv_filtered, 2, function(col) calculate_f_value_svm(col, y_train_factor))
      f_values[is.na(f_values)] <- 0

      if (all(f_values == 0)) {
        warning("SVM Linear: All ANOVA F-values are zero. Using all features post-NZV.")
        selected_feature_names_anova <- colnames(X_train_nzv_filtered)
      } else if (sum(f_values > 0) == 0) {
         warning("SVM Linear: No features with non-zero F-value. Using all features post-NZV.")
         selected_feature_names_anova <- colnames(X_train_nzv_filtered)
      } else {
        # Select top features based on F-values
        ordered_indices <- order(f_values, decreasing = TRUE)
        num_to_actually_select <- min(num_features_to_select_anova, sum(f_values > 0))
        selected_feature_names_anova <- colnames(X_train_nzv_filtered)[ordered_indices[1:num_to_actually_select]]
      }
  } else { 
      selected_feature_names_anova <- colnames(X_train_nzv_filtered) 
  }
  
  # Apply feature selection to both training and test sets
  X_train_selected_df <- as.data.frame(X_train_nzv_filtered[, selected_feature_names_anova, drop = FALSE])
  X_test_selected_df <- as.data.frame(X_test_nzv_filtered[, selected_feature_names_anova, drop = FALSE])

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

# Function to train an SVM Linear model
train_svmlinear_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  # Preprocess data
  processed_data <- preprocess_data_svm(X_train_df, X_test_df, y_train_df, y_test_df)
  
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

  # Train SVM model with linear kernel (scaling handled by SVM)
  svm_model_linear <- svm(X_train_selected_mat, y_train_factor, kernel = "linear", probability = TRUE, scale = TRUE)
  
  # Convert test data to matrix format
  X_test_selected_mat <- as.matrix(X_test_selected_df)
  
  # Create prediction probabilities
  y_pred_svm_factor <- predict(svm_model_linear, X_test_selected_mat, probability = TRUE)
  y_prob_values_matrix <- attr(y_pred_svm_factor, "probabilities")
  y_pred_prob_vector <- y_prob_values_matrix[, positive_class_label_factor]
  y_pred_numeric <- as.numeric(y_pred_svm_factor == positive_class_label_factor)

  # Create predictions dataframe
  predictions_output_df <- data.frame(
    y_test = y_test_numeric_output,
    y_pred = y_pred_numeric,
    y_pred_prob = as.numeric(y_pred_prob_vector)
  )

  # Calculate feature importance from linear SVM coefficients
  w <- t(svm_model_linear$coefs) %*% svm_model_linear$SV
  
  feature_importance_df <- data.frame(
    Feature = colnames(X_train_selected_mat),
    Value = as.numeric(w)
  )

  # Add any features not used in modeling with zero importance
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

  return(list(
    model_object = svm_model_linear,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
