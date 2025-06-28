library(tidyverse)
library(caret)
library(glmnet)

# ElasticNet specific preprocessing
preprocess_data_en <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  # Get name of target variable
  y_train_col_name <- names(y_train_df)[1]
  y_test_col_name <- names(y_test_df)[1]

  # Define expected levels consistently
  expected_levels <- c("Normal", "Tumor")
  
  # Convert to factor with explicit levels
  y_train_factor <- factor(y_train_df[[y_train_col_name]], levels = expected_levels)
  y_test_factor <- factor(y_test_df[[y_test_col_name]], levels = expected_levels)
  
  # Set the positive class label
  positive_class_label <- "Tumor"
  
  # Convert target to numeric where tumor is 1 and normal is 0
  y_train_numeric <- as.numeric(y_train_factor == positive_class_label)
  y_test_numeric <- as.numeric(y_test_factor == positive_class_label)
  
  # Store original feature names
  original_features <- colnames(X_train_df)

  # Convert to matrix for glmnet
  X_train_matrix <- as.matrix(X_train_df)
  X_test_matrix <- as.matrix(X_test_df)

  return(list(
    X_train_matrix = X_train_matrix,
    X_test_matrix = X_test_matrix,
    y_train_numeric = y_train_numeric,
    y_test_numeric = y_test_numeric,
    kept_feature_names = original_features,
    positive_class_label = positive_class_label
  ))
}

# Train ElasticNet model
train_elasticnet_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  # Preprocess the data
  processed_data <- preprocess_data_en(X_train_df, X_test_df, y_train_df, y_test_df)
  
  # Set variables from processed data
  X_train_matrix <- processed_data$X_train_matrix
  X_test_matrix <- processed_data$X_test_matrix
  y_train_numeric <- processed_data$y_train_numeric
  positive_class_label <- processed_data$positive_class_label

  # Define alpha values for tuning
  alpha_values <- c(0.01, 0.05, 0.1, 0.3, 0.5)
  
  # Init variables
  best_cvm <- Inf 
  best_alpha <- NA
  best_lambda <- NA
  
  # Assign each sample to one of 3 folds manually
  set.seed(42)
  fold_ids <- sample(1:3, size = length(y_train_numeric), replace = TRUE)

  # Tune alpha parameter
  for (alpha_val in alpha_values) {
    cv_fit <- cv.glmnet(
      x = X_train_matrix, 
      y = y_train_numeric, 
      family = "binomial", 
      alpha = alpha_val,
      foldid = fold_ids,
      type.measure = "deviance",
      standardize = TRUE
    )
    
    # Retrieve best loss
    current_best_cvm_for_alpha <- min(cv_fit$cvm, na.rm = TRUE)

    # Overwrite best values
    if (current_best_cvm_for_alpha < best_cvm) {
      best_cvm <- current_best_cvm_for_alpha
      best_alpha <- alpha_val
      best_lambda_idx <- which.min(cv_fit$cvm)
      best_lambda <- cv_fit$lambda[best_lambda_idx]
    }
  }

  # Train the final model
  final_model <- glmnet(
    x = X_train_matrix, 
    y = y_train_numeric, 
    family = "binomial", 
    alpha = best_alpha, 
    lambda = best_lambda,
    standardize = TRUE
  )
  
  # Create prediction probabilities vector on test set
  raw_preds_vector <- as.vector(predict(
    final_model, 
    newx = X_test_matrix, 
    type = "response"
  ))

  # Convert probabilities to class predictions
  y_pred_class <- ifelse(raw_preds_vector > 0.5, 1, 0)

  # Create dataframe for predictions output
  predictions_output_df <- data.frame(
    y_test = processed_data$y_test_numeric,
    y_pred = y_pred_class,
    y_pred_prob = raw_preds_vector
  )
  
  # Calculate feature importance from coefficients
  coef_matrix <- as.matrix(coef(final_model))
  
  # Create feature importance dataframe
  feature_importance_df <- data.frame(
    Feature = rownames(coef_matrix),
    Value = coef_matrix[,1]
  ) %>%
    filter(Feature != "(Intercept)")

  # Add any features not used in modeling with zero importance
  features_not_in_model <- setdiff(processed_data$kept_feature_names, feature_importance_df$Feature)
  if (length(features_not_in_model) > 0) {
    zero_importance_df <- data.frame(
      Feature = features_not_in_model,
      Value = 0.0
    )
    feature_importance_df <- rbind(feature_importance_df, zero_importance_df)
  }
  
  # Sort feature importance by absolute value
  feature_importance_df <- feature_importance_df %>%
    arrange(desc(abs(Value)))

  return(list(
    model_object = final_model,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
