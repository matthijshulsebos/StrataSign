library(tidyverse)
library(caret)
library(xgboost)

# Function to preprocess data for XGBoost
preprocess_data_xgb <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
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
  
  # Convert target variables to factors with explicit levels
  y_train_factor <- factor(y_train_df[[y_train_col_name]], levels = expected_levels)
  y_test_factor <- factor(y_test_df[[y_test_col_name]], levels = expected_levels)

  # Set the positive class label
  positive_class_label <- "Tumor"
  
  # Convert target to numeric where tumor is 1 and normal is 0
  y_train_numeric <- as.numeric(y_train_factor == positive_class_label)
  y_test_numeric <- as.numeric(y_test_factor == positive_class_label)

  # No near zero variance filtering for XGBoost
  X_train_matrix <- as.matrix(X_train_df)
  X_test_matrix <- as.matrix(X_test_df)

  # Create DMatrices for xgb.train
  dtrain <- xgb.DMatrix(data = X_train_matrix, label = y_train_numeric)
  dtest <- xgb.DMatrix(data = X_test_matrix, label = y_test_numeric)

  return(list(
    X_train_matrix_for_caret = X_train_matrix,
    y_train_factor_for_caret = y_train_factor,
    dtrain_for_xgb = dtrain,
    dtest_for_xgb = dtest,
    y_test_numeric_for_output = y_test_numeric,
    all_original_features = all_original_features_train,
    kept_feature_names = all_original_features_train,
    positive_class_label = positive_class_label
  ))
}

# Function to train an XGBoost model
train_xgboost_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  # Preprocess data
  processed_data <- preprocess_data_xgb(X_train_df, X_test_df, y_train_df, y_test_df)
  
  # Set preprocessed variables
  X_train_matrix_for_caret <- processed_data$X_train_matrix_for_caret
  y_train_factor_for_caret <- processed_data$y_train_factor_for_caret
  dtrain_for_xgb <- processed_data$dtrain_for_xgb
  dtest_for_xgb <- processed_data$dtest_for_xgb
  y_test_numeric_for_output <- processed_data$y_test_numeric_for_output
  all_original_features <- processed_data$all_original_features
  kept_feature_names <- processed_data$kept_feature_names
  positive_class_label <- processed_data$positive_class_label

  message("XGBoost: Starting hyperparameter tuning...")
  
  # Set seed
  set.seed(42)

  # Define hyperparameter grid
  tuning_grid_xgb <- expand.grid(
    nrounds = c(200, 300),
    max_depth = c(6, 8, 10),
    eta = c(0.05, 0.1),
    gamma = c(0),
    colsample_bytree = c(1.0),
    min_child_weight = c(1),
    subsample = c(1.0),
    reg_alpha = c(0),
    reg_lambda = c(0, 0.01)
  )

  best_auc <- -Inf
  best_params <- NULL
  best_nrounds <- NULL
  
  for (i in seq_len(nrow(tuning_grid_xgb))) {
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = tuning_grid_xgb$max_depth[i],
      eta = tuning_grid_xgb$eta[i],
      gamma = tuning_grid_xgb$gamma[i],
      colsample_bytree = tuning_grid_xgb$colsample_bytree[i],
      min_child_weight = tuning_grid_xgb$min_child_weight[i],
      subsample = tuning_grid_xgb$subsample[i]
    )
    
    # Cross validation
    cv_results <- xgb.cv(
      params = params,
      data = dtrain_for_xgb,
      nrounds = tuning_grid_xgb$nrounds[i],
      nfold = 3,
      early_stopping_rounds = 15,
      verbose = FALSE,
      showsd = FALSE
    )
    
    # Get best iteration and score
    best_iter <- cv_results$best_iteration
    best_score <- cv_results$evaluation_log$test_auc_mean[best_iter]
    
    if (best_score > best_auc) {
      best_auc <- best_score
      best_params <- params
      best_nrounds <- best_iter
    }
  }
  
  # Train final XGBoost model with best parameters
  final_xgb_model <- xgb.train(
    params = best_params,
    data = dtrain_for_xgb,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  # Create prediction probabilities
  y_pred_prob_vector <- predict(final_xgb_model, newdata = dtest_for_xgb)
  y_pred_numeric <- ifelse(y_pred_prob_vector > 0.5, 1, 0)
  
  # Create predictions dataframe
  predictions_output_df <- data.frame(
    y_test = y_test_numeric_for_output,
    y_pred = y_pred_numeric,
    y_pred_prob = as.numeric(y_pred_prob_vector)
  )
  
  # Calculate feature importance using gain
  importance_matrix <- xgb.importance(feature_names = kept_feature_names, model = final_xgb_model)
  
  # Create feature importance dataframe
  feature_importance_df <- data.frame(
      Feature = importance_matrix$Feature,
      Value = importance_matrix$Gain
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
  feature_importance_df <- feature_importance_df %>%
    arrange(desc(abs(Value)))

  return(list(
    model_object = final_xgb_model,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
