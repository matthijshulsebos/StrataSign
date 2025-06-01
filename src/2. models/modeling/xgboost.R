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
  
  # Convert target variables to factors
  y_train_factor_caret <- factor(y_train_df[[y_train_col_name]], levels = expected_levels)
  y_test_factor_eval <- factor(y_test_df[[y_test_col_name]], levels = expected_levels)
  
  # Set the positive class label
  positive_class_label <- "Tumor"
  
  # Convert target to numeric where tumor is 1 and normal is 0
  y_train_numeric_xgb <- as.numeric(y_train_factor_caret == positive_class_label)
  y_test_numeric_output <- as.numeric(y_test_factor_eval == positive_class_label)

  # Near zero variance filter
  nzv_features <- nearZeroVar(X_train_df, saveMetrics = FALSE, names = TRUE)
  cols_to_keep_names <- setdiff(colnames(X_train_df), nzv_features)
  
  # Apply NZV filter to both sets
  X_train_filtered <- X_train_df[, cols_to_keep_names, drop = FALSE]
  X_test_filtered <- X_test_df[, cols_to_keep_names, drop = FALSE]

  # Convert to matrix for XGBoost
  X_train_matrix <- as.matrix(X_train_filtered)
  X_test_matrix <- as.matrix(X_test_filtered)
  
  # Create DMatrices for xgb.train
  dtrain <- xgb.DMatrix(data = X_train_matrix, label = y_train_numeric_xgb)
  dtest <- xgb.DMatrix(data = X_test_matrix, label = y_test_numeric_output)

  return(list(
    X_train_matrix_for_caret = X_train_matrix,
    y_train_factor_for_caret = y_train_factor_caret,
    dtrain_for_xgb = dtrain,
    dtest_for_xgb = dtest,
    y_test_numeric_for_output = y_test_numeric_output,
    all_original_features = all_original_features_train,
    kept_feature_names = cols_to_keep_names,
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
    nrounds = c(100),
    max_depth = c(3, 6),
    eta = c(0.1, 0.3),
    gamma = c(0),
    colsample_bytree = c(0.8),
    min_child_weight = c(1),
    subsample = c(0.8)
  )

  # Configure cross-validation
  control <- trainControl(
    method = "cv",
    number = 3,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE,
    allowParallel = TRUE
  )

  # Tune hyperparameters using caret
  xgb_tune <- train(
    x = X_train_matrix_for_caret,
    y = y_train_factor_for_caret,
    method = "xgbTree",
    trControl = control,
    tuneGrid = tuning_grid_xgb,
    metric = "ROC",
    verbose = FALSE
  )

  # Set best parameters for final model
  best_params_list <- list(
    max_depth = xgb_tune$bestTune$max_depth,
    eta = xgb_tune$bestTune$eta,
    gamma = xgb_tune$bestTune$gamma,
    colsample_bytree = xgb_tune$bestTune$colsample_bytree,
    min_child_weight = xgb_tune$bestTune$min_child_weight,
    subsample = xgb_tune$bestTune$subsample,
    objective = "binary:logistic",
    eval_metric = "logloss"
  )
  
  # Train final XGBoost model
  final_xgb_model <- xgb.train(
    params = best_params_list,
    data = dtrain_for_xgb,
    nrounds = xgb_tune$bestTune$nrounds,
    watchlist = list(train = dtrain_for_xgb, test = dtest_for_xgb),
    early_stopping_rounds = 10,
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
  
  # Calculate feature importance using Gain
  importance_matrix <- xgb.importance(feature_names = kept_feature_names, model = final_xgb_model)
  
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
