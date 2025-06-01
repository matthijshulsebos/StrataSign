library(tidyverse)
library(caret)
library(lightgbm)

# Function to preprocess data for LightGBM
preprocess_data_lgbm <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
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

  # Near zero variance filter
  nzv_features <- nearZeroVar(X_train_df, saveMetrics = FALSE, names = TRUE)
  cols_to_keep_after_nzv <- setdiff(colnames(X_train_df), nzv_features)
  
  # Apply nzv filter to both sets  
  X_train_filtered <- X_train_df[, cols_to_keep_after_nzv, drop = FALSE]
  X_test_filtered <- X_test_df[, cols_to_keep_after_nzv, drop = FALSE]

  # Convert to matrix
  X_train_matrix <- as.matrix(X_train_filtered)
  X_test_matrix <- as.matrix(X_test_filtered)
  
  # Create LightGBM format
  lgb_train_data <- lgb.Dataset(data = X_train_matrix, label = y_train_numeric)

  return(list(
    lgb_train_data = lgb_train_data,
    X_test_matrix_for_predict = X_test_matrix,
    y_test_numeric_for_output = y_test_numeric,
    all_original_features = all_original_features_train,
    kept_feature_names = cols_to_keep_after_nzv
  ))
}

# Function to train a LightGBM model
train_lightgbm_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  # Preprocess data
  processed_data <- preprocess_data_lgbm(X_train_df, X_test_df, y_train_df, y_test_df)
  
  # Set preprocessed variables
  lgb_train_data <- processed_data$lgb_train_data
  X_test_matrix_for_predict <- processed_data$X_test_matrix_for_predict
  y_test_numeric_for_output <- processed_data$y_test_numeric_for_output
  all_original_features <- processed_data$all_original_features
  kept_feature_names <- processed_data$kept_feature_names

  message("LightGBM: Starting hyperparameter tuning with lgb.cv...")

  # Set seed
  set.seed(42)

  # Define hyperparameter grid
  tuning_grid_lgbm_cv <- expand.grid(
    learning_rate = c(0.1),              # Reduced from 2 to 1
    num_leaves = c(31),                   # Reduced from 2 to 1
    max_depth = c(6),                     # Reduced from 2 to 1
    feature_fraction = c(0.8),            # Reduced from 2 to 1
    bagging_fraction = c(0.8),            # Reduced from 2 to 1
    min_child_samples = c(10),            # Reduced from 2 to 1
    lambda_l1 = c(0),                     # Reduced from 2 to 1
    lambda_l2 = c(0),                     # Reduced from 2 to 1
    stringsAsFactors = FALSE
  )

  # Init variables for best parameters
  best_metric_val <- Inf
  best_params_cv <- NULL
  best_nrounds_cv <- -1

  # Maximum rounds for cv
  nrounds_cv_search <- 50  
  nfold_cv <- 3

  for (i in 1:nrow(tuning_grid_lgbm_cv)) {
    # Set training parameters
    current_params <- as.list(tuning_grid_lgbm_cv[i, ])
    params_for_lgb_cv <- list(
      objective = "binary",
      metric = "binary_logloss",
      learning_rate = current_params$learning_rate,
      num_leaves = current_params$num_leaves,
      max_depth = current_params$max_depth,
      feature_fraction = current_params$feature_fraction,
      bagging_fraction = current_params$bagging_fraction,
      min_child_samples = current_params$min_child_samples,
      lambda_l1 = current_params$lambda_l1,
      lambda_l2 = current_params$lambda_l2,
      bagging_freq = 5, 
      verbosity = -1,
      seed = 42,
      feature_pre_filter = FALSE
    )

    # Perform cv
    cv_results <- NULL
    tryCatch({
      cv_results <- lgb.cv(
        params = params_for_lgb_cv,
        data = lgb_train_data,
        nrounds = nrounds_cv_search,
        nfold = nfold_cv,
        early_stopping_rounds = 10,
        verbose = -1
      )
      
      # Force garbage collection after each CV run
      gc()
      
    }, error = function(e) {
      warning(paste("lgb.cv failed for param set", i, ":", e$message))
      # Force garbage collection even on error
      gc()
    })

    # Save best parameters
    if (!is.null(cv_results) && !is.null(cv_results$best_score)) {
      if (best_metric_val == Inf || cv_results$best_score < best_metric_val) {
        best_metric_val <- cv_results$best_score
        best_params_cv <- params_for_lgb_cv
        best_nrounds_cv <- cv_results$best_iter
      }
    }
  }
  
  # Train with best parameters
  final_lgbm_model <- lgb.train(
    params = best_params_cv,
    data = lgb_train_data,
    nrounds = best_nrounds_cv,
    verbose = -1
  )
  
  # Create prediction probabilities
  y_pred_prob_vector <- predict(final_lgbm_model, newdata = X_test_matrix_for_predict)

  # Convert predictions to classes
  y_pred_numeric <- ifelse(y_pred_prob_vector > 0.5, 1, 0)
  
  # Create predictions dataframe
  predictions_output_df <- data.frame(
    y_test = y_test_numeric_for_output,
    y_pred = y_pred_numeric,
    y_pred_prob = as.numeric(y_pred_prob_vector)
  )
  
  # Create feature importance object
  importance_matrix_lgb <- lgb.importance(model = final_lgbm_model)
  
  # Create feature importance dataframe with model importances
  if (nrow(importance_matrix_lgb) > 0) {
    feature_importance_df <- data.frame(
        Feature = importance_matrix_lgb$Feature,
        Value = importance_matrix_lgb$Gain,
        stringsAsFactors = FALSE
    )
  } else {
    feature_importance_df <- data.frame(Feature = character(0), Value = numeric(0), stringsAsFactors = FALSE)
  }
  
  # Add any features not used in modeling with zero importance  
  features_not_in_model <- setdiff(all_original_features, feature_importance_df$Feature)
  if (length(features_not_in_model) > 0) {
    zero_importance_df <- data.frame(
        Feature = features_not_in_model,
        Value = 0.0
    )
    feature_importance_df <- rbind(feature_importance_df, zero_importance_df)
  }
  
  feature_importance_df <- feature_importance_df %>%
    arrange(desc(abs(Value)))

  return(list(
    model_object = final_lgbm_model,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
