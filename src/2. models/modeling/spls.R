library(tidyverse)
library(caret)
library(spls)

# Function to preprocess data for spls
preprocess_data_spls <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
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
  
  # Apply nzv filter to training set
  X_train_processed <- X_train_df[, cols_to_keep_after_nzv, drop = FALSE]
  
  # Align test set to match training features
  X_test_aligned <- data.frame(matrix(0.0, nrow = nrow(X_test_df), ncol = length(cols_to_keep_after_nzv)))
  colnames(X_test_aligned) <- cols_to_keep_after_nzv
  
  # Filter test set to keep only columns that are also in training set
  common_cols_in_test <- intersect(cols_to_keep_after_nzv, colnames(X_test_df))
  if(length(common_cols_in_test) > 0) {
      X_test_aligned[, common_cols_in_test] <- X_test_df[, common_cols_in_test, drop = FALSE]
  }
  X_test_processed <- X_test_aligned
  
  # Set selected feature names
  selected_feature_names <- cols_to_keep_after_nzv

  return(list(
    X_train_processed = X_train_processed,
    X_test_processed = X_test_processed,
    y_train_factor = y_train_factor,
    y_test_numeric = y_test_numeric_output,
    selected_feature_names = selected_feature_names,
    all_original_features = all_original_features_train,
    positive_class_label = positive_class_label_factor
  ))
}

# Function to train an spls model
train_spls_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  # Preprocess data
  processed_data <- preprocess_data_spls(X_train_df, X_test_df, y_train_df, y_test_df)
  
  # Set preprocessed variables
  X_train_processed <- processed_data$X_train_processed
  X_test_processed <- processed_data$X_test_processed 
  y_train_factor <- processed_data$y_train_factor
  y_test_numeric_output <- processed_data$y_test_numeric 
  features_for_model <- processed_data$selected_feature_names
  all_original_features <- processed_data$all_original_features
  positive_class_label_factor <- processed_data$positive_class_label

  # Convert to matrix format for sPLS
  X_train_matrix <- as.matrix(X_train_processed)
  X_test_matrix <- as.matrix(X_test_processed)
  y_train_numeric <- as.numeric(y_train_factor == positive_class_label_factor)
  
  # Set seed
  set.seed(42)

  # Define hyperparameter grid
  tuning_grid_spls <- expand.grid(
    K = c(2, 3, 5),
    eta = c(0.3, 0.5, 0.7),
    stringsAsFactors = FALSE
  )

  # Manual cv hyperparameter tuning
  best_error <- Inf
  best_params <- NULL
  nfolds <- 3
  
  # Set fold indices
  fold_indices <- sample(rep(1:nfolds, length.out = nrow(X_train_matrix)))
  
  for (i in 1:nrow(tuning_grid_spls)) {
    current_params <- tuning_grid_spls[i, ]
    cv_errors <- numeric(nfolds)
    
    for (fold in 1:nfolds) {
      # Create folds
      train_idx <- which(fold_indices != fold)
      val_idx <- which(fold_indices == fold)
      
      X_train_cv <- X_train_matrix[train_idx, , drop = FALSE]
      y_train_cv <- y_train_numeric[train_idx]
      X_val_cv <- X_train_matrix[val_idx, , drop = FALSE]
      y_val_cv <- y_train_numeric[val_idx]
      
      # Train model with current parameters
      tryCatch({
        cv_model <- spls::spls(
          X_train_cv, y_train_cv,
          K = current_params$K,
          eta = current_params$eta,
          scale.x = TRUE,
          scale.y = TRUE,
          select = "pls2",
          fit = "kernelpls"
        )
        
        # Make predictions and calculate error
        cv_pred <- predict(cv_model, X_val_cv)
        cv_pred_binary <- ifelse(cv_pred > 0.5, 1, 0)
        cv_errors[fold] <- mean(cv_pred_binary != y_val_cv)
        
      }, error = function(e) {
        cv_errors[fold] <- 1
      })
    }
    
    # Calculate average cv error
    avg_cv_error <- mean(cv_errors)
    
    # Update best parameters if this is better
    if (avg_cv_error < best_error) {
      best_error <- avg_cv_error
      best_params <- current_params
    }
  }

  # Train final model with best parameters
  final_model <- NULL
  tryCatch({
    final_model <- spls::spls(
      X_train_matrix,
      y_train_numeric,
      K = best_params$K,
      eta = best_params$eta,
      scale.x = TRUE,
      scale.y = TRUE,
      select = "pls2", 
      fit = "simpls"
    )
  }, error = function(e) {
    stop("Error training sPLS model.")
  })

  # Create prediction probabilities
  y_pred_raw_scores <- predict(final_model, X_test_matrix)
  y_pred_numeric <- ifelse(y_pred_raw_scores > 0.5, 1, 0)

  # Values can be outside 0 and 1 range so set them within them
  y_pred_prob <- pmin(pmax(as.numeric(y_pred_raw_scores), 0), 1)
  
  # Create predictions dataframe
  predictions_output_df <- data.frame(
    y_test = y_test_numeric_output,
    y_pred = y_pred_numeric,
    y_pred_prob = as.numeric(y_pred_prob)
  )

  # Calculate feature importance using regression coefficients ie betahat
  feature_importance_df <- data.frame(Feature = all_original_features, Value = 0.0, stringsAsFactors = FALSE)
  
  # Extract feature importance from final regression coefficients
  if (!is.null(final_model$betahat) && is.matrix(final_model$betahat)) {
    importance_values <- final_model$betahat[, 1]
    names(importance_values) <- colnames(X_train_matrix)
    
    # Overwrite 0 feature importance for modeled features
    for (feature in names(importance_values)) {
      if (feature %in% feature_importance_df$Feature) {
        feature_importance_df$Value[feature_importance_df$Feature == feature] <- importance_values[feature]
      }
    }
  } else {
    stop("Betahat not available.")
  }
  
  # Sort by absolute value of importance
  feature_importance_df <- feature_importance_df %>%
    arrange(desc(abs(Value)))

  return(list(
    model_object = final_model,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
