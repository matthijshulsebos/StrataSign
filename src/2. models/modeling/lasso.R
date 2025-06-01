library(tidyverse)
library(caret)
library(glmnet)

# LASSO specific preprocessing function
preprocess_data_lasso <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  # Get name of target variable
  y_train_col_name <- names(y_train_df)[1]
  y_test_col_name <- names(y_test_df)[1]

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
  
  # Store original feature names
  original_features_train <- colnames(X_train_df)

  # Near zero variance filtering
  nzv_to_remove_names <- nearZeroVar(X_train_df, saveMetrics = FALSE, names = TRUE)

  # Feature cols to keep after removing nzv features
  features_to_keep <- setdiff(original_features_train, nzv_to_remove_names)

  # Apply nzv filter to both sets
  X_train_filtered <- X_train_df[, features_to_keep, drop = FALSE]
  X_test_filtered <- X_test_df[, features_to_keep, drop = FALSE]

  # Convert feature data to matrix format for glmnet
  X_train_matrix <- as.matrix(X_train_filtered)
  X_test_matrix <- as.matrix(X_test_filtered)

  return(list(
    X_train_matrix = X_train_matrix,
    X_test_matrix = X_test_matrix,
    y_train_numeric = y_train_numeric,
    y_test_numeric = y_test_numeric,
    kept_feature_names = features_to_keep,
    all_original_feature_names = original_features_train
  ))
}


# Train LASSO model with bootstrapping
train_lasso_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  message("Preprocessing data for LASSO...")
  processed_data <- preprocess_data_lasso(X_train_df, X_test_df, y_train_df, y_test_df)
  
  # Set preprocessed variables
  X_train_matrix <- processed_data$X_train_matrix
  X_test_matrix <- processed_data$X_test_matrix
  y_train_numeric <- processed_data$y_train_numeric
  y_test_numeric_output <- processed_data$y_test_numeric
  kept_feature_names <- processed_data$kept_feature_names
  all_original_feature_names <- processed_data$all_original_feature_names

  message("Training LASSO model with bootstrapping...")

  # Set bootstrapping parameters
  n_bootstraps <- 100
  n_samples_train <- nrow(X_train_matrix)
  
  # Init empty variables for results
  all_coefs_list <- list() 
  
  # Track best model across bootstraps
  best_deviance <- Inf
  best_model <- NULL
  best_predictions <- NULL
  
  # Create coefficient names which are features and intercept
  all_possible_coef_names <- c("(Intercept)", colnames(X_train_matrix))

  for (i in 1:n_bootstraps) {
    # Bootstrapping takes ages to run so give some progress indication
    if (i %% 10 == 0 || i == n_bootstraps) {
        message(paste("  Bootstrap iteration", i, "/", n_bootstraps))
    }

    # Shift sampling for each bootstrap iteration
    set.seed(42 + i)
    
    # Perform bootstrap sampling
    boot_indices <- sample(1:n_samples_train, n_samples_train, replace = TRUE)
    X_train_boot_matrix <- X_train_matrix[boot_indices, , drop = FALSE]
    y_train_boot_numeric <- y_train_numeric[boot_indices]
    
    # Fit LASSO model with cv on the bootstrap sample
    cv_fit_boot <- cv.glmnet(X_train_boot_matrix, y_train_boot_numeric, 
                             alpha = 1,
                             family = "binomial",
                             nfolds = 3,  # Changed from 5 to 3
                             type.measure = "deviance",
                             standardize = TRUE
    )
    
    # Check if this is the best model (lower deviance is better)
    current_deviance <- min(cv_fit_boot$cvm, na.rm = TRUE)
    if (current_deviance < best_deviance) {
      best_deviance <- current_deviance
      best_model <- cv_fit_boot
      best_predictions <- predict(cv_fit_boot, X_test_matrix, s = "lambda.min", type = "response")[,1]
    }
    
    # Get all non-zero coefficients
    coef_vector_boot_sparse <- coef(cv_fit_boot, s = "lambda.min")

    # Create a vector for all coefficients and set all values to 0
    coef_vector_boot_dense <- rep(0, length(all_possible_coef_names))

    # Assign names to the coefficient vector
    names(coef_vector_boot_dense) <- all_possible_coef_names
    
    # Get common coefficients which is the same as sparse coefficients
    common_coef_names <- intersect(names(coef_vector_boot_dense), rownames(coef_vector_boot_sparse))

    # Assign the non zero coefficients to the 0 vector
    coef_vector_boot_dense[common_coef_names] <- coef_vector_boot_sparse[common_coef_names, 1]

    # Save coefficients for this bootstrap run
    all_coefs_list[[i]] <- coef_vector_boot_dense
  }
  message("Bootstrap loop finished.")
  
  # Use predictions from the best model
  y_pred_prob_vector <- best_predictions
  y_pred_numeric <- ifelse(y_pred_prob_vector > 0.5, 1, 0)
  
  # Create a dataframe for predictions
  predictions_output_df <- data.frame(
    y_test = y_test_numeric_output,
    y_pred = y_pred_numeric,
    y_pred_prob = as.numeric(y_pred_prob_vector)
  )
  
  message("Calculating feature importance for LASSO...")

  # Transpose the list of coefficient vectors into a matrix
  coef_matrix_all_bootstraps <- do.call(rbind, all_coefs_list)

  # Calculate the average of non zero coefficients across all bootstraps
  avg_non_zero_coeffs <- apply(coef_matrix_all_bootstraps, 2, function(feature_coeffs_across_bootstraps) {
    non_zero_vals <- feature_coeffs_across_bootstraps[feature_coeffs_across_bootstraps != 0]
    if (length(non_zero_vals) > 0) {
      mean(non_zero_vals)
    } else {
      0
    }
  })
  
  # Create feature importance dataframe
  feature_importance_df <- data.frame(
    Feature = names(avg_non_zero_coeffs),
    Value = avg_non_zero_coeffs
  ) %>% filter(Feature != "(Intercept)")

  # Add any features not used in modeling with zero importance
  features_not_in_model <- setdiff(all_original_feature_names, feature_importance_df$Feature)
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


  # Return the best model from bootstrapping
  return(list(
    model_object = best_model,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
