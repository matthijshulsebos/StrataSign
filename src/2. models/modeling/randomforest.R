library(tidyverse)
library(caret)
library(randomForest)

# Function to preprocess data for RandomForest
preprocess_data_rf <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  # Get name of target variable
  y_train_col_name <- names(y_train_df)[1]
  y_test_col_name <- names(y_test_df)[1]

  # Define expected levels
  expected_levels <- c("Normal", "Tumor")
  
  # Convert target variables to factors with explicit levels
  y_train_factor <- factor(y_train_df[[y_train_col_name]], levels = expected_levels)
  y_test_factor <- factor(y_test_df[[y_test_col_name]], levels = expected_levels)

  # Set the positive class label
  positive_class_label <- "Tumor"
  
  # Convert target to numeric where tumor is 1 and normal is 0
  y_train_numeric <- as.numeric(y_train_factor == positive_class_label)
  y_test_numeric <- as.numeric(y_test_factor == positive_class_label)

  # No near zero variance filtering for RandomForest
  X_train_processed <- X_train_df
  X_test_processed <- X_test_df

  return(list(
    X_train_processed = X_train_processed,
    X_test_processed = X_test_processed,
    y_train_factor = y_train_factor,
    y_test_numeric = y_test_numeric,
    positive_class_label = positive_class_label
  ))
}

# Train a Random Forest model
train_randomforest_model <- function(X_train_df, X_test_df, y_train_df, y_test_df) {
  
  # Preprocess data
  processed_data <- preprocess_data_rf(X_train_df, X_test_df, y_train_df, y_test_df)
  
  # Set preprocessed variables
  X_train_processed <- processed_data$X_train_processed
  X_test_processed <- processed_data$X_test_processed
  y_train_factor <- processed_data$y_train_factor
  y_test_numeric_for_output <- processed_data$y_test_numeric
  positive_class_label <- processed_data$positive_class_label
  
  # Set seed
  set.seed(42)

  # Define tuning grid for ranger
  num_cols_train <- ncol(X_train_processed)
  
  # Use standard mtry values for classification
  sqrt_features <- floor(sqrt(num_cols_train))
  mtry_values <- c(sqrt_features)

  # Create tuning grid
  tuning_grid_rf <- expand.grid(
    mtry = mtry_values,
    splitrule = "gini",
    min.node.size = c(1, 3, 5)
  )

  # Set up cv
  control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE
  )

  # Perform hyperparameter tuning
  rf_tune <- train(
    x = X_train_processed,
    y = y_train_factor,
    method = "ranger",
    num.trees = 200,
    importance = "permutation",
    tuneGrid = tuning_grid_rf,
    trControl = control,
    metric = "ROC"
  )

  # Train final model with best parameters
  final_rf_model <- randomForest(
    x = X_train_processed,
    y = y_train_factor,
    ntree = 200,
    mtry = rf_tune$bestTune$mtry,
    nodesize = rf_tune$bestTune$min.node.size,
    importance = TRUE
  )
  
  # Create prediction probabilities
  y_pred_prob_vector <- predict(final_rf_model, newdata = X_test_processed, type = "prob")[, positive_class_label]
  
  # Convert predictions to classes
  y_pred_factor <- ifelse(y_pred_prob_vector > 0.5, positive_class_label, levels(y_train_factor)[1])
  y_pred_factor <- factor(y_pred_factor, levels = levels(y_train_factor))
  y_pred_numeric <- as.numeric(y_pred_factor == positive_class_label)
  
  # Create predictions dataframe
  predictions_output_df <- data.frame(
    y_test = y_test_numeric_for_output,
    y_pred = y_pred_numeric,
    y_pred_prob = as.numeric(y_pred_prob_vector)
  )
  
  # Calculate feature importance
  feature_imp_matrix <- importance(final_rf_model, type = 1, scale = FALSE) 
  
  # Create feature importance dataframe
  feature_importance_df <- data.frame(
    Feature = rownames(feature_imp_matrix),
    Value = feature_imp_matrix[, "MeanDecreaseAccuracy"]
  ) %>%
    arrange(desc(Value))

  # Sort feature importance by absolute value
  feature_importance_df <- feature_importance_df %>%
    arrange(desc(abs(Value)))

  return(list(
    model_object = final_rf_model,
    predictions_df = predictions_output_df,
    raw_feature_importance = feature_importance_df
  ))
}
