library(tidyverse)
library(pROC)
library(caret)

# Calculate performance metrics for binary classification
calculate_performance_metrics <- function(y_true, y_pred, y_pred_prob) {
  
  # Create confusion matrix of caret that has a bunch of useful metrics
  cm <- confusionMatrix(factor(y_pred, levels = c(0, 1)), 
                       factor(y_true, levels = c(0, 1)), 
                       positive = "1")
  
  # ROC AUC
  roc_obj <- roc(y_true, y_pred_prob, quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  
  # Precision, recall, f1 extraction from confusion matrix
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"] 
  specificity <- cm$byClass["Specificity"]
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  # Replace missing values with zeros
  precision <- ifelse(is.na(precision), 0, precision)
  recall <- ifelse(is.na(recall), 0, recall)
  specificity <- ifelse(is.na(specificity), 0, specificity)
  f1_score <- ifelse(is.na(f1_score), 0, f1_score)
  
  # Additional metrics
  accuracy <- cm$overall["Accuracy"]
  balanced_accuracy <- (recall + specificity) / 2
  
  # Prediction class counts
  tn <- cm$table[1,1]
  fp <- cm$table[1,2] 
  fn <- cm$table[2,1]
  tp <- cm$table[2,2]
  
  # Return a dataframe with the metrics
  return(data.frame(
    Accuracy = round(accuracy, 4),
    Precision = round(precision, 4),
    Recall = round(recall, 4),
    Specificity = round(specificity, 4),
    F1_Score = round(f1_score, 4),
    AUC = round(auc_value, 4),
    Balanced_Accuracy = round(balanced_accuracy, 4),
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn
  ))
}


# Retrieves model metadata such as subsetting information
extract_model_info <- function(file_path) {
  # Extract relative path components
  path_parts <- strsplit(file_path, "[/\\\\]")[[1]]
  
  # Find the model based on path
  model_idx <- which(path_parts %in% c("xgboost", "svmrbf", "svmlinear", "spls", 
                                       "randomforest", "lightgbm", "lasso", "elasticnet"))
  
  # If path contains a model name then extract its data and metadata
  if (length(model_idx) > 0) {
    algorithm <- path_parts[model_idx[1]]
    
    # Adjust SVM names for consistency
    if (algorithm == "svmrbf") algorithm <- "svm_rbf"
    if (algorithm == "svmlinear") algorithm <- "svm_linear"
    
    # Find base models directory
    models_idx <- which(path_parts == "2. models")
    
    # Given that the model is in a subdirectory of 2. models
    if (length(models_idx) > 0 && model_idx[1] > models_idx[1] + 1) {
      # Read carefully it is modelS for the top directory where we want to start and model for the low level model directory
      start_idx <- models_idx[1] + 1
      end_idx <- model_idx[1] - 1
      
      # Make sure the end index is not before the start index
      if (end_idx >= start_idx) {
        # Get all path parts between models and model directory
        components <- path_parts[start_idx:end_idx]
        
        # Structure is dataset/cell_type/gene_set
        dataset <- if (length(components) >= 1) components[1] else "unknown"
        cell_type <- if (length(components) >= 2) components[2] else "unknown"
        gene_set <- if (length(components) >= 3) components[3] else "unknown"
        
        return(list(
          algorithm = algorithm, 
          dataset = dataset,
          cell_type = cell_type,
          gene_set = gene_set
        ))
      }
    }
  }
  
  # If no model information is found then return unknowns
  return(list(
    algorithm = "unknown", 
    dataset = "unknown",
    cell_type = "unknown", 
    gene_set = "unknown"
  ))
}


# Main function to generate performance summary
generate_model_performance_summary <- function() {
  # Output directory path should not be hardcoded but oopsie
  output_dir <- "c:/Users/mchul/Documents/StrataSign/output/2. models"
  
  # Find all prediction csv files
  prediction_files <- list.files(output_dir, 
                                pattern = "predictions.*\\.csv$", 
                                recursive = TRUE, 
                                full.names = TRUE)
  

  # Throw a hard error if no prediction files found
  if (length(prediction_files) == 0) {
    stop("No prediction files found in the output directory.")
  }
  
  # Init results dataframe
  performance_results <- data.frame()
  
  # Process each prediction file
  for (pred_file in prediction_files) {
    tryCatch({
      # Read predictions
      pred_data <- read.csv(pred_file, stringsAsFactors = FALSE)
      
      # Extract model information
      model_info <- extract_model_info(pred_file)
      
      # Calculate performance metrics
      metrics <- calculate_performance_metrics(pred_data$y_test, 
                                             pred_data$y_pred, 
                                             pred_data$y_pred_prob)
      
      # Add model information
      result_row <- data.frame(
        Model = model_info$algorithm,
        Dataset = model_info$dataset,
        Cell_Type = model_info$cell_type,
        Gene_Set = model_info$gene_set,
        N_Samples = nrow(pred_data),
        stringsAsFactors = FALSE
      )
      
      # Combine with metrics by adding columns
      result_row <- cbind(result_row, metrics)
      
      # Add to results by adding as rows
      performance_results <- rbind(performance_results, result_row)
      
    }, error = function(e) {
      warning(paste("Error processing", pred_file, ":", e$message))
    })
  }
  
  # Sort results by descending AUC and F1 Score
  performance_results <- performance_results %>%
    arrange(desc(AUC), desc(F1_Score))
  
  # Save results
  output_file <- file.path(output_dir, "model_performance.csv")
  write.csv(performance_results, output_file, row.names = FALSE)
  
  message(paste("\nResults saved to:", output_file))
  
  # Return results for further analysis
  return(performance_results)
}


# Run the performance summary generation
generate_model_performance_summary()
