library(ggplot2)
library(data.table)
library(caret)
library(dplyr)
library(pROC)
library(tidyr)

# Define base path
base_path <- "src/modeling/ablation_study"

# Define datasets structure
datasets <- list(
  all_clusters = list(),
  lcam_hi = list(),
  lcam_lo = list(),
  lcam_both = list()
)

versions <- c("default", "random1", "random2", "random3")

plot_confusion_matrix <- function(predictions_path, title, output_path) {
  print(paste("Reading predictions from:", predictions_path))
  predictions <- fread(predictions_path)
  predictions$y_pred <- factor(predictions$y_pred, levels = c(0, 1))
  predictions$y_test <- factor(predictions$y_test, levels = c(0, 1))
  print("Generating confusion matrix...")
  confusion <- confusionMatrix(predictions$y_pred, predictions$y_test)
  confusion_matrix_plot <- as.data.frame(confusion$table)
  p <- ggplot(confusion_matrix_plot, aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(title) +
    theme_minimal()
  print(paste("Saving confusion matrix plot to:", output_path))
  ggsave(output_path, plot = p)
}

plot_feature_importance <- function(feature_importance_path, dataset_name, version, max_features = 50) {
  print(paste("Reading feature importance from:", feature_importance_path))
  feature_importance <- fread(feature_importance_path)
  
  # Select only non-zero features, sorted by absolute importance
  important_features <- feature_importance %>%
    filter(Coefficient != 0) %>%           # Remove zero coefficients
    arrange(desc(Abs_Coefficient)) %>%     # Sort by absolute value
    head(max_features)                     # Take top features (up to 50)
  
  # Create plot if there are any important features
  if (nrow(important_features) > 0) {
    num_features <- nrow(important_features)
    p <- ggplot(important_features, 
                aes(x = reorder(Feature, Abs_Coefficient), 
                    y = Coefficient, 
                    fill = Direction)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("Negative" = "red", "Positive" = "blue")) +
      labs(
        title = paste0("Most Important Features (n=", num_features, ")\n",
                      dataset_name, " - ", version),
        x = "Feature",
        y = "LASSO Coefficient"
      ) +
      theme_minimal()
    
    # Save plot
    output_dir <- paste0("src/modeling/ablation_study/plotting_output_lasso/", 
                        dataset_name, "/", version)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_path <- paste0(output_dir, "/feature_importance_", 
                         dataset_name, "_", version, ".png")
    print(paste("Saving plot with", num_features, "features to:", output_path))
    ggsave(output_path, plot = p, width = 10, height = min(12, max(6, num_features * 0.25)))
  } else {
    warning(paste("No non-zero coefficients found for", dataset_name, version))
  }
}

# Ensure output directory exists
output_dir <- "src/modeling/ablation_study/plotting_output_lasso"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate plots for each dataset and version
for (dataset_name in names(datasets)) {
  for (version in versions) {
    print(paste("Processing dataset:", dataset_name, "version:", version))
    
    # Plot confusion matrix
    plot_confusion_matrix(
      paste0("src/modeling/ablation_study/lasso_output/", dataset_name, "/", version, 
             "/predictions_", dataset_name, "_", version, ".csv"),
      paste("Confusion Matrix -", dataset_name, "-", version),
      paste0(output_dir, "/", dataset_name, "/", version, 
             "/confusion_matrix_", dataset_name, "_", version, ".png")
    )
    
    # Plot feature importance
    plot_feature_importance(
      paste0("src/modeling/ablation_study/lasso_output/", dataset_name, "/", version, 
             "/feature_importance_", dataset_name, "_", version, ".csv"),
      dataset_name,
      version,
      max_features = 50
    )
  }
}

print("LASSO plotting completed successfully.")