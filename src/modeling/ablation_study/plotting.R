library(ggplot2)
library(data.table)
library(caret)
library(dplyr)
library(pROC)
library(tidyr)

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
    ggtitle(title)
  print(paste("Saving confusion matrix plot to:", output_path))
  ggsave(output_path, plot = p)
}

plot_shap_values <- function(shapley_path, title, output_path) {
  print(paste("Reading SHAP values from:", shapley_path))
  shapley <- readRDS(shapley_path)
  shapley_plot <- shapley$plot() +
    ggtitle(title)
  print(paste("Saving SHAP values plot to:", output_path))
  ggsave(output_path, plot = shapley_plot)
}

plot_top_genes_per_cluster <- function(original_feature_contributions_path, dataset_name, kernel_type, n = 10) {
  print(paste("Reading original feature contributions from:", original_feature_contributions_path))
  original_feature_contributions <- fread(original_feature_contributions_path)
  print("Finished reading original feature contributions.")
  
  # Extract cluster name from feature names
  original_feature_contributions <- original_feature_contributions %>%
    gather(key = "Feature", value = "Contribution") %>%
    mutate(ClusterName = sub(".*\\|(.*)\\|.*$", "\\1", Feature))
  
  # Sanitize cluster names for file paths
  original_feature_contributions$ClusterName <- gsub("/", "_", original_feature_contributions$ClusterName)
  
  # Group by cluster name and plot top genes
  clusters <- unique(original_feature_contributions$ClusterName)
  for (cluster in clusters) {
    print(paste("Processing cluster:", cluster))
    cluster_contributions <- original_feature_contributions %>%
      filter(ClusterName == cluster) %>%
      arrange(desc(abs(Contribution))) %>%
      head(n)
    
    # Separate positive and negative contributions
    positive_contributions <- cluster_contributions %>%
      filter(Contribution > 0) %>%
      arrange(desc(Contribution)) %>%
      head(n)
    
    negative_contributions <- cluster_contributions %>%
      filter(Contribution < 0) %>%
      arrange(Contribution) %>%
      head(n)
    
    # Combine positive and negative contributions
    combined_contributions <- bind_rows(positive_contributions, negative_contributions)
    
    p <- ggplot(combined_contributions, aes(x = reorder(Feature, Contribution), y = Contribution, fill = Contribution > 0)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("red", "blue"), labels = c("Negative", "Positive")) +
      ggtitle(paste("Top Genes for Cluster", cluster, "-", dataset_name, "-", kernel_type)) +
      labs(fill = "Contribution")
    
    output_path <- paste0("src/modeling/ablation_study/plotting_output/top_genes/", dataset_name, "_", kernel_type, "_cluster_", cluster, ".png")
    print(paste("Saving top genes plot to:", output_path))
    ggsave(output_path, plot = p)
  }
}

plot_overview_top_genes <- function(original_feature_contributions_path, dataset_name, kernel_type, n = 10) {
  print(paste("Reading original feature contributions from:", original_feature_contributions_path))
  original_feature_contributions <- fread(original_feature_contributions_path)
  print("Finished reading original feature contributions.")
  
  # Extract cluster name from feature names
  original_feature_contributions <- original_feature_contributions %>%
    gather(key = "Feature", value = "Contribution") %>%
    mutate(ClusterName = sub(".*\\|(.*)\\|.*$", "\\1", Feature))
  
  # Separate positive and negative contributions
  positive_contributions <- original_feature_contributions %>%
    filter(Contribution > 0) %>%
    arrange(desc(Contribution)) %>%
    head(n)
  
  negative_contributions <- original_feature_contributions %>%
    filter(Contribution < 0) %>%
    arrange(Contribution) %>%
    head(n)
  
  # Combine positive and negative contributions
  combined_contributions <- bind_rows(positive_contributions, negative_contributions)
  
  p <- ggplot(combined_contributions, aes(x = reorder(Feature, Contribution), y = Contribution, fill = Contribution > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("red", "blue"), labels = c("Negative", "Positive")) +
    ggtitle(paste("Top Contributing Genes Across All Clusters -", dataset_name, "-", kernel_type)) +
    labs(fill = "Contribution")
  
  output_path <- paste0("src/modeling/ablation_study/plotting_output/top_genes_overview_", dataset_name, "_", kernel_type, ".png")
  print(paste("Saving top genes overview plot to:", output_path))
  ggsave(output_path, plot = p)
}

plot_combined_roc_curve <- function(datasets, output_path) {
  print("Generating combined ROC curve...")
  roc_curves <- list()
  
  for (dataset_name in datasets) {
    print(paste("Reading predictions for dataset:", dataset_name))
    predictions_linear <- fread(paste0("src/modeling/ablation_study/model_output/predictions_linear_", dataset_name, ".csv"))
    predictions_rbf <- fread(paste0("src/modeling/ablation_study/model_output/predictions_rbf_", dataset_name, ".csv"))
    
    roc_curve_linear <- roc(predictions_linear$y_test, predictions_linear$y_pred)
    roc_curve_rbf <- roc(predictions_rbf$y_test, predictions_rbf$y_pred)
    
    roc_curves[[paste(dataset_name, "Linear")]] <- roc_curve_linear
    roc_curves[[paste(dataset_name, "RBF")]] <- roc_curve_rbf
  }
  
  p <- ggroc(roc_curves) +
    ggtitle("Combined ROC Curves for All Datasets")
  
  print(paste("Saving combined ROC curve plot to:", output_path))
  ggsave(output_path, plot = p)
}

# Ensure the output directory exists
output_dir <- "src/modeling/ablation_study/plotting_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Example usage for plotting (you can call these functions as needed)
datasets <- c("all_clusters", "lcam_hi", "lcam_lo", "lcam_both")

for (dataset_name in datasets) {
  print(paste("Processing dataset:", dataset_name))
  
  # Confusion matrices
  print("Plotting confusion matrix for linear kernel...")
  plot_confusion_matrix(
    paste0("src/modeling/ablation_study/model_output/predictions_linear_", dataset_name, ".csv"),
    paste("Confusion Matrix for Linear Kernel -", dataset_name),
    paste0(output_dir, "/confusion_matrix_linear_", dataset_name, ".png")
  )
  
  print("Plotting confusion matrix for RBF kernel...")
  plot_confusion_matrix(
    paste0("src/modeling/ablation_study/model_output/predictions_rbf_", dataset_name, ".csv"),
    paste("Confusion Matrix for RBF Kernel -", dataset_name),
    paste0(output_dir, "/confusion_matrix_rbf_", dataset_name, ".png")
  )
  
  # SHAP values
  print("Plotting SHAP values for linear kernel...")
  plot_shap_values(
    paste0("src/modeling/ablation_study/model_output/shapley_linear_", dataset_name, ".rds"),
    paste("SHAP Values for Linear Kernel -", dataset_name),
    paste0(output_dir, "/shap_linear_", dataset_name, ".png")
  )
  
  print("Plotting SHAP values for RBF kernel...")
  plot_shap_values(
    paste0("src/modeling/ablation_study/model_output/shapley_rbf_", dataset_name, ".rds"),
    paste("SHAP Values for RBF Kernel -", dataset_name),
    paste0(output_dir, "/shap_rbf_", dataset_name, ".png")
  )
  
  # Top genes per cluster
  print("Plotting top genes per cluster for linear kernel...")
  plot_top_genes_per_cluster(
    paste0("src/modeling/ablation_study/model_output/original_feature_contributions_linear_", dataset_name, ".csv"),
    dataset_name,
    "linear",
    n = 10
  )
  
  print("Plotting top genes per cluster for RBF kernel...")
  plot_top_genes_per_cluster(
    paste0("src/modeling/ablation_study/model_output/original_feature_contributions_rbf_", dataset_name, ".csv"),
    dataset_name,
    "rbf",
    n = 10
  )
  
  # Overview of top contributing genes across all clusters
  print("Plotting overview of top contributing genes across all clusters for linear kernel...")
  plot_overview_top_genes(
    paste0("src/modeling/ablation_study/model_output/original_feature_contributions_linear_", dataset_name, ".csv"),
    dataset_name,
    "linear",
    n = 10
  )
  
  print("Plotting overview of top contributing genes across all clusters for RBF kernel...")
  plot_overview_top_genes(
    paste0("src/modeling/ablation_study/model_output/original_feature_contributions_rbf_", dataset_name, ".csv"),
    dataset_name,
    "rbf",
    n = 10
  )
}

# Combined ROC curves for all datasets
plot_combined_roc_curve(datasets, paste0(output_dir, "/combined_roc_curves.png"))

print("Script completed successfully.")
