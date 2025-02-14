library(ggplot2)
library(data.table)
library(caret)
library(dplyr)
library(pROC)
library(tidyr)
library(sva)

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

plot_shap_values <- function(shapley_path, title, output_path) {
  print(paste("Reading SHAP values from:", shapley_path))
  shapley <- readRDS(shapley_path)
  shapley_plot <- shapley$plot() +
    ggtitle(title) +
    theme_minimal()
  print(paste("Saving SHAP values plot to:", output_path))
  ggsave(output_path, plot = shapley_plot)
}

plot_top_genes_per_cluster <- function(original_feature_contributions_path, dataset_name, version, kernel_type, n = 50) {
  print(paste("Reading original feature contributions from:", original_feature_contributions_path))
  original_feature_contributions <- fread(original_feature_contributions_path)
  print("Finished reading original feature contributions.")
  
  # Extract cluster name from feature names
  original_feature_contributions <- original_feature_contributions %>%
    gather(key = "Feature", value = "Contribution") %>%
    mutate(ClusterName = sub(".*\\|(.*)\\|.*$", "\\1", Feature))
  
  # Sanitize cluster names for file paths
  original_feature_contributions$ClusterName <- gsub("/", "_", original_feature_contributions$ClusterName)
  
  # Create output directory for top genes with kernel-specific subdirectory
  output_dir <- paste0("src/modeling/ablation_study/plotting_output/", dataset_name, "/", version, "/top_genes/", kernel_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Group by cluster name and plot top genes
  clusters <- unique(original_feature_contributions$ClusterName)
  for (cluster in clusters) {
    print(paste("Processing cluster:", cluster))
    cluster_contributions <- original_feature_contributions %>%
      filter(ClusterName == cluster) %>%
      arrange(desc(abs(Contribution))) %>%
      head(n)
    
    p <- ggplot(cluster_contributions, aes(x = reorder(Feature, Contribution), y = Contribution, fill = Contribution > 0)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("red", "blue"), labels = c("Negative", "Positive")) +
      ggtitle(paste("Top Genes for Cluster", cluster, "-", dataset_name, "-", version, "-", kernel_type)) +
      labs(fill = "Contribution") +
      theme_minimal()
    
    ggsave(paste0(output_dir, "/top_genes_cluster_", cluster, "_", kernel_type, ".png"), plot = p)
  }
}

plot_overview_top_genes <- function(original_feature_contributions_path, dataset_name, version, kernel_type, n = 50) {
  print(paste("Reading original feature contributions from:", original_feature_contributions_path))
  original_feature_contributions <- fread(original_feature_contributions_path)
  print("Finished reading original feature contributions.")
  
  # Create output directory for top genes overview with kernel-specific subdirectory
  output_dir <- paste0("src/modeling/ablation_study/plotting_output/", dataset_name, "/", version, "/top_genes/", kernel_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract cluster name from feature names
  original_feature_contributions <- original_feature_contributions %>%
    gather(key = "Feature", value = "Contribution") %>%
    mutate(ClusterName = sub(".*\\|(.*)\\|.*$", "\\1", Feature))
  
  # Select top n most impactful genes by absolute contribution
  top_contributions <- original_feature_contributions %>%
    arrange(desc(abs(Contribution))) %>% 
    head(n)
  
  p <- ggplot(top_contributions, aes(x = reorder(Feature, Contribution), y = Contribution, fill = Contribution > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("red", "blue"), labels = c("Negative", "Positive")) +
    ggtitle(paste("Top Contributing Genes Across All Clusters -", dataset_name, "-", version, "-", kernel_type)) +
    labs(fill = "Contribution") +
    theme_minimal()
  
  ggsave(paste0(output_dir, "/top_genes_overview_", kernel_type, ".png"), plot = p)
}

plot_feature_contributions_distribution <- function(original_feature_contributions_path, dataset_name, version, kernel_type) {
  print(paste("Reading original feature contributions from:", original_feature_contributions_path))
  original_feature_contributions <- fread(original_feature_contributions_path)
  print("Finished reading original feature contributions.")
  
  # Create output directory for feature contributions
  output_dir <- paste0("src/modeling/ablation_study/plotting_output/", dataset_name, "/", version, "/feature_contributions")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Gather the data for plotting
  original_feature_contributions <- original_feature_contributions %>%
    gather(key = "PrincipalComponent", value = "Contribution")
  
  p <- ggplot(original_feature_contributions, aes(x = Contribution)) +
    geom_density(fill = "blue", alpha = 0.5) +
    ggtitle(paste("Distribution of Original Feature Contributions -", dataset_name, "-", version, "-", kernel_type)) +
    xlab("Contribution") +
    ylab("Density") +
    theme_minimal()
  
  ggsave(paste0(output_dir, "/feature_contributions_distribution_", kernel_type, ".png"), plot = p)
}

plot_combined_roc_curve <- function(datasets, output_path) {
  print("Generating combined ROC curve...")
  roc_curves <- list()
  
  for (dataset_name in datasets) {
    print(paste("Reading predictions for dataset:", dataset_name))
    predictions_linear <- fread(paste0("src/modeling/ablation_study/model_output/", dataset_name, "/default/predictions_linear_", dataset_name, "_default.csv"))
    predictions_rbf <- fread(paste0("src/modeling/ablation_study/model_output/", dataset_name, "/default/predictions_rbf_", dataset_name, "_default.csv"))
    
    roc_curve_linear <- roc(predictions_linear$y_test, predictions_linear$y_pred)
    roc_curve_rbf <- roc(predictions_rbf$y_test, predictions_rbf$y_pred)
    
    roc_curves[[paste(dataset_name, "Linear")]] <- roc_curve_linear
    roc_curves[[paste(dataset_name, "RBF")]] <- roc_curve_rbf
  }
  
  p <- ggroc(roc_curves) +
    ggtitle("Combined ROC Curves for Default Models") +
    theme_minimal()
  
  print(paste("Saving combined ROC curve plot to:", output_path))
  ggsave(output_path, plot = p)
}

# Function to plot PCA diagnostics
plot_pca_diagnostics <- function(pca_scores, output_dir, plot_type = "boxplot") {
  # Define metadata columns
  metadata_columns <- c("amp_batch_ID", "tissue", "disease", "library_chemistry", "prep")
  
  # Ensure amp_batch_ID is treated as a factor
  pca_scores$amp_batch_ID <- as.factor(pca_scores$amp_batch_ID)
  
  # Create PCA diagnostics directory
  pca_dir <- paste0(output_dir, "/pca_diagnostics")
  dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Plot PC1 vs metadata columns
  for (col in metadata_columns) {
    if (col == "amp_batch_ID") {
      p <- ggplot(pca_scores, aes_string(x = col, y = "PC1", color = "tissue", shape = "tissue")) +
        geom_point() +
        ggtitle(paste("PC1 vs", col)) +
        xlab(col) +
        ylab("PC1") +
        theme_minimal()
    } else {
      if (plot_type == "violin") {
        p <- ggplot(pca_scores, aes_string(x = col, y = "PC1")) +
          geom_violin() +
          ggtitle(paste("PC1 vs", col)) +
          xlab(col) +
          ylab("PC1") +
          theme_minimal()
      } else {
        p <- ggplot(pca_scores, aes_string(x = col, y = "PC1")) +
          geom_boxplot() +
          ggtitle(paste("PC1 vs", col)) +
          xlab(col) +
          ylab("PC1") +
          theme_minimal()
      }
    }
    ggsave(paste0(pca_dir, "/pc1_vs_", col, ".png"), plot = p)
  }
  
  # Scatterplot of PC1 vs. PC2, colored by Use.in.Clustering.Model.
  p2 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Use.in.Clustering.Model.)) +
    geom_point() +
    ggtitle("PC1 vs. PC2 Colored by Use in Clustering Model") +
    xlab("PC1") +
    ylab("PC2") +
    theme_minimal()
  ggsave(paste0(pca_dir, "/pc1_vs_pc2_by_use_in_clustering_model.png"), plot = p2)
}

# Example usage for plotting
datasets <- list(
  all_clusters = c("default", "random1", "random2", "random3"),
  lcam_hi = c("default", "random1", "random2", "random3"),
  lcam_lo = c("default", "random1", "random2", "random3"),
  lcam_both = c("default", "random1", "random2", "random3")
)

# Ensure the output directory exists
output_dir <- "src/modeling/ablation_study/plotting_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (dataset_name in names(datasets)) {
  for (version in datasets[[dataset_name]]) {
    print(paste("Processing dataset:", dataset_name, "version:", version))
    
    # Create output directory for this dataset and version
    dir.create(paste0(output_dir, "/", dataset_name, "/", version), recursive = TRUE, showWarnings = FALSE)
    
    # Confusion matrices
    print("Plotting confusion matrix for linear kernel...")
    plot_confusion_matrix(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/predictions_linear_", dataset_name, "_", version, ".csv"),
      paste("Confusion Matrix for Linear Kernel -", dataset_name, "-", version),
      paste0(output_dir, "/", dataset_name, "/", version, "/confusion_matrix_linear_", dataset_name, "_", version, ".png")
    )
    
    print("Plotting confusion matrix for RBF kernel...")
    plot_confusion_matrix(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/predictions_rbf_", dataset_name, "_", version, ".csv"),
      paste("Confusion Matrix for RBF Kernel -", dataset_name, "-", version),
      paste0(output_dir, "/", dataset_name, "/", version, "/confusion_matrix_rbf_", dataset_name, "_", version, ".png")
    )
    
    # SHAP values
    print("Plotting SHAP values for linear kernel...")
    plot_shap_values(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/shapley_linear_", dataset_name, "_", version, ".rds"),
      paste("SHAP Values for Linear Kernel -", dataset_name, "-", version),
      paste0(output_dir, "/", dataset_name, "/", version, "/shap_linear_", dataset_name, "_", version, ".png")
    )
    
    print("Plotting SHAP values for RBF kernel...")
    plot_shap_values(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/shapley_rbf_", dataset_name, "_", version, ".rds"),
      paste("SHAP Values for RBF Kernel -", dataset_name, "-", version),
      paste0(output_dir, "/", dataset_name, "/", version, "/shap_rbf_", dataset_name, "_", version, ".png")
    )
    
    # Top genes per cluster
    print("Plotting top genes per cluster for linear kernel...")
    plot_top_genes_per_cluster(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/original_feature_contributions_linear_", dataset_name, "_", version, ".csv"),
      dataset_name,  # Pass dataset_name and version separately
      version,
      "linear",
      n = 50
    )
    
    print("Plotting top genes per cluster for RBF kernel...")
    plot_top_genes_per_cluster(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/original_feature_contributions_rbf_", dataset_name, "_", version, ".csv"),
      dataset_name,  # Pass dataset_name and version separately
      version,
      "rbf",
      n = 50
    )
    
    # Overview of top contributing genes across all clusters
    print("Plotting overview of top contributing genes...")
    plot_overview_top_genes(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/original_feature_contributions_linear_", dataset_name, "_", version, ".csv"),
      dataset_name,  # Pass dataset_name and version separately
      version,
      "linear",
      n = 50
    )
    
    print("Plotting overview of top contributing genes for RBF kernel...")
    plot_overview_top_genes(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/original_feature_contributions_rbf_", dataset_name, "_", version, ".csv"),
      dataset_name,  # Pass dataset_name and version separately
      version,
      "rbf",
      n = 50
    )
    
    # Distribution of original feature contributions
    print("Plotting distribution of original feature contributions...")
    plot_feature_contributions_distribution(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/original_feature_contributions_linear_", dataset_name, "_", version, ".csv"),
      dataset_name,  # Pass dataset_name and version separately
      version,
      "linear"
    )
    
    print("Plotting distribution of original feature contributions for RBF kernel...")
    plot_feature_contributions_distribution(
      paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/original_feature_contributions_rbf_", dataset_name, "_", version, ".csv"),
      dataset_name,  # Pass dataset_name and version separately
      version,
      "rbf"
    )
    
    # PCA diagnostics
    print("Plotting PCA diagnostics...")
    pca_scores <- fread(paste0("src/modeling/ablation_study/model_output/", dataset_name, "/", version, "/pca_scores_", dataset_name, "_", version, ".csv"))
    plot_pca_diagnostics(pca_scores, paste0(output_dir, "/", dataset_name, "/", version), plot_type = "boxplot")
  }
}

# Combined ROC curves for all datasets and versions
plot_combined_roc_curve(names(datasets), paste0(output_dir, "/combined_roc_curves.png"))

print("Script completed successfully.")
