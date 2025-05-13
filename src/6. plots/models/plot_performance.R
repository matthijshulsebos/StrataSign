library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(viridis)
library(stringr)
library(patchwork)

# Helper function to safely save plots without fs package
safe_save_plot <- function(plot, filename, width = 8, height = 6) {
  full_path <- normalizePath(filename, mustWork = FALSE)
  
  # Use base R directory creation instead of fs::dir_create
  dir.create(dirname(full_path), recursive = TRUE, showWarnings = FALSE)
  
  message("Saving plot to: ", full_path)
  
  # Save the plot
  png(full_path, width = width * 300, height = height * 300, res = 300)
  print(plot)
  dev.off()
}

compare_model_performances <- function() {
  # Use absolute paths for everything
  base_dir <- normalizePath(".", winslash = "/")
  output_dir <- file.path(base_dir, "output/figures/performance")
  
  # Use base R for directory creation
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Find all performance summary files (using list.files instead of dir_ls)
  models_dir <- file.path(base_dir, "output/models")
  perf_pattern <- "*_performance_summary.csv"
  
  # Use base R to find files
  all_files <- list.files(models_dir, 
                         pattern = perf_pattern, 
                         recursive = TRUE, 
                         full.names = TRUE)
  
  message(paste("Found", length(all_files), "performance summary files"))
  
  # Load and combine all performance data
  all_performance <- data.frame()
  for (file in all_files) {
    # Extract model name from path
    model_dir <- dirname(file)
    model_name <- basename(model_dir)
    
    # Read data
    perf_data <- read_csv(file, show_col_types = FALSE)
    perf_data$Model <- model_name
    
    # Append to combined dataset
    all_performance <- bind_rows(all_performance, perf_data)
  }
  
  # Process data
  all_performance <- all_performance %>%
    mutate(
      Gene_Set = case_when(
        endsWith(Dataset, "_metabolic") ~ "metabolic",
        endsWith(Dataset, "_nonmetabolic") ~ "nonmetabolic",
        endsWith(Dataset, "_random") ~ "random",
        TRUE ~ "unknown"
      ),
      Cluster = case_when(
        startsWith(Dataset, "lcam_hi_") ~ "lcam_hi",
        startsWith(Dataset, "lcam_lo_") ~ "lcam_lo",
        startsWith(Dataset, "lcam_both_") ~ "lcam_both",
        startsWith(Dataset, "all_clusters_") ~ "all_clusters",
        TRUE ~ "other"
      )
    )
  
  # 1. RESTORED: Create bar chart comparison of key metrics
  key_metrics_plot <- all_performance %>%
    pivot_longer(cols = c("Accuracy", "F1_score", "ROC_AUC"), 
                 names_to = "Metric", values_to = "Value") %>%
    ggplot(aes(x = Model, y = Value, fill = Model)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(Metric ~ Dataset) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    labs(title = "Model Performance Comparison", 
         subtitle = "Showing Accuracy, F1 Score, and ROC AUC metrics",
         y = "Score Value (higher is better)", 
         x = "Model") +
    scale_fill_viridis_d()
  
  # Save key metrics plot
  safe_save_plot(
    key_metrics_plot, 
    file.path(output_dir, "all_metrics_by_model.png"),
    width = 12,
    height = 10
  )
  
  # 2. RESTORED: Precision vs Recall plot
  precision_recall_plot <- all_performance %>%
    filter(!is.na(Precision) & !is.na(Recall)) %>%
    ggplot(aes(x = Recall, y = Precision, color = Model)) +
    geom_point(size = 3) + 
    facet_wrap(~ Gene_Set) +
    theme_minimal() +
    labs(title = "Precision vs Recall by Model", 
         x = "Recall", 
         y = "Precision") +
    scale_color_viridis_d()
  
  # Save precision-recall plot
  safe_save_plot(
    precision_recall_plot, 
    file.path(output_dir, "precision_vs_recall.png"),
    width = 10,
    height = 8
  )
  
  # 3. Metabolic gene accuracy heatmap
  metabolic_accuracy_heatmap <- all_performance %>%
    filter(Gene_Set == "metabolic") %>%
    select(Dataset, Model, Accuracy, Cluster) %>%
    ggplot(aes(x = Model, y = Cluster, fill = Accuracy)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Accuracy, 2)), color = "black", size = 3.5) +
    scale_fill_gradient2(low = "white", mid = "#56B4E9", high = "#0072B2", 
                        midpoint = 0.85, limits = c(0.7, 1), name = "Accuracy") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = "Model Accuracy for Metabolic Genes",
      x = NULL,
      y = NULL
    )
  
  # Save metabolic heatmap
  safe_save_plot(
    metabolic_accuracy_heatmap, 
    file.path(output_dir, "metabolic_accuracy_heatmap.png"),
    width = 8,
    height = 6
  )
  
  # 4. Nonmetabolic gene accuracy heatmap
  nonmetabolic_accuracy_heatmap <- all_performance %>%
    filter(Gene_Set == "nonmetabolic") %>%
    select(Dataset, Model, Accuracy, Cluster) %>%
    ggplot(aes(x = Model, y = Cluster, fill = Accuracy)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Accuracy, 2)), color = "black", size = 3.5) +
    scale_fill_gradient2(low = "white", mid = "#56B4E9", high = "#0072B2", 
                        midpoint = 0.85, limits = c(0.7, 1), name = "Accuracy") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = "Model Accuracy for Non-metabolic Genes",
      x = NULL,
      y = NULL
    )
  
  # Save nonmetabolic heatmap
  safe_save_plot(
    nonmetabolic_accuracy_heatmap, 
    file.path(output_dir, "nonmetabolic_accuracy_heatmap.png"),
    width = 8,
    height = 6
  )
  
  # 5. Random gene accuracy heatmap
  random_accuracy_heatmap <- all_performance %>%
    filter(Gene_Set == "random") %>%
    select(Dataset, Model, Accuracy, Cluster) %>%
    ggplot(aes(x = Model, y = Cluster, fill = Accuracy)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Accuracy, 2)), color = "black", size = 3.5) +
    scale_fill_gradient2(low = "white", mid = "#56B4E9", high = "#0072B2", 
                        midpoint = 0.85, limits = c(0.7, 1), name = "Accuracy") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = "Model Accuracy for Random Genes",
      x = NULL,
      y = NULL
    )
  
  # Save random heatmap
  safe_save_plot(
    random_accuracy_heatmap, 
    file.path(output_dir, "random_accuracy_heatmap.png"),
    width = 8,
    height = 6
  )
  
  # 6. RESTORED: Cluster-specific bar charts comparing gene sets
  clusters <- unique(all_performance$Cluster)
  for(cluster in clusters) {
    # Skip invalid clusters
    if(is.na(cluster) || cluster == "other") next
    
    # Filter for current cluster
    cluster_data <- all_performance %>% 
      filter(Cluster == cluster)
    
    # Create bar chart for this cluster
    cluster_plot <- cluster_data %>%
      pivot_longer(cols = c("Accuracy", "F1_score", "ROC_AUC"), 
                   names_to = "Metric", values_to = "Value") %>%
      ggplot(aes(x = Model, y = Value, fill = Gene_Set)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ Metric, ncol = 1) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank(),
        legend.position = "top"
      ) +
      labs(
        title = paste(cluster, "Performance"),
        y = "Score", 
        x = NULL
      ) +
      scale_fill_viridis_d()
    
    # Save cluster plot
    safe_save_plot(
      cluster_plot,
      file.path(output_dir, paste0("cluster_performance_", cluster, ".png")),
      width = 10,
      height = 8
    )
  }
  
  # 7. RESTORED: Gene set-specific bar charts (one file per gene set)
  gene_sets <- unique(all_performance$Gene_Set)
  for(gene_set in gene_sets) {
    # Skip invalid gene sets
    if(is.na(gene_set) || gene_set == "unknown") next
    
    # Filter for current gene set
    gene_set_data <- all_performance %>% 
      filter(Gene_Set == gene_set)
    
    # Create simplified bar chart for this gene set
    gene_set_plot <- gene_set_data %>%
      pivot_longer(cols = c("Accuracy", "F1_score", "ROC_AUC"), 
                  names_to = "Metric", values_to = "Value") %>%
      ggplot(aes(x = Model, y = Value, fill = Cluster)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ Metric, ncol = 1) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank(),
        legend.position = "top"
      ) +
      labs(
        title = paste(gene_set, "Genes"),
        y = "Score", 
        x = NULL
      ) +
      scale_fill_viridis_d()
    
    # Save gene set plot
    safe_save_plot(
      gene_set_plot,
      file.path(output_dir, paste0("geneset_performance_", gene_set, ".png")),
      width = 10,
      height = 8
    )
  }
  
  # 8. Create best models summary
  best_models <- all_performance %>%
    group_by(Dataset) %>%
    summarize(
      Best_Accuracy_Model = Model[which.max(Accuracy)],
      Best_Accuracy = max(Accuracy),
      Best_F1_Model = Model[which.max(F1_score)],
      Best_F1 = max(F1_score),
      Best_AUC_Model = Model[which.max(ROC_AUC)],
      Best_AUC = max(ROC_AUC)
    )
  
  # Save summary
  write_csv(best_models, file.path(output_dir, "best_models_summary.csv"))
  
  message("Performance visualization complete. Results saved to ", output_dir)
}

# Run the analysis
compare_model_performances()