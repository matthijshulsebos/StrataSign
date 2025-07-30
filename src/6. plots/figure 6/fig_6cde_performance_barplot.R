library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(viridis)
library(stringr)
library(patchwork)

# Load and clean performance data
load_model_performance <- function(file_path = "output/2. models/model_performance.csv") {
  perf_data <- read_csv(file_path, show_col_types = FALSE)
  
  # Clean and filter the data
  metric_cols <- intersect(c("AUC", "Accuracy", "F1_Score", "Precision", "Recall", "Balanced_Accuracy", "Specificity"), names(perf_data))
  clean_data <- perf_data %>%
    mutate(
      # Clean model names
      Model = case_when(
        str_to_lower(Model) == "lasso" ~ "LASSO",
        str_to_lower(Model) == "lightgbm" ~ "LightGBM",
        str_to_lower(Model) == "xgboost" ~ "XGBoost",
        str_to_lower(Model) == "randomforest" ~ "Random Forest",
        str_to_lower(Model) == "elasticnet" ~ "Elastic Net",
        str_to_lower(Model) == "spls" ~ "sPLS",
        str_to_lower(Model) == "svm_rbf" ~ "SVM (RBF)",
        str_to_lower(Model) == "svm_linear" ~ "SVM (Linear)",
        TRUE ~ str_to_title(Model)
      )
    )
  if (length(metric_cols) > 0) {
    clean_data <- clean_data %>% mutate(across(all_of(metric_cols), as.numeric))
  }
  
  return(clean_data)
}


# Save barplots for all combinations
save_all_performance_barplots <- function(perf_data, 
                                          metrics_to_plot = c("AUC", "Accuracy", "F1_Score", "Precision", "Recall", "Balanced_Accuracy", "Specificity", "MCC")) {
  # Use raw values for iteration
  raw_datasets <- unique(perf_data$Dataset)
  raw_cell_types <- unique(perf_data$Cell_Type)
  raw_gene_sets <- unique(perf_data$Gene_Set)

  # Output root for figure 6
  performance_root <- file.path("output/6. plots/figure 6/performance")
  dir.create(performance_root, recursive = TRUE, showWarnings = FALSE)

  for (dataset in raw_datasets) {
    for (cell_type in raw_cell_types) {
      for (gene_set in raw_gene_sets) {
        filtered_data <- perf_data %>%
          filter(Dataset == dataset, Cell_Type == cell_type, Gene_Set == gene_set)
        if (nrow(filtered_data) == 0) next

        # Order models alphabetically
        model_order <- sort(unique(filtered_data$Model), decreasing = TRUE)
        filtered_data <- filtered_data %>% mutate(Model = factor(Model, levels = model_order))

        # Set the metrics to plot
        existing_metrics <- metrics_to_plot[metrics_to_plot %in% names(filtered_data)]
        if (length(existing_metrics) == 0) next

        long_data <- filtered_data %>%
          dplyr::select(Model, all_of(existing_metrics)) %>%
          pivot_longer(cols = all_of(existing_metrics), names_to = "Metric", values_to = "Value") %>%
          mutate(Metric = factor(Metric, levels = existing_metrics))

        for (metric in existing_metrics) {
          metric_data <- long_data %>% filter(Metric == metric)
          p_metric <- ggplot(metric_data, aes(x = Model, y = Value, fill = Model)) +
            geom_col(alpha = 0.85, position = position_dodge(width = 0.8), width = 0.7) +
            geom_text(aes(label = round(Value, 3)), 
                      position = position_dodge(width = 0.8), 
                      hjust = -0.1, size = 4.5, family = "sans") +
            scale_fill_brewer(palette = "Set1") +
            scale_y_continuous(expand = expansion(mult = c(0.01, 0.15))) +
            coord_flip() +
            theme_minimal(base_size = 18, base_family = "sans") +
            theme(
              axis.text.x = element_text(size = 15, color = "black"),
              axis.text.y = element_text(size = 15, color = "black"),
              axis.title.x = element_text(size = 17, color = "black", face = "plain"),
              axis.title.y = element_text(size = 17, color = "black", face = "plain"),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(color = "grey80"),
              plot.title = element_blank(),
              plot.margin = margin(10, 24, 10, 10),
              strip.background = element_blank(),
              strip.text = element_text(face = "plain", size = 15, color = "black"),
              axis.line = element_line(color = "black", linewidth = 0.8)
            ) +
            labs(
              x = NULL,
              y = metric,
              title = sprintf("%s\n%s | %s | %s", metric, dataset, cell_type, gene_set)
            )

          # Directory for this combination
          safe_dataset <- gsub("[^a-zA-Z0-9_]+", "_", dataset)
          safe_celltype <- gsub("[^a-zA-Z0-9_]+", "_", cell_type)
          safe_geneset <- gsub("[^a-zA-Z0-9_]+", "_", gene_set)
          combo_dir <- file.path(performance_root, safe_dataset, safe_celltype, safe_geneset)
          dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
          metric_file <- file.path(combo_dir, sprintf("barplot_%s.png", metric))
          ggsave(metric_file, p_metric, width = 8, height = 5, dpi = 300, bg = "transparent")
        }
      }
    }
  }
  message("Completed writing figure 6CDE to file.")
}


# Main function
generate_performance_plots <- function(performance_file = "output/2. models/model_performance.csv") {
  # Use relative paths
  if (!file.exists(performance_file) && !startsWith(performance_file, "/") && !grepl("^[A-Za-z]:", performance_file)) {
    performance_file_path <- file.path(performance_file)
  } else {
    performance_file_path <- performance_file
  }

  perf_data <- load_model_performance(performance_file_path)

  # Save all barplots for all combinations (all datasets)
  save_all_performance_barplots(perf_data)

  return(list(data = perf_data))
}


# Run plotting function
performance_results <- generate_performance_plots()
