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
      ),
      # Clean normalization names
      Dataset = case_when(
        Dataset == "ctnorm_global" ~ "CT-norm (global)",
        Dataset == "ctnorm_relative" ~ "CT-norm (relative)",
        Dataset == "read_depth" ~ "Read depth",
        TRUE ~ Dataset
      ),
      # Clean cell type names
      Cell_Type = case_when(
        Cell_Type == "all_clusters" ~ "All cell types",
        Cell_Type == "lcam_both" ~ "LCAM both",
        Cell_Type == "lcam_lo" ~ "LCAM low",
        Cell_Type == "lcam_hi" ~ "LCAM high",
        Cell_Type == "macrophages" ~ "Macrophages",
        TRUE ~ str_to_title(Cell_Type)
      ),
      # Clean gene type names
      Gene_Set = case_when(
        Gene_Set == "metabolic" ~ "Metabolic",
        Gene_Set == "nonmetabolic" ~ "Non-metabolic",
        Gene_Set == "random" ~ "Random",
        TRUE ~ str_to_title(Gene_Set)
      )
    )
  
  return(clean_data)
}


# Save barplots for all combinations
save_all_performance_barplots <- function(perf_data, 
                                          datasets = NULL, 
                                          cell_types = NULL, 
                                          gene_sets = NULL, 
                                          metrics_to_plot = c("AUC", "Accuracy", "F1_Score", "Precision", "Recall", "Balanced_Accuracy", "Specificity", "MCC")) {

  datasets <- unique(perf_data$Dataset)
  cell_types <- unique(perf_data$Cell_Type)
  gene_sets <- unique(perf_data$Gene_Set)

  # Output path
  s1_metrics_root <- file.path("output/6. plots/figure s1/metrics")
  dir.create(s1_metrics_root, recursive = TRUE, showWarnings = FALSE)

  # Output for figure 4 (accuracy barplot)
  fig4_dir <- file.path("output/6. plots/figure 4")
  dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)

  for (dataset in datasets) {
    for (cell_type in cell_types) {
      for (gene_set in gene_sets) {
        filtered_data <- perf_data %>%
          filter(Dataset == dataset, Cell_Type == cell_type, Gene_Set == gene_set)
        if (nrow(filtered_data) == 0) next

        # Order models alphabetically (descending)
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
              axis.title.y = element_blank(),
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
              x = "Model",
              y = NULL,
              title = sprintf("%s\n%s | %s | %s", metric, dataset, cell_type, gene_set)
            )

          # Directory for this metric
          metric_dir <- file.path(s1_metrics_root, metric)
          dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
          # File name of combination but remove special characters
          safe_dataset <- gsub("[^a-zA-Z0-9_]", "_", dataset)
          safe_celltype <- gsub("[^a-zA-Z0-9_]", "_", cell_type)
          safe_geneset <- gsub("[^a-zA-Z0-9_]", "_", gene_set)
          metric_file <- file.path(metric_dir, sprintf("barplot_%s_%s_%s.png", safe_dataset, safe_celltype, safe_geneset))
          ggsave(metric_file, p_metric, width = 8, height = 5, dpi = 300)

          # If this is the main barplot also save to figure 4
          if (tolower(metric) == "accuracy" &&
              dataset == "CT-norm (global)" &&
              cell_type == "All cell types" &&
              gene_set == "Metabolic") {
            fig4_file <- file.path(fig4_dir, "fig_4b_accuracy_barplot.png")
            ggsave(fig4_file, p_metric, width = 8, height = 5, dpi = 300)
          }
        }
      }
    }
  }
  message("Completed writing figure 4B to file.")
}


# Main function
generate_performance_plots <- function(performance_file = "output/2. models/model_performance.csv",
                                       focused_dataset = "ctnorm_global") {
  # Use relative paths
  if (!file.exists(performance_file) && !startsWith(performance_file, "/") && !grepl("^[A-Za-z]:", performance_file)) {
    performance_file_path <- file.path(performance_file)
  } else {
    performance_file_path <- performance_file
  }

  perf_data <- load_model_performance(performance_file_path)

  # Save all barplots for all combinations
  save_all_performance_barplots(perf_data)

  return(list(data = perf_data))
}


# Run plotting function
performance_results <- generate_performance_plots(focused_dataset = "ctnorm_global")
