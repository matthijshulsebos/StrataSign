library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(ggrepel)

# Source utility functions
source("src/0. utils/feature_name_utils.R")

# Load detection rates data once
detection_rates <- read_csv("output/6. plots/data/detection/detection_rates.csv", show_col_types = FALSE) %>%
  mutate(
    det_diff = (n_detected_tumor / n_tumor) - (n_detected_normal / n_normal)
  ) %>%
  select(feature_id, det_diff)

# Define a reusable plotting function
create_and_save_volcano <- function(data, x_var, x_label, cutoff_val, file_path) {
  # Create significance categories based on the x-variable
  data <- data %>%
    mutate(
      significance = case_when(
        meta_score >= meta_score_cutoff & abs(.data[[x_var]]) >= cutoff_val ~ "High importance & High Effect",
        meta_score >= meta_score_cutoff & abs(.data[[x_var]]) < cutoff_val ~ "High importance & Low Effect",
        meta_score < meta_score_cutoff & abs(.data[[x_var]]) >= cutoff_val ~ "Low importance & High Effect",
        TRUE ~ "Low importance & Low Effect"
      ),
      significance = factor(significance, levels = c(
        "High importance & High Effect",
        "High importance & Low Effect", 
        "Low importance & High Effect",
        "Low importance & Low Effect"
      ))
    )

  # Create volcano plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = meta_score)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 1.5) + 
    geom_hline(yintercept = meta_score_cutoff, linetype = "dashed", color = "red", alpha = 0.7) +
    annotate(
      "text",
      x = -Inf, y = meta_score_cutoff, hjust = -0.1, vjust = -0.5, size = 4, color = "black",
      label = sprintf("100th highest meta-score = %.2f", meta_score_cutoff)
    ) +
    geom_vline(xintercept = c(-cutoff_val, cutoff_val), linetype = "dashed", color = "grey40") +
    annotate("text", x = cutoff_val, y = 0, label = as.character(cutoff_val), vjust = -0.5, hjust = -0.2, size = 3.5, color = "grey40") +
    annotate("text", x = -cutoff_val, y = 0, label = as.character(-cutoff_val), vjust = -0.5, hjust = 1.2, size = 3.5, color = "grey40") +
    geom_text_repel(
      # Corrected data filtering: include zero values by removing the lower bound check
      data = data %>% 
               filter(abs(.data[[x_var]]) < cutoff_val) %>% 
               arrange(desc(meta_score)) %>% 
               slice_head(n = 5),
      aes(label = paste(gene, sublineage, sep = ", ")), 
      size = 3, 
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = 'black', 
      segment.size = 0.6,      
      min.segment.length = 0,  
      max.overlaps = Inf 
    ) +
    scale_color_manual(
      name = "Significance",
      values = c(
        "High importance & High Effect" = "#e74c3c",
        "High importance & Low Effect" = "#f39c12", 
        "Low importance & High Effect" = "#3498db",
        "Low importance & Low Effect" = "#95a5a6"
      )
    ) +
    labs(
      x = x_label,
      y = "Meta Importance Score"
    ) +
    theme_bw(base_size = 16, base_family = "sans") +
    theme(
      legend.position = "right",
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.margin = margin(8, 8, 8, 8),
      panel.grid.minor = element_blank(),
      plot.title = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA)
    ) +
    guides(
      color = guide_legend(title = "Significance", override.aes = list(size = 3, alpha = 1))
    )

  # Ensure directory exists and save the plot
  if (!dir.exists(dirname(file_path))) dir.create(dirname(file_path), recursive = TRUE)
  ggsave(file_path, p, width = 12, height = 8, dpi = 300, bg = "transparent")
}


# Define all combinations to plot
norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both")
gene_sets <- c("metabolic", "nonmetabolic", "random")

# Load metabolic genes for classification
hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", show_col_types = FALSE)

# Iterate over all combinations
for (norm in norms) {
  for (cell_type in cell_types) {
    for (gene_set in gene_sets) {
      meta_scores_file <- file.path("output/3. intersector", norm, cell_type, gene_set, "meta_scores.csv")
      fc_file <- file.path("output/4. fold changes", norm, cell_type, gene_set, "fold_changes.csv")
      if (!file.exists(meta_scores_file) || !file.exists(fc_file)) {
        message(sprintf("Missing file for: norm=%s, cell_type=%s, gene_set=%s", norm, cell_type, gene_set))
        next
      }
      meta_scores <- read_csv(meta_scores_file, show_col_types = FALSE)
      fc <- read_csv(fc_file, show_col_types = FALSE)
      combined_data <- meta_scores %>%
        left_join(fc %>% dplyr::rename(log2_fold_change = Value), by = c("feature_id" = "Feature")) %>%
        left_join(detection_rates, by = "feature_id") %>% # Join detection difference
        mutate(
          gene = get_gene_from_feature(feature_id),
          sublineage = get_simplified_sublineage(feature_id)
        )

      # Only filter to metabolic genes for the metabolic gene set
      if (gene_set == "metabolic") {
        combined_data <- combined_data %>%
          filter(gene %in% hsa01100_genes$Symbol) %>%
          mutate(gene_type = "Metabolic")
      }

      # Filter out extreme fold changes and detection differences
      combined_data <- combined_data %>% 
        filter(abs(log2_fold_change) < 10) %>%
        filter(abs(det_diff) <= 1.0)

      # Remove rows with NAs in columns critical for plotting
      combined_data <- combined_data %>% 
        filter(!is.na(log2_fold_change) & !is.na(meta_score) & !is.na(sublineage) & !is.na(det_diff))

      if (nrow(combined_data) == 0) next

      # Calculate the 100th highest meta score as cutoff (or max if <100 rows)
      meta_score_vec <- combined_data %>% arrange(desc(meta_score)) %>% pull(meta_score)
      meta_score_cutoff <- if (length(meta_score_vec) >= 100) meta_score_vec[100] else utils::tail(meta_score_vec, 1)
      
      # Generate and save the fold change volcano plot
      fc_volcano_dir <- file.path("output/6. plots/figure 9/fold_change_volcano", norm, cell_type, gene_set)
      fc_volcano_path <- file.path(fc_volcano_dir, sprintf("volcano_fc_%s_%s_%s.png", norm, cell_type, gene_set))
      create_and_save_volcano(
        data = combined_data,
        x_var = "log2_fold_change",
        x_label = "Log2 Fold Change (Tumor vs Normal)",
        cutoff_val = 1,
        file_path = fc_volcano_path
      )

      # Generate and save the detection difference volcano plot
      dd_volcano_dir <- file.path("output/6. plots/figure 9/detection_diff_volcano", norm, cell_type, gene_set)
      dd_volcano_path <- file.path(dd_volcano_dir, sprintf("volcano_dd_%s_%s_%s.png", norm, cell_type, gene_set))
      create_and_save_volcano(
        data = combined_data,
        x_var = "det_diff",
        x_label = "Detection Difference (tumor - normal)",
         # Using 25% as the significance cutoff for detection difference
        cutoff_val = 0.25,
        file_path = dd_volcano_path
      )
    }
  }
}
