# Figure 6F/G: Distribution comparison and statistical test for meta score bias
library(dplyr)
library(ggplot2)
library(readr)



# Helper to format p-values for plotting
format_p_value <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  } else {
    return(paste("p =", round(p, 3)))
  }
}

# Iterate over all combinations with nested for loops
norms <- c("ctnorm_global", "ctnorm_global_zscaled", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both")
gene_sets <- c("metabolic", "nonmetabolic", "random")

detection_rates_all <- read_csv("output/6. plots/data/detection/detection_rates.csv", show_col_types = FALSE)

for (norm in norms) {
  for (cell_type in cell_types) {
    for (gene_set in gene_sets) {
      meta_scores_file <- file.path("output/3. intersector", norm, cell_type, gene_set, "meta_scores.csv")
      fold_changes_file <- file.path("output/4. fold changes", norm, cell_type, gene_set, "fold_changes.csv")
      if (!file.exists(meta_scores_file) || !file.exists(fold_changes_file)) {
        message(sprintf("[SKIP] Missing input for: %s %s %s", norm, cell_type, gene_set))
        next
      }
      meta_scores <- read_csv(meta_scores_file, show_col_types = FALSE)
      fold_changes <- read_csv(fold_changes_file, show_col_types = FALSE)
      detection_rates <- detection_rates_all

      # Detection difference = (n_detected_tumor / n_tumor) - (n_detected_normal / n_normal)
      if (all(c("feature_id", "n_detected_tumor", "n_tumor", "n_detected_normal", "n_normal") %in% colnames(detection_rates))) {
        detection_rates <- detection_rates %>%
          mutate(
            det_rate_tumor = n_detected_tumor / n_tumor,
            det_rate_normal = n_detected_normal / n_normal,
            det_diff = det_rate_tumor - det_rate_normal
          )
      } else if (all(c("feature_id", "det_rate_tumor", "det_rate_normal") %in% colnames(detection_rates))) {
        detection_rates <- detection_rates %>%
          mutate(det_diff = det_rate_tumor - det_rate_normal)
      } else {
        stop("Detection rates file does not have expected columns.")
      }

      plot_data <- fold_changes %>%
        rename(feature_id = Feature, fold_change = Value) %>%
        left_join(meta_scores %>% select(feature_id, meta_score), by = "feature_id") %>%
        left_join(detection_rates %>% select(feature_id, det_diff, n_detected_normal, n_detected_tumor), by = "feature_id") %>%
        filter(!(n_detected_normal == 0 & n_detected_tumor == 0))

      plot_data <- plot_data %>%
        mutate(has_meta_score = !is.na(meta_score)) %>%
        filter(
          is.finite(det_diff),
          is.finite(fold_change),
          !is.na(has_meta_score)
        )

      plot_data$has_meta_score <- factor(plot_data$has_meta_score, levels = c(FALSE, TRUE), labels = c("No meta score", "Has meta score"))

      # Only plot if both groups are present
      if (length(unique(plot_data$has_meta_score)) < 2) {
        message(sprintf("[SKIP] Not enough groups for: %s %s %s", norm, cell_type, gene_set))
        next
      }

      # 1. Perform statistical tests first
      fc_test <- wilcox.test(fold_change ~ has_meta_score, data = plot_data)
      dd_test <- wilcox.test(det_diff ~ has_meta_score, data = plot_data)

      fc_p_value <- format_p_value(fc_test$p.value)
      dd_p_value <- format_p_value(dd_test$p.value)

      # 2. Fold change distribution plot
      p_fc <- ggplot(plot_data, aes(x = fold_change, fill = has_meta_score)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = c("grey70", "#1f77b4"), labels = c("No meta score", "Has meta score")) +
        annotate("text", x = Inf, y = Inf, label = fc_p_value, hjust = 1.1, vjust = 1.5, size = 5) +
        labs(x = "log2(Fold Change)", y = "Density", fill = "Meta score") +
        theme_bw(base_size = 16)

      # 3. Detection difference distribution plot
      p_dr <- ggplot(plot_data, aes(x = det_diff, fill = has_meta_score)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = c("grey70", "#1f77b4"), labels = c("No meta score", "Has meta score")) +
        annotate("text", x = Inf, y = Inf, label = dd_p_value, hjust = 1.1, vjust = 1.5, size = 5) +
        labs(x = "Detection difference (T - N)", y = "Density", fill = "Meta score") +
        theme_bw(base_size = 16)

      # Output directory structure
      outdir <- file.path("output/6. plots/figure s7/dist_test", norm, cell_type, gene_set)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      ggsave(file.path(outdir, "fold_change_dist.png"), p_fc, width = 7, height = 5, dpi = 300)
      ggsave(file.path(outdir, "detection_difference_dist.png"), p_dr, width = 7, height = 5, dpi = 300)

      # Print test results to console
      cat(sprintf("\n[%s | %s | %s] Fold change distribution test (Wilcoxon):\n", norm, cell_type, gene_set))
      print(fc_test)
      cat(sprintf("\n[%s | %s | %s] Detection difference distribution test (Wilcoxon):\n", norm, cell_type, gene_set))
      print(dd_test)
    }
  }
}
