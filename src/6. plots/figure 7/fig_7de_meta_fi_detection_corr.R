library(dplyr)
library(stringr)
library(RColorBrewer)
library(fs)
library(readr)
library(ggplot2)
library(viridis)


# Source utils
source("src/0. utils/feature_name_utils.R")

# === Load detection rates from CSV (precomputed) ===
detection_rates_df <- readr::read_csv("output/6. plots/data/detection/detection_rates.csv", show_col_types = FALSE)


# === For all dataset combinations, create and save per-sublineage heatmaps ===
run_all_sublineage_heatmaps <- function(
    intersector_parent_input_dir = "output/3. intersector",
    dataset_types_to_process = c("ctnorm_global"),
    cell_type_groups_to_process = c("all_clusters"),
    gene_sets_to_process = c("metabolic")
) {
  for (type in dataset_types_to_process) {
    for (cell_type_group in cell_type_groups_to_process) {
      for (gene_set in gene_sets_to_process) {
        # Construct input directory and check if files exist
        current_input_dir <- file.path(intersector_parent_input_dir, type, cell_type_group, gene_set)
        meta_scores_file <- file.path(current_input_dir, "meta_scores.csv")
        fold_changes_file <- file.path("output", "4. fold changes", type, cell_type_group, gene_set, "fold_changes.csv")
        if (!file.exists(meta_scores_file) || !file.exists(fold_changes_file)) {
          message(sprintf("[SKIP] Missing input for: %s %s %s", type, cell_type_group, gene_set))
          next
        }
        meta_scores <- read_csv(meta_scores_file, show_col_types = FALSE)
        fold_changes <- read_csv(fold_changes_file, show_col_types = FALSE)
        create_and_save_sublineage_heatmaps_full(meta_scores, fold_changes, detection_rates_df)
      }
    }
  }
}


create_and_save_sublineage_heatmaps_full <- function(meta_scores, fold_changes, detection_rates_df) {
  # Always extract sublineage and cluster_ID from feature_id for all datasets
  extract_sublineage <- function(fid) sub("_[0-9]+$", "", sapply(strsplit(fid, "@"), function(x) x[[length(x)]]))
  extract_cluster <- function(fid) sub(".*_([0-9]+)$", "\\1", sapply(strsplit(fid, "@"), function(x) x[[length(x)]]))

  meta_scores$sublineage <- extract_sublineage(meta_scores$feature_id)
  meta_scores$cluster_ID <- extract_cluster(meta_scores$feature_id)
  detection_rates_df$sublineage <- extract_sublineage(detection_rates_df$feature_id)
  detection_rates_df$cluster_ID <- extract_cluster(detection_rates_df$feature_id)
  fold_changes$sublineage <- extract_sublineage(fold_changes$Feature)
  fold_changes$cluster_ID <- extract_cluster(fold_changes$Feature)

  combos <- unique(meta_scores[, c("sublineage", "cluster_ID")])
  combos <- combos[!is.na(combos$sublineage) & !is.na(combos$cluster_ID), ]
  # Use a subdirectory structure for each combination
  outdir_s7 <- file.path(
    "output/6. plots/figure s7/detection_correlation",
    norm, cell_type, gene_set
  )
  dir.create(outdir_s7, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(combos))) {

    sub <- combos$sublineage[i]
    clust <- combos$cluster_ID[i]
    ms_sub <- meta_scores %>% filter(sublineage == sub, cluster_ID == clust)
    dr_sub <- detection_rates_df %>% filter(sublineage == sub, cluster_ID == clust)
    fc_sub <- fold_changes %>% filter(sublineage == sub, cluster_ID == clust)

    # Ensure join columns are character
    ms_sub$feature_id <- trimws(as.character(ms_sub$feature_id))
    fc_sub$Feature <- trimws(as.character(fc_sub$Feature))

    # Remove fold_change column from dr_sub if present to avoid duplicate columns in join
    if ("fold_change" %in% names(dr_sub)) {
      dr_sub <- dr_sub[, setdiff(names(dr_sub), "fold_change"), drop = FALSE]
    }

    # Remove fold_change from ms_sub if present (shouldn't be, but for robustness)
    if ("fold_change" %in% names(ms_sub)) {
      ms_sub <- ms_sub[, setdiff(names(ms_sub), "fold_change"), drop = FALSE]
    }

    # Use base R merge for fold_change join, then left_join for detection rates
    plot_data <- merge(ms_sub, fc_sub[, c("Feature", "Value")], by.x = "feature_id", by.y = "Feature", all.x = TRUE)
    names(plot_data)[names(plot_data) == "Value"] <- "fold_change"

    # Remove fold_change from plot_data if present (should only be one, but for robustness)
    if ("fold_change" %in% names(plot_data) && sum(names(plot_data) == "fold_change") > 1) {
      # Keep only the first occurrence
      keep_cols <- which(names(plot_data) == "fold_change")
      plot_data <- plot_data[, -keep_cols[-1], drop = FALSE]
    }

    plot_data <- left_join(plot_data, dr_sub, by = "feature_id")

    if (nrow(plot_data) < 2) next
    plot_data <- plot_data %>% mutate(
      detection_rate = rowMeans(cbind(det_rate_normal, det_rate_tumor), na.rm = TRUE),
      det_ratio = (det_rate_tumor + 1e-4) / (det_rate_normal + 1e-4),
      det_ratio_normal = det_rate_normal,
      det_ratio_tumor = det_rate_tumor
    )
    cor_vars <- c("fold_change", "meta_score", "det_ratio", "det_ratio_normal", "det_ratio_tumor")
    available_vars <- intersect(cor_vars, colnames(plot_data))
    keep_vars <- sapply(plot_data[, available_vars, drop = FALSE], function(x) length(unique(na.omit(x))) > 1)
    use_vars <- available_vars[keep_vars]
    if (length(use_vars) < 2) next
    cor_matrix <- plot_data %>% select(any_of(use_vars)) %>% cor(use = "pairwise.complete.obs", method = "pearson")
    cor_df <- as.data.frame(as.table(cor_matrix))
    colnames(cor_df) <- c("Var1", "Var2", "Correlation")
    pretty_names <- c(
      fold_change = "Fold Change",
      meta_score = "Meta Score",
      det_ratio = "Detection Ratio (T/N)",
      det_ratio_normal = "Detection Rate (N)",
      det_ratio_tumor = "Detection Rate (T)"
    )
    cor_df$Var1 <- pretty_names[as.character(cor_df$Var1)]
    cor_df$Var2 <- pretty_names[as.character(cor_df$Var2)]
    sub_title <- paste0("Sublineage: ", sub, " (", clust, ")")
    p_cor <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
      geom_tile(color = "white", width = 0.95, height = 0.95) +
      scale_fill_viridis_c(option = "C", limits = c(-1, 1)) +
      geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 4, family = "sans") +
      labs(x = NULL, y = NULL, fill = "Pearson r", title = sub_title) +
      theme_bw(base_size = 16, base_family = "sans") +
      theme(
        plot.title = element_text(color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.grid = element_blank(),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(0.7, "cm"),
        strip.text = element_text(color = "black", face = "bold")
      )
    safe_sub <- gsub("[^a-zA-Z0-9]+", "_", paste0(sub, "_", clust))
    ggsave(file.path(outdir_s7, paste0("detection_correlation_", safe_sub, ".png")), p_cor, width = 11, height = 8, dpi = 300)
  }
}




# === ITERATE OVER ALL COMBINATIONS WITH NESTED FOR LOOPS ===
norms <- c("ctnorm_global", "ctnorm_global_zscaled", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both")
gene_sets <- c("metabolic", "nonmetabolic", "random")

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
      create_and_save_sublineage_heatmaps_full(meta_scores, fold_changes, detection_rates_df)
    }
  }
}
