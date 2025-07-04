library(dplyr)
library(stringr)
library(RColorBrewer)
library(fs)
library(readr)
library(ggplot2)
library(viridis)


# Source utils
source("src/0. utils/feature_name_utils.R")


# Compute detection rates per feature
compute_detection_rates_per_feature <- function(umitab_filtered, cell_metadata_final, table_s1) {
  # Add tissue info to cell metadata
  cell_metadata_final <- cell_metadata_final %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")

  all_cell_types <- unique(cell_metadata_final$cell_type)
  detection_list <- list()

  for (ct in all_cell_types) {
    # Subset metadata for this cell type
    meta_ct <- cell_metadata_final %>% filter(cell_type == ct)
    if (nrow(meta_ct) == 0) next

    # Ensure ids are character
    meta_ct$sample_ID <- as.character(meta_ct$sample_ID)
    meta_ct$cell_ID <- as.character(meta_ct$cell_ID)

    # Get cell columns for this cell type
    cell_idx <- which(colnames(umitab_filtered) %in% meta_ct$cell_ID)
    if (length(cell_idx) == 0) next
    umi_mat_ct <- umitab_filtered[, cell_idx, drop = FALSE]

    # Map cell_ID to sample_ID for this cell type
    cell_to_sample <- meta_ct$sample_ID
    names(cell_to_sample) <- meta_ct$cell_ID
    sample_vec <- cell_to_sample[colnames(umi_mat_ct)]
    unique_samples <- unique(sample_vec)

    # Get tissue type for each sample
    sample_tissue <- table_s1 %>% mutate(sample_ID = as.character(sample_ID)) %>% select(sample_ID, tissue)
    sample_tissue_sub <- sample_tissue[match(unique_samples, sample_tissue$sample_ID), ]

    # For each sample, get column indices
    sample_col_list <- split(seq_along(sample_vec), sample_vec)

    # Aggregate counts by sample (columns)
    gene_sample_matrix <- Matrix::Matrix(0, nrow = nrow(umi_mat_ct), ncol = length(unique_samples), sparse = TRUE)
    for (j in seq_along(unique_samples)) {
      cols <- sample_col_list[[unique_samples[j]]]
      if (length(cols) == 0) next
      gene_sample_matrix[, j] <- Matrix::rowSums(umi_mat_ct[, cols, drop = FALSE])
    }
    rownames(gene_sample_matrix) <- rownames(umi_mat_ct)
    colnames(gene_sample_matrix) <- unique_samples

    # Calculate detection rates for normal and tumor
    is_normal <- sample_tissue_sub$tissue == "Normal"
    is_tumor <- sample_tissue_sub$tissue == "Tumor"
    n_normal <- sum(is_normal, na.rm = TRUE)
    n_tumor <- sum(is_tumor, na.rm = TRUE)
    det_rate_normal <- if (any(is_normal)) Matrix::rowSums(gene_sample_matrix[, is_normal, drop = FALSE] > 0) / n_normal else rep(NA, nrow(gene_sample_matrix))
    det_rate_tumor  <- if (any(is_tumor)) Matrix::rowSums(gene_sample_matrix[, is_tumor, drop = FALSE] > 0) / n_tumor else rep(NA, nrow(gene_sample_matrix))

    # Map samples to clusters
    sample_to_cluster <- sapply(unique_samples, function(sid) {
      cl <- meta_ct$cluster_ID[meta_ct$sample_ID == sid]
      if (length(cl) == 0) return(NA)
      as.character(cl[1])
    })

    # Aggregate detection rates by cluster
    detection_df_list <- list()
    for (clust in unique(sample_to_cluster)) {
      if (is.na(clust)) next
      cols_in_clust <- which(sample_to_cluster == clust)
      if (length(cols_in_clust) == 0) next
      gsm_sub <- gene_sample_matrix[, cols_in_clust, drop = FALSE]
      sts_sub <- sample_tissue_sub[cols_in_clust, , drop = FALSE]
      is_normal_sub <- sts_sub$tissue == "Normal"
      is_tumor_sub <- sts_sub$tissue == "Tumor"
      n_normal_sub <- sum(is_normal_sub, na.rm = TRUE)
      n_tumor_sub <- sum(is_tumor_sub, na.rm = TRUE)
      det_rate_normal_sub <- if (any(is_normal_sub)) Matrix::rowSums(gsm_sub[, is_normal_sub, drop = FALSE] > 0) / n_normal_sub else rep(NA, nrow(gsm_sub))
      det_rate_tumor_sub  <- if (any(is_tumor_sub)) Matrix::rowSums(gsm_sub[, is_tumor_sub, drop = FALSE] > 0) / n_tumor_sub else rep(NA, nrow(gsm_sub))
      n_detected_normal <- if (any(is_normal_sub)) Matrix::rowSums(gsm_sub[, is_normal_sub, drop = FALSE] > 0) else rep(NA, nrow(gsm_sub))
      n_detected_tumor  <- if (any(is_tumor_sub)) Matrix::rowSums(gsm_sub[, is_tumor_sub, drop = FALSE] > 0) else rep(NA, nrow(gsm_sub))

      # Create feature_id for each gene/cluster/celltype
      feature_id_vec <- create_feature_identifier(
        gene = rownames(gsm_sub),
        cluster_name_part = ct,
        cluster_id_numeric = as.numeric(clust)
      )

      detection_df_list[[clust]] <- data.frame(
        gene = rownames(gsm_sub),
        cell_type = ct,
        cluster_ID = clust,
        feature_id = feature_id_vec,
        det_rate_normal = det_rate_normal_sub,
        det_rate_tumor = det_rate_tumor_sub,
        n_detected_normal = n_detected_normal,
        n_detected_tumor = n_detected_tumor,
        n_normal = n_normal_sub,
        n_tumor = n_tumor_sub,
        stringsAsFactors = FALSE
      )
    }
    detection_df <- do.call(rbind, detection_df_list)
    detection_list[[ct]] <- detection_df
  }

  do.call(rbind, detection_list)
}


# === Load CP10K normalized counts and metadata ===
umitab_filtered <- readRDS("output/6. plots/data/cp10k/cp10k_normalized_umitab.rds")
cell_metadata_final <- readr::read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
annots_list <- readr::read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
table_s1 <- readr::read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE)

# Add cell_type column to metadata
cell_metadata_final <- cell_metadata_final %>%
  left_join(annots_list %>% select(cluster = cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage))

# Compute detection rates for all features
detection_rates_df <- compute_detection_rates_per_feature(umitab_filtered, cell_metadata_final, table_s1)


# Create all correlation heatmaps
run_all_meta_fold_detection_correlation_heatmaps <- function(
    intersector_parent_input_dir = "output/3. intersector",
    figures_parent_output_dir = "output/6. plots/figure s3/biascorrelation", 
    dataset_types_to_process = c("ctnorm_global"),
    cell_type_groups_to_process = c("all_clusters"),
    gene_sets_to_process = c("metabolic")
) {
  total_batches <- length(dataset_types_to_process) * length(cell_type_groups_to_process) * length(gene_sets_to_process)
  batch_num <- 1

  for (type in dataset_types_to_process) {
    for (cell_type_group in cell_type_groups_to_process) {
      for (gene_set in gene_sets_to_process) {
        # Set up paths
        current_input_dir <- file.path(intersector_parent_input_dir, type, cell_type_group, gene_set)
        current_figures_output_dir <- file.path(figures_parent_output_dir, type, cell_type_group, gene_set)

        # Create output directory if it doesn't exist
        if (!dir.exists(current_figures_output_dir)) {
          dir.create(current_figures_output_dir, recursive = TRUE)
        }

        # Input files
        meta_scores_file <- file.path(current_input_dir, "meta_scores.csv")
        fold_changes_file <- file.path("output", "4. fold changes", type, cell_type_group, gene_set, "fold_changes.csv")
        if (!file.exists(meta_scores_file) || !file.exists(fold_changes_file)) {
          message(sprintf("[SKIP] Missing input for: %s %s %s", type, cell_type_group, gene_set))
          batch_num <- batch_num + 1
          next
        }

        # Read input data
        meta_scores <- read_csv(meta_scores_file, show_col_types = FALSE)
        fold_changes <- read_csv(fold_changes_file, show_col_types = FALSE)

        # Join meta_scores and fold_changes on feature_id/Feature (per feature)
        meta_fc_merged <- meta_scores %>%
          left_join(fold_changes %>% select(Feature, fold_change = Value), by = c("feature_id" = "Feature"))

        # Join detection rates by feature_id (gene@celltype)
        plot_data <- meta_fc_merged %>% left_join(detection_rates_df, by = "feature_id") %>%
          filter(feature_id %in% meta_scores$feature_id)

        # Optionally load annotation file (not used in plot, but available)
        annots_file <- file.path("base", "input_tables", "annots_list.csv")
        if (file.exists(annots_file)) {
          annots <- read_csv(annots_file, show_col_types = FALSE)
        } else {
          annots <- data.frame(cluster = character(), sub_lineage = character())
        }

        # Only plot if all data is available
        if (nrow(meta_scores) > 0 && nrow(fold_changes) > 0 && nrow(detection_rates_df) > 0) {
          # Unify fold_change column
          plot_data <- plot_data %>% mutate(fold_change = coalesce(fold_change.y, fold_change.x))
          # Add detection_rate column as mean of det_rate_normal and det_rate_tumor
          plot_data <- plot_data %>% mutate(detection_rate = rowMeans(cbind(det_rate_normal, det_rate_tumor), na.rm = TRUE))
          # Calculate detection rate ratio
          epsilon <- 1e-4
          plot_data <- plot_data %>% mutate(det_ratio = (det_rate_tumor + epsilon) / (det_rate_normal + epsilon))

          # Add tissue-specific detection ratios
          plot_data <- plot_data %>%
            mutate(
              det_ratio_normal = det_rate_normal, # for clarity
              det_ratio_tumor = det_rate_tumor
            )

          # Calculate correlations between all variables of interest
          cor_vars <- c("fold_change", "meta_score", "det_ratio", "det_ratio_normal", "det_ratio_tumor")
          cor_matrix <- plot_data %>%
            select(any_of(cor_vars)) %>%
            cor(use = "pairwise.complete.obs", method = "pearson")

          # Plot correlation matrix as a heatmap with pretty variable names and no title
          cor_df <- as.data.frame(as.table(cor_matrix))
          colnames(cor_df) <- c("Var1", "Var2", "Correlation")

          # Map variable names to pretty labels
          pretty_names <- c(
            fold_change = "Fold Change",
            meta_score = "Meta Score",
            det_ratio = "Detection Ratio (T/N)",
            det_ratio_normal = "Detection Rate (N)",
            det_ratio_tumor = "Detection Rate (T)"
          )
          cor_df$Var1 <- pretty_names[as.character(cor_df$Var1)]
          cor_df$Var2 <- pretty_names[as.character(cor_df$Var2)]

          p_cor <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
            geom_tile(color = "white", width = 0.95, height = 0.95) +
            scale_fill_viridis_c(option = "C", limits = c(-1, 1)) +
            geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 6) +
            labs(x = NULL, y = NULL, fill = "Pearson r") +
            theme_minimal(base_size = 28) +
            theme(
              axis.text = element_text(size = 18),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              plot.margin = margin(10, 10, 10, 10, "pt"),
              panel.grid = element_blank(),
              legend.key.height = unit(1.5, "cm"),
              legend.key.width = unit(0.7, "cm")
            )

          # --- Additional: Unfiltered (no meta_score) correlation heatmap ---
          cor_vars_unfiltered <- c("fold_change", "det_ratio", "det_ratio_normal", "det_ratio_tumor")
          cor_matrix_unfiltered <- plot_data %>%
            select(any_of(cor_vars_unfiltered)) %>%
            cor(use = "pairwise.complete.obs", method = "pearson")

          cor_df_unfiltered <- as.data.frame(as.table(cor_matrix_unfiltered))
          colnames(cor_df_unfiltered) <- c("Var1", "Var2", "Correlation")

          pretty_names_unfiltered <- c(
            fold_change = "Fold Change",
            det_ratio = "Detection Ratio (T/N)",
            det_ratio_normal = "Detection Rate (N)",
            det_ratio_tumor = "Detection Rate (T)"
          )
          cor_df_unfiltered$Var1 <- pretty_names_unfiltered[as.character(cor_df_unfiltered$Var1)]
          cor_df_unfiltered$Var2 <- pretty_names_unfiltered[as.character(cor_df_unfiltered$Var2)]

          p_cor_unfiltered <- ggplot(cor_df_unfiltered, aes(x = Var1, y = Var2, fill = Correlation)) +
            geom_tile(color = "white", width = 0.95, height = 0.95) +
            scale_fill_viridis_c(option = "C", limits = c(-1, 1)) +
            geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 6) +
            labs(x = NULL, y = NULL, fill = "Pearson r") +
            theme_minimal(base_size = 28) +
            theme(
              axis.text = element_text(size = 18),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              plot.margin = margin(10, 10, 10, 10, "pt"),
              panel.grid = element_blank(),
              legend.key.height = unit(1.5, "cm"),
              legend.key.width = unit(0.7, "cm")
            )

          # Save correlation heatmap to biascorrelation directory
          ggsave(sprintf("output/6. plots/figure s3/biascorrelation/correlation_heatmap_%s_%s_%s.png", type, cell_type_group, gene_set), p_cor, width = 11, height = 8, dpi = 300)

          # If this is the main version for figure 6 panels d and e
          if (type == "ctnorm_global" && cell_type_group == "all_clusters" && gene_set == "metabolic") {
            fig6_dir <- "output/6. plots/figure 6"
            if (!dir.exists(fig6_dir)) dir.create(fig6_dir, recursive = TRUE)
            ggsave(file.path(fig6_dir, "fig_6f_correlation_heatmap.png"), p_cor, width = 13, height = 10, dpi = 300)
            ggsave(file.path(fig6_dir, "fig_6e_correlation_heatmap_unfiltered.png"), p_cor_unfiltered, width = 13, height = 10, dpi = 300)
          }
        } else {
          message(sprintf("[SKIP] No data to plot for: %s %s %s", type, cell_type_group, gene_set))
        }
        batch_num <- batch_num + 1
      }
    }
  }
}


# === EXECUTE CORRELATION HEATMAP PLOTTING ===
run_all_meta_fold_detection_correlation_heatmaps()
