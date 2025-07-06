# Load metabolic gene symbols
metabolic_genes <- readr::read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", show_col_types = FALSE)$Symbol
# Scatter plot of cp10k normalized total counts vs. number of detected genes per cell type per sample
library(ggplot2)
library(dplyr)
library(Matrix)
library(readr)

# --- Load data ---
# Use the same data sources as detection.R
umitab_filtered <- readRDS("output/6. plots/data/cp10k/cp10k_normalized_umitab.rds")
cell_metadata_final <- readr::read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
annots_list <- readr::read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
table_s1 <- readr::read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE)

# Add cell_type column to metadata (using sub_lineage if available, else lineage)
cell_metadata_final <- cell_metadata_final %>%
  left_join(annots_list %>% select(cluster = cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage))

# Output root for figure s6
ds6_root <- file.path("output/6. plots/figure s6/detection_counts_scatter")
dir.create(ds6_root, recursive = TRUE, showWarnings = FALSE)

all_cell_types <- unique(cell_metadata_final$cell_type)

scatter_df_list <- list()

for (ct in all_cell_types) {
  # Subset metadata for this cell type
  meta_ct <- cell_metadata_final %>% filter(cell_type == ct)
  if (nrow(meta_ct) == 0) next

  # Ensure IDs are character
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

  # For each sample, get column indices
  sample_col_list <- split(seq_along(sample_vec), sample_vec)

  # Aggregate counts by sample (columns)
  sample_names <- names(sample_col_list)
  gene_sample_matrix <- Matrix::Matrix(0, nrow = nrow(umi_mat_ct), ncol = length(sample_names), sparse = TRUE)
  for (j in seq_along(sample_names)) {
    cols <- sample_col_list[[sample_names[j]]]
    if (length(cols) == 0) next
    gene_sample_matrix[, j] <- Matrix::rowSums(umi_mat_ct[, cols, drop = FALSE])
  }
  rownames(gene_sample_matrix) <- rownames(umi_mat_ct)
  colnames(gene_sample_matrix) <- sample_names

  # For each sample, count number of detected genes (nonzero rows)
  n_detected_genes <- Matrix::colSums(gene_sample_matrix > 0)
  # For each sample, sum total cp10k counts
  total_counts_cp10k <- Matrix::colSums(gene_sample_matrix)

  scatter_df <- data.frame(
    cell_type = ct,
    sample_ID = unique_samples,
    total_counts_cp10k = as.numeric(total_counts_cp10k),
    n_detected_genes = as.numeric(n_detected_genes),
    stringsAsFactors = FALSE
  )
  scatter_df_list[[ct]] <- scatter_df
}

plot_df <- do.call(rbind, scatter_df_list)

# Add tissue type (Normal/Tumor) to plot_df using table_s1
plot_df <- plot_df %>%
  left_join(table_s1 %>% mutate(sample_ID = as.character(sample_ID)) %>% select(sample_ID, tissue), by = "sample_ID")

# Remove outliers based on 1st and 99th percentiles for both axes
x_q <- quantile(plot_df$total_counts_cp10k, probs = c(0.01, 0.99), na.rm = TRUE)
y_q <- quantile(plot_df$n_detected_genes, probs = c(0.01, 0.99), na.rm = TRUE)
plot_df_filt <- plot_df %>%
  filter(
    total_counts_cp10k >= x_q[1], total_counts_cp10k <= x_q[2],
    n_detected_genes >= y_q[1], n_detected_genes <= y_q[2]
  )

# Define color mapping for tissue: Tumor = red, Normal = blue
tissue_colors <- c("Tumor" = "#D62728", "Normal" = "#377EB8")

# Define p_tissue for all-cell-types colored by tissue, using specified colors
p_tissue <- ggplot(plot_df_filt, aes(x = total_counts_cp10k, y = n_detected_genes, color = tissue)) +
  geom_point(alpha = 0.7, size = 2.8) +
  scale_color_manual(values = tissue_colors, na.value = "grey50") +
  theme_bw(base_size = 18) +
  labs(
    x = "Total counts (cp10k, per cell type/sample)",
    y = "Number of detected genes",
    title = NULL,
    color = "Tissue"
  )

# Plot all cell types together, color by cell type
p <- ggplot(plot_df_filt, aes(x = total_counts_cp10k, y = n_detected_genes, color = cell_type)) +
  geom_point(alpha = 0.7, size = 2.8) +
  theme_bw(base_size = 18) +
  labs(
    x = "Total counts (cp10k, per cell type/sample)",
    y = "Number of detected genes",
    title = NULL,
    color = "Cell type"
  )

# Save to figure s6/detection_counts_scatter/all_cell_types_scatter.png
s6_file <- file.path(ds6_root, "all_cell_types_scatter.png")
ggsave(s6_file, p, width = 12, height = 8)

# Save per-cell-type figures in figure s6/scatter (no outlier filtering)
s6_scatter_dir <- file.path("output/6. plots/figure s6/scatter")
dir.create(s6_scatter_dir, recursive = TRUE, showWarnings = FALSE)
for (ct in unique(plot_df$cell_type)) {
  df_ct <- plot_df %>% filter(cell_type == ct)
  # Compute axis limits for both all-genes and metabolic plots before filtering
  # (use all values for this cell type, not just filtered)
  # All-genes axis limits
  xlim_all <- range(df_ct$total_counts_cp10k, na.rm = TRUE)
  ylim_all <- range(df_ct$n_detected_genes, na.rm = TRUE)

  # Prepare metabolic data for this cell type (for axis limits)
  meta_ct <- cell_metadata_final %>% filter(cell_type == ct)
  cell_idx <- which(colnames(umitab_filtered) %in% meta_ct$cell_ID)
  xlim_metabolic <- NULL
  ylim_metabolic <- NULL
  scatter_df_metabolic <- NULL
  if (length(cell_idx) > 0) {
    umi_mat_ct <- umitab_filtered[, cell_idx, drop = FALSE]
    metabolic_idx <- which(rownames(umi_mat_ct) %in% metabolic_genes)
    if (length(metabolic_idx) > 0) {
      umi_mat_metabolic <- umi_mat_ct[metabolic_idx, , drop = FALSE]
      cell_to_sample <- meta_ct$sample_ID
      names(cell_to_sample) <- meta_ct$cell_ID
      sample_vec <- cell_to_sample[colnames(umi_mat_metabolic)]
      sample_col_list <- split(seq_along(sample_vec), sample_vec)
      sample_names <- names(sample_col_list)
      gene_sample_matrix <- Matrix::Matrix(0, nrow = nrow(umi_mat_metabolic), ncol = length(sample_names), sparse = TRUE)
      for (j in seq_along(sample_names)) {
        cols <- sample_col_list[[sample_names[j]]]
        if (length(cols) == 0) next
        gene_sample_matrix[, j] <- Matrix::rowSums(umi_mat_metabolic[, cols, drop = FALSE])
      }
      rownames(gene_sample_matrix) <- rownames(umi_mat_metabolic)
      colnames(gene_sample_matrix) <- sample_names
      n_detected_genes <- Matrix::colSums(gene_sample_matrix > 0)
      total_counts_cp10k <- Matrix::colSums(gene_sample_matrix)
      scatter_df_metabolic <- data.frame(
        cell_type = ct,
        sample_ID = sample_names,
        total_counts_cp10k = as.numeric(total_counts_cp10k),
        n_detected_genes = as.numeric(n_detected_genes),
        stringsAsFactors = FALSE
      )
      scatter_df_metabolic$sample_ID <- as.character(scatter_df_metabolic$sample_ID)
      scatter_df_metabolic <- scatter_df_metabolic %>%
        left_join(table_s1 %>% mutate(sample_ID = as.character(sample_ID)) %>% select(sample_ID, tissue), by = "sample_ID")
      xlim_metabolic <- range(scatter_df_metabolic$total_counts_cp10k, na.rm = TRUE)
      ylim_metabolic <- range(scatter_df_metabolic$n_detected_genes, na.rm = TRUE)
    }
  }
  # Use the union of the ranges for both plots
  xlim_both <- range(c(xlim_all, xlim_metabolic), na.rm = TRUE)
  ylim_both <- range(c(ylim_all, ylim_metabolic), na.rm = TRUE)
  if (all(is.na(df_ct$tissue))) {
    warning(sprintf("No tissue annotation for cell type: %s. Plotting with default color.", ct))
    p_ct <- ggplot(df_ct, aes(x = total_counts_cp10k, y = n_detected_genes)) +
      geom_point(alpha = 0.7, size = 2.8, color = "#377EB8") +
      theme_bw(base_size = 18) +
      labs(
        x = "Total counts (cp10k, per cell type/sample)",
        y = "Number of detected genes",
        title = ct,
        color = NULL
      ) +
      xlim(xlim_both) + ylim(ylim_both)
  } else {
    df_ct$tissue <- as.factor(df_ct$tissue)
    p_ct <- ggplot(df_ct, aes(x = total_counts_cp10k, y = n_detected_genes, color = tissue)) +
      geom_point(alpha = 0.7, size = 2.8) +
      scale_color_manual(values = tissue_colors, na.value = "grey50") +
      theme_bw(base_size = 18) +
      labs(
        x = "Total counts (cp10k, per cell type/sample)",
        y = "Number of detected genes",
        title = ct,
        color = "Tissue"
      ) +
      xlim(xlim_both) + ylim(ylim_both)
  }
  safe_ct <- gsub("[^a-zA-Z0-9_]+", "_", ct)
  ggsave(file.path(s6_scatter_dir, paste0("scatter_", safe_ct, ".png")), p_ct, width = 12, height = 8)

  # Metabolic genes only: use precomputed scatter_df_metabolic (if available)
  if (!is.null(scatter_df_metabolic)) {
    if (all(is.na(scatter_df_metabolic$tissue))) {
      p_ct_metabolic <- ggplot(scatter_df_metabolic, aes(x = total_counts_cp10k, y = n_detected_genes)) +
        geom_point(alpha = 0.7, size = 2.8, color = "#377EB8") +
        theme_bw(base_size = 18) +
        labs(
          x = "Total counts (cp10k, per cell type/sample)",
          y = "Number of detected metabolic genes",
          title = paste0(ct, " (metabolic)"),
          color = NULL
        ) +
        xlim(xlim_both) + ylim(ylim_both)
    } else {
      scatter_df_metabolic$tissue <- as.factor(scatter_df_metabolic$tissue)
      p_ct_metabolic <- ggplot(scatter_df_metabolic, aes(x = total_counts_cp10k, y = n_detected_genes, color = tissue)) +
        geom_point(alpha = 0.7, size = 2.8) +
        scale_color_manual(values = tissue_colors, na.value = "grey50") +
        theme_bw(base_size = 18) +
        labs(
          x = "Total counts (cp10k, per cell type/sample)",
          y = "Number of detected metabolic genes",
          title = paste0(ct, " (metabolic)"),
          color = "Tissue"
        ) +
        xlim(xlim_both) + ylim(ylim_both)
    }
    ggsave(file.path(s6_scatter_dir, paste0("scatter_", safe_ct, "_metabolic.png")), p_ct_metabolic, width = 12, height = 8)
  }
}

# Save all-cell-types and tissue plots in s6/scatter as well
all_types_file <- file.path(s6_scatter_dir, "all_cell_types_scatter.png")
ggsave(all_types_file, p, width = 12, height = 8)
all_types_tissue_file <- file.path(s6_scatter_dir, "all_cell_types_scatter_by_tissue.png")
ggsave(all_types_tissue_file, p_tissue, width = 12, height = 8)

# Save the main paper figures in figure 6 directory
fig6_dir <- "output/6. plots/figure 6"
dir.create(fig6_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(fig6_dir, "fig_6h_detection_counts_scatter.png"), p, width = 12, height = 8)
ggsave(file.path(fig6_dir, "fig_6i_detection_counts_scatter.png"), p_tissue, width = 12, height = 8)
