# Scatter plots of gene detection in normal vs tumor for selected cell types
library(ggplot2)
library(Matrix)
library(dplyr)
library(readr)
library(tidyr)

# Load CP10K normalized counts and metadata
umitab_filtered <- readRDS("output/6. plots/data/cp10k/cp10k_normalized_umitab.rds")
cell_metadata_final <- readr::read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
annots_list <- readr::read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
table_s1 <- readr::read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE)

# Add (sub)lineage and tissue annotations to each cell
cell_metadata_final <- cell_metadata_final %>%
  left_join(annots_list %>% select(cluster = cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage)) %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")


# Set combination variables
plot_config <- data.frame(
  sublineage = c("MoMac-II", "T_activated", "Treg", "AM"),
  meta_celltype_set = c("all_clusters"),
  gene_type = c("metabolic"),
  normalization_method = c("ctnorm_global"),
  stringsAsFactors = FALSE
)

cell_types_focus <- plot_config$sublineage
all_sublineages <- sort(unique(cell_metadata_final$cell_type))


# Check how many unique samples per tissue type are present in the cp10k-filtered data
cp10k_sample_ids <- unique(cell_metadata_final$sample_ID)
cp10k_sample_tissue <- table_s1 %>% filter(sample_ID %in% cp10k_sample_ids)


# Get total number of unique normal and tumor samples in cp10k-filtered data
cp10k_n_normal <- cp10k_sample_tissue %>% filter(tissue == "Normal") %>% pull(sample_ID) %>% unique() %>% length()
cp10k_n_tumor  <- cp10k_sample_tissue %>% filter(tissue == "Tumor") %>% pull(sample_ID) %>% unique() %>% length()


# First, determine a global color scale limit (99th percentile of all bin counts across all cell types)
all_bin_counts <- c()
for (ct in cell_types_focus) {
  meta_ct <- cell_metadata_final %>% filter(cell_type == ct)
  if (nrow(meta_ct) == 0) next
  meta_ct$sample_ID <- as.character(meta_ct$sample_ID)
  sample_cell_list <- split(meta_ct$cell_ID, meta_ct$sample_ID)
  gene_sample_matrix <- sapply(sample_cell_list, function(cells) {
    if (length(cells) == 0) return(rep(0, nrow(umitab_filtered)))
    Matrix::rowSums(umitab_filtered[, cells, drop = FALSE])
  })
  colnames(gene_sample_matrix) <- as.character(names(sample_cell_list))
  sample_tissue <- table_s1 %>% mutate(sample_ID = as.character(sample_ID)) %>% select(sample_ID, tissue)
  shared_samples <- intersect(colnames(gene_sample_matrix), sample_tissue$sample_ID)
  if (length(shared_samples) == 0) next
  gene_sample_matrix <- gene_sample_matrix[, shared_samples, drop = FALSE]
  sample_tissue_sub <- sample_tissue[match(shared_samples, sample_tissue$sample_ID), ]
  is_normal <- sample_tissue_sub$tissue == "Normal"
  is_tumor <- sample_tissue_sub$tissue == "Tumor"
  detection_df <- data.frame(
    gene = rownames(gene_sample_matrix),
    detected_in_normal = if (any(is_normal)) Matrix::rowSums(gene_sample_matrix[, is_normal, drop = FALSE] > 0) else 0,
    detected_in_tumor = if (any(is_tumor)) Matrix::rowSums(gene_sample_matrix[, is_tumor, drop = FALSE] > 0) else 0
  )
  # Add feature_id column for robust feature-level joins (gene@celltype)
  detection_df$feature_id <- paste0(detection_df$gene, "@", ct)
  # Bin counts for this cell type
  bin_table <- as.data.frame(table(detection_df$detected_in_normal, detection_df$detected_in_tumor))
  all_bin_counts <- c(all_bin_counts, bin_table$Freq)
}
color_limit <- unname(quantile(all_bin_counts, 0.95, na.rm = TRUE))


# ---- PAIRED PLOTTING LOOP FOR ALL SUBLINEAGES ----
all_sublineages <- sort(unique(cell_metadata_final$cell_type))

for (ct in all_sublineages) {
  meta_ct <- cell_metadata_final %>% filter(cell_type == ct)
  if (nrow(meta_ct) == 0) next
  meta_ct$sample_ID <- as.character(meta_ct$sample_ID)
  sample_cell_list <- split(meta_ct$cell_ID, meta_ct$sample_ID)
  gene_sample_matrix <- sapply(sample_cell_list, function(cells) {
    if (length(cells) == 0) return(rep(0, nrow(umitab_filtered)))
    Matrix::rowSums(umitab_filtered[, cells, drop = FALSE])
  })
  colnames(gene_sample_matrix) <- as.character(names(sample_cell_list))
  sample_tissue <- table_s1 %>% mutate(sample_ID = as.character(sample_ID)) %>% select(sample_ID, tissue)
  shared_samples <- intersect(colnames(gene_sample_matrix), sample_tissue$sample_ID)
  if (length(shared_samples) == 0) {
    message(sprintf("No shared samples for cell type %s", ct))
    next
  }
  gene_sample_matrix <- gene_sample_matrix[, shared_samples, drop = FALSE]
  sample_tissue_sub <- sample_tissue[match(shared_samples, sample_tissue$sample_ID), ]
  is_normal <- sample_tissue_sub$tissue == "Normal"
  is_tumor <- sample_tissue_sub$tissue == "Tumor"
  detection_df <- data.frame(
    gene = rownames(gene_sample_matrix),
    detected_in_normal = if (any(is_normal)) Matrix::rowSums(gene_sample_matrix[, is_normal, drop = FALSE] > 0) else 0,
    detected_in_tumor = if (any(is_tumor)) Matrix::rowSums(gene_sample_matrix[, is_tumor, drop = FALSE] > 0) else 0
  )

  # Determine if this sublineage is in the user's focus list (figure 6) or not (figure s3/detection)
  in_focus <- ct %in% cell_types_focus
  # Get meta config for this sublineage if present, else use first row as fallback (for non-focus sublineages)
  if (in_focus) {
    meta_row <- plot_config[plot_config$sublineage == ct, ]
  } else {
    meta_row <- plot_config[1, ]
  }

  # Output directory and filename logic
  if (in_focus) {
    outdir <- "output/6. plots/figure 6/"
    # Assign 6a, 6b, 6c, ... for each focus cell type
    prefix <- paste0("fig_6", letters[which(cell_types_focus == ct)])
  } else {
    outdir <- "output/6. plots/figure s3/detection/"
    prefix <- paste0("fig_s3detection_", gsub("[ /]", "_", tolower(ct)))
  }
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # 1. Unfiltered plot (all genes)
  p_all <- ggplot(detection_df, aes(x = detected_in_normal, y = detected_in_tumor)) +
    geom_bin2d(binwidth = c(1, 1), color = NA) +
    {
      bin_freqs <- as.numeric(table(detection_df$detected_in_normal, detection_df$detected_in_tumor))
      bin_freqs <- bin_freqs[bin_freqs > 0]
      if (length(bin_freqs) > 1) {
        max_bin <- max(bin_freqs)
        bin_freqs_wo_max <- bin_freqs[bin_freqs < max_bin]
        if (length(bin_freqs_wo_max) > 0) {
          quants <- quantile(bin_freqs_wo_max, probs = c(0, 0.25, 0.5, 0.75))
          breaks <- unique(c(round(quants), max_bin))
        } else {
          breaks <- unique(bin_freqs)
        }
      } else {
        breaks <- unique(bin_freqs)
      }
      scale_fill_viridis_c(
        option = "C",
        trans = "log1p",
        name = "# genes",
        guide = guide_colorbar(
          barwidth = 1,
          barheight = 20,
          ticks = TRUE,
          draw.ulim = TRUE,
          draw.llim = TRUE,
          frame.colour = "black",
          frame.linewidth = 0.5,
          label.theme = element_text(size = 10)
        ),
        breaks = breaks,
        labels = breaks
      )
    } +
    scale_x_continuous(
      limits = c(0, cp10k_n_normal),
      breaks = seq(0, cp10k_n_normal, by = 5),
      minor_breaks = seq(0, cp10k_n_normal, by = 1)
    ) +
    scale_y_continuous(
      limits = c(0, cp10k_n_tumor),
      breaks = seq(0, cp10k_n_tumor, by = 5),
      minor_breaks = seq(0, cp10k_n_tumor, by = 1)
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    coord_fixed() +
    theme_bw(base_size = 14) +
    labs(
      x = paste0("# normal samples (", cp10k_n_normal, ")"),
      y = paste0("# tumor samples (", cp10k_n_tumor, ")"),
      subtitle = paste0("Sublineage: ", ct),
      fill = "# genes"
    )
  outname_all <- file.path(outdir, paste0(prefix, "_gene_detection_bin2d_", gsub("[ /]", "_", tolower(ct)), "_all.png"))
  ggsave(outname_all, p_all, width = 6, height = 6, dpi = 300)

  # 2. Filtered plot (important genes from meta_scores)
  meta_scores_file <- file.path(
    "output/3. intersector",
    meta_row$normalization_method,
    meta_row$meta_celltype_set,
    meta_row$gene_type,
    "meta_scores.csv"
  )
  if (!is.na(meta_scores_file) && file.exists(meta_scores_file)) {
    meta_scores <- suppressMessages(readr::read_csv(meta_scores_file, show_col_types = FALSE))
    important_genes <- meta_scores %>% filter(sublineage == ct) %>% pull(gene) %>% unique()
    detection_df_filt <- detection_df %>% filter(gene %in% important_genes)
    if (nrow(detection_df_filt) > 0) {
      p_imp <- ggplot(detection_df_filt, aes(x = detected_in_normal, y = detected_in_tumor)) +
        geom_bin2d(binwidth = c(1, 1), color = NA) +
        {
          bin_freqs <- as.numeric(table(detection_df_filt$detected_in_normal, detection_df_filt$detected_in_tumor))
          bin_freqs <- bin_freqs[bin_freqs > 0]
          if (length(bin_freqs) > 1) {
            max_bin <- max(bin_freqs)
            bin_freqs_wo_max <- bin_freqs[bin_freqs < max_bin]
            if (length(bin_freqs_wo_max) > 0) {
              quants <- quantile(bin_freqs_wo_max, probs = c(0, 0.25, 0.5, 0.75))
              breaks <- unique(c(round(quants), max_bin))
            } else {
              breaks <- unique(bin_freqs)
            }
          } else {
            breaks <- unique(bin_freqs)
          }
          scale_fill_viridis_c(
            option = "C",
            trans = "log1p",
            name = "# genes",
            guide = guide_colorbar(
              barwidth = 1,
              barheight = 20,
              ticks = TRUE,
              draw.ulim = TRUE,
              draw.llim = TRUE,
              frame.colour = "black",
              frame.linewidth = 0.5,
              label.theme = element_text(size = 10)
            ),
            breaks = breaks,
            labels = breaks
          )
        } +
        scale_x_continuous(
          limits = c(0, cp10k_n_normal),
          breaks = seq(0, cp10k_n_normal, by = 5),
          minor_breaks = seq(0, cp10k_n_normal, by = 1)
        ) +
        scale_y_continuous(
          limits = c(0, cp10k_n_tumor),
          breaks = seq(0, cp10k_n_tumor, by = 5),
          minor_breaks = seq(0, cp10k_n_tumor, by = 1)
        ) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        coord_fixed() +
        theme_bw(base_size = 14) +
        labs(
          x = paste0("# normal samples (", cp10k_n_normal, ")"),
          y = paste0("# tumor samples (", cp10k_n_tumor, ")"),
          subtitle = paste0("Sublineage: ", ct),
          fill = "# genes"
        )
      outname_imp <- file.path(outdir, paste0(prefix, "_gene_detection_bin2d_", gsub("[ /]", "_", tolower(ct)), "_important.png"))
      ggsave(outname_imp, p_imp, width = 6, height = 6, dpi = 300)
    } else {
      message(sprintf("No important genes found for sublineage %s in %s", ct, meta_scores_file))
    }
  } else if (!is.na(meta_scores_file)) {
    message(sprintf("meta_scores.csv not found for sublineage %s at %s", ct, meta_scores_file))
  }
}

message("Completed writing figure 6bcde to file.")
