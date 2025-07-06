# Calculate detection rates per feature and save to output/6. plots/data/detection

library(dplyr)
library(Matrix)
library(readr)

# --- Load data ---

umitab_filtered <- readRDS("output/6. plots/data/cp10k/cp10k_normalized_umitab.rds")
cell_metadata_final <- readr::read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
annots_list <- readr::read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
# Ensure sample_ID is consistently treated as a character string
table_s1 <- readr::read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>%
  mutate(sample_ID = as.character(sample_ID))

# --- Function to compute detection rates per feature ---

compute_detection_rates_per_feature <- function(umitab_filtered, cell_metadata_final, table_s1, annots_list) {
  # --- 1. Pre-join metadata ---
  # Add tissue and cell type information to the cell metadata table for easier access.
  cell_metadata_final <- cell_metadata_final %>%
    mutate(sample_ID = as.character(sample_ID)) %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
    left_join(annots_list %>% select(cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
    mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage))

  # Get all unique cluster IDs to iterate over
  all_cluster_ids <- unique(cell_metadata_final$cluster_ID)
  detection_list <- list()

  # --- 2. Iterate over each cluster to calculate detection rates ---
  for (current_cluster_id in all_cluster_ids) {
    # Filter metadata for the current cluster
    meta_cluster <- cell_metadata_final %>% filter(cluster_ID == current_cluster_id)
    if (nrow(meta_cluster) == 0) next

    # Get the sublineage/cell_type for this cluster
    current_cell_type <- meta_cluster$cell_type[1]

    # Get the columns (cells) from the UMI matrix that belong to this cluster
    cell_ids_in_cluster <- meta_cluster$cell_ID
    umi_mat_cluster <- umitab_filtered[, colnames(umitab_filtered) %in% cell_ids_in_cluster, drop = FALSE]
    if (ncol(umi_mat_cluster) == 0) next

    # --- 3. Aggregate cell counts into sample-level counts for this cluster ---
    # Create a mapping from cell ID to sample ID
    cell_to_sample_map <- setNames(as.character(meta_cluster$sample_ID), meta_cluster$cell_ID)
    sample_vector <- cell_to_sample_map[colnames(umi_mat_cluster)]
    unique_samples_in_cluster <- unique(sample_vector)

    # Aggregate UMI counts by sample
    # This creates a matrix where rows are genes and columns are samples,
    # and values are the sum of counts for all cells from that sample in this cluster.
    gene_sample_matrix <- sapply(unique_samples_in_cluster, function(sid) {
      cells_in_sample <- names(sample_vector)[sample_vector == sid]
      Matrix::rowSums(umi_mat_cluster[, cells_in_sample, drop = FALSE])
    })

    # Coerce the result to a matrix to handle the edge case of a single sample,
    # where sapply would otherwise return a vector.
    if (is.vector(gene_sample_matrix)) {
      gene_sample_matrix <- matrix(gene_sample_matrix, ncol = 1, dimnames = list(rownames(umi_mat_cluster), unique_samples_in_cluster))
    }

    # Ensure it's a sparse matrix for efficiency
    gene_sample_matrix <- as(gene_sample_matrix, "sparseMatrix")
    
    # --- 4. Calculate detection counts per gene for Normal and Tumor tissues ---
    # Get tissue information for the samples present in this cluster
    tissue_info <- table_s1 %>%
      filter(sample_ID %in% colnames(gene_sample_matrix)) %>%
      select(sample_ID, tissue)

    normal_samples <- tissue_info$sample_ID[tissue_info$tissue == "Normal"]
    tumor_samples <- tissue_info$sample_ID[tissue_info$tissue == "Tumor"]

    # Count how many normal/tumor samples have a non-zero count for each gene
    n_detected_normal <- if (length(normal_samples) > 0) {
      Matrix::rowSums(gene_sample_matrix[, normal_samples, drop = FALSE] > 0)
    } else {
      rep(0, nrow(gene_sample_matrix))
    }

    n_detected_tumor <- if (length(tumor_samples) > 0) {
      Matrix::rowSums(gene_sample_matrix[, tumor_samples, drop = FALSE] > 0)
    } else {
      rep(0, nrow(gene_sample_matrix))
    }

    # Calculate detection rates (fraction of samples with non-zero counts)
    det_rate_normal <- if (length(normal_samples) > 0) {
      n_detected_normal / length(normal_samples)
    } else {
      rep(NA, nrow(gene_sample_matrix))
    }

    det_rate_tumor <- if (length(tumor_samples) > 0) {
      n_detected_tumor / length(tumor_samples)
    } else {
      rep(NA, nrow(gene_sample_matrix))
    }


    # --- 5. Create feature identifiers and format output ---
    source("src/0. utils/feature_name_utils.R")
    # Ensure cluster_id is a character before splitting
    cluster_id_numeric <- suppressWarnings(as.numeric(tail(strsplit(as.character(current_cluster_id), "_")[[1]], 1)))
    feature_id_vec <- create_feature_identifier(
      gene = rownames(gene_sample_matrix),
      cluster_name_part = rep(current_cell_type, nrow(gene_sample_matrix)),
      cluster_id_numeric = rep(cluster_id_numeric, nrow(gene_sample_matrix))
    )

    # Store results for this cluster in a data frame
    detection_list[[current_cluster_id]] <- data.frame(
      gene = rownames(gene_sample_matrix),
      cell_type = current_cell_type,
      cluster_ID = current_cluster_id,
      feature_id = feature_id_vec,
      det_rate_normal = det_rate_normal,
      det_rate_tumor = det_rate_tumor,
      n_detected_normal = n_detected_normal,
      n_detected_tumor = n_detected_tumor,
      n_normal = length(normal_samples),
      n_tumor = length(tumor_samples),
      stringsAsFactors = FALSE
    )
  }

  # --- 6. Combine results from all clusters ---
  # This single data frame contains the detection rates for every feature (gene@cluster)
  final_detection_df <- do.call(rbind, detection_list)
  rownames(final_detection_df) <- NULL # Clean up row names
  return(final_detection_df)
}

# --- Run detection rate calculation ---
detection_rates_df <- compute_detection_rates_per_feature(umitab_filtered, cell_metadata_final, table_s1, annots_list)

# --- Save to output/6. plots/data/detection ---
outdir <- "output/6. plots/data/detection"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
write_csv(detection_rates_df, file.path(outdir, "detection_rates.csv"))
cat("Detection rates written to", file.path(outdir, "detection_rates.csv"), "
")
