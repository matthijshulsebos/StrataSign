library(dplyr)
library(Matrix)
library(readr)

source("src/0. utils/feature_name_utils.R")

# Load UMI table and cell metadata
umitab_filtered <- readRDS("output/6. plots/data/cp10k/cp10k_normalized_umitab.rds")
cell_metadata_final <- read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>%
  mutate(sample_ID = as.character(sample_ID))

# Function to compute detection rates per feature
compute_detection_rates_per_feature <- function(umitab_filtered, cell_metadata_final, table_s1, annots_list) {
  cell_metadata_final <- cell_metadata_final %>%
    mutate(sample_ID = as.character(sample_ID)) %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
    left_join(annots_list %>% select(cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
    mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage))

  # Get unique cluster IDs from cell metadata
  all_cluster_ids <- unique(cell_metadata_final$cluster_ID)
  detection_list <- list()

  # Iterate over each cluster
  for (current_cluster_id in all_cluster_ids) {
    # Filter metadata for the current cluster
    meta_cluster <- cell_metadata_final %>% filter(cluster_ID == current_cluster_id)
    if (nrow(meta_cluster) == 0) next

    # Get the current cell type name
    current_cell_type <- meta_cluster$cell_type[1]

    # Get cells in the cluster
    cell_ids_in_cluster <- meta_cluster$cell_ID
    umi_mat_cluster <- umitab_filtered[, colnames(umitab_filtered) %in% cell_ids_in_cluster, drop = FALSE]
    if (ncol(umi_mat_cluster) == 0) next

    # Map cell IDs to sample IDs
    cell_to_sample_map <- setNames(as.character(meta_cluster$sample_ID), meta_cluster$cell_ID)
    sample_vector <- cell_to_sample_map[colnames(umi_mat_cluster)]
    unique_samples_in_cluster <- unique(sample_vector)

    # Aggregate UMI counts by sample
    gene_sample_matrix <- sapply(unique_samples_in_cluster, function(sid) {
      cells_in_sample <- names(sample_vector)[sample_vector == sid]
      rowSums(umi_mat_cluster[, cells_in_sample, drop = FALSE])
    })

    # If there is a single sample then gene sample matrix is a vector and then convert it to a matrix to be sure
    if (is.vector(gene_sample_matrix)) {
      gene_sample_matrix <- matrix(gene_sample_matrix, ncol = 1, dimnames = list(rownames(umi_mat_cluster), unique_samples_in_cluster))
    }

    # Convert to sparse matrix
    gene_sample_matrix <- as(gene_sample_matrix, "sparseMatrix")
    
    # Get tissue information for the samples
    tissue_info <- table_s1 %>%
      filter(sample_ID %in% colnames(gene_sample_matrix)) %>%
      select(sample_ID, tissue)

    # Get detection rates for normal and tumor samples
    normal_samples <- tissue_info$sample_ID[tissue_info$tissue == "Normal"]
    tumor_samples <- tissue_info$sample_ID[tissue_info$tissue == "Tumor"]

    # Calculate detection rates
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
    
    # Create feature identifiers
    cluster_id_numeric <- as.numeric(tail(strsplit(as.character(current_cluster_id), "_")[[1]], 1))
    feature_id_vec <- create_feature_identifier(
      gene = rownames(gene_sample_matrix),
      cluster_name_part = rep(current_cell_type, nrow(gene_sample_matrix)),
      cluster_id_numeric = rep(cluster_id_numeric, nrow(gene_sample_matrix))
    )

    # Create a data frame for the current cluster detection rates
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

  # Row bind all the detection data frames
  final_detection_df <- do.call(rbind, detection_list)
  rownames(final_detection_df) <- NULL

  return(final_detection_df)
}

# Compute detection rates per feature
detection_rates_df <- compute_detection_rates_per_feature(umitab_filtered, cell_metadata_final, table_s1, annots_list)

# Save the detection rates to file
outdir <- "output/6. plots/data/detection"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
write_csv(detection_rates_df, file.path(outdir, "detection_rates.csv"))
message("Finished computing detection rates.")
