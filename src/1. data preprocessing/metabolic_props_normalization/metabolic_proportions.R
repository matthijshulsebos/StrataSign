calculate_metabolic_proportions <- function(ds_matrix, cell_metadata, hsa01100_genes) {
  message("Calculating metabolic proportions for raw, global, and relative normalization stages.")

  metabolic_genes <- intersect(hsa01100_genes$SYMBOL, rownames(ds_matrix))
  results <- list()
  sample_ids <- unique(cell_metadata$sample_ID)

  # Load annotation and sample tables for cell type and tissue lookup
  annots_list <- readr::read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
  table_s1 <- readr::read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE)

  # Raw (aggregate to long format)
  message("Aggregating raw counts to long format.")
  raw_long <- aggregate_umitab_to_long(ds_matrix, cell_metadata)

  # Global normalization
  message("Applying global cell type normalization.")
  global_norm_long <- apply_celltype_normalization(raw_long)

  # Relative normalization
  message("Filter on metabolic genes.")
  ds_metabolic <- ds_matrix[metabolic_genes, , drop=FALSE]
  message("Applying CP-median normalization.")
  rel_cpmed_long <- apply_cp_median_normalization(ds_metabolic, cell_metadata)
  message("Applying cell type normalization.")
  rel_norm_long <- apply_celltype_normalization(rel_cpmed_long)

  total_steps <- length(sample_ids)
  step <- 0
  for (sample in sample_ids) {
    step <- step + 1
    message(sprintf("[Progress: %d/%d] Processing sample_ID: %s", step, total_steps, sample))
    sample_clusters <- unique(cell_metadata$cluster_ID[cell_metadata$sample_ID == sample])
    for (cluster in sample_clusters) {
      # Read depth
      raw_sub <- raw_long %>% filter(sample_ID == sample, cluster_ID == cluster)
      total_raw <- sum(raw_sub$count)
      metabolic_raw <- sum(raw_sub$count[raw_sub$gene %in% metabolic_genes])
      prop_raw <- ifelse(total_raw > 0, metabolic_raw / total_raw, NA)

      # Global normalization
      global_sub <- global_norm_long %>% filter(sample_ID == sample, cluster_ID == cluster)
      total_global <- sum(global_sub$normalized_count)
      metabolic_global <- sum(global_sub$normalized_count[global_sub$gene %in% metabolic_genes])
      prop_global <- ifelse(total_global > 0, metabolic_global / total_global, NA)

      # Relative normalization
      rel_sub <- rel_norm_long %>% filter(sample_ID == sample, cluster_ID == cluster)
      total_relative <- sum(rel_sub$normalized_count)
      metabolic_relative <- total_relative
      prop_relative <- ifelse(total_relative > 0, metabolic_relative / total_relative, NA)

      # Lookup cell type and tissue
      cell_type_row <- annots_list[annots_list$cluster == cluster, ]
      if (nrow(cell_type_row) == 0) {
        cell_type <- NA
      } else if (!is.na(cell_type_row$sub_lineage) && cell_type_row$sub_lineage != "") {
        cell_type <- cell_type_row$sub_lineage
      } else {
        cell_type <- cell_type_row$lineage
      }
      tissue <- table_s1$tissue[table_s1$sample_ID == sample]
      if (length(tissue) == 0) tissue <- NA
      # Combine cell type and cluster ID
      celltype_cluster <- if (!is.na(cell_type) && !is.na(cluster)) paste(cell_type, cluster, sep = "_") else NA

      results[[length(results) + 1]] <- data.frame(
        sample_ID = sample,
        cluster_ID = cluster,
        cell_type = cell_type,
        celltype_cluster = celltype_cluster,
        tissue = tissue,
        metabolic_prop_raw = prop_raw,
        metabolic_prop_global = prop_global,
        metabolic_prop_relative = prop_relative,
        metabolic_counts_raw = metabolic_raw,
        total_counts_raw = total_raw,
        metabolic_counts_global = metabolic_global,
        total_counts_global = total_global,
        metabolic_counts_relative = metabolic_relative,
        total_counts_relative = total_relative
      )
    }
  }
  out <- do.call(rbind, results)
  output_path <- file.path("output", "0. intermediates", "metabolic_proportions_by_normalization.csv")
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(out, output_path)
  message("Metabolic proportions written to: ", output_path)
  return(out)
}
