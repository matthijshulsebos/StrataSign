library(dplyr)
library(tidyr)


# Process a single sample-cluster pair for CP10K normalization
process_sample_cluster_pair <- function(sample_id, cluster_id, cell_metadata_filtered_df, umitab_filtered, all_relevant_genes) {
  # Select all relevant cells using cell metadata
  cells_in_pair <- cell_metadata_filtered_df %>%
    filter(sample_ID == sample_id, cluster_ID == cluster_id) %>%
    pull(cell_ID)
  
  # Create empty data frame for results
  result_df <- data.frame(
    sample_ID = sample_id,
    cluster_ID = cluster_id,
    gene = all_relevant_genes,
    normalized_count = 0,
    stringsAsFactors = FALSE
  )
  
  if (length(cells_in_pair) > 0) {
    # Filter umitab for relevant cells
    sub_umitab <- umitab_filtered[, cells_in_pair, drop = FALSE]

    # Total counts for each cell
    total_counts_sub_cells <- colSums(sub_umitab)

    # If total counts of cell is 0 then set to 1 to avoid division by 0
    total_counts_sub_cells[total_counts_sub_cells == 0] <- 1
    
    # Divide each gene count by the total counts of the cell and multiply by 10K
    sub_umitab_cp10k <- sweep(sub_umitab, 2, total_counts_sub_cells, FUN = "/") * 10000

    # If any values are empty or infinite set them to 0
    sub_umitab_cp10k[is.na(sub_umitab_cp10k) | is.infinite(sub_umitab_cp10k)] <- 0
    
    # Aggregate counts for each gene across all cells
    gene_sums_cp10k <- rowSums(sub_umitab_cp10k)
    result_df$normalized_count <- gene_sums_cp10k
    
    # Clear intermediate variables to free memory
    rm(sub_umitab, sub_umitab_cp10k)
  }
  
  return(result_df)
}


# CP10K normalization
normalize_cp10k <- function(cell_metadata_filtered_df, umitab_filtered, annots_list, DOUBLETS) {
  message("Starting CP10K normalization...")

  # Keep track of all samples, genes, and clusters
  all_relevant_samples <- unique(cell_metadata_filtered_df$sample_ID)
  all_relevant_genes <- rownames(umitab_filtered)
  all_relevant_clusters <- unique(cell_metadata_filtered_df$cluster_ID)
  
  # Create all combinations of sample and cluster IDs
  sample_cluster_pairs <- expand_grid(
    sample_ID = all_relevant_samples,
    cluster_ID = all_relevant_clusters
  )
  
  # Set batch size
  batch_size <- 500
  total_batches <- ceiling(nrow(sample_cluster_pairs) / batch_size)
  
  # Aggregated results list
  list_of_aggregated_dfs <- list()
  
  # Loop over batches
  for (batch_idx in 1:total_batches) {
    # Set start and end pairs for current batch
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, nrow(sample_cluster_pairs))
    
    message(sprintf("Processing batch %d/%d (pairs %d-%d)", 
                    batch_idx, total_batches, start_idx, end_idx))
    
    # Select current batch combinations
    batch_pairs <- sample_cluster_pairs[start_idx:end_idx, ]
    
    batch_results <- lapply(1:nrow(batch_pairs), function(pair_idx) {
      process_sample_cluster_pair(
        sample_id = batch_pairs$sample_ID[pair_idx],
        cluster_id = batch_pairs$cluster_ID[pair_idx],
        cell_metadata_filtered_df = cell_metadata_filtered_df,
        umitab_filtered = umitab_filtered,
        all_relevant_genes = all_relevant_genes
      )
    })
    
    # Combine separate pairs from batch processing
    batch_combined <- dplyr::bind_rows(batch_results)

    # Combines batch result with other batches
    list_of_aggregated_dfs[[batch_idx]] <- batch_combined
    
    # Clear batch data and force garbage collection
    rm(batch_results, batch_combined)
    gc()
    
  }
  
  message("Combining all batches...")
  aggregated_counts_long <- dplyr::bind_rows(list_of_aggregated_dfs)
  
  # Clear memory
  rm(list_of_aggregated_dfs)
  gc()
  
  message("CP10K normalization complete.")
  return(aggregated_counts_long)
}


# Combination function for CP10K and cell type proportion normalization
normalize_cp10k_celltype <- function(cell_metadata_filtered_df, umitab_filtered, annots_list, DOUBLETS) {
  message("Performing CP10K with cell type proportion normalization.")
  
  # Perform CP10K normalization
  aggregated_counts_long <- normalize_cp10k(cell_metadata_filtered_df, umitab_filtered, annots_list, DOUBLETS)
  
  # Calculate total counts per cell type and sample
  celltype_totals <- aggregated_counts_long %>%
    group_by(sample_ID, cluster_ID) %>%
    summarise(celltype_total = sum(normalized_count), .groups = 'drop')
  
  # Calculate global average cell type total
  global_avg_celltype_total <- celltype_totals %>%
    summarise(mean_total = mean(celltype_total)) %>%
    pull(mean_total)
  
  # Calculate gene fraction and multiply by global average
  normalized_counts <- aggregated_counts_long %>%
    left_join(celltype_totals, by = c("sample_ID", "cluster_ID")) %>%
    mutate(
      gene_fraction = ifelse(celltype_total > 0, normalized_count / celltype_total, 0),
      normalized_count = gene_fraction * global_avg_celltype_total
    ) %>%
    select(sample_ID, gene, cluster_ID, normalized_count)
  
  message("CP10K with cell type proportion normalization complete.")
  return(normalized_counts)
}


# Helper function to apply cell type normalization
apply_celltype_normalization <- function(normalized_counts) {
  # Calculate total counts per cell type and sample
  celltype_totals <- normalized_counts %>%
    group_by(sample_ID, cluster_ID) %>%
    summarise(celltype_total = sum(normalized_count), .groups = 'drop')
  
  # Calculate global average cell type total
  global_avg_celltype_total <- celltype_totals %>%
    summarise(mean_total = mean(celltype_total)) %>%
    pull(mean_total)
  
  # Apply cell type proportion normalization
  counts_celltype_normalized <- normalized_counts %>%
    left_join(celltype_totals, by = c("sample_ID", "cluster_ID")) %>%
    mutate(
      gene_fraction = ifelse(celltype_total > 0, normalized_count / celltype_total, 0),
      normalized_count = gene_fraction * global_avg_celltype_total
    ) %>%
    select(sample_ID, gene, cluster_ID, normalized_count)
  
  return(counts_celltype_normalized)
}


# Helper function to apply proportion cell type normalization
apply_proportion_celltype_normalization <- function(counts_long) {
  message("Starting proportion normalization with cell type normalization...")
  
  # Calculate total counts per cell type and sample
  total_counts_celltype_sample <- counts_long %>%
    group_by(sample_ID, cluster_ID) %>%
    summarise(total_counts_in_celltype = sum(count), .groups = 'drop')
  
  # Calculate global average total counts across all cell types
  global_avg_celltype_total <- total_counts_celltype_sample %>%
    summarise(mean_total = mean(total_counts_in_celltype)) %>%
    pull(mean_total)
  
  # Apply proportion with cell type normalization
  counts_normalized <- counts_long %>%
    left_join(total_counts_celltype_sample, by = c("sample_ID", "cluster_ID")) %>%
    mutate(
      gene_fraction_in_celltype = ifelse(total_counts_in_celltype > 0, count / total_counts_in_celltype, 0),
      normalized_count = gene_fraction_in_celltype * global_avg_celltype_total
    ) %>%
    select(sample_ID, gene, cluster_ID, normalized_count)
  
  message("Proportion normalization with cell type normalization complete.")
  return(counts_normalized)
}


# Sample depth normalization function
normalize_sample_depth <- function(counts_long) {
  message("Starting sample read depth normalization...")

  # Calculate total counts per sample
  total_counts_per_sample <- counts_long %>%
    group_by(sample_ID) %>%
    summarise(total_counts = sum(count), .groups = 'drop')
  
  # Calculate target count as mean of total counts across samples
  target_count <- mean(total_counts_per_sample$total_counts)
  message("Target count for sample depth normalization: ", target_count)
  
  # Normalize counts to be proportional to mean sample count
  counts_normalized <- counts_long %>%
    left_join(total_counts_per_sample, by = "sample_ID") %>%
    mutate(
      scaling_factor = target_count / total_counts,
      normalized_count = count * scaling_factor
    ) %>%
    select(sample_ID, gene, cluster_ID, normalized_count)
  
  message("Sample read depth normalization complete.")
  return(counts_normalized)
}
