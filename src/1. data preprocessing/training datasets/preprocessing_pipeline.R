library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(tibble)

# Load utilities
source("src/1. data preprocessing/training datasets/data_loader.R")

# Configuration
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
GENE_TYPE_SETS <- c("metabolic", "nonmetabolic", "random")
CELL_TYPE_SETS <- list(
  lcam_hi = c(44, 9, 17, 28, 46, 11, 42),
  lcam_lo = c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16),
  lcam_both = c(c(44, 9, 17, 28, 46, 11, 42), c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16)),
  macrophages = c(5, 8, 10, 11, 25, 32, 33, 35, 38, 42, 47, 54, 55, 57),
  has_10_in_50pct = c(1, 10, 11, 13, 14, 17, 18, 19, 2, 20, 23, 25, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 42, 45, 46, 47, 48, 49, 5, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 7, 8, 9),
  has_10_in_70pct = c(1, 10, 14, 18, 2, 25, 29, 30, 32, 33, 34, 35, 36, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59),
  all_clusters = NULL
)

# Load datasets
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
table_s1 <- select_representative_samples(datasets$table_s1)
annots_list <- datasets$annots_list
hsa01100_genes <- datasets$hsa01100_genes

# Datasets from Leader et al. 2021
ds_matrix <- lung_ldm$dataset$ds[[1]]
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)
umitab_filtered <- filter_umitab(ds_matrix, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))


# Transform sparse UMI matrix to long format aggregated at cell type level
aggregate_umitab_to_long <- function(umitab, cell_metadata) {
  message("Converting sparse UMI matrix to aggregated long format.")
  
  # Convert sparse matrix to triplet format
  umitab_triplet <- summary(umitab)
  colnames(umitab_triplet) <- c("gene_idx", "cell_idx", "count")
  
  # Map indices to names
  umitab_triplet$gene <- rownames(umitab)[umitab_triplet$gene_idx]
  umitab_triplet$cell_ID <- colnames(umitab)[umitab_triplet$cell_idx]
  
  # Join with metadata and aggregate
  message("Aggregating counts by sample, cell type, gene.")
  aggregated_data <- umitab_triplet %>%
    select(gene, cell_ID, count) %>%
    inner_join(cell_metadata %>% select(cell_ID, sample_ID, cluster_ID), by = "cell_ID") %>%
    group_by(sample_ID, cluster_ID, gene) %>%
    summarise(count = sum(count), .groups = 'drop') %>%
    mutate(
      sample_ID = as.character(sample_ID),
      cluster_ID = as.numeric(cluster_ID),
      gene = as.character(gene)
    )
  
  return(aggregated_data)
}


# Make sure all cell types for each sample have equal total counts
apply_celltype_normalization <- function(long_data) {
  message("Applying cell type normalization.")

  # Calculate total counts for each cell type in each sample
  celltype_totals <- long_data %>%
    group_by(sample_ID, cluster_ID) %>%
    summarise(celltype_total = sum(count), .groups = 'drop')

  # Calculate global average cell type total
  global_avg_celltype_total <- mean(celltype_totals$celltype_total)

  # Normalize counts by cell type totals and global average
  normalized_data <- long_data %>%
    left_join(celltype_totals, by = c("sample_ID", "cluster_ID")) %>%
    mutate(
      gene_fraction = ifelse(celltype_total > 0, count / celltype_total, 0),
      normalized_count = gene_fraction * global_avg_celltype_total
    ) %>%
    select(sample_ID, gene, cluster_ID, normalized_count)

  return(normalized_data)
}


# Filter genes by type
filter_genes_by_type <- function(all_genes, gene_type) {
  metabolic_genes <- hsa01100_genes$SYMBOL[hsa01100_genes$SYMBOL %in% all_genes]
  
  if (gene_type == "metabolic") {
    return(metabolic_genes)
  } else if (gene_type == "nonmetabolic") {
    nonmetabolic_genes <- setdiff(all_genes, metabolic_genes)
    set.seed(42)
    return(sample(nonmetabolic_genes, length(metabolic_genes)))
  } else if (gene_type == "random") {
    set.seed(43)
    return(sample(all_genes, length(metabolic_genes)))
  }
}


# Create training datasets and write to file
save_results <- function(long_data, method_name, cell_type, gene_type) {
  # Pivot from long to wide format
  data_wide <- long_data %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster, 
                values_from = log_normalized_count, 
                values_fill = 0)
  
  # Join with metadata
  data_with_meta <- data_wide %>%
    left_join(table_s1 %>% select(sample_ID, patient_ID, tissue), by = "sample_ID") %>%
    filter(!is.na(tissue))
  
  if (nrow(data_with_meta) > 0) {
    # Construct output path
    suffix <- switch(gene_type, "metabolic" = "_metabolic", "nonmetabolic" = "_nonmetabolic", "random" = "_random")
    output_dir <- file.path("output", "1. data preprocessing", "training datasets", method_name, cell_type, gene_type)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save complete dataset
    X_complete <- data_with_meta %>% select(-sample_ID, -patient_ID, -tissue)
    y_complete <- data_with_meta$tissue
    metadata_complete <- table_s1 %>% filter(sample_ID %in% data_with_meta$sample_ID)
    
    write_csv(X_complete, file.path(output_dir, paste0("X_complete_", cell_type, suffix, ".csv")))
    write_csv(data.frame(y = y_complete), file.path(output_dir, paste0("y_complete_", cell_type, suffix, ".csv")))
    write_csv(metadata_complete, file.path(output_dir, paste0("metadata_complete_", cell_type, suffix, ".csv")))
    
    # Create train/test split
    set.seed(42)
    patients <- unique(data_with_meta$patient_ID)
    train_patients <- sample(patients, size = round(length(patients) * 0.7))
    
    train_data <- data_with_meta %>% filter(patient_ID %in% train_patients)
    test_data <- data_with_meta %>% filter(!patient_ID %in% train_patients)
    
    if (nrow(train_data) > 0 && nrow(test_data) > 0) {
      # For training input leave out target variables
      X_train <- train_data %>% select(-sample_ID, -patient_ID, -tissue)
      X_test <- test_data %>% select(-sample_ID, -patient_ID, -tissue)

      # Target variable is the tissue type
      y_train <- train_data$tissue
      y_test <- test_data$tissue
      
      metadata_train <- table_s1 %>% filter(sample_ID %in% train_data$sample_ID)
      metadata_test <- table_s1 %>% filter(sample_ID %in% test_data$sample_ID)
      
      # Write datasets to csv files
      write_csv(X_train, file.path(output_dir, paste0("X_train_", cell_type, suffix, ".csv")))
      write_csv(X_test, file.path(output_dir, paste0("X_test_", cell_type, suffix, ".csv")))
      write_csv(data.frame(y = y_train), file.path(output_dir, paste0("y_train_", cell_type, suffix, ".csv")))
      write_csv(data.frame(y = y_test), file.path(output_dir, paste0("y_test_", cell_type, suffix, ".csv")))
      write_csv(metadata_train, file.path(output_dir, paste0("metadata_train_", cell_type, suffix, ".csv")))
      write_csv(metadata_test, file.path(output_dir, paste0("metadata_test_", cell_type, suffix, ".csv")))
      
      message(sprintf("Saved results for %s normalization, %s cell types, %s genes.", method_name, cell_type, gene_type))
    }
  }
}


# CP-median normalization and aggregation function
apply_cp_median_normalization <- function(umitab_subset, cell_metadata_subset) {
  # Calculate total counts per cell for the subset
  total_counts_per_cell <- colSums(umitab_subset)

  # Calculate median total counts across all cells (excluding zeros)
  median_counts <- median(total_counts_per_cell[total_counts_per_cell > 0])
  
  # Avoid division by zero
  total_counts_per_cell <- ifelse(total_counts_per_cell == 0, 1, total_counts_per_cell)

  # Per column applies division of umitabs by total counts per cell, scaled to median
  cp_median_cell_level_matrix <- sweep(umitab_subset, 2, total_counts_per_cell, FUN = "/") * median_counts

  # Ensure cells with zero total counts have zero normalized counts
  cols_to_zero_indices <- which(colSums(umitab_subset) == 0)
  if (length(cols_to_zero_indices) > 0) {
    # Ensure that cp_median_cell_level_matrix has columns before attempting to subset
    if (ncol(cp_median_cell_level_matrix) > 0) {
      cp_median_cell_level_matrix[, cols_to_zero_indices] <- 0
    }
  }
  
  # Set invalid values to zero
  cp_median_cell_level_matrix[is.na(cp_median_cell_level_matrix) | is.infinite(cp_median_cell_level_matrix)] <- 0
  
  # Convert sparse matrix to triplet format
  cp_median_triplet <- summary(cp_median_cell_level_matrix)
  
  # Rename columns for clarity
  colnames(cp_median_triplet) <- c("gene_idx", "cell_idx", "normalized_value") 

  cp_median_triplet_df <- data.frame(
    gene = rownames(cp_median_cell_level_matrix)[cp_median_triplet$gene_idx],
    cell_ID = colnames(cp_median_cell_level_matrix)[cp_median_triplet$cell_idx],
    normalized_value = cp_median_triplet$normalized_value
  )

  # Join with cell_metadata_subset to get sample_ID and cluster_ID
  cp_median_long_df <- cp_median_triplet_df %>%
    inner_join(cell_metadata_subset %>% select(cell_ID, sample_ID, cluster_ID), by = "cell_ID")

  # Aggregate the cell-level CP-median normalized counts
  aggregated_cp_median_data <- cp_median_long_df %>%
    group_by(sample_ID, cluster_ID, gene) %>%
    summarise(count = sum(normalized_value), .groups = 'drop') %>%
    mutate(
      sample_ID = as.character(sample_ID),
      cluster_ID = as.numeric(cluster_ID),
      gene = as.character(gene)
    )
  
  return(aggregated_cp_median_data)
}


# Global normalization approach
run_global_approach <- function(aggregated_data) {
  message("\n GLOBAL NORMALIZATION APPROACH")
  
  # Apply cell type normalization on the entire dataset
  global_ctnorm_data <- apply_celltype_normalization(aggregated_data)

  # Create each ablation subset for cluster and gene sets
  for (cell_type in names(CELL_TYPE_SETS)) {
    cluster_ids <- CELL_TYPE_SETS[[cell_type]]
    
    # Filter by clusters
    if (!is.null(cluster_ids)) {
      subset_data <- global_ctnorm_data %>% filter(cluster_ID %in% cluster_ids)
    } else {
      subset_data <- global_ctnorm_data
    }
    
    for (gene_type in GENE_TYPE_SETS) {
      # Retrieve the filtered gene set
      all_genes_in_data <- unique(subset_data$gene)
      filtered_genes <- filter_genes_by_type(all_genes_in_data, gene_type)

      # Filter the subset_data based on the filtered gene set
      gene_filtered_data <- subset_data %>% filter(gene %in% filtered_genes)
      
      # Log1p transformation
      log_data <- gene_filtered_data %>% mutate(log_normalized_count = log1p(normalized_count))
      
      # Add feature names
      final_data <- create_features(log_data, annots_list)
      
      # Save results
      save_results(final_data, "ctnorm_global", cell_type, gene_type)
    }
  }
  
  # Write sum per cell type per sample
  global_sum <- global_ctnorm_data %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
    group_by(tissue, sample_ID, cluster_ID) %>%
    summarise(total_count = sum(normalized_count), .groups = 'drop')
  write_csv(global_sum, file.path("output", "0. intermediates", "sum_counts_per_celltype_per_sample_global.csv"))
}


# Relative normalization approach
run_relative_approach <- function(umitab_filtered, cell_metadata_final) {
  message("\n RELATIVE NORMALIZATION APPROACH")

  # Iterate over all ablation subsets
  for (cell_type in names(CELL_TYPE_SETS)) {
    cluster_ids <- CELL_TYPE_SETS[[cell_type]]
    
    for (gene_type in GENE_TYPE_SETS) {
      message(sprintf("Processing %s cell type and %s genes for relative normalization.", cell_type, gene_type))
      
      # Filter cell metadata by clusters
      subset_cell_metadata <- if (!is.null(cluster_ids)) {
        cell_metadata_final %>% filter(cluster_ID %in% cluster_ids)
      } else {
        # For all_clusters
        cell_metadata_final
      }
      
      # Filter genes by type using helper function
      all_genes_in_umitab <- rownames(umitab_filtered)
      filtered_genes <- filter_genes_by_type(all_genes_in_umitab, gene_type)

      # Ensure cells are present in both cell metadata and umitab
      cells_to_keep_in_subset <- intersect(subset_cell_metadata$cell_ID, colnames(umitab_filtered))

      # Apply the subset filter
      subset_umitab <- umitab_filtered[filtered_genes, cells_to_keep_in_subset, drop = FALSE]
      
      # Apply Custom CP10K normalization 
      cp_median_data <- apply_cp_median_normalization(subset_umitab, subset_cell_metadata)

      # Apply cell type normalization
      relative_ctnorm_data <- apply_celltype_normalization(cp_median_data)
      
      # Log1p transformation
      log_data <- relative_ctnorm_data %>% mutate(log_normalized_count = log1p(normalized_count))
      
      # Add feature names
      final_data <- create_features(log_data, annots_list)
      
      # Save results
      save_results(final_data, "ctnorm_relative", cell_type, gene_type)
    }
  }
  
  # Write sum per cell type per sample
  rel_sum <- relative_ctnorm_data %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
    group_by(tissue, sample_ID, cluster_ID) %>%
    summarise(total_count = sum(normalized_count), .groups = 'drop')
  write_csv(rel_sum, file.path("output", "0. intermediates", "sum_counts_per_celltype_per_sample_relative.csv"))
}


# Read depth normalization approach
run_read_depth_approach <- function(aggregated_data) {
  message("\n READ DEPTH NORMALIZATION")
  for (cell_type in names(CELL_TYPE_SETS)) {
    cluster_ids <- CELL_TYPE_SETS[[cell_type]]
    
    # Filter by cell types first
    if (!is.null(cluster_ids)) {
      subset_data <- aggregated_data %>% filter(cluster_ID %in% cluster_ids)
    } else {
      subset_data <- aggregated_data
    }
    
    for (gene_type in GENE_TYPE_SETS) {
      # Retrieve the filtered gene set
      all_genes_in_data <- unique(subset_data$gene)
      filtered_genes <- filter_genes_by_type(all_genes_in_data, gene_type)

      # Filter the subset_data based on the filtered gene set
      gene_filtered_data <- subset_data %>% filter(gene %in% filtered_genes)
      
      # Rename count column to normalized_count for consistency with downstream functions
      read_depth_data <- gene_filtered_data %>% 
        rename(normalized_count = count)
        
      # Log1p transformation
      log_data <- read_depth_data %>% mutate(log_normalized_count = log1p(normalized_count))
        
      # Add feature names
      final_data <- create_features(log_data, annots_list)
        
      # Save results
      save_results(final_data, "read_depth", cell_type, gene_type)
    }
  }
  
  # Write sum per cell type per sample
  read_depth_sum <- aggregated_data %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
    group_by(tissue, sample_ID, cluster_ID) %>%
    summarise(total_count = sum(count), .groups = 'drop')
  write_csv(read_depth_sum, file.path("output", "0. intermediates", "sum_counts_per_celltype_per_sample_read_depth.csv"))
}


# Convert downsampled data to aggregated format for global and read_depth normalization
message("Converting downsampled data to aggregated format.")
aggregated_data <- aggregate_umitab_to_long(umitab_filtered, cell_metadata_final)

# Run global approach
run_global_approach(aggregated_data)

# Run relative approach
run_relative_approach(umitab_filtered, cell_metadata_final)

# Run read depth normalization approach
run_read_depth_approach(aggregated_data)

message("\n Preprocessing pipeline completed.")


# THIS IS OPTIONAL FOR FIGURE 1CDE

# Calculate metabolic gene proportions per sample and cell type at three normalization stages
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


# Calculate metabolic proportions efficiently
#calculate_metabolic_proportions(ds_matrix, cell_metadata, hsa01100_genes)
