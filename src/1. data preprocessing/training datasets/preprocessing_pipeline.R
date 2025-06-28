library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(tibble)

# Load utilities
source("src/1. data preprocessing/training datasets/data_loader.R")
source("src/4. fold changes/fold_changes.R")

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
raw_umitab <- lung_ldm$dataset$umitab
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)
umitab_filtered <- filter_umitab(raw_umitab, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# Remove large object from memory after extraction
rm(lung_ldm)
gc()

# Remove raw_umitab if not needed later
rm(raw_umitab)
gc()

# Remove cell_metadata if not needed later
rm(cell_metadata)
gc()


# Transform sparse UMI matrix to long format aggregated at cell type level
aggregate_umitab_to_long <- function(umitab, cell_metadata) {
  message("Converting sparse UMI matrix to aggregated long format.")
  
  # Convert sparse matrix to triplet format
  umitab_triplet <- summary(umitab)
  colnames(umitab_triplet) <- c("gene_idx", "cell_idx", "count")
  
  # Map indices to names
  umitab_triplet$gene <- rownames(umitab)[umitab_triplet$gene_idx]
  umitab_triplet$cell_ID <- colnames(umitab)[umitab_triplet$cell_idx]
  
  # Convert to tibble for memory efficiency
  umitab_triplet <- as_tibble(umitab_triplet)
  
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
  metabolic_genes <- hsa01100_genes$Symbol[hsa01100_genes$Symbol %in% all_genes]
  
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



# CP10K normalization function for raw umitab data (sparse matrix aware)
apply_cp10k_normalization <- function(umitab_matrix) {
  message("Applying CP10K normalization to raw umitab data.")
  
  # Calculate total counts per cell
  total_counts_per_cell <- Matrix::colSums(umitab_matrix)
  
  # Set scaling factor to 10,000 (CP10K)
  scaling_factor <- 10000
  
  # Avoid division by zero
  total_counts_per_cell[total_counts_per_cell == 0] <- 1
  
  # Create a diagonal matrix for scaling
  scaling_diag <- Matrix::Diagonal(x = scaling_factor / total_counts_per_cell)
  
  # Multiply umitab_matrix by the diagonal scaling matrix
  cp10k_matrix <- umitab_matrix %*% scaling_diag
  
  # Set cells with originally zero total counts back to zero
  cols_to_zero <- which(Matrix::colSums(umitab_matrix) == 0)
  if (length(cols_to_zero) > 0) {
    cp10k_matrix[, cols_to_zero] <- 0
  }
  
  # Handle any remaining invalid values
  cp10k_matrix[is.na(cp10k_matrix) | is.infinite(cp10k_matrix)] <- 0
  
  # Ensure row and column names are preserved
  rownames(cp10k_matrix) <- rownames(umitab_matrix)
  colnames(cp10k_matrix) <- colnames(umitab_matrix)
  
  return(cp10k_matrix)
}


# Calculate fold changes for normalized data



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
      
      # Calculate fold changes before log transformation
      calculate_fold_changes_for_normalization(gene_filtered_data, table_s1, annots_list, "ctnorm_global", cell_type, gene_type)
      
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
  global_sum_path <- file.path("output", "0. intermediates", "sum_counts_per_celltype_per_sample_global.csv")
  dir.create(dirname(global_sum_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(global_sum, global_sum_path)
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
      
      # Calculate fold changes before log transformation
      calculate_fold_changes_for_normalization(relative_ctnorm_data, table_s1, annots_list, "ctnorm_relative", cell_type, gene_type)
      
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
  rel_sum_path <- file.path("output", "0. intermediates", "sum_counts_per_celltype_per_sample_relative.csv")
  dir.create(dirname(rel_sum_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(rel_sum, rel_sum_path)
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
      
      # Calculate fold changes before log transformation
      calculate_fold_changes_for_normalization(read_depth_data, table_s1, annots_list, "read_depth", cell_type, gene_type)
        
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
  read_depth_sum_path <- file.path("output", "0. intermediates", "sum_counts_per_celltype_per_sample_read_depth.csv")
  dir.create(dirname(read_depth_sum_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(read_depth_sum, read_depth_sum_path)
}


# Apply CP10K normalization to the filtered umitab data for global and read depth approaches
cp10k_normalized_umitab <- apply_cp10k_normalization(umitab_filtered)

# Convert CP10K normalized data to aggregated format
aggregated_data <- aggregate_umitab_to_long(cp10k_normalized_umitab, cell_metadata_final)

# Run global approach
run_global_approach(aggregated_data)

# Run relative approach
run_relative_approach(cp10k_normalized_umitab, cell_metadata_final)

# Run read depth normalization approach
run_read_depth_approach(aggregated_data)

message("\n Preprocessing pipeline completed.")


# THIS IS OPTIONAL FOR FIGURE 1CDE
# source("src/1. data preprocessing/metabolic_props_normalization/metabolic_proportions.R")

# Calculate metabolic proportions
# calculate_metabolic_proportions(ds_matrix, cell_metadata, hsa01100_genes)
