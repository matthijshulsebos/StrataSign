###############################################################################
# preprocessing_pipeline.R — Fully annotated
#
# What this script does (high level):
# 1) Loads core data objects (UMI matrix, per-cell metadata, sample metadata,
#    cluster annotations, KEGG metabolic genes) via data_loader.R.
# 2) Defines subsets of clusters (cell-type sets) and gene-type sets to create
#    multiple training datasets (metabolic / nonmetabolic / random).
# 3) Implements three normalization “approaches”:
#    - Global: rescale each (sample × cluster) to a single global-average total.
#    - Relative: CP10K per cell → log1p → aggregate to sample×cluster → rescale
#                each cell type to the single global-average total.
#    - Read depth: treat raw counts as normalized (baseline), log1p, and save.
# 4) Produces wide “training matrices” (features = gene_cluster columns) and
#    train/test splits, writes metadata, and calls fold-change helpers.
#
###############################################################################

# --- Libraries ---------------------------------------------------------------
# Tidy manipulation, reading/writing CSVs, and sparse matrix handling
library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(tibble)

# --- Project helpers ---------------------------------------------------------
# Data loading (UMIs, metadata, annotations, KEGG genes)
source("src/1. data preprocessing/training datasets/data_loader.R")
# Fold-change utilities (used per normalization/cell/gene set)
source("src/4. fold changes/fold_changes.R")


# --- Configuration: subsets & labels ----------------------------------------
# Clusters suspected to be doublets (to flag/exclude earlier in pipeline)
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Which gene subsets to generate datasets for
GENE_TYPE_SETS <- c("metabolic", "nonmetabolic", "random")

# Cluster subsets for various experimental conditions/plots
CELL_TYPE_SETS <- list(
  lcam_hi   = c(44, 9, 17, 28, 46, 11, 42),
  lcam_lo   = c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16),
  lcam_both = c(c(44, 9, 17, 28, 46, 11, 42),
                c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16)),
  macrophages = c(5, 8, 10, 11, 25, 32, 33, 35, 38, 42, 47, 54, 55, 57),
  all_clusters = NULL  # NULL means "no filtering" later (include all clusters)
)

# --- Dataset loading & preparation -------------------------------------------
# Pull in all objects the rest of the pipeline depends on. The loader returns a
# list; we unpack it and (importantly) select representative samples (Table S1).
datasets   <- load_all_datasets()
lung_ldm   <- datasets$lung_ldm
table_s1   <- select_representative_samples(datasets$table_s1)
annots_list <- datasets$annots_list
hsa01100_genes <- datasets$hsa01100_genesa01100_genes <- datasets$hsa01100_genes

# Select the raw UMI matrix variant from the `lung_ldm` object
raw_umitab <- lung_ldm$dataset$umitab

# Build comprehensive per-cell metadata and remove doublets (uses DOUBLETS)
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)

# Filter the UMI matrix to the cells that remain after metadata filtering
umitab_filtered <- filter_umitab(raw_umitab, cell_metadata)

# Keep only metadata rows for cells that are in the filtered UMI matrix
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# Free large objects to reduce memory pressure (common with scRNA-seq scale data)
# (We deliberately keep `umitab_filtered` and `cell_metadata_final` for downstream work)
# Clear memory to avoid issues with large datasets
rm(lung_ldm)
rm(raw_umitab)
rm(cell_metadata)
gc()


# Transform sparse UMI matrix to long format aggregated at cell type level
# --- Helper: aggregate_umitab_to_long ----------------------------------------
# Purpose:
#   Convert a sparse genes×cells matrix (UMI counts) into a long table with rows
#   (sample_ID, cluster_ID, gene, count), summing across all cells belonging to
#   the same (sample_ID × cluster_ID). 
# Also appends a column `gene_cluster` = "<gene>_<cluster_ID>" for wide reshaping.
#     Steps:
#     1) summary(umitab) → triplet nonzero entries (gene_idx, cell_idx, count).
#     2) Map indices to row/column names to get 'gene' and 'cell_ID'.
#     3) Keep only needed columns (gene, cell_ID, count).
#     4) inner_join() with per-cell metadata to add sample_ID and cluster_ID.
#     5) group_by(sample_ID, cluster_ID, gene) and sum counts.
#     6) Coerce key columns to appropriate types for consistent joins later.
aggregate_umitab_to_long <- function(umitab, cell_metadata) {
  message("Converting sparse UMI matrix to aggregated long format.")
  
  # Convert sparseMatrix to triplet form (row index, col index, value) for nonzeros
  # This avoids expanding zeros and keeps memory/time manageable.
  umitab_triplet <- summary(umitab)
  colnames(umitab_triplet) <- c("gene_idx", "cell_idx", "count")
  
  # Map indices back to gene IDs and cell barcodes
  umitab_triplet$gene    <- rownames(umitab)[umitab_triplet$gene_idx]
  umitab_triplet$cell_ID <- colnames(umitab)[umitab_triplet$cell_idx]
  
  # Use tibble for clearer dplyr pipelines
  umitab_triplet <- as_tibble(umitab_triplet)
  
  # Join with cell metadata and aggregate counts
  aggregated_data <- umitab_triplet %>%
    # Keep only the columns needed for joining/aggregation; drop index cols (gene_idx, cell_idx).
    select(gene, cell_ID, count) %>%
    # Join with cell metadata to get sample and cluster IDs for the groupby
    # Exact match join: only cells present in both the counts and metadata are retained.
    inner_join(cell_metadata %>% select(cell_ID, sample_ID, cluster_ID), by = "cell_ID") %>%
    # Create transcriptional profiles per sample, cluster, and gene
    # Aggregate counts within each (sample × cluster × gene).
    group_by(sample_ID, cluster_ID, gene) %>%
    # Sum counts across all cells belonging to the same group.
    summarise(count = sum(count), .groups = 'drop') %>%
    # Make sure they are characters for joining later
    mutate(
      sample_ID = as.character(sample_ID),
      cluster_ID = as.numeric(cluster_ID),
      gene = as.character(gene)
    )
  
  return(aggregated_data)
}


# Scales total counts per cell type to the global average
# --- Helper: apply_celltype_normalization ------------------------------------
# Purpose:
#   Rescale each (sample × cluster) gene profile to a single *global* average
#   magnitude so that totals are comparable across samples and clusters.
#   Implementation:
#     - Compute per-(sample, cluster) totals.
#     - Take global_avg_celltype_total = mean of those totals across ALL groups.
#     - Within each group, counts → fractions, multiply by global average.
#   Output: long table with normalized_count values per (sample, cluster, gene).
apply_celltype_normalization <- function(long_data) {
  message("Applying cell type normalization.")
  

  # Calculate global average cell type total
  celltype_totals <- long_data %>%
    group_by(sample_ID, cluster_ID) %>%
    # Sum counts across all cells belonging to the same group.
    summarise(celltype_total = sum(count), .groups = 'drop')
  
  global_avg_celltype_total <- mean(celltype_totals$celltype_total)
  
  # Scale cell type counts to the global average
  normalized_data <- long_data %>%
    left_join(celltype_totals, by = c("sample_ID", "cluster_ID")) %>%
    mutate(
      gene_fraction = ifelse(celltype_total > 0, count / celltype_total, 0),
      normalized_count = gene_fraction * global_avg_celltype_total
    ) %>%
    select(sample_ID, gene, cluster_ID, normalized_count)
  
  return(normalized_data)
}


# Takes a list of genes and spits out the filtered genes based on the gene type
# --- Helper: filter_genes_by_type --------------------------------------------
# Purpose:
#   Select a gene set of a given type, keeping the same size as the metabolic
#   set present in the current data. This allows fair comparisons across gene
#   types when training models or computing fold changes.
filter_genes_by_type <- function(all_genes, gene_type) {
  # Filter on metabolic genes that occur in the data
  metabolic_genes <- hsa01100_genes$Symbol[hsa01100_genes$Symbol %in% all_genes]
  
  # Sample genes depending on type parameter and make sure they are the same gene set size
  if (gene_type == "metabolic") {
    return(metabolic_genes)
  } else if (gene_type == "nonmetabolic") {
    nonmetabolic_genes <- setdiff(all_genes, metabolic_genes)
    # Reproducibility: fix RNG seed so splits/samples are stable across runs.
    set.seed(42)
    return(sample(nonmetabolic_genes, length(metabolic_genes)))
  } else if (gene_type == "random") {
    # Reproducibility: fix RNG seed so splits/samples are stable across runs.
    set.seed(43)
    return(sample(all_genes, length(metabolic_genes)))
  }
}


# Creates training datasets and saves them but also calculates fold changes with sourced function
# --- Helper: save_results -----------------------------------------------------
# Purpose:
#   Create and persist datasets for modeling:
#     - pivot longer→wider (features = gene_cluster columns),
#     - join patient/tissue metadata,
#     - split by patient into 70/30 train/test,
#     - save X_train/X_test, y_train/y_test, and metadata CSVs.
#   Notes:
#     - Splitting by patient avoids leakage across samples from the same patient.
#     - Output directory encodes method/celltypes/genetype for traceability.
save_results <- function(long_data, method_name, cell_type, gene_type) {
  # Convert from long format with only non-zero values to wide format with zeros
  data_wide <- long_data %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    # Long → wide: each gene_cluster becomes a feature column; missing values filled with 0.
    pivot_wider(names_from = gene_cluster, 
                values_from = log_normalized_count, 
                values_fill = 0)
  
  # Join with metadata
  data_with_meta <- data_wide %>%
    left_join(
      table_s1 %>%
        select(sample_ID, patient_ID, tissue) %>%
        mutate(sample_ID = as.character(sample_ID)),
      by = "sample_ID"
    ) %>%
    filter(!is.na(tissue))
  
  # Save results if there are any rows
  if (nrow(data_with_meta) > 0) {
    # The suffix is added to the file name
    suffix <- switch(gene_type, "metabolic" = "_metabolic", "nonmetabolic" = "_nonmetabolic", "random" = "_random")
    
    # Create output directory for the subset
    output_dir <- file.path("output", "1. data preprocessing", "training datasets", method_name, cell_type, gene_type)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Remove non transcriptional data for the complete version of the dataset
    X_complete <- data_with_meta %>% select(-sample_ID, -patient_ID, -tissue)
    y_complete <- data_with_meta$tissue
    
    # Filter table s1 which contains sample metadata on the training data
    metadata_complete <- table_s1 %>% 
      mutate(sample_ID = as.character(sample_ID)) %>%
      filter(sample_ID %in% data_with_meta$sample_ID)
    
    # Write the complete datasets to files
    write_csv(X_complete, file.path(output_dir, paste0("X_complete_", cell_type, suffix, ".csv")))
    write_csv(data.frame(y = y_complete), file.path(output_dir, paste0("y_complete_", cell_type, suffix, ".csv")))
    write_csv(metadata_complete, file.path(output_dir, paste0("metadata_complete_", cell_type, suffix, ".csv")))
    
    # Seed for reproducibility of patient split
    # Reproducibility: fix RNG seed so splits/samples are stable across runs.
    set.seed(42)

    # Get list of all patients and set which are used for training
    patients <- unique(data_with_meta$patient_ID)
    train_patients <- sample(patients, size = round(length(patients) * 0.7))
    
    # Create the train and test datasets
    train_data <- data_with_meta %>% filter(patient_ID %in% train_patients)
    test_data <- data_with_meta %>% filter(!patient_ID %in% train_patients)

    # If the split is successful
    if (nrow(train_data) > 0 && nrow(test_data) > 0) {
      # Make sure to remove non transcriptional data
      X_train <- train_data %>% select(-sample_ID, -patient_ID, -tissue)
      X_test <- test_data %>% select(-sample_ID, -patient_ID, -tissue)
      
      # Set the target variables
      y_train <- train_data$tissue
      y_test <- test_data$tissue
      
      # Filter train metadata on actual train samples
      metadata_train <- table_s1 %>% 
        mutate(sample_ID = as.character(sample_ID)) %>%
        filter(sample_ID %in% train_data$sample_ID)
      
      # Filter test metadata on actual test samples
      metadata_test <- table_s1 %>% 
        mutate(sample_ID = as.character(sample_ID)) %>%
        filter(sample_ID %in% test_data$sample_ID)
      
      # Write all the files to the output directory
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


# Scales counts per cell to the median counts per cell to compensate for gene type filtering
# --- Helper: apply_cp_median_normalization -----------------------------------
# Purpose:
#   Normalize counts per cell by that cell's median count (counts-per-median),
#   then aggregate normalized values to (sample, cluster, gene).
# Rationale:
#   Median-based scaling reduces sensitivity to highly expressed genes and
#   per-cell depth differences vs total-sum scaling.
apply_cp_median_normalization <- function(umitab_subset, cell_metadata_subset) {
  # Each column is a cell and rows are genes so colsums are total counts per cell
  total_counts_per_cell <- colSums(umitab_subset)
  
  # Calculate median counts per cell
  median_counts <- median(total_counts_per_cell[total_counts_per_cell > 0])

  # Set zero values to 1 to avoid division by zero
  total_counts_per_cell <- ifelse(total_counts_per_cell == 0, 1, total_counts_per_cell)
  
  # Calculate expression proportions and scale by median count
  cp_median_cell_level_matrix <- sweep(umitab_subset, 2, total_counts_per_cell, FUN = "/") * median_counts
  
  # Get cells that dont contain expression values
  cols_to_zero_indices <- which(colSums(umitab_subset) == 0)
  
  # Set cells with zero counts to zero for safe measure
  if (length(cols_to_zero_indices) > 0) {
    if (ncol(cp_median_cell_level_matrix) > 0) {
      cp_median_cell_level_matrix[, cols_to_zero_indices] <- 0
    }
  }
  
  # Set invalid values to zero
  cp_median_cell_level_matrix[is.na(cp_median_cell_level_matrix) | is.infinite(cp_median_cell_level_matrix)] <- 0
  
  # Convert to triplet format which is more efficient and removes zero values
  cp_median_triplet <- summary(cp_median_cell_level_matrix)
  colnames(cp_median_triplet) <- c("gene_idx", "cell_idx", "normalized_value") 
  
  # Convert indices to genes and cell IDs
  cp_median_triplet_df <- data.frame(
    gene = rownames(cp_median_cell_level_matrix)[cp_median_triplet$gene_idx],
    cell_ID = colnames(cp_median_cell_level_matrix)[cp_median_triplet$cell_idx],
    normalized_value = cp_median_triplet$normalized_value
  )

  # Add cell type metadata to the cell data
  cp_median_long_df <- cp_median_triplet_df %>%
    # Exact match join: only cells present in both the counts and metadata are retained.
    inner_join(cell_metadata_subset %>% select(cell_ID, sample_ID, cluster_ID), by = "cell_ID")
  
  # Aggregate the data to get the total counts per gene per cell type per sample
  aggregated_cp_median_data <- cp_median_long_df %>%
    # Aggregate counts within each (sample × cluster × gene).
    group_by(sample_ID, cluster_ID, gene) %>%
    # Sum counts across all cells belonging to the same group.
    summarise(count = sum(normalized_value), .groups = 'drop') %>%
    mutate(
      sample_ID = as.character(sample_ID),
      cluster_ID = as.numeric(cluster_ID),
      gene = as.character(gene)
    )
  
  return(aggregated_cp_median_data)
}


# Scales individual cells to total sum of 10K counts
# --- Helper: apply_cp10k_normalization ---------------------------------------
# Purpose:
#   Scale each cell so its total counts sum to 10,000 (CP10K), drop empty cells,
#   and log1p-transform to stabilize variance. Returns a gene×cell matrix.
apply_cp10k_normalization <- function(umitab_matrix) {
  message("Applying CP10K normalization to raw umitab data.")
  
  # Each column is a cell so colsums is the total counts per cell
  total_counts_per_cell <- colSums(umitab_matrix)
  
  # Scale to 10K counts in total
  scaling_factor <- 10000
  
  # Set zero values to 1 to avoid division by zero
  total_counts_per_cell[total_counts_per_cell == 0] <- 1

  # Matrix multiplication with a diagonal matrix for scaling to 10K
  scaling_diag <- Diagonal(x = scaling_factor / total_counts_per_cell)
  cp10k_matrix <- umitab_matrix %*% scaling_diag
  
  # Get cells that dont contain expression values
  cols_to_zero <- which(colSums(umitab_matrix) == 0)
  if (length(cols_to_zero) > 0) {
    cp10k_matrix[, cols_to_zero] <- 0
  }
  
  # Set invalid values to zero
  cp10k_matrix[is.na(cp10k_matrix) | is.infinite(cp10k_matrix)] <- 0
  
  # Set row and column names to match the original matrix
  rownames(cp10k_matrix) <- rownames(umitab_matrix)
  colnames(cp10k_matrix) <- colnames(umitab_matrix)
  
  return(cp10k_matrix)
}


# Saves cell type total counts after normalization for figure 4
# --- Helper: save_celltype_total_counts --------------------------------------
# Purpose:
#   Aggregate normalized counts to totals per (sample, cluster), join tissue and
#   cluster annotations, and save to disk for QC/plotting.
save_celltype_total_counts <- function(normalized_data, normalization_method, cell_type_set, gene_type) {
  # Aggregate total counts per cell type and sample
  total_counts <- normalized_data %>%
    group_by(sample_ID, cluster_ID) %>%
    # Sum counts across all cells belonging to the same group.
    summarise(total_count = sum(normalized_count), .groups = 'drop') %>%
    left_join(table_s1 %>% select(sample_ID, tissue) %>% mutate(sample_ID = as.character(sample_ID)), by = "sample_ID") %>%
    left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
    mutate(
      normalization_method = normalization_method,
      cell_type_set = cell_type_set,
      gene_type = gene_type,
      cell_type_label = paste0(sub_lineage, "_", cluster_ID)
    ) %>%
    select(sample_ID, cluster_ID, cell_type_label, tissue, normalization_method, 
           cell_type_set, gene_type, total_count)
  
  # Set and create the output directory
  output_dir <- file.path("output", "6. plots", "data", "cell_type_norm_counts", gene_type, cell_type_set)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create the filename based on normalization strategy
  filename <- paste0(normalization_method, "_total_counts.csv")
  output_path <- file.path(output_dir, filename)
  
  # If the file already exists then read it and combine with new data because it contains data from multiple normalization methods
  if (file.exists(output_path)) {
    existing_data <- read_csv(output_path, show_col_types = FALSE) %>%
      mutate(sample_ID = as.character(sample_ID))
    combined_data <- bind_rows(existing_data, total_counts)
    write_csv(combined_data, output_path)
  } else {
    # If the file doesnt exist then create it
    write_csv(total_counts, output_path)
  }
}


# Applies the global normalization strategy
# === Driver: Global normalization ============================================
# Flow:
#   - Apply single-global-mean rescaling to all (sample×cluster) profiles.
#   - For each cell-type set and gene set:
#       • save totals, • compute fold changes, • log1p & build gene_cluster,
#       • save training matrices and splits.
run_global_approach <- function(aggregated_data) {
  message("\n GLOBAL NORMALIZATION APPROACH")
  
  # This is the directory where the results will be saved
  norm_method_name <- "ctnorm_global"
  
  # Global normalizes the entire dataset first and then does filtering
  global_ctnorm_data <- apply_celltype_normalization(aggregated_data)
  
  # Filter on cell types and gene types
  for (cell_type in names(CELL_TYPE_SETS)) {
    cluster_ids <- CELL_TYPE_SETS[[cell_type]]
    
    # Filter on cell types
    if (!is.null(cluster_ids)) {
      subset_data <- global_ctnorm_data %>% filter(cluster_ID %in% cluster_ids)
    } else {
      # Null indicates all cell types as we dont specify a subset
      subset_data <- global_ctnorm_data
    }
    
    for (gene_type in GENE_TYPE_SETS) {
      # Filters on metabolic, nonmetabolic or random genes
      all_genes_in_data <- unique(subset_data$gene)
      filtered_genes <- filter_genes_by_type(all_genes_in_data, gene_type)
      gene_filtered_data <- subset_data %>% filter(gene %in% filtered_genes)
      
      # Save cell type total counts after normalization for plotting
      save_celltype_total_counts(gene_filtered_data, norm_method_name, cell_type, gene_type)
      
      # Calculate fold changes for this dataset
      calculate_fold_changes_for_normalization(gene_filtered_data, table_s1, annots_list, norm_method_name, cell_type, gene_type)
      
      # Log1p normalization
      log_data <- gene_filtered_data %>% mutate(log_normalized_count = log1p(normalized_count))
      
      # Create feature identifiers
      final_data <- create_features(log_data, annots_list)
      
      # Create the training splits and save the results
      save_results(final_data, norm_method_name, cell_type, gene_type)
    }
  }
}


# Applies the relative normalization strategy
# === Driver: Relative normalization (CP10K-based) ============================
# Flow:
#   - Subset cells by cluster set, apply CP10K+log1p per cell,
#   - Aggregate to (sample, cluster, gene),
#   - Rescale via apply_celltype_normalization to single global average,
#   - Select gene sets, save totals & fold changes, and persist matrices.
run_relative_approach <- function(umitab_filtered, cell_metadata_final) {
  message("\n RELATIVE NORMALIZATION APPROACH")
  
  # Iterate over cell types and gene types
  for (cell_type in names(CELL_TYPE_SETS)) {
    cluster_ids <- CELL_TYPE_SETS[[cell_type]]
    
    for (gene_type in GENE_TYPE_SETS) {
      message(sprintf("Processing %s cell type and %s genes for relative normalization.", cell_type, gene_type))
      
      # Filter cell metadata based on cluster IDs
      subset_cell_metadata <- if (!is.null(cluster_ids)) {
        cell_metadata_final %>% filter(cluster_ID %in% cluster_ids)
      } else {
        cell_metadata_final
      }
      
      # List of genes we want to retain
      all_genes_in_umitab <- rownames(umitab_filtered)
      filtered_genes <- filter_genes_by_type(all_genes_in_umitab, gene_type)
      
      # List of cells we want to retain
      cells_to_keep_in_subset <- intersect(subset_cell_metadata$cell_ID, colnames(umitab_filtered))
      
      # Omit all data outside of what we want to retain
      subset_umitab <- umitab_filtered[filtered_genes, cells_to_keep_in_subset, drop = FALSE]
      
      # Scale all cell counts to the median total counts per cell
      cp_median_data <- apply_cp_median_normalization(subset_umitab, subset_cell_metadata)
      
      # Apply cell type normalization to the median normalized data
      relative_ctnorm_data <- apply_celltype_normalization(cp_median_data)
      
      # Save cell type total counts after normalization for plotting
      save_celltype_total_counts(relative_ctnorm_data, "ctnorm_relative", cell_type, gene_type)
      
      # Calculate fold changes for this dataset
      calculate_fold_changes_for_normalization(relative_ctnorm_data, table_s1, annots_list, "ctnorm_relative", cell_type, gene_type)
      
      # Log1p normalization
      log_data <- relative_ctnorm_data %>% mutate(log_normalized_count = log1p(normalized_count))
      
      # Create feature identifiers
      final_data <- create_features(log_data, annots_list)
      
      # Create the training splits and save the results
      save_results(final_data, "ctnorm_relative", cell_type, gene_type)
    }
  }
}

# Runs the read depth normalization approach as baseline
# === Driver: Read-depth baseline =============================================
# Flow:
#   - Use aggregated raw counts (no rescaling) as normalized_count,
#   - Select gene sets, save totals & fold changes, log1p, and persist matrices.
run_read_depth_approach <- function(aggregated_data) {
  message("\n READ DEPTH NORMALIZATION")
  
  # Iterate over cell types and gene types
  for (cell_type in names(CELL_TYPE_SETS)) {
    cluster_ids <- CELL_TYPE_SETS[[cell_type]]
    
    # Filter on cell types
    if (!is.null(cluster_ids)) {
      subset_data <- aggregated_data %>% filter(cluster_ID %in% cluster_ids)
    } else {
      # Again null indicates we dont filter on cell types so it indicates all clusters
      subset_data <- aggregated_data
    }
    
    for (gene_type in GENE_TYPE_SETS) {
      # Filters on metabolic, nonmetabolic or random genes
      all_genes_in_data <- unique(subset_data$gene)
      filtered_genes <- filter_genes_by_type(all_genes_in_data, gene_type)
      gene_filtered_data <- subset_data %>% filter(gene %in% filtered_genes)
      
      # Rename the count column to normalized_count for consistency (as we dont actually cell type normalize)
      read_depth_data <- gene_filtered_data %>% 
        rename(normalized_count = count)
      
      # Save cell type total counts after normalization for plotting
      save_celltype_total_counts(read_depth_data, "read_depth", cell_type, gene_type)
      
      # Calculate fold changes for this dataset
      calculate_fold_changes_for_normalization(read_depth_data, table_s1, annots_list, "read_depth", cell_type, gene_type)
      
      # Log1p normalization
      log_data <- read_depth_data %>% mutate(log_normalized_count = log1p(normalized_count))
      
      # Create feature identifiers
      final_data <- create_features(log_data, annots_list)
      
      # Create the training splits and save the results
      save_results(final_data, "read_depth", cell_type, gene_type)
    }
  }
}




# Apply and save CP10K normalization to the filtered umitab data
cp10k_normalized_umitab <- apply_cp10k_normalization(umitab_filtered)

# Save the CP10K normalized data and cell metadata
output_dir_cp10k <- "output/6. plots/data/cp10k"
dir.create(output_dir_cp10k, recursive = TRUE, showWarnings = FALSE)
saveRDS(cp10k_normalized_umitab, file.path(output_dir_cp10k, "cp10k_normalized_umitab.rds"))
write_csv(cell_metadata_final, file.path(output_dir_cp10k, "cell_metadata_final.csv"))

# Aggregated version of CP10K normalized data
aggregated_data <- aggregate_umitab_to_long(cp10k_normalized_umitab, cell_metadata_final)

# Run normalization strategies
run_global_approach(aggregated_data)
run_relative_approach(cp10k_normalized_umitab, cell_metadata_final)
run_read_depth_approach(aggregated_data)

message("\n Preprocessing pipeline completed.")
