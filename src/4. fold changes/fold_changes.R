library(dplyr)
library(tidyverse)
library(readr)
library(Matrix)

# Source utility functions
source("src/0. utils/feature_name_utils.R")
source("src/1. data preprocessing/training datasets/data_loader.R")

# Create output directory
dir.create("output/4. fold changes", recursive = TRUE, showWarnings = FALSE)

# Load Leader et al. data and prepare cell metadata
if (!exists("lung_ldm")) {
  load("base/data/lung_ldm.rd")
}

# Load metadata
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv")
annots_list <- read_csv("base/input_tables/annots_list.csv")

# Define doublet clusters to exclude
clusters_to_exclude <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Extract downsampled UMI matrix and prepare cell metadata
ds_matrix <- lung_ldm$dataset$ds[[1]]
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, clusters_to_exclude)
umitab_filtered <- filter_umitab(ds_matrix, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))


# === FUNCTION DEFINITIONS ===

# Aggregate single-cell data to cell type level
aggregate_to_celltype_level <- function(umitab, cell_metadata) {
  message("Converting sparse UMI matrix to aggregated long format.")
  
  # Ensure input is sparse
  if (!inherits(umitab, "Matrix")) {
    umitab <- as(umitab, "sparseMatrix")
  }
  
  # Convert sparse matrix to triplet format
  umitab_triplet <- summary(umitab)
  colnames(umitab_triplet) <- c("gene_idx", "cell_idx", "count")
  
  # Map indices to names
  umitab_triplet$gene <- rownames(umitab)[umitab_triplet$gene_idx]
  umitab_triplet$cell_ID <- colnames(umitab)[umitab_triplet$cell_idx]
  
  # Join with metadata and aggregate
  message("Aggregating counts by sample-cluster-gene.")
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


# Apply cell type normalization
apply_celltype_normalization <- function(data) {
  message("Applying cell type normalization.")
  
  celltype_totals <- data %>%
    group_by(sample_ID, cluster_ID) %>%
    summarise(celltype_total = sum(count), .groups = 'drop')
  
  global_avg_celltype_total <- mean(celltype_totals$celltype_total)
  
  normalized_data <- data %>%
    left_join(celltype_totals, by = c("sample_ID", "cluster_ID")) %>%
    mutate(
      gene_fraction = ifelse(celltype_total > 0, count / celltype_total, 0),
      normalized_count = gene_fraction * global_avg_celltype_total
    ) %>%
    select(sample_ID, gene, cluster_ID, normalized_count)
  
  message("Cell type normalization complete.")
  return(normalized_data)
}


# Function to calculate fold changes for normalized cell type data
calculate_fold_changes <- function(normalized_data, sample_metadata) {
  message("Calculating fold changes.")
  
  # Add condition information
  data_with_condition <- normalized_data %>%
    filter(sample_ID %in% rownames(sample_metadata)) %>%
    mutate(condition = sample_metadata[sample_ID, "condition"])
  
  # Calculate mean expression by condition and cluster
  condition_means <- data_with_condition %>%
    group_by(gene, cluster_ID, condition) %>%
    summarize(
      mean_expr = mean(normalized_count),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      names_from = condition, 
      values_from = c(mean_expr, n_samples)
    )
  
  # Calculate log2 fold changes
  result_df <- condition_means %>%
    mutate(
      mean_expr_normal = get(paste0("mean_expr_", ref_level)),
      mean_expr_tumor = get(paste0("mean_expr_", comparison_level)),
      log2FoldChange = log2((mean_expr_tumor + 0.1) / (mean_expr_normal + 0.1)),
      n_normal = get(paste0("n_samples_", ref_level)),
      n_tumor = get(paste0("n_samples_", comparison_level))
    ) %>%
    filter(!is.na(log2FoldChange) & is.finite(log2FoldChange))
  
  return(result_df)
}


# === MAIN EXECUTION ===

# Aggregate and normalize the data
message("Aggregating single-cell data to cell type level.")
aggregated_data <- aggregate_to_celltype_level(umitab_filtered, cell_metadata_final)
normalized_data <- apply_celltype_normalization(aggregated_data)

# Set reference and comparison levels
ref_level <- "Normal" 
comparison_level <- "Tumor"

# Preprocess sample metadata
sample_metadata <- table_s1 %>%
  select(sample_ID, tissue) %>%
  filter(sample_ID %in% unique(normalized_data$sample_ID)) %>%
  mutate(
    condition = tissue,
    condition = factor(condition, levels = c(ref_level, comparison_level))
  ) %>%
  column_to_rownames("sample_ID")

# Calculate fold changes
fold_change_results <- calculate_fold_changes(normalized_data, sample_metadata)

# Create feature identifiers using utility functions
fold_changes_feature_value <- fold_change_results %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    Feature = create_feature_identifier(gene, cluster_name, cluster_ID)
  ) %>%
  select(Feature, Value = log2FoldChange)

# Write results
write_csv(fold_changes_feature_value, "output/4. fold changes/feature_fold_changes.csv")

# Write a more detailed results file with additional information
detailed_results <- fold_change_results %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    Feature = create_feature_identifier(gene, cluster_name, cluster_ID)
  )

write_csv(detailed_results, "output/4. fold changes/detailed_feature_fold_changes.csv")

# Print summary statistics
message("Completed fold change calculation.")
