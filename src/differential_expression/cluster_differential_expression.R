# Load required libraries
library(dplyr)
library(tidyverse)
library(readr)
library(Matrix)

# Create output directory
dir.create("output/differential_expression", recursive = TRUE, showWarnings = FALSE)

# Load Leader et al. data
if (!exists("lung_ldm")) {
  load("base/data/lung_ldm.rd")
}

# Load metadata
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv")
annots_list <- read_csv("base/input_tables/annots_list.csv")
cell_metadata <- read_csv("base/input_tables/cell_metadata.csv") # Added cell metadata

# Filter samples to those used in the clustering model 
table_s1 <- table_s1 %>% 
  filter(Use.in.Clustering.Model. == "Yes") %>%
  filter(!is.na(tissue)) %>%
  mutate(tissue = toupper(tissue))

# Define doublet clusters to exclude
clusters_to_exclude <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Prepare count data
counts <- lung_ldm$dataset$counts
 
# Filter out doublets
all_clusters <- as.numeric(dimnames(counts)[[3]])
clusters_to_keep <- all_clusters[!all_clusters %in% clusters_to_exclude]
counts <- counts[, , as.character(clusters_to_keep), drop = FALSE]

# Filter samples to match table_s1
sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
counts <- counts[dimnames(counts)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]

# Calculate cluster proportions
message("Calculating cell type proportions for normalization...")
cluster_proportions <- cell_metadata %>%
  filter(!cluster_ID %in% clusters_to_exclude) %>%
  mutate(
    sample_ID = as.character(sample_ID),
    cluster_ID = as.character(cluster_ID)
  ) %>%
  group_by(sample_ID, cluster_ID) %>%
  summarize(cluster_cell_count = n(), .groups = 'drop') %>%
  group_by(sample_ID) %>%
  mutate(total_cells = sum(cluster_cell_count),
         cluster_proportion = cluster_cell_count / total_cells) %>%
  ungroup()

# Set reference (normal) and comparison (tumor) levels
unique_tissues <- sort(unique(table_s1$tissue))
ref_level <- unique_tissues[1]
comparison_level <- unique_tissues[2]

# Create sample metadata
sample_metadata <- table_s1 %>%
  select(sample_ID, tissue) %>%
  filter(sample_ID %in% rownames(counts)) %>%
  mutate(
    condition = tissue,
    condition = factor(condition, levels = unique_tissues)
  ) %>%
  column_to_rownames("sample_ID")

# Function to process counts for a cluster-specific dataset
process_counts <- function(count_matrix, sample_metadata, dataset_name, cluster_id) {
  # Create a long-format data frame
  counts_long <- data.frame(
    sample_ID = rep(rownames(count_matrix), each = ncol(count_matrix)),
    gene = rep(colnames(count_matrix), times = nrow(count_matrix)),
    count = as.vector(count_matrix),
    stringsAsFactors = FALSE
  )
  
  # Add condition information
  counts_long$condition <- sample_metadata[counts_long$sample_ID, "condition"]
  
  # Filter genes with low counts
  gene_totals <- counts_long %>%
    group_by(gene) %>%
    summarize(total_count = sum(count), .groups = 'drop')
  
  genes_to_keep <- gene_totals %>%
    filter(total_count >= 10) %>%
    pull(gene)
  
  if (length(genes_to_keep) < 10) {
    message(paste("  Skipping", dataset_name, "- too few genes after filtering"))
    return(NULL)
  }
  
  counts_long <- counts_long %>% filter(gene %in% genes_to_keep)
  
  # Calculate mean sample size for normalization
  mean_sample_size <- counts_long %>%
    group_by(sample_ID) %>%
    summarize(sample_size = sum(count), .groups = 'drop') %>%
    summarize(mean_size = mean(sample_size)) %>%
    pull(mean_size)
  
  # Normalize counts by total sample counts
  normalized_counts <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(
      total_sample_counts = sum(count),
      sample_normalized_count = (count / total_sample_counts) * mean_sample_size,
      normalized_count = sample_normalized_count
    ) %>%
    ungroup()
  
  # Apply cell type proportion normalization for cluster-specific analysis
  cluster_id_str <- as.character(cluster_id)
  
  normalized_counts <- normalized_counts %>%
    # Add cluster_ID column for joining with cluster_proportions
    mutate(cluster_ID = cluster_id_str) %>%
    left_join(
      cluster_proportions %>% 
        filter(cluster_ID == cluster_id_str),
      by = c("sample_ID", "cluster_ID")
    ) %>%
    mutate(
      # Default to 1 if cluster_proportion is NA (no matching data)
      cluster_proportion = ifelse(is.na(cluster_proportion), 1, cluster_proportion),
      # Correct for over/under-representation by dividing by cluster proportion
      normalized_count = sample_normalized_count / cluster_proportion
    ) %>%
    # Re-normalize to maintain the same total counts per sample
    group_by(sample_ID) %>%
    mutate(
      total_after_norm = sum(normalized_count),
      normalized_count = normalized_count * (mean_sample_size / total_after_norm)
    ) %>%
    ungroup()
  
  # Calculate log-normalized counts
  normalized_counts <- normalized_counts %>%
    mutate(log_normalized_count = log1p(normalized_count))
  
  # Calculate mean expression by condition
  condition_means <- normalized_counts %>%
    group_by(gene, condition) %>%
    summarize(
      mean_log_expr = mean(log_normalized_count),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      names_from = condition, 
      values_from = mean_log_expr
    )
  
  # Calculate log2 fold changes
  result_df <- condition_means %>%
    mutate(
      log2FoldChange = log2((exp(.data[[comparison_level]]) + 1) / (exp(.data[[ref_level]]) + 1))
    )
  
  return(result_df)
}

# Process individual clusters
results_by_cluster <- list()
message("Processing clusters with normalized fold change approach...")

for (cluster_id in clusters_to_keep) {
  message(paste("Processing cluster", cluster_id))
  
  # Get counts for this cluster
  cluster_counts <- counts[, , as.character(cluster_id)]
  common_samples <- intersect(rownames(cluster_counts), rownames(sample_metadata))
  
  # Skip clusters with too few samples
  if (length(common_samples) < 10) {
    message(paste("  Skipping cluster", cluster_id, "- too few samples"))
    next
  }
  
  # Prepare count data for this cluster
  cluster_counts_ordered <- cluster_counts[common_samples, ]
  sample_metadata_ordered <- sample_metadata[common_samples, , drop=FALSE]
  
  # Check if we have samples from both conditions
  condition_counts <- table(sample_metadata_ordered$condition)
  if (length(condition_counts) < 2 || any(condition_counts < 3)) {
    message(paste("  Skipping cluster", cluster_id, "- insufficient samples per condition"))
    next
  }
  
  tryCatch({
    # Process this cluster's counts - passing cluster_id for proportion normalization
    cluster_results <- process_counts(cluster_counts_ordered, sample_metadata_ordered, 
                                     paste("cluster", cluster_id), cluster_id)
    
    if (!is.null(cluster_results)) {
      # Add cluster ID
      cluster_results$cluster_ID <- cluster_id
      
      # Store results
      results_by_cluster[[as.character(cluster_id)]] <- cluster_results
      
      message(paste("  Successfully processed cluster", cluster_id))
    }
  }, error = function(e) {
    message(paste("  Error processing cluster", cluster_id, ":", e$message))
  })
}

# Combine cluster results
combined_cluster_results <- bind_rows(results_by_cluster) %>%
  mutate(cluster = as.numeric(cluster_ID))

# Create output with only cluster specific features
fold_changes_feature_value <- combined_cluster_results %>%
  left_join(annots_list, by = "cluster") %>%
  mutate(
    # Match exactly the data preprocessing script logic
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    feature_cluster = paste(cluster_name, cluster_ID, sep = "_"),
    Feature = paste(gene, feature_cluster, sep = "@")
  ) %>%
  select(Feature, Value = log2FoldChange)

# Write results
write_csv(fold_changes_feature_value, "output/differential_expression/feature_fold_changes.csv")
message("Created feature_fold_changes.csv with cluster-specific normalized log2 fold changes")
