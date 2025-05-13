library(dplyr)
library(tidyverse)
library(readr)
library(Matrix)

# Create output directory
dir.create("output/differential_expression", recursive = TRUE, showWarnings = FALSE)

# Check for required files
required_files <- c(
  "base/data/lung_ldm.rd",
  "base/input_tables/table_s1_sample_table.csv",
  "base/input_tables/annots_list.csv",
  "base/input_tables/cell_metadata.csv"
)

for (file_path in required_files) {
  if (!file.exists(file_path)) {
    stop(paste("Required input file missing:", file_path))
  }
}

# Load Leader et al. data
if (!exists("lung_ldm")) {
  load("base/data/lung_ldm.rd")
}

# Load metadata
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv")
annots_list <- read_csv("base/input_tables/annots_list.csv")
cell_metadata <- read_csv("base/input_tables/cell_metadata.csv") # Added cell metadata

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

# Validate remaining data
if (length(dim(counts)[1]) == 0 || dim(counts)[1] == 0) {
  stop("No samples remain after filtering. Check sample IDs in metadata and count data.")
}
message(paste("Processing", dim(counts)[1], "samples across", dim(counts)[3], "clusters"))

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

# Set reference and comparison levels
ref_level <- "Normal" 
comparison_level <- "Tumor"

# Create sample metadata, filtering out samples with NA tissue
sample_metadata <- table_s1 %>%
  select(sample_ID, tissue) %>%
  filter(sample_ID %in% rownames(counts)) %>%
  mutate(
    condition = tissue,
    condition = factor(condition, levels = c(ref_level, comparison_level))
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
  
  # Calculate mean sample size for normalization
  mean_sample_size <- counts_long %>%
    group_by(sample_ID) %>%
    summarize(sample_size = sum(count), .groups = 'drop') %>%
    summarize(mean_size = mean(sample_size)) %>%
    pull(mean_size)
  
  # Get cluster-specific proportion information for this cluster
  cluster_props <- cluster_proportions %>%
    filter(cluster_ID == as.character(cluster_id)) %>%
    select(sample_ID, cluster_proportion)
  
  # Normalize counts by total sample counts first (CPM normalization)
  normalized_counts <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(
      total_sample_counts = sum(count),
      normalized_count = (count / total_sample_counts) * mean_sample_size
    ) %>%
    ungroup()
  
  # Apply cluster proportion correction
  normalized_counts <- normalized_counts %>%
    left_join(cluster_props, by = "sample_ID") %>%
    mutate(
      # Default to 1 if cluster_proportion is NA or too small
      cluster_proportion = ifelse(is.na(cluster_proportion) | cluster_proportion < 0.001, 
                                 1, cluster_proportion),
      # Correct for over/under-representation
      normalized_count = normalized_count / cluster_proportion
    ) %>%
    # Re-normalize to maintain the same total counts
    group_by(sample_ID) %>%
    mutate(
      total_after_norm = sum(normalized_count),
      normalized_count = normalized_count * (mean_sample_size / total_after_norm)
    ) %>%
    ungroup()
  
  # Calculate mean expression by condition
  condition_means <- normalized_counts %>%
    group_by(gene, condition) %>%
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
      # Using explicit column names matching the exact case
      mean_expr_normal = get(paste0("mean_expr_", ref_level)), # Will look for mean_expr_Normal
      mean_expr_tumor = get(paste0("mean_expr_", comparison_level)), # Will look for mean_expr_Tumor
      log2FoldChange = log2((mean_expr_tumor + 0.1) / (mean_expr_normal + 0.1)),
      n_normal = get(paste0("n_samples_", ref_level)), # Will look for n_samples_Normal
      n_tumor = get(paste0("n_samples_", comparison_level)) # Will look for n_samples_Tumor
    )
  
  return(result_df)
}

# Process individual clusters
results_by_cluster <- list()
message("Processing clusters with normalized fold change approach...")
clusters_processed <- 0
clusters_skipped <- 0

for (cluster_id in clusters_to_keep) {
  message(paste("Processing cluster:", cluster_id))
  
  # Get counts for this cluster
  cluster_counts <- counts[, , as.character(cluster_id)]
  common_samples <- intersect(rownames(cluster_counts), rownames(sample_metadata))
  
  # Skip clusters with too few samples
  if (length(common_samples) < 10) {
    message(paste("  Skipping cluster:", cluster_id, "- too few samples"))
    clusters_skipped <- clusters_skipped + 1
    next
  }
  
  # Prepare count data for this cluster
  cluster_counts_ordered <- cluster_counts[common_samples, ]
  sample_metadata_ordered <- sample_metadata[common_samples, , drop=FALSE]
  
  # Check if we have samples from both conditions
  condition_counts <- table(sample_metadata_ordered$condition)
  if (length(condition_counts) < 2 || any(condition_counts < 3)) {
    message(paste("  Skipping cluster", cluster_id, "- insufficient samples per condition"))
    clusters_skipped <- clusters_skipped + 1
    next
  }
  
  tryCatch({
    # Process this cluster's counts
    cluster_results <- process_counts(cluster_counts_ordered, sample_metadata_ordered, 
                                     paste("cluster", cluster_id), cluster_id)
    
    if (!is.null(cluster_results)) {
      # Add cluster ID
      cluster_results$cluster_ID <- cluster_id
      
      # Store results
      results_by_cluster[[as.character(cluster_id)]] <- cluster_results
      
      message(paste("Successfully processed cluster", cluster_id))
      clusters_processed <- clusters_processed + 1
    } else {
      clusters_skipped <- clusters_skipped + 1
    }
  }, error = function(e) {
    message(paste("  Error processing cluster", cluster_id, ":", e$message))
    clusters_skipped <- clusters_skipped + 1
  })
}

# Check if we processed any clusters
if (length(results_by_cluster) == 0) {
  stop("No clusters were successfully processed. Check filtering criteria and input data.")
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

# Write a more detailed results file with additional information
detailed_results <- combined_cluster_results %>%
  left_join(annots_list, by = "cluster") %>%
  mutate(
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    feature_cluster = paste(cluster_name, cluster_ID, sep = "_"),
    Feature = paste(gene, feature_cluster, sep = "@")
  )

write_csv(detailed_results, "output/differential_expression/detailed_feature_fold_changes.csv")

# Print summary statistics
message(paste("Processed", clusters_processed, "clusters successfully"))
message(paste("Skipped", clusters_skipped, "clusters due to filtering criteria"))
message(paste("Total features generated:", nrow(fold_changes_feature_value)))
message("Created feature_fold_changes.csv with cluster-specific normalized log2 fold changes")
message("Created detailed_feature_fold_changes.csv with additional information per feature")
