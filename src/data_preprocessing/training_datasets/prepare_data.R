library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(scales)

# Load data
if (!exists("lung_ldm")) {
  message("Loading Mt. Sinai scRNAseq data into R")
  load("base/data/lung_ldm.rd")
}

table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv")
annots_list <- read_csv("base/input_tables/annots_list.csv")
hsa01100_genes <- read_csv("output/kegg/hsa01100_genes.csv")

# Extract metadata
table_s1 <- table_s1 %>% filter(Use.in.Clustering.Model. == "Yes")

# Convert 'sample_ID' column to character
table_s1 <- table_s1 %>% mutate(sample_ID = as.character(sample_ID))

# Define clusters to exclude
clusters_to_exclude <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
message(paste("Excluding doublets:", paste(clusters_to_exclude, collapse = ", ")))

# Define cluster groups
lcam_hi_clusters <- c(44, 9, 17, 28, 46, 11, 42)
lcam_lo_clusters <- c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16)
lcam_both_clusters <- c(lcam_hi_clusters, lcam_lo_clusters)

# Function to preprocess count data for a given cluster subset
preprocess_counts <- function(cluster_subset, subset_name) {
  counts <- lung_ldm$dataset$counts  # 3D matrix: sample × gene × cluster

  # First, filter out the unwanted clusters
  all_clusters <- as.numeric(dimnames(counts)[[3]])
  clusters_to_keep <- all_clusters[!all_clusters %in% clusters_to_exclude]
  counts <- counts[, , as.character(clusters_to_keep), drop = FALSE]

  # Filter counts to match table_s1
  sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
  counts <- counts[dimnames(counts)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]
  
  # Filter genes based on metabolic genes
  gene_symbols <- hsa01100_genes$SYMBOL
  existing_genes <- gene_symbols[gene_symbols %in% dimnames(counts)[[2]]]
  counts <- counts[, existing_genes, , drop = FALSE]
  
  # Subset clusters
  if (!is.null(cluster_subset)) {
    counts <- counts[, , dimnames(counts)[[3]] %in% cluster_subset, drop = FALSE]
  }
  
  # Convert 3D to 2D (sample × gene-cluster)
  counts_long <- as.data.frame(as.table(counts))
  colnames(counts_long) <- c("sample_ID", "gene", "cluster", "count")

  # Convert cluster column from factor to numeric
  counts_long$cluster <- as.numeric(as.character(counts_long$cluster))

  # Convert cluster numbers to sublineage or lineage names
  counts_long <- counts_long %>%
    left_join(annots_list, by = c("cluster" = "cluster")) %>%
    mutate(
      cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
      cluster_name = gsub("/", "-", cluster_name),  # Sanitize cluster names
      cluster = paste(cluster_name, cluster, sep = "_"),  # Create cluster identifier
      gene_cluster = paste(gene, cluster, sep = "@")  # Create feature name
    )

  # Ensure count column is numeric
  counts_long$count <- as.numeric(counts_long$count)

  # Step 1: Normalize by sample total to account for sequencing depth differences
  # Calculate mean sample size for sample normalization
  mean_sample_size <- counts_long %>%
    group_by(sample_ID) %>%
    summarize(sample_size = sum(count)) %>%
    ungroup() %>%
    summarize(mean_size = mean(sample_size)) %>%
    pull(mean_size)

  # Normalize counts by total sample counts
  counts_long <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(total_sample_counts = sum(count),
           sample_normalized_count = (count / total_sample_counts) * mean_sample_size) %>%
    ungroup()

  # Step 2: Now normalize each cluster within each sample to have equal total counts
  # First determine the target count per cluster per sample
  target_cluster_count <- counts_long %>%
    # Count the number of clusters per sample
    group_by(sample_ID) %>%
    summarize(num_clusters = n_distinct(cluster),
              # Divide the sample total by number of clusters to get target per cluster
              target_count = mean_sample_size / num_clusters) %>%
    ungroup()

  # Then normalize each cluster to have the target count
  counts_long <- counts_long %>%
    # Calculate total counts per cluster per sample
    group_by(sample_ID, cluster) %>%
    mutate(cluster_total = sum(sample_normalized_count)) %>%
    ungroup() %>%
    # Join with target counts
    left_join(target_cluster_count, by = "sample_ID") %>%
    # Scale each cluster to the target count - with protection against division by zero
    mutate(
      # Avoid division by zero by adding small epsilon
      cluster_total = ifelse(cluster_total == 0, 1e-10, cluster_total),
      scaling_factor = target_count / cluster_total,
      normalized_count = sample_normalized_count * scaling_factor
    ) %>%
    # Clean up intermediate columns
    select(-total_sample_counts, -cluster_total, -num_clusters, 
           -target_count, -scaling_factor)

  # Apply log transformation
  counts_long <- counts_long %>%
    mutate(log_normalized_count = log1p(normalized_count))

  # Pivot so each gene-cluster combo becomes its own column
  counts_wide <- counts_long %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster, values_from = log_normalized_count, values_fill = list(log_normalized_count = 0))

  # Merge with metadata - but only keep essential columns
  counts_wide <- counts_wide %>% 
    mutate(sample_ID = as.character(sample_ID)) %>%
    left_join(
      table_s1 %>% select(sample_ID, patient_ID, tissue),
      by = "sample_ID"
    )

  # Split train/test by patient_ID
  set.seed(42)
  patients <- unique(counts_wide$patient_ID)
  train_patients <- sample(patients, size = round(length(patients) * 0.7))
  train <- counts_wide %>% filter(patient_ID %in% train_patients)
  test <- counts_wide %>% filter(!patient_ID %in% train_patients)
  
  # Extract and remove the target variable 'tissue'
  y_train <- train$tissue
  y_test <- test$tissue
  train <- train %>% select(-tissue)
  test <- test %>% select(-tissue)
  
  # Add right before saving - check for NAs
  if (anyNA(train) || anyNA(test)) {
    stop("Final training/test data still has NA values - aborting save")
  }
  
  # Create directories if they don't exist - use new name directly
  output_dir <- paste0("output/training_data/raw/", subset_name, "/metabolic")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save metadata in the order of the training set - add _metabolic suffix
  metadata_train <- table_s1 %>% filter(sample_ID %in% train$sample_ID) %>% arrange(match(sample_ID, train$sample_ID))
  write_csv(metadata_train, paste0(output_dir, "/metadata_", subset_name, "_metabolic.csv"))

  # Remove 'sample_ID' and 'patient_ID' from train and test datasets
  train <- train %>% select(-sample_ID, -patient_ID)
  test <- test %>% select(-sample_ID, -patient_ID)

  # Save outputs with consistent naming (add _metabolic suffix)
  write_csv(train, paste0(output_dir, "/X_train_", subset_name, "_metabolic.csv"))
  write_csv(test, paste0(output_dir, "/X_test_", subset_name, "_metabolic.csv"))
  write_csv(data.frame(x = y_train), paste0(output_dir, "/y_train_", subset_name, "_metabolic.csv"))
  write_csv(data.frame(x = y_test), paste0(output_dir, "/y_test_", subset_name, "_metabolic.csv"))
}

# Function to preprocess count data with random genes for a given cluster subset
preprocess_counts_random_genes <- function(cluster_subset, subset_name, random_version, seed) {
  counts <- lung_ldm$dataset$counts  # 3D matrix: sample × gene × cluster
  
  # First, filter out the unwanted clusters
  all_clusters <- as.numeric(dimnames(counts)[[3]])
  clusters_to_keep <- all_clusters[!all_clusters %in% clusters_to_exclude]
  counts <- counts[, , as.character(clusters_to_keep), drop = FALSE]

  # Filter counts to match table_s1
  sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
  counts <- counts[dimnames(counts)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]
  
  # Get all genes and metabolic genes
  all_genes <- dimnames(counts)[[2]]
  gene_symbols <- hsa01100_genes$SYMBOL
  existing_metabolic_genes <- gene_symbols[gene_symbols %in% all_genes]
  
  # Select random genes based on version
  set.seed(seed)
  if (random_version == "matched_nonmetabolic") {
    # For matched_nonmetabolic, exclude metabolic genes
    non_metabolic_genes <- setdiff(all_genes, existing_metabolic_genes)
    random_genes <- sample(non_metabolic_genes, 1609)
    # Get overlap (should be 0)
    overlap <- intersect(random_genes, existing_metabolic_genes)
  } else {
    # For matched_random_a and matched_random_b, sample from all genes
    random_genes <- sample(all_genes, 1609)
    # Get overlap for information
    overlap <- intersect(random_genes, existing_metabolic_genes)
  }
  
  # Create directories if they don't exist - use new names directly
  output_dir <- paste0("output/training_data/raw/", subset_name, "/", random_version)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save overlap information but limit gene output to console
  overlap_info <- data.frame(
    dataset = subset_name,
    version = random_version,
    total_genes = length(random_genes),
    metabolic_genes = length(overlap),
    overlap_percentage = round(length(overlap) / length(random_genes) * 100, 2)
  )

  # Include overlapping genes in the CSV file
  if (length(overlap) > 0) {
    overlap_info$overlapping_genes <- paste(overlap, collapse = ", ")
  } else {
    overlap_info$overlapping_genes <- "None"
  }
  
  counts <- counts[, random_genes, , drop = FALSE]
  
  # Subset clusters
  if (!is.null(cluster_subset)) {
    counts <- counts[, , dimnames(counts)[[3]] %in% cluster_subset, drop = FALSE]
  }
  
  # Convert 3D to 2D (sample × gene-cluster)
  counts_long <- as.data.frame(as.table(counts))
  colnames(counts_long) <- c("sample_ID", "gene", "cluster", "count")

  # Convert cluster column from factor to numeric
  counts_long$cluster <- as.numeric(as.character(counts_long$cluster))

  # Convert cluster numbers to sublineage or lineage names
  counts_long <- counts_long %>%
    left_join(annots_list, by = c("cluster" = "cluster")) %>%
    mutate(
      cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
      cluster_name = gsub("/", "-", cluster_name),  # Sanitize cluster names
      cluster = paste(cluster_name, cluster, sep = "_"),  # Create cluster identifier
      gene_cluster = paste(gene, cluster, sep = "@")  # Create feature name
    )

  # Ensure count column is numeric
  counts_long$count <- as.numeric(counts_long$count)

  # Step 1: Normalize by sample total to account for sequencing depth differences
  # Calculate mean sample size for sample normalization
  mean_sample_size <- counts_long %>%
    group_by(sample_ID) %>%
    summarize(sample_size = sum(count)) %>%
    ungroup() %>%
    summarize(mean_size = mean(sample_size)) %>%
    pull(mean_size)

  # Normalize counts by total sample counts
  counts_long <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(total_sample_counts = sum(count),
           sample_normalized_count = (count / total_sample_counts) * mean_sample_size) %>%
    ungroup()

  # Step 2: Now normalize each cluster within each sample to have equal total counts
  # First determine the target count per cluster per sample
  target_cluster_count <- counts_long %>%
    # Count the number of clusters per sample
    group_by(sample_ID) %>%
    summarize(num_clusters = n_distinct(cluster),
              # Divide the sample total by number of clusters to get target per cluster
              target_count = mean_sample_size / num_clusters) %>%
    ungroup()

  # Then normalize each cluster to have the target count
  counts_long <- counts_long %>%
    # Calculate total counts per cluster per sample
    group_by(sample_ID, cluster) %>%
    mutate(cluster_total = sum(sample_normalized_count)) %>%
    ungroup() %>%
    # Join with target counts
    left_join(target_cluster_count, by = "sample_ID") %>%
    # Scale each cluster to the target count - with protection against division by zero
    mutate(
      # Avoid division by zero by adding small epsilon
      cluster_total = ifelse(cluster_total == 0, 1e-10, cluster_total),
      scaling_factor = target_count / cluster_total,
      normalized_count = sample_normalized_count * scaling_factor
    ) %>%
    # Clean up intermediate columns
    select(-total_sample_counts, -cluster_total, -num_clusters, 
           -target_count, -scaling_factor)

  # Apply log transformation
  counts_long <- counts_long %>%
    mutate(log_normalized_count = log1p(normalized_count))

  # Pivot so each gene-cluster combo becomes its own column
  counts_wide <- counts_long %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster, values_from = log_normalized_count, values_fill = list(log_normalized_count = 0))

  # Merge with metadata - but only keep essential columns
  counts_wide <- counts_wide %>% 
    mutate(sample_ID = as.character(sample_ID)) %>%
    left_join(
      table_s1 %>% select(sample_ID, patient_ID, tissue),
      by = "sample_ID"
    )

  # Split train/test by patient_ID
  set.seed(42)
  patients <- unique(counts_wide$patient_ID)
  train_patients <- sample(patients, size = round(length(patients) * 0.7))
  train <- counts_wide %>% filter(patient_ID %in% train_patients)
  test <- counts_wide %>% filter(!patient_ID %in% train_patients)
  
  # Extract and remove the target variable 'tissue'
  y_train <- train$tissue
  y_test <- test$tissue
  train <- train %>% select(-tissue)
  test <- test %>% select(-tissue)
  
  # Add right before saving - check for NAs
  if (anyNA(train) || anyNA(test)) {
    stop("Final training/test data still has NA values - aborting save")
  }
  
  # Save metadata in the order of the training set
  metadata_train <- table_s1 %>% filter(sample_ID %in% train$sample_ID) %>% arrange(match(sample_ID, train$sample_ID))
  write_csv(metadata_train, paste0(output_dir, "/metadata_", subset_name, "_", random_version, ".csv"))

  # Remove 'sample_ID' and 'patient_ID' from train and test datasets
  train <- train %>% select(-sample_ID, -patient_ID)
  test <- test %>% select(-sample_ID, -patient_ID)

  # Save outputs with new naming convention
  write_csv(train, paste0(output_dir, "/X_train_", subset_name, "_", random_version, ".csv"))
  write_csv(test, paste0(output_dir, "/X_test_", subset_name, "_", random_version, ".csv"))
  write_csv(data.frame(x = y_train), paste0(output_dir, "/y_train_", subset_name, "_", random_version, ".csv"))
  write_csv(data.frame(x = y_test), paste0(output_dir, "/y_test_", subset_name, "_", random_version, ".csv"))
  write_csv(overlap_info, paste0(output_dir, "/gene_overlap_info.csv"))
}

# Process all cluster subsets with metabolic genes
preprocess_counts(lcam_hi_clusters, "lcam_hi")
preprocess_counts(lcam_lo_clusters, "lcam_lo")
preprocess_counts(lcam_both_clusters, "lcam_both")
preprocess_counts(NULL, "all_clusters")

# Process all cluster subsets with random genes using CORRECT naming
preprocess_counts_random_genes(lcam_hi_clusters, "lcam_hi", "matched_nonmetabolic", seed = 101)
preprocess_counts_random_genes(lcam_lo_clusters, "lcam_lo", "matched_nonmetabolic", seed = 101)
preprocess_counts_random_genes(lcam_both_clusters, "lcam_both", "matched_nonmetabolic", seed = 101)
preprocess_counts_random_genes(NULL, "all_clusters", "matched_nonmetabolic", seed = 101)

preprocess_counts_random_genes(lcam_hi_clusters, "lcam_hi", "matched_random_a", seed = 102)
preprocess_counts_random_genes(lcam_lo_clusters, "lcam_lo", "matched_random_a", seed = 102)
preprocess_counts_random_genes(lcam_both_clusters, "lcam_both", "matched_random_a", seed = 102)
preprocess_counts_random_genes(NULL, "all_clusters", "matched_random_a", seed = 102)

preprocess_counts_random_genes(lcam_hi_clusters, "lcam_hi", "matched_random_b", seed = 103)
preprocess_counts_random_genes(lcam_lo_clusters, "lcam_lo", "matched_random_b", seed = 103)
preprocess_counts_random_genes(lcam_both_clusters, "lcam_both", "matched_random_b", seed = 103)
preprocess_counts_random_genes(NULL, "all_clusters", "matched_random_b", seed = 103)

print("Data preprocessing completed. Files saved to output/training_data/raw directory.")
