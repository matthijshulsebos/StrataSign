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
cell_metadata <- read_csv("base/input_tables/cell_metadata.csv")

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

# Function to preprocess count data for a given cluster subset and gene set type
preprocess_counts <- function(cluster_subset, subset_name, gene_set_type = "metabolic") {
  counts <- lung_ldm$dataset$counts  # 3D matrix: sample × gene × cluster

  # First, filter out the unwanted clusters
  all_clusters <- as.numeric(dimnames(counts)[[3]])
  clusters_to_keep <- all_clusters[!all_clusters %in% clusters_to_exclude]
  counts <- counts[, , as.character(clusters_to_keep), drop = FALSE]

  # Filter counts to match table_s1
  sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
  counts <- counts[dimnames(counts)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]
  
  # Filter genes based on gene set type
  all_genes <- dimnames(counts)[[2]]
  metabolic_genes <- hsa01100_genes$SYMBOL[hsa01100_genes$SYMBOL %in% all_genes]
  
  # Filter genes based on the gene set type
  if (gene_set_type == "metabolic") {
    # Use metabolic genes
    genes_to_use <- metabolic_genes
    suffix <- "_metabolic"
  } else if (gene_set_type == "nonmetabolic") {
    # Use non-metabolic genes
    genes_to_use <- setdiff(all_genes, metabolic_genes)
    # Limit to the same number of genes as metabolic for fair comparison
    if (length(genes_to_use) > length(metabolic_genes)) {
      set.seed(42)  # For reproducibility
      genes_to_use <- sample(genes_to_use, length(metabolic_genes))
    }
    suffix <- "_nonmetabolic"
  } else if (gene_set_type == "random") {
    # Use random genes
    set.seed(43) 
    genes_to_use <- sample(all_genes, length(metabolic_genes))
    suffix <- "_random"
  } else {
    stop("Invalid gene set type. Must be one of: metabolic, nonmetabolic, random")
  }
  
  counts <- counts[, genes_to_use, , drop = FALSE]
  
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
      cluster_ID = cluster,  # Preserve the numeric cluster ID for joining
      cluster = paste(cluster_name, cluster, sep = "_"),  # Create cluster identifier
      gene_cluster = paste(gene, cluster, sep = "@")  # Create feature name
    )

  # Ensure count column is numeric
  counts_long$count <- as.numeric(counts_long$count)

  # Convert 'sample_ID' to character
  counts_long$sample_ID <- as.character(counts_long$sample_ID)

  # Step 1: Normalize by sample total to account for sequencing depth differences
  # Calculate mean sample size for sample normalization
  mean_sample_size <- counts_long %>%
    group_by(sample_ID) %>%
    summarize(sample_size = sum(count), .groups = 'drop') %>%
    summarize(mean_size = mean(sample_size)) %>%
    pull(mean_size)

  # Normalize counts by total sample counts
  counts_long <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(total_sample_counts = sum(count),
           sample_normalized_count = (count / total_sample_counts) * mean_sample_size) %>%
    ungroup()

  # Step 2: Normalize each cluster within each sample based on real cluster proportions
  # Calculate the number of cells per cluster per sample
  cluster_proportions <- cell_metadata %>%
    group_by(sample_ID, cluster_ID) %>%
    summarize(cluster_cell_count = n(), .groups = 'drop') %>%
    group_by(sample_ID) %>%
    mutate(total_cells = sum(cluster_cell_count),
           cluster_proportion = cluster_cell_count / total_cells) %>%
    ungroup()

  # Convert 'sample_ID' to character in cluster_proportions
  cluster_proportions$sample_ID <- as.character(cluster_proportions$sample_ID)

  # Join cluster proportions with counts_long
  counts_long <- counts_long %>%
    left_join(cluster_proportions, by = c("sample_ID", "cluster_ID")) %>%
    # Add a check for missing cluster proportions
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
    ungroup() %>%
    # Clean up intermediate columns
    select(-total_sample_counts, -cluster_cell_count, -total_cells, -cluster_proportion, -cluster_ID, -total_after_norm)

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
  
  # Check for and handle NAs
  if (anyNA(train) || anyNA(test)) {
    # Identify columns with NA values
    na_cols_train <- colnames(train)[colSums(is.na(train)) > 0]
    na_cols_test <- colnames(test)[colSums(is.na(test)) > 0]
    
    if (length(na_cols_train) > 0 || length(na_cols_test) > 0) {
      message("Replacing NA values with 0 for gene set type: ", gene_set_type)
      # Replace NA values with 0
      train[is.na(train)] <- 0
      test[is.na(test)] <- 0
    }
  }
  
  # Create directories if they don't exist
  output_dir <- paste0("output/training_data/raw/", subset_name, "/", gene_set_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save metadata in the order of the training set
  metadata_train <- table_s1 %>% filter(sample_ID %in% train$sample_ID) %>% arrange(match(sample_ID, train$sample_ID))
  write_csv(metadata_train, paste0(output_dir, "/metadata_", subset_name, suffix, ".csv"))

  # Remove 'sample_ID' and 'patient_ID' from train and test datasets
  train <- train %>% select(-sample_ID, -patient_ID)
  test <- test %>% select(-sample_ID, -patient_ID)

  # Save outputs with consistent naming
  write_csv(train, paste0(output_dir, "/X_train_", subset_name, suffix, ".csv"))
  write_csv(test, paste0(output_dir, "/X_test_", subset_name, suffix, ".csv"))
  write_csv(data.frame(x = y_train), paste0(output_dir, "/y_train_", subset_name, suffix, ".csv"))
  write_csv(data.frame(x = y_test), paste0(output_dir, "/y_test_", subset_name, suffix, ".csv"))
  
  message("Finished processing ", gene_set_type, " genes for ", subset_name)
}

# Process with metabolic genes
message("Processing metabolic gene sets...")
preprocess_counts(lcam_hi_clusters, "lcam_hi", "metabolic")
preprocess_counts(lcam_lo_clusters, "lcam_lo", "metabolic")
preprocess_counts(lcam_both_clusters, "lcam_both", "metabolic")
preprocess_counts(NULL, "all_clusters", "metabolic")

# Process with non-metabolic genes 
message("Processing non-metabolic gene sets...")
preprocess_counts(lcam_hi_clusters, "lcam_hi", "nonmetabolic")
preprocess_counts(lcam_lo_clusters, "lcam_lo", "nonmetabolic")
preprocess_counts(lcam_both_clusters, "lcam_both", "nonmetabolic")
preprocess_counts(NULL, "all_clusters", "nonmetabolic")

# Process with random genes
message("Processing random gene sets...")
preprocess_counts(lcam_hi_clusters, "lcam_hi", "random")
preprocess_counts(lcam_lo_clusters, "lcam_lo", "random")
preprocess_counts(lcam_both_clusters, "lcam_both", "random")
preprocess_counts(NULL, "all_clusters", "random")

message("All preprocessing complete!")
