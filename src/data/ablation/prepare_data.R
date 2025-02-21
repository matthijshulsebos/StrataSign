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
hsa01100_genes <- read_csv("data/kegg/hsa01100_genes.csv")

# Extract metadata
table_s1 <- table_s1 %>% filter(Use.in.Clustering.Model. == "Yes")

# Convert 'sample_ID' column to character
table_s1 <- table_s1 %>% mutate(sample_ID = as.character(sample_ID))

# Define cluster groups
lcam_hi_clusters <- c(44, 9, 17, 28, 46, 11, 42)
lcam_lo_clusters <- c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16)
lcam_both_clusters <- c(lcam_hi_clusters, lcam_lo_clusters)

# Function to preprocess count data for a given cluster subset
preprocess_counts <- function(cluster_subset, subset_name) {
  counts <- lung_ldm$dataset$counts  # 3D matrix: sample × gene × cluster

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
    mutate(cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage)) %>%
    mutate(cluster_name = gsub("/", "-", cluster_name)) %>%  # Sanitize cluster names
    mutate(gene_cluster = paste(gene, cluster_name, cluster, sep = "|"))

  # Ensure count column is numeric
  counts_long$count <- as.numeric(counts_long$count)

  # Normalize counts per sample using CP10K
  counts_long <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(total_counts = sum(count)) %>%
    ungroup() %>%
    mutate(normalized_count = (count / total_counts) * 10000)

  # Apply log transformation
  counts_long <- counts_long %>%
    mutate(log_normalized_count = log1p(normalized_count))

  # Pivot so each gene-cluster combo becomes its own column
  counts_wide <- counts_long %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster, values_from = log_normalized_count, values_fill = list(log_normalized_count = 0))

  # Convert 'sample_ID' column to character
  counts_wide <- counts_wide %>% mutate(sample_ID = as.character(sample_ID))

  # Merge with metadata
  counts_wide <- counts_wide %>% inner_join(table_s1, by = "sample_ID")

  # Remove metadata columns before splitting, but keep 'tissue' as the target variable
  counts_wide <- counts_wide %>% select(-c(amp_batch_ID, old_lib_name, HTO, disease, `Use.in.Clustering.Model.`, library_chemistry, prime, vdj_kit, prep, metadata_indicator, biopsy_site))

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
  
  # Create directories if they don't exist
  dir.create(paste0("data/ablation/datasets/", subset_name, "/default"), recursive = TRUE, showWarnings = FALSE)
  
  # Save metadata in the order of the training set
  metadata_train <- table_s1 %>% filter(sample_ID %in% train$sample_ID) %>% arrange(match(sample_ID, train$sample_ID))
  write_csv(metadata_train, paste0("data/ablation/datasets/", subset_name, "/default/metadata_", subset_name, ".csv"))

  # Remove 'sample_ID' and 'patient_ID' from train and test datasets
  train <- train %>% select(-sample_ID, -patient_ID)
  test <- test %>% select(-sample_ID, -patient_ID)

  # Check for NA values in train and test dataframes
  if (anyNA(train)) {
    warning(paste("Train dataframe contains NA values for subset:", subset_name))
    na_cols <- colnames(train)[sapply(train, anyNA)]
    print(paste("Columns with NA in train:", paste(na_cols, collapse = ", ")))
  }
  if (anyNA(test)) {
    warning(paste("Test dataframe contains NA values for subset:", subset_name))
    na_cols <- colnames(test)[sapply(test, anyNA)]
    print(paste("Columns with NA in test:", paste(na_cols, collapse = ", ")))
  }

  # Save outputs
  write_csv(train, paste0("data/ablation/datasets/", subset_name, "/default/X_train_", subset_name, ".csv"))
  write_csv(test, paste0("data/ablation/datasets/", subset_name, "/default/X_test_", subset_name, ".csv"))
  write_csv(data.frame(x = y_train), paste0("data/ablation/datasets/", subset_name, "/default/y_train_", subset_name, ".csv"))
  write_csv(data.frame(x = y_test), paste0("data/ablation/datasets/", subset_name, "/default/y_test_", subset_name, ".csv"))
}

# Function to preprocess count data with random genes for a given cluster subset
preprocess_counts_random_genes <- function(cluster_subset, subset_name, random_version, seed) {
  counts <- lung_ldm$dataset$counts  # 3D matrix: sample × gene × cluster

  # Filter counts to match table_s1
  sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
  counts <- counts[dimnames(counts)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]
  
  # Get all genes and metabolic genes
  all_genes <- dimnames(counts)[[2]]
  gene_symbols <- hsa01100_genes$SYMBOL
  existing_metabolic_genes <- gene_symbols[gene_symbols %in% all_genes]
  
  # Select random genes based on version
  set.seed(seed)
  if (random_version == "random1") {
    # For random1, exclude metabolic genes
    non_metabolic_genes <- setdiff(all_genes, existing_metabolic_genes)
    random_genes <- sample(non_metabolic_genes, 1609)
    # Get overlap (should be 0)
    overlap <- intersect(random_genes, existing_metabolic_genes)
  } else {
    # For random2 and random3, sample from all genes
    random_genes <- sample(all_genes, 1609)
    # Get overlap for information
    overlap <- intersect(random_genes, existing_metabolic_genes)
  }
  
  # Create directories if they don't exist
  dir.create(paste0("data/ablation/datasets/", subset_name, "/", random_version), recursive = TRUE, showWarnings = FALSE)
  
  # Save overlap information
  overlap_info <- data.frame(
    dataset = subset_name,
    version = random_version,
    total_genes = length(random_genes),
    metabolic_genes = length(overlap),
    overlap_percentage = round(length(overlap) / length(random_genes) * 100, 2),
    overlapping_genes = paste(overlap, collapse = ", ")
  )
  write_csv(overlap_info, paste0("data/ablation/datasets/", subset_name, "/", random_version, "/gene_overlap_info.csv"))
  
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
    mutate(cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage)) %>%
    mutate(cluster_name = gsub("/", "-", cluster_name)) %>%  # Sanitize cluster names
    mutate(gene_cluster = paste(gene, cluster_name, cluster, sep = "|"))

  # Ensure count column is numeric
  counts_long$count <- as.numeric(counts_long$count)

  # Normalize counts per sample using CP10K
  counts_long <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(total_counts = sum(count)) %>%
    ungroup() %>%
    mutate(normalized_count = (count / total_counts) * 10000)

  # Apply log transformation
  counts_long <- counts_long %>%
    mutate(log_normalized_count = log1p(normalized_count))

  # Pivot so each gene-cluster combo becomes its own column
  counts_wide <- counts_long %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster, values_from = log_normalized_count, values_fill = list(log_normalized_count = 0))

  # Convert 'sample_ID' column to character
  counts_wide <- counts_wide %>% mutate(sample_ID = as.character(sample_ID))

  # Merge with metadata
  counts_wide <- counts_wide %>% inner_join(table_s1, by = "sample_ID")

  # Remove metadata columns before splitting, but keep 'tissue' as the target variable
  counts_wide <- counts_wide %>% select(-c(amp_batch_ID, old_lib_name, HTO, disease, `Use.in.Clustering.Model.`, library_chemistry, prime, vdj_kit, prep, metadata_indicator, biopsy_site))

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
  
  # Create directories if they don't exist
  dir.create(paste0("data/ablation/datasets/", subset_name, "/", random_version), recursive = TRUE, showWarnings = FALSE)
  
  # Save metadata in the order of the training set
  metadata_train <- table_s1 %>% filter(sample_ID %in% train$sample_ID) %>% arrange(match(sample_ID, train$sample_ID))
  write_csv(metadata_train, paste0("data/ablation/datasets/", subset_name, "/", random_version, "/metadata_", subset_name, "_", random_version, ".csv"))

  # Remove 'sample_ID' and 'patient_ID' from train and test datasets
  train <- train %>% select(-sample_ID, -patient_ID)
  test <- test %>% select(-sample_ID, -patient_ID)

  # Check for NA values in train and test dataframes
  if (anyNA(train)) {
    warning(paste("Train dataframe contains NA values for subset:", subset_name, random_version))
    na_cols <- colnames(train)[sapply(train, anyNA)]
    print(paste("Columns with NA in train:", paste(na_cols, collapse = ", ")))
  }
  if (anyNA(test)) {
    warning(paste("Test dataframe contains NA values for subset:", subset_name, random_version))
    na_cols <- colnames(test)[sapply(test, anyNA)]
    print(paste("Columns with NA in test:", paste(na_cols, collapse = ", ")))
  }

  # Save outputs
  write_csv(train, paste0("data/ablation/datasets/", subset_name, "/", random_version, "/X_train_", subset_name, "_", random_version, ".csv"))
  write_csv(test, paste0("data/ablation/datasets/", subset_name, "/", random_version, "/X_test_", subset_name, "_", random_version, ".csv"))
  write_csv(data.frame(x = y_train), paste0("data/ablation/datasets/", subset_name, "/", random_version, "/y_train_", subset_name, "_", random_version, ".csv"))
  write_csv(data.frame(x = y_test), paste0("data/ablation/datasets/", subset_name, "/", random_version, "/y_test_", subset_name, "_", random_version, ".csv"))
  write_csv(overlap_info, paste0("data/ablation/datasets/", subset_name, "/", random_version, "/gene_overlap_info.csv"))
}

# Process all cluster subsets with metabolic genes
preprocess_counts(lcam_hi_clusters, "lcam_hi")
preprocess_counts(lcam_lo_clusters, "lcam_lo")
preprocess_counts(lcam_both_clusters, "lcam_both")
preprocess_counts(NULL, "all_clusters")

# Process all cluster subsets with random genes
preprocess_counts_random_genes(lcam_hi_clusters, "lcam_hi", "random1", seed = 101)
preprocess_counts_random_genes(lcam_lo_clusters, "lcam_lo", "random1", seed = 101)
preprocess_counts_random_genes(lcam_both_clusters, "lcam_both", "random1", seed = 101)
preprocess_counts_random_genes(NULL, "all_clusters", "random1", seed = 101)

preprocess_counts_random_genes(lcam_hi_clusters, "lcam_hi", "random2", seed = 102)
preprocess_counts_random_genes(lcam_lo_clusters, "lcam_lo", "random2", seed = 102)
preprocess_counts_random_genes(lcam_both_clusters, "lcam_both", "random2", seed = 102)
preprocess_counts_random_genes(NULL, "all_clusters", "random2", seed = 102)

preprocess_counts_random_genes(lcam_hi_clusters, "lcam_hi", "random3", seed = 103)
preprocess_counts_random_genes(lcam_lo_clusters, "lcam_lo", "random3", seed = 103)
preprocess_counts_random_genes(lcam_both_clusters, "lcam_both", "random3", seed = 103)
preprocess_counts_random_genes(NULL, "all_clusters", "random3", seed = 103)
