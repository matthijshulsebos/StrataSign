library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(scales)

DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
GENE_SET_TYPES <- c("metabolic", "nonmetabolic", "random")
BASE_OUTPUT_DIR <- file.path("output", "1. data preprocessing", "training datasets", "relative")
CLUSTER_DEFINITIONS <- list(
  lcam_hi = c(44, 9, 17, 28, 46, 11, 42),
  lcam_lo = c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16),
  lcam_both = c(c(44, 9, 17, 28, 46, 11, 42), c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16)),
  all_clusters = NULL
)

# Load the Leader dataset
message("Loading Leader et al. data...")
if (!exists("lung_ldm")) {
  load("base/data/lung_ldm.rd")
}

table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv") %>% 
  mutate(sample_ID = as.character(sample_ID))
cell_metadata <- read_csv("base/input_tables/cell_metadata.csv") %>%
  mutate(sample_ID = as.character(sample_ID))
annots_list <- read_csv("base/input_tables/annots_list.csv")
hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv")

# 3D matrix: sample × gene × cluster
counts_initial <- lung_ldm$dataset$counts  

# Filter out doublets
all_clusters_initial <- as.numeric(dimnames(counts_initial)[[3]])
clusters_to_keep_initial <- all_clusters_initial[!all_clusters_initial %in% CLUSTERS_TO_EXCLUDE]
counts_filtered <- counts_initial[, , as.character(clusters_to_keep_initial), drop = FALSE]

# Filter counts to match samples present in table_s1
sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
counts_filtered <- counts_filtered[dimnames(counts_filtered)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]

# Process a specific subset defined by cluster and gene type
preprocess_subset <- function(counts_input_3d, # Use the doublet filtered 3D matrix
                              cluster_subset_ids, 
                              subset_name, 
                              gene_set_type,
                              hsa01100_genes_df,
                              annots_list_df,
                              cell_metadata_df,
                              table_s1_df,
                              base_output_path) {
  
  message(paste("Processing Subset:", subset_name, "| Gene Set:", gene_set_type))
  
  # Filter genes based on type
  all_genes <- dimnames(counts_input_3d)[[2]]
  metabolic_genes <- hsa01100_genes_df$SYMBOL[hsa01100_genes_df$SYMBOL %in% all_genes]
  
  if (gene_set_type == "metabolic") {
    genes_to_use <- metabolic_genes
    suffix <- "_metabolic"
  } else if (gene_set_type == "nonmetabolic") {
    nonmetabolic_genes <- setdiff(all_genes, metabolic_genes)
    if (length(nonmetabolic_genes) > length(metabolic_genes)) {
      set.seed(42)
      genes_to_use <- sample(nonmetabolic_genes, length(metabolic_genes))
    } else {
       genes_to_use <- nonmetabolic_genes
    }
    suffix <- "_nonmetabolic"
  } else if (gene_set_type == "random") {
    set.seed(43) 
    sample_size <- min(length(metabolic_genes), length(all_genes))
    genes_to_use <- sample(all_genes, sample_size)
    suffix <- "_random"
  } else {
    stop("Invalid gene set type.")
  }
  
  counts_genes_filtered <- counts_input_3d[, genes_to_use, , drop = FALSE]
  
  # Filter clusters for the specific subset
  if (!is.null(cluster_subset_ids)) {
    counts_subset <- counts_genes_filtered[, , dimnames(counts_genes_filtered)[[3]] %in% cluster_subset_ids, drop = FALSE]
  } else {
    counts_subset <- counts_genes_filtered
  }
  
  # Check dimensions after filtering
  if (any(dim(counts_subset) == 0)) {
      warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero dimensions after filtering. Skipping."))
      return() 
  }

  # Convert to long format and cast types
  counts_long <- as.data.frame(as.table(counts_subset))
  colnames(counts_long) <- c("sample_ID", "gene", "cluster_ID", "count")
  
  counts_long <- counts_long %>%
    mutate(
      sample_ID = as.character(sample_ID),
      cluster_ID = as.numeric(as.character(cluster_ID)),
      gene = as.character(gene)
    ) %>%
  filter(count > 0) # Filter zero counts for efficiency

  # Check if data remains after filtering
  if (nrow(counts_long) == 0) {
      warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after filtering."))
      return()
  }

  # Create feature names
  counts_long <- counts_long %>%
    left_join(annots_list_df, by = c("cluster_ID" = "cluster")) %>%
    mutate(
      cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
      cluster_name = gsub("/", "-", cluster_name), 
      cluster_identifier = paste(cluster_name, cluster_ID, sep = "_"),
      gene_cluster = paste(gene, cluster_identifier, sep = "@") 
    ) %>%
    select(sample_ID, gene, cluster_ID, gene_cluster, count)

  # Check for data
  if (nrow(counts_long) == 0) {
      warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after annotation step."))
      return()
  }

  # Calculate mean total counts per sample for this specific subset (using the modified counts_long)
  mean_sample_size <- counts_long %>%
    group_by(sample_ID) %>%
    summarize(sample_size = sum(count), .groups = 'drop') %>%
    summarize(mean_size = mean(sample_size)) %>%
    pull(mean_size)

  # Normalize counts by total sample counts
  counts_sample_norm <- counts_long %>%
    group_by(sample_ID) %>%
    mutate(total_sample_counts = sum(count),
           sample_normalized_count = ifelse(total_sample_counts > 0, (count / total_sample_counts) * mean_sample_size, 0)) %>%
    ungroup()
  
  # Calculate cluster proportions using cell_metadata
  cluster_proportions <- cell_metadata_df %>%
    filter(sample_ID %in% unique(counts_sample_norm$sample_ID)) %>%
    group_by(sample_ID, cluster_ID) %>%
    summarize(cluster_cell_count = n(), .groups = 'drop') %>%
    group_by(sample_ID) %>%
    mutate(total_cells = sum(cluster_cell_count),
           cluster_proportion = ifelse(total_cells > 0, cluster_cell_count / total_cells, 0)) %>%
    ungroup() %>%
    select(sample_ID, cluster_ID, cluster_proportion)

  # Join proportions and normalize by proportion then renormalize 
  counts_relative_norm <- counts_sample_norm %>%
    left_join(cluster_proportions, by = c("sample_ID", "cluster_ID")) %>%
    mutate(
      # Handle NA/zero proportions
      cluster_proportion_safe = ifelse(is.na(cluster_proportion) | cluster_proportion == 0, 1, cluster_proportion),
      normalized_count_prop = sample_normalized_count / cluster_proportion_safe
    ) %>%
    group_by(sample_ID) %>%
    mutate(
      total_after_prop_norm = sum(normalized_count_prop),
      # Handle potential division by zero during renormalization
      normalized_count = ifelse(total_after_prop_norm > 0, normalized_count_prop * (mean_sample_size / total_after_prop_norm), 0)
    ) %>%
    ungroup() %>%
    select(sample_ID, gene_cluster, normalized_count) # Keep only needed columns

  # Filter samples with NaN/zero total normalized counts after all steps
  final_norm_check <- counts_relative_norm %>%
    group_by(sample_ID) %>%
    summarize(total = sum(normalized_count), .groups = 'drop')
  
  samples_to_filter <- final_norm_check %>%
    filter(is.nan(total) | total < 1e-9) %>% 
    pull(sample_ID)
  
  if (length(samples_to_filter) > 0) {
    message(paste("Filtering", length(samples_to_filter), 
                  "samples with NaN/zero total normalized counts for:", subset_name, "|", gene_set_type))
    counts_relative_norm <- counts_relative_norm %>%
      filter(!sample_ID %in% samples_to_filter)
  }
  
  if (nrow(counts_relative_norm) == 0) {
      warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after normalization. Skipping."))
      return() 
  }

  # Log transform and pivot to wide format
  counts_wide <- counts_relative_norm %>%
    mutate(log_normalized_count = log1p(normalized_count)) %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster, 
                values_from = log_normalized_count, 
                values_fill = 0)

  # Merge metadata and split train/test
  counts_wide_meta <- counts_wide %>% 
    left_join(
      table_s1_df %>% select(sample_ID, patient_ID, tissue),
      by = "sample_ID"
    ) %>%
    filter(!is.na(tissue)) 

  if (nrow(counts_wide_meta) == 0) {
    warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after merging metadata/filtering NAs."))
    return()
  }

  set.seed(42)
  patients <- unique(counts_wide_meta$patient_ID)
  train_patients <- sample(patients, size = round(length(patients) * 0.7))
  
  train_df <- counts_wide_meta %>% filter(patient_ID %in% train_patients)
  test_df <- counts_wide_meta %>% filter(!patient_ID %in% train_patients)

  if (nrow(train_df) == 0 || nrow(test_df) == 0) {
     warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in empty train/test set. Skipping."))
     return()
  }

  y_train <- train_df$tissue
  y_test <- test_df$tissue
  
  # Prepare metadata for saving
  metadata_train <- table_s1_df %>% 
    filter(sample_ID %in% train_df$sample_ID) %>% 
    arrange(match(sample_ID, train_df$sample_ID))

  # Prepare feature matrices (remove metadata)
  X_train <- train_df %>% select(-sample_ID, -patient_ID, -tissue)
  X_test <- test_df %>% select(-sample_ID, -patient_ID, -tissue)

  # Ensure numeric and handle NAs (should be 0 from pivot)
  X_train <- X_train %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))
  X_test <- X_test %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))

  # Save outputs
  output_dir <- file.path(base_output_path, subset_name, gene_set_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  write_csv(metadata_train, file.path(output_dir, paste0("metadata_", subset_name, suffix, ".csv")))
  write_csv(X_train, file.path(output_dir, paste0("X_train_", subset_name, suffix, ".csv")))
  write_csv(X_test, file.path(output_dir, paste0("X_test_", subset_name, suffix, ".csv")))
  write_csv(data.frame(x = y_train), file.path(output_dir, paste0("y_train_", subset_name, suffix, ".csv")))
  write_csv(data.frame(x = y_test), file.path(output_dir, paste0("y_test_", subset_name, suffix, ".csv")))
  
  message(paste("Finished:", subset_name, "|", gene_set_type))
}

message("\n Preprocessing all defined subsets...")

for (subset_name in names(CLUSTER_DEFINITIONS)) {
  cluster_ids <- CLUSTER_DEFINITIONS[[subset_name]]
  
  for (gene_type in GENE_SET_TYPES) {
    
    # Call the processing function for the current combination
    preprocess_subset(
      counts_input_3d = counts_filtered,
      cluster_subset_ids = cluster_ids,
      subset_name = subset_name,
      gene_set_type = gene_type,
      hsa01100_genes_df = hsa01100_genes,
      annots_list_df = annots_list,
      cell_metadata_df = cell_metadata,
      table_s1_df = table_s1,
      base_output_path = BASE_OUTPUT_DIR
    )
    
  } 
}

message("\n All preprocessing complete!")
