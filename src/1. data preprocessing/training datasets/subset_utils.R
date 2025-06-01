library(dplyr)
library(tidyr)
library(readr)

# Filter genes by type
filter_genes_by_type <- function(all_genes, gene_set_type, hsa01100_genes_df) {
  # Select metabolic genes that occur in gene set
  metabolic_genes <- hsa01100_genes_df$SYMBOL[hsa01100_genes_df$SYMBOL %in% all_genes]
  
  if (gene_set_type == "metabolic") {
    genes_to_keep <- metabolic_genes
    suffix <- "_metabolic"
  } else if (gene_set_type == "nonmetabolic") {
    nonmetabolic_genes <- setdiff(all_genes, metabolic_genes)
    if (length(nonmetabolic_genes) > length(metabolic_genes)) {
      set.seed(42)
      genes_to_keep <- sample(nonmetabolic_genes, length(metabolic_genes))
    } else {
      genes_to_keep <- nonmetabolic_genes
    }
    suffix <- "_nonmetabolic"
  } else if (gene_set_type == "random") {
    set.seed(43)
    sample_size <- min(length(metabolic_genes), length(all_genes))
    genes_to_keep <- sample(all_genes, sample_size)
    suffix <- "_random"
  } else {
    stop("Invalid gene set type.")
  }
  
  return(list(genes = genes_to_keep, suffix = suffix))
}

# Split data into train test sets
split_train_test <- function(counts_wide_meta, train_ratio = 0.7) {
  set.seed(42)

  # Get patient IDs
  patients <- unique(counts_wide_meta$patient_ID)

  # Sample patients for training set
  train_patients <- sample(patients, size = round(length(patients) * train_ratio))
  
  # Split data into train and test based on patient IDs
  train_df <- counts_wide_meta %>% filter(patient_ID %in% train_patients)
  test_df <- counts_wide_meta %>% filter(!patient_ID %in% train_patients)
  
  return(list(train = train_df, test = test_df))
}

# Prepare final outputs for saving
prepare_outputs <- function(train_df, test_df, table_s1_df) {
  y_train <- train_df$tissue
  y_test <- test_df$tissue
  
  # Prepare metadata for saving
  metadata_train <- table_s1_df %>%
    filter(sample_ID %in% train_df$sample_ID) %>%
    arrange(match(sample_ID, train_df$sample_ID))
  
  # Prepare feature datasets by removing metadata
  X_train <- train_df %>% select(-sample_ID, -patient_ID, -tissue)
  X_test <- test_df %>% select(-sample_ID, -patient_ID, -tissue)
  
  # Ensure numeric and replace empty values with 0
  X_train <- X_train %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))
  X_test <- X_test %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))
  
  return(list(
    metadata_train = metadata_train,
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test
  ))
}

# Save all outputs to files
save_outputs <- function(outputs, output_dir, subset_name, suffix) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_csv(outputs$metadata_train, 
            file.path(output_dir, paste0("metadata_", subset_name, suffix, ".csv")))
  write_csv(outputs$X_train, 
            file.path(output_dir, paste0("X_train_", subset_name, suffix, ".csv")))
  write_csv(outputs$X_test, 
            file.path(output_dir, paste0("X_test_", subset_name, suffix, ".csv")))
  write_csv(data.frame(x = outputs$y_train), 
            file.path(output_dir, paste0("y_train_", subset_name, suffix, ".csv")))
  write_csv(data.frame(x = outputs$y_test), 
            file.path(output_dir, paste0("y_test_", subset_name, suffix, ".csv")))
}

# Main subset processing function
process_and_save_subset <- function(globally_processed_data_long, 
                                    cluster_subset_ids,
                                    subset_name,
                                    gene_set_type,
                                    hsa01100_genes_df,
                                    table_s1_df,
                                    base_output_path) {

  message(paste("Processing Subset:", subset_name, "| Gene Set:", gene_set_type))

  # Filter by cluster subset
  if (!is.null(cluster_subset_ids)) {
    subset_data_long <- globally_processed_data_long %>% filter(cluster_ID %in% cluster_subset_ids)
  } else {
    subset_data_long <- globally_processed_data_long
  }

  # Filter by gene set type
  all_genes_in_subset <- unique(subset_data_long$gene)
  gene_filter_result <- filter_genes_by_type(all_genes_in_subset, gene_set_type, hsa01100_genes_df)
  genes_to_keep <- gene_filter_result$genes
  suffix <- gene_filter_result$suffix

  # Apply gene filter
  subset_data_long <- subset_data_long %>% filter(gene %in% genes_to_keep)

  # Check if filtering resulted in empty data
  if (nrow(subset_data_long) == 0) {
    warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after filtering."))
    return()
  }

  # Pivot to wide format
  message("Pivoting data to wide format...")
  counts_wide <- subset_data_long %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster,
                values_from = log_normalized_count,
                values_fill = 0)

  # Merge metadata and filter
  counts_wide_meta <- counts_wide %>%
    left_join(
      table_s1_df %>% select(sample_ID, patient_ID, tissue),
      by = "sample_ID"
    ) %>%
    filter(!is.na(tissue))

  # Check if merge resulted in empty data
  if (nrow(counts_wide_meta) == 0) {
    warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after merging metadata/filtering NAs."))
    return()
  }

  # Split into train test
  split_data <- split_train_test(counts_wide_meta)
  
  # Check if split failed
  if (nrow(split_data$train) == 0 || nrow(split_data$test) == 0) {
     warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in empty train/test set."))
     return()
  }

  # Prepare outputs
  outputs <- prepare_outputs(split_data$train, split_data$test, table_s1_df)

  # Save outputs
  output_dir <- file.path(base_output_path, subset_name, gene_set_type)
  save_outputs(outputs, output_dir, subset_name, suffix)

  message(paste("Finished:", subset_name, "|", gene_set_type))
}

# Process all subsets for given definitions
process_all_subsets <- function(counts_processed_long, 
                                cluster_definitions,
                                gene_set_types,
                                hsa01100_genes,
                                table_s1,
                                base_output_path) {
  
  message("\nPreprocessing all defined subsets...")
  
  for (subset_name in names(cluster_definitions)) {
    cluster_ids <- cluster_definitions[[subset_name]]
    
    for (gene_type in gene_set_types) {
      process_and_save_subset(
        globally_processed_data_long = counts_processed_long, 
        cluster_subset_ids = cluster_ids,
        subset_name = subset_name,
        gene_set_type = gene_type,
        hsa01100_genes_df = hsa01100_genes,
        table_s1_df = table_s1,
        base_output_path = base_output_path
      )
    }
  }
  
  message("\nAll preprocessing complete!")
}
