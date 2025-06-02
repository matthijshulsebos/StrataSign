library(dplyr)
library(tidyr)
library(readr)

# Filter genes by type
filter_genes_by_type <- function(all_genes, gene_set_type, hsa01100_genes_df) {
  # Select metabolic genes
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
  
  # Prepare metadata for both training and test sets
  metadata_train <- table_s1_df %>%
    filter(sample_ID %in% train_df$sample_ID) %>%
    arrange(match(sample_ID, train_df$sample_ID))
  
  metadata_test <- table_s1_df %>%
    filter(sample_ID %in% test_df$sample_ID) %>%
    arrange(match(sample_ID, test_df$sample_ID))
  
  # Prepare feature datasets by removing metadata
  X_train <- train_df %>% select(-sample_ID, -patient_ID, -tissue)
  X_test <- test_df %>% select(-sample_ID, -patient_ID, -tissue)
  
  # Ensure numeric and replace empty values with 0
  X_train <- X_train %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))
  X_test <- X_test %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))
  
  return(list(
    metadata_train = metadata_train,
    metadata_test = metadata_test,
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
            file.path(output_dir, paste0("metadata_train_", subset_name, suffix, ".csv")))
  write_csv(outputs$metadata_test, 
            file.path(output_dir, paste0("metadata_test_", subset_name, suffix, ".csv")))
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

# Process both global and relative normalization
process_all_subsets_unified <- function(input_data, cluster_definitions, gene_set_types, 
                                       hsa01100_genes, table_s1, base_output_path, 
                                       is_relative = FALSE, normalization_method = NULL, 
                                       annots_list = NULL) {
  
  message("\nProcessing all subsets with unified approach...")
  message(sprintf("Mode: %s", if(is_relative) "Relative" else "Global"))
  
  for (subset_name in names(cluster_definitions)) {
    cluster_ids <- cluster_definitions[[subset_name]]
    
    for (gene_type in gene_set_types) {
      message(sprintf("Processing %s - %s", subset_name, gene_type))
      
      if (is_relative) {
        # Relative: Filter, normalize, post-process, save
        process_relative_subset(input_data, cluster_ids, subset_name, gene_type, 
                               hsa01100_genes, table_s1, base_output_path, 
                               normalization_method, annots_list)
      } else {
        # Filter already normalized data, save
        process_global_subset(input_data, cluster_ids, subset_name, gene_type, 
                             hsa01100_genes, table_s1, base_output_path)
      }
    }
  }
  
  message("All subsets processed!")
}

# Process a single subset for relative normalization
process_relative_subset <- function(base_data, cluster_ids, subset_name, gene_type, 
                                   hsa01100_genes, table_s1, base_output_path, 
                                   normalization_method, annots_list) {
  
  # Apply subset filters to base data
  filtered_data <- apply_subset_filters(base_data, cluster_ids, gene_type, hsa01100_genes)
  
  if (is.null(filtered_data)) {
    warning(sprintf("Subset %s - %s resulted in no data after filtering", subset_name, gene_type))
    return()
  }
  
  # Apply normalization to filtered data
  base_normalization_method <- gsub("_relative$", "", normalization_method)
  
  if (!exists("DOUBLETS")) {
    DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
  }
  
  # Use data type to determine pipeline path (raw umitab vs counts)
  if (filtered_data$type == "cell_level") {
    # Use CP10K normalization functions that work with umitab
    if (base_normalization_method == "cp10k") {
      counts_normalized <- normalize_cp10k(filtered_data$cell_metadata, filtered_data$umitab, NULL, DOUBLETS)
    } else if (base_normalization_method == "cp10k_ctnorm") {
      subset_normalized <- normalize_cp10k(filtered_data$cell_metadata, filtered_data$umitab, NULL, DOUBLETS)
      counts_normalized <- apply_celltype_normalization(subset_normalized)
    }
    
  } else if (filtered_data$type == "aggregated") {
    # Use proportion/sample depth functions that work with aggregated counts
    if (base_normalization_method == "proportion_ctnorm") {
      counts_normalized <- apply_proportion_celltype_normalization(filtered_data$counts_long)
    } else if (base_normalization_method == "sample_depth") {
      counts_normalized <- normalize_sample_depth(filtered_data$counts_long)
    }
  }
  
  # Apply post processing
  counts_processed <- apply_subset_postprocessing(counts_normalized, annots_list)
  
  # Save results
  gene_filter_result <- filter_genes_by_type(unique(counts_processed$gene), gene_type, hsa01100_genes)
  suffix <- gene_filter_result$suffix
  save_subset_results(counts_processed, subset_name, gene_type, suffix, table_s1, base_output_path)
}

# Process a single subset for global normalization (existing functionality)
process_global_subset <- function(globally_processed_data_long, cluster_ids, subset_name, gene_type, 
                                 hsa01100_genes, table_s1, base_output_path) {
  
  # Use existing process_and_save_subset function
  process_and_save_subset(globally_processed_data_long, cluster_ids, subset_name, gene_type, 
                         hsa01100_genes, table_s1, base_output_path)
}

# Apply filters to base data for a specific subset
apply_subset_filters <- function(base_data, cluster_ids, gene_type, hsa01100_genes) {
  
  if (base_data$type == "cell_level") {
    # Filter cells by cluster
    if (!is.null(cluster_ids)) {
      filtered_cell_metadata <- base_data$cell_metadata %>% filter(cluster_ID %in% cluster_ids)
    } else {
      filtered_cell_metadata <- base_data$cell_metadata
    }
    
    # Filter genes by type
    all_genes <- rownames(base_data$umitab)
    gene_filter_result <- filter_genes_by_type(all_genes, gene_type, hsa01100_genes)
    genes_to_keep <- gene_filter_result$genes
    
    if (length(genes_to_keep) == 0 || nrow(filtered_cell_metadata) == 0) {
      return(NULL)
    }
    
    # Filter umitabs by genes
    subset_umitab <- base_data$umitab[genes_to_keep, , drop = FALSE]
    
    return(list(
      type = "cell_level",
      cell_metadata = filtered_cell_metadata,
      umitab = subset_umitab
    ))
    
  } else if (base_data$type == "aggregated") {
    # Filter by cluster subset
    if (!is.null(cluster_ids)) {
      subset_data_long <- base_data$counts_long %>% filter(cluster_ID %in% cluster_ids)
    } else {
      subset_data_long <- base_data$counts_long
    }
    
    # Filter genes by type
    all_genes_in_subset <- unique(subset_data_long$gene)
    gene_filter_result <- filter_genes_by_type(all_genes_in_subset, gene_type, hsa01100_genes)
    genes_to_keep <- gene_filter_result$genes
    
    # Apply gene filter
    subset_data_long <- subset_data_long %>% filter(gene %in% genes_to_keep)
    
    if (nrow(subset_data_long) == 0) {
      return(NULL)
    }
    
    return(list(
      type = "aggregated",
      counts_long = subset_data_long
    ))
  }
  
  return(NULL)
}

# Apply post processing to normalized subset data
apply_subset_postprocessing <- function(counts_normalized, annots_list) {
  counts_log_transformed <- apply_log_transform(counts_normalized)
  counts_processed_long <- create_features(counts_log_transformed, annots_list)
  return(counts_processed_long)
}

# Save subset results
save_subset_results <- function(counts_processed_long, subset_name, gene_type, suffix, table_s1, base_output_path) {
  
  # Convert to wide format
  counts_wide <- counts_processed_long %>%
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster, values_from = log_normalized_count, values_fill = 0)
  
  # Merge with metadata
  counts_wide_meta <- counts_wide %>%
    left_join(table_s1 %>% select(sample_ID, patient_ID, tissue), by = "sample_ID") %>%
    filter(!is.na(tissue))
  
  # Save if data exists
  if (nrow(counts_wide_meta) > 0) {
    split_data <- split_train_test(counts_wide_meta)
    if (nrow(split_data$train) > 0 && nrow(split_data$test) > 0) {
      outputs <- prepare_outputs(split_data$train, split_data$test, table_s1)
      output_dir <- file.path(base_output_path, subset_name, gene_type)
      save_outputs(outputs, output_dir, subset_name, suffix)
    }
  }
}
