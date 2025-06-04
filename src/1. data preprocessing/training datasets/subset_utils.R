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


# Unified function to save subset results
save_subset_results <- function(counts_processed_long, subset_name, gene_type, table_s1, base_output_path) {
  
  # Determine suffix based on gene type
  suffix <- switch(gene_type,
                   "metabolic" = "_metabolic",
                   "nonmetabolic" = "_nonmetabolic", 
                   "random" = "_random")
  
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
      write_outputs(outputs, output_dir, subset_name, suffix)
    }
  }
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
write_outputs <- function(outputs, output_dir, subset_name, suffix) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_csv(outputs$metadata_train, 
            file.path(output_dir, paste0("metadata_train_", subset_name, suffix, ".csv")))
  write_csv(outputs$metadata_test, 
            file.path(output_dir, paste0("metadata_test_", subset_name, suffix, ".csv")))
  
  write_csv(outputs$X_train, 
            file.path(output_dir, paste0("X_train_", subset_name, suffix, ".csv")))
  write_csv(outputs$X_test, 
            file.path(output_dir, paste0("X_test_", subset_name, suffix, ".csv")))
  
  write_csv(data.frame(y = outputs$y_train), 
            file.path(output_dir, paste0("y_train_", subset_name, suffix, ".csv")))
  write_csv(data.frame(y = outputs$y_test), 
            file.path(output_dir, paste0("y_test_", subset_name, suffix, ".csv")))
}


# Main subset processing function
create_and_save_subset <- function(data_to_filter_long,
                                    cluster_subset_ids,
                                    subset_name,
                                    gene_set_type,
                                    hsa01100_genes_df,
                                    table_s1_df,
                                    base_output_path,
                                    annots_list_param) {

  message(paste("Processing Subset:", subset_name, "| Gene Set:", gene_set_type))

  # Filter by cluster subset first
  if (!is.null(cluster_subset_ids)) {
    subset_data_long_filtered <- data_to_filter_long %>% filter(cluster_ID %in% cluster_subset_ids)
  } else {
    subset_data_long_filtered <- data_to_filter_long
  }

  # Filter by gene set type
  all_genes_in_subset <- unique(subset_data_long_filtered$gene)
  gene_filter_result <- filter_genes_by_type(all_genes_in_subset, gene_set_type, hsa01100_genes_df)
  genes_to_keep <- gene_filter_result$genes

  # Apply gene filter
  subset_data_long_filtered <- subset_data_long_filtered %>% filter(gene %in% genes_to_keep)

  # Check if filtering resulted in empty data
  if (nrow(subset_data_long_filtered) == 0) {
    warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after filtering."))
    return()
  }

  # Create features after filtering
  counts_processed_long_with_features <- create_features(subset_data_long_filtered, annots_list_param)

  # Use the unified save function
  save_subset_results(counts_processed_long_with_features, subset_name, gene_set_type, table_s1_df, base_output_path)

  message(paste("Finished:", subset_name, "|", gene_set_type))
}


# Process both global and relative normalization
process_all_subsets <- function(input_data, cluster_definitions, gene_set_types, 
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
        # Filter, normalize, post process, save
        process_relative_subset(input_data, cluster_ids, subset_name, gene_type, 
                               hsa01100_genes, table_s1, base_output_path, 
                               normalization_method, annots_list)
      } else {
        # Filter already normalized data, save
        create_and_save_subset(input_data, cluster_ids, subset_name, gene_type, 
                               hsa01100_genes, table_s1, base_output_path, 
                               annots_list_param = annots_list)
      }
    }
  }
  
  message("All subsets processed!")
}


# Process a single subset for relative normalization
process_relative_subset <- function(prepared_data, cluster_ids, subset_name, gene_type, 
                                   hsa01100_genes, table_s1, base_output_path, 
                                   normalization_method, annots_list) {
  
  # Apply subset filters to base data
  filtered_data <- apply_subset_filters(prepared_data, cluster_ids, gene_type, hsa01100_genes)
  
  if (is.null(filtered_data)) {
    warning(sprintf("Subset %s - %s resulted in no data after filtering", subset_name, gene_type))
    return()
  }
  
  # Apply normalization to filtered data using the same function as global processing
  base_normalization_method <- gsub("_relative$", "", normalization_method)
  counts_normalized <- apply_normalization_method(base_normalization_method, filtered_data, NULL)
  
  # Apply log1p transformation
  counts_log_transformed <- apply_log_transform(counts_normalized)

  # Create feature identifiers
  counts_processed <- create_features(counts_log_transformed, annots_list)
  
  # Save results
  save_subset_results(counts_processed, subset_name, gene_type, table_s1, base_output_path)
}


# Apply filters to base data for a specific subset
apply_subset_filters <- function(prepared_data, cluster_ids, gene_type, hsa01100_genes) {
  
  if (prepared_data$type == "cell_level") {
    # Filter cell metadata by clusters
    if (!is.null(cluster_ids)) {
      filtered_cell_metadata <- prepared_data$cell_metadata %>% filter(cluster_ID %in% cluster_ids)
    } else {
      filtered_cell_metadata <- prepared_data$cell_metadata
    }
    
    # Get all available genes before filtering by cells
    all_genes <- rownames(prepared_data$umitab)
    gene_filter_result <- filter_genes_by_type(all_genes, gene_type, hsa01100_genes)
    genes_to_keep <- gene_filter_result$genes
    
    if (length(genes_to_keep) == 0 || nrow(filtered_cell_metadata) == 0) {
      return(NULL)
    }
    
    # Use filter_umitab function to filter by cell metadata
    subset_umitab <- filter_umitab(prepared_data$umitab, filtered_cell_metadata)
    
    # Then filter by genes
    subset_umitab <- subset_umitab[genes_to_keep, , drop = FALSE]
    
    return(list(
      type = "cell_level",
      cell_metadata = filtered_cell_metadata,
      umitab = subset_umitab
    ))
    
  } else if (prepared_data$type == "cell_type_level") {
    # Filter by cluster subset
    if (!is.null(cluster_ids)) {
      subset_data_long <- prepared_data$counts_long %>% filter(cluster_ID %in% cluster_ids)
    } else {
      subset_data_long <- prepared_data$counts_long
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
      type = "cell_type_level",
      counts_long = subset_data_long
    ))
  }
  
  return(NULL)
}
