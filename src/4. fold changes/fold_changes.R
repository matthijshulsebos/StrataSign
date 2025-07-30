
# Calculate fold changes for normalized data
calculate_fold_changes_for_normalization <- function(normalized_data, table_s1, annots_list, method_name, cell_type, gene_type) {
  # Load feature name utility functions inside this function because otherwise it wont load in data preprocessing
  source("src/0. utils/feature_name_utils.R")

  # Add tissue information from table_s1
  data_with_tissue <- normalized_data %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
    filter(!is.na(tissue))

  # Get all unique sample IDs for each tissue type
  all_tissues <- unique(data_with_tissue$tissue)

  # All samples stratified by tissue type
  all_samples_by_tissue <- lapply(all_tissues, function(t) unique(data_with_tissue$sample_ID[data_with_tissue$tissue == t]))
  
  # Set names for the list of samples by tissue
  names(all_samples_by_tissue) <- all_tissues

  # Get all unique cluster_IDs and genes
  all_clusters <- unique(data_with_tissue$cluster_ID)
  all_genes <- unique(data_with_tissue$gene)

  # Create a complete grid of gene, cluster_ID, and sample_ID for each tissue
  complete_data <- bind_rows(lapply(all_tissues, function(tissue_type) {
    expand.grid(
      gene = all_genes,
      cluster_ID = all_clusters,
      sample_ID = all_samples_by_tissue[[tissue_type]],
      tissue = tissue_type,
      stringsAsFactors = FALSE
    )
  }))

  # Fill in zeros for missing combinations just to be sure
  data_with_zeros <- complete_data %>%
    left_join(data_with_tissue %>% select(gene, cluster_ID, sample_ID, tissue, normalized_count),
              by = c("gene", "cluster_ID", "sample_ID", "tissue")) %>%
    mutate(normalized_count = ifelse(is.na(normalized_count), 0, normalized_count))

  # Calculate mean expression and nonzero counts by tissue and cluster
  tissue_means <- data_with_zeros %>%
    group_by(gene, cluster_ID, tissue) %>%
    summarise(
      mean_expr = mean(normalized_count),
      n_samples = n(),
      n_nonzero = sum(normalized_count > 0),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      # This creates names like mean_expr_Normal
      names_from = tissue,
      values_from = c(mean_expr, n_samples, n_nonzero),
      values_fill = 0
    )

  # Set tissue column names
  normal_col <- "Normal"
  tumor_col <- "Tumor"

  # Calculate log2 fold changes with pseudocount of 0.0001
  fold_changes <- tissue_means %>%
    mutate(
      # .data can help get columns dynamically from strings
      normal_mean = .data[[paste0("mean_expr_", normal_col)]],
      tumor_mean = .data[[paste0("mean_expr_", tumor_col)]],
      log2FoldChange = log2((tumor_mean + 0.0001) / (normal_mean + 0.0001)),
      n_nonzero_normal = .data[[paste0("n_nonzero_", normal_col)]],
      n_nonzero_tumor = .data[[paste0("n_nonzero_", tumor_col)]],
      n_samples_normal = .data[[paste0("n_samples_", normal_col)]],
      n_samples_tumor = .data[[paste0("n_samples_", tumor_col)]]
    ) %>%
    filter(!is.na(log2FoldChange) & is.finite(log2FoldChange)) %>%
    select(gene, cluster_ID, log2FoldChange, normal_mean, tumor_mean, n_nonzero_normal, n_nonzero_tumor, n_samples_normal, n_samples_tumor)

  # Add cell type annotations and create feature identifiers
  fold_changes_with_features <- fold_changes %>%
    left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
    mutate(
      cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
      cluster_name = gsub("/", "-", cluster_name),
      Feature = create_feature_identifier(gene, cluster_name, cluster_ID)
    ) %>%
    select(Feature, Value = log2FoldChange, gene, cluster_ID, cluster_name, normal_mean, tumor_mean, n_nonzero_normal, n_nonzero_tumor, n_samples_normal, n_samples_tumor)

  # Create output directory and save
  output_dir <- file.path("output", "4. fold changes", method_name, cell_type, gene_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Save main results with Feature and Value columns
  main_results <- fold_changes_with_features %>%
    select(Feature, Value)
  output_path <- file.path(output_dir, "fold_changes.csv")
  write_csv(main_results, output_path)

  # Also save detailed results with all information
  detailed_output_path <- file.path(output_dir, "detailed_fold_changes.csv")
  detailed_results_lower <- fold_changes_with_features
  names(detailed_results_lower) <- tolower(names(detailed_results_lower))
  write_csv(detailed_results_lower, detailed_output_path)

  message(sprintf("Saved fold changes: %s", output_path))

  return(fold_changes_with_features)
}
