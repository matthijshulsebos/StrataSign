library(dplyr)
library(readr)


# Load and cache the main dataset
load_lung_data <- function() {
  if (!exists("lung_ldm", envir = .GlobalEnv)) {
    message("Loading lung dataset...")
    
    load("base/data/lung_ldm.rd", envir = .GlobalEnv)
    message("Lung dataset loaded successfully.")
  } else {
    message("Using cached lung dataset.")
  }
  return(get("lung_ldm", envir = .GlobalEnv))
}


# Load sample metadata
load_table_s1 <- function() {
  message("Loading table_s1...")

  table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", 
                       show_col_types = FALSE) %>% 
    mutate(sample_ID = as.character(sample_ID))
  return(table_s1)
}


# Load annotations list
load_annotations <- function() {
  message("Loading annotations...")
  annots_list <- read_csv("base/input_tables/annots_list.csv", 
                         show_col_types = FALSE)
  return(annots_list)
}


# Load metabolic genes
load_metabolic_genes <- function() {
  message("Loading metabolic genes...")
  hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", 
                            show_col_types = FALSE)
  return(hsa01100_genes)
}

# Load all datasets
load_all_datasets <- function() {
  message("Loading all datasets...")

  lung_ldm <- load_lung_data()
  table_s1 <- load_table_s1()
  annots_list <- load_annotations()
  hsa01100_genes <- load_metabolic_genes()
  
  return(list(
    lung_ldm = lung_ldm,
    table_s1 = table_s1,
    annots_list = annots_list,
    hsa01100_genes = hsa01100_genes
  ))
}


# Prepare cell metadata for CP10K methods
prepare_cell_metadata <- function(lung_ldm, table_s1, doublets) {
  message("Preparing cell metadata...")
  
  # Create cell to sample dataframe
  cell_to_sample_df <- data.frame(
    cell_ID = names(lung_ldm$dataset$cell_to_sample),
    sample_ID = unname(lung_ldm$dataset$cell_to_sample),
    stringsAsFactors = FALSE
  )
  
  # Create cell to cluster dataframe
  cell_to_cluster_df <- data.frame(
    cell_ID = names(lung_ldm$dataset$cell_to_cluster),
    cluster_ID = as.numeric(unname(lung_ldm$dataset$cell_to_cluster)), 
    stringsAsFactors = FALSE
  )
  
  # Join cell metadata with sample and cluster information
  cell_metadata_df <- inner_join(cell_to_sample_df, cell_to_cluster_df, by = "cell_ID")
  
  # Select samples from table_s1
  sample_ids_to_keep <- table_s1 %>% pull(sample_ID) %>% unique()
  
  # Remove doublets and filter samples
  cell_metadata_filtered_df <- cell_metadata_df %>%
    filter(!cluster_ID %in% doublets) %>%
    filter(sample_ID %in% sample_ids_to_keep)
  
  return(cell_metadata_filtered_df)
}


# Prepare 3D counts array for aggregated methods
prepare_counts_array <- function(lung_ldm, table_s1, doublets) {
  message("Preparing counts array for aggregated methods...")
  
  # Extract counts from leader dataset
  counts <- lung_ldm$dataset$counts
  
  # Filter out doublets
  all_clusters_initial <- as.numeric(dimnames(counts)[[3]])
  clusters_to_keep_initial <- all_clusters_initial[!all_clusters_initial %in% doublets]
  counts_filtered <- counts[, , as.character(clusters_to_keep_initial), drop = FALSE]
  
  # Filter samples
  sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
  counts_filtered <- counts_filtered[dimnames(counts_filtered)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]
  
  return(counts_filtered)
}


# Filter umitab for CP10K methods
filter_umitab <- function(umitab, cell_metadata_filtered_df) {
  message("Filtering UMI table using cell metadata...")
  
  # Select cell IDs from metadata
  cells_to_keep <- cell_metadata_filtered_df$cell_ID
  message(sprintf("Keeping %d cells out of %d total", 
                  length(cells_to_keep), ncol(umitab)))
  
  # Filter umitab columns based on cell ids and reorder to match metadata
  umitab_filtered <- umitab[, colnames(umitab) %in% cells_to_keep, drop = FALSE]
  umitab_filtered <- umitab_filtered[, cell_metadata_filtered_df$cell_ID, drop = FALSE]
  
  message(sprintf("Filtered UMI table dimensions: %d genes x %d cells", 
                  nrow(umitab_filtered), ncol(umitab_filtered)))
  
  return(umitab_filtered)
}


# Convert 3D array to long format for aggregated methods
array_to_long <- function(counts_filtered) {
  message("Converting cell type aggregated array to long format...")
  
  # Convert 3D array to long format
  counts_long <- as.data.frame(as.table(counts_filtered))

  # Set column names
  colnames(counts_long) <- c("sample_ID", "gene", "cluster_ID", "count")
  
  # Set data types
  counts_long <- counts_long %>%
    mutate(
      sample_ID = as.character(sample_ID),
      cluster_ID = as.numeric(as.character(cluster_ID)),
      gene = as.character(gene)
    )
  
  return(counts_long)
}


# Unified data preparation function for any normalization method
prepare_data_for_normalization <- function(normalization_method, lung_ldm, table_s1, doublets) {
  base_normalization_method <- gsub("_relative$", "", normalization_method)
  
  # If cp10K then use single cell data
  if (grepl("^cp10k", base_normalization_method)) {
    # Filter out samples and doublets
    cell_metadata_filtered_df <- prepare_cell_metadata(lung_ldm, table_s1, doublets)

    # Filter umitab based on cell metadata df
    umitab_filtered <- filter_umitab(lung_ldm$dataset$umitab, cell_metadata_filtered_df)

    result <- list(
      type = "cell_level",
      cell_metadata = cell_metadata_filtered_df,
      umitab = umitab_filtered
    )
    
  } else {
    # Use cell type aggregated counts data
    counts_filtered <- prepare_counts_array(lung_ldm, table_s1, doublets)

    # Converts the 3D counts dataset to long 2D format
    counts_long <- array_to_long(counts_filtered)

    result <- list(
      type = "cell_type_level",
      counts_long = counts_long
    )
    rm(counts_filtered)
  }
  
  return(result)
}


# Apply log transformation
apply_log_transform <- function(counts_normalized) {
  message("Applying log transformation...")

  counts_normalized %>%
    mutate(log_normalized_count = log1p(normalized_count))
}


# Create feature names
create_features <- function(counts_log_transformed, annots_list) {
  # Create feature identifiers
  counts_log_transformed %>%
    left_join(annots_list, by = c("cluster_ID" = "cluster")) %>% 
    mutate(
      cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
      cluster_name = gsub("/", "-", cluster_name),
      cluster_identifier = paste(cluster_name, cluster_ID, sep = "_"),
      gene_cluster = paste(gene, cluster_identifier, sep = "@") 
    ) %>%
    select(sample_ID, gene, cluster_ID, gene_cluster, log_normalized_count)
}
