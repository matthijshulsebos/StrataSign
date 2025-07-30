library(dplyr)
library(readr)


# Loads the lung data of Leader et al. (2021)
load_lung_data <- function() {
  # Check if we have lung_ldm in the env to save time
  if (!exists("lung_ldm", envir = .GlobalEnv)) {
    # Print because this is very large and takes a while to load
    message("Loading lung dataset.")
    load("base/data/lung_ldm.rd", envir = .GlobalEnv)
  } else {
    # Data is cached in env so use that instead
    message("Using cached lung dataset.")
  }

  # Return the object but you may as well use from env directly
  return(get("lung_ldm", envir = .GlobalEnv))
}


# Loads sample metadata which contains tissue type and histological subtype
load_table_s1 <- function() {
  # We join on sample ID so make it a character
  table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", 
                       show_col_types = FALSE) %>% 
    mutate(sample_ID = as.character(sample_ID))

  return(table_s1)
}


# Loads annotations list containing cluster information on (sub)lineage
load_annotations <- function() {
  annots_list <- read_csv("base/input_tables/annots_list.csv", 
                         show_col_types = FALSE)

  return(annots_list)
}


# Loads the metabolic genes we pulled from KEGG
load_metabolic_genes <- function() {
  hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", 
                            show_col_types = FALSE)

  return(hsa01100_genes)
}


# Loads all datasets for the data preprocessing pipeline
load_all_datasets <- function() {
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


# Joins each cell ID with sample metadata and cluster metadata
prepare_cell_metadata <- function(lung_ldm, table_s1, doublets) {
  # Cell to sample mapping with metadata
  cell_to_sample_df <- data.frame(
    cell_ID = names(lung_ldm$dataset$cell_to_sample),
    sample_ID = unname(lung_ldm$dataset$cell_to_sample),
    stringsAsFactors = FALSE
  )
  
  # Cell to cluster and cell type mapping
  cell_to_cluster_df <- data.frame(
    cell_ID = names(lung_ldm$dataset$cell_to_cluster),
    cluster_ID = as.numeric(unname(lung_ldm$dataset$cell_to_cluster)), 
    stringsAsFactors = FALSE
  )
  
  # Create complete cell metadata df
  cell_metadata_df <- inner_join(cell_to_sample_df, cell_to_cluster_df, by = "cell_ID")
  
  # Filters cells on filtered table_s1
  sample_ids_to_keep <- table_s1 %>% pull(sample_ID) %>% unique()
  
  # Filter out doublets and check if sample IDs are in the filtered table_s1
  cell_metadata_filtered_df <- cell_metadata_df %>%
    filter(!cluster_ID %in% doublets) %>%
    filter(sample_ID %in% sample_ids_to_keep)
  
  return(cell_metadata_filtered_df)
}

# Filters umitab using the cell metadata which is already filtered
filter_umitab <- function(umitab, cell_metadata_filtered_df) {
  # All filtered cells in the metadata
  cells_to_keep <- cell_metadata_filtered_df$cell_ID

  # All cells in the umitab to retain
  available_cells <- intersect(colnames(umitab), cells_to_keep)
  
  # If metadata filters out everything then throw an error
  if (length(available_cells) == 0) {
    stop("No common cells found between umitab and cell metadata.")
  }
  
  # Apply the filter to the umitab matrix
  umitab_filtered <- umitab[, available_cells, drop = FALSE]
  
  return(umitab_filtered)
}


# For the processed data create feature identifiers
create_features <- function(counts_log_transformed, annots_list) {
  # Create feature identifiers
  counts_log_transformed %>%
    # Get cluster and sublineage from annotations list
    left_join(annots_list, by = c("cluster_ID" = "cluster")) %>% 
    # Manually create the feature identifiers but you can also use the feature utils
    mutate(
      cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
      cluster_name = gsub("/", "-", cluster_name),
      cluster_identifier = paste(cluster_name, cluster_ID, sep = "_"),
      gene_cluster = paste(gene, cluster_identifier, sep = "@")
    ) %>%
    select(sample_ID, gene, cluster_ID, gene_cluster, log_normalized_count)
}


# Avoids overrepresentation of patients and selects one sample per patient per tissue type
select_representative_samples <- function(table_s1) {
  # Prefer the clustering model samples
  table_s1 %>%
    group_by(patient_ID, tissue) %>%
    arrange(desc(Use.in.Clustering.Model. == "Yes"), sample_ID) %>%
    slice(1) %>%
    ungroup()
}
