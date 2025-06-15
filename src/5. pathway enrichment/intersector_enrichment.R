library(dplyr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fs)
library(stringr)

# Utils for mapping gene symbols to entrez ids
source("src/0. utils/format_utils.R") 

# Set paths and parameters
INTERSECTOR_PARENT_INPUT_DIR <- "output/3. intersector"
PATHWAY_MAPPING_OUTPUT_PARENT_DIR <- "output/5. pathway enrichment"
DATASET_TYPES_TO_PROCESS <- c("ctnorm_global", "ctnorm_relative", "read_depth")
KEGG_ORGANISM_CODE <- "hsa"
METABOLIC_PATHWAYS_FILE <- "output/1. data preprocessing/kegg/hsa01100_rel_pathways.csv"

# Load metabolic pathways from the CSV file
load_metabolic_pathways <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("Metabolic pathways file not found:", file_path))
  }
  metabolic_pathways <- read.csv(file_path)$PathwayID
  return(metabolic_pathways)
}

# Filter KEGG results to include only metabolic pathways
filter_metabolic_pathways <- function(kegg_result_df, metabolic_pathways) {
  kegg_result_df <- kegg_result_df %>%
    dplyr::filter(PathwayID %in% metabolic_pathways)
  return(kegg_result_df)
}

# Parse KEGG output format to get gene symbols
map_kegg_geneid_to_symbols <- function(gene_id_str, entrez_to_symbol_map) {
  entrez_ids <- stringr::str_split(gene_id_str, "/")[[1]]
  symbols <- entrez_to_symbol_map[entrez_ids]
  symbols <- ifelse(is.na(symbols), entrez_ids, symbols)
  paste(na.omit(symbols), collapse="/")
}

# Map gene symbols to entrez ids
process_kegg_results <- function(kegg_result_df, entrez_to_symbol_map, meta_scores_dt = NULL) {
  kegg_result_df$InputGeneSymbolsInPathway <- sapply(
    kegg_result_df$geneID, map_kegg_geneid_to_symbols, entrez_to_symbol_map
  )

  kegg_result_df$GeneCountFromInput <- sapply(
    stringr::str_split(kegg_result_df$geneID, "/"), length
  )

  if (!is.null(meta_scores_dt)) {
    kegg_result_df$ContributingCellTypes <- sapply(
      kegg_result_df$InputGeneSymbolsInPathway, function(gs_string) {
        current_pathway_gene_symbols <- stringr::str_split(gs_string, "/")[[1]]
        paste(
          sort(unique(meta_scores_dt[gene %in% current_pathway_gene_symbols, sublineage])),
          collapse="/"
        )
      }
    )

    kegg_result_df$NumberOfContributingCellTypes <- sapply(
      kegg_result_df$ContributingCellTypes, function(ct_string)
        length(stringr::str_split(ct_string, "/")[[1]])
    )

    kegg_result_df <- kegg_result_df %>%
      dplyr::select(
        PathwayID = ID,
        PathwayName = Description,
        GeneCountFromInput,
        NumberOfContributingCellTypes,
        InputGeneSymbolsInPathway,
        ContributingCellTypes
      ) %>%
      dplyr::arrange(desc(GeneCountFromInput), desc(NumberOfContributingCellTypes))

  } else {
    kegg_result_df <- kegg_result_df %>%
      dplyr::select(
        PathwayID = ID,
        PathwayName = Description,
        GeneCountFromInput,
        InputGeneSymbolsInPathway
      ) %>%
      dplyr::arrange(desc(GeneCountFromInput))
  }

  kegg_result_df
}

# Retrieve all enriched pathways from KEGG
perform_kegg_mapping <- function(entrez_vector) {
  enrichKEGG(
    gene         = entrez_vector,
    organism     = KEGG_ORGANISM_CODE,
    pvalueCutoff = 1.0,
    qvalueCutoff = 1.0,
    minGSSize    = 1
  )
}

# Load metabolic pathways
metabolic_pathways <- load_metabolic_pathways(METABOLIC_PATHWAYS_FILE)

# Set output dir
dir_create(PATHWAY_MAPPING_OUTPUT_PARENT_DIR, recurse = TRUE)

# Main loop to process each dataset type (normalization method)
for (current_dataset_type in DATASET_TYPES_TO_PROCESS) {
  message(paste("Processing dataset type:", current_dataset_type))
  
  dataset_type_dir <- file.path(INTERSECTOR_PARENT_INPUT_DIR, current_dataset_type)
  if (!dir.exists(dataset_type_dir)) {
    warning(paste("Dataset type directory not found:", dataset_type_dir))
    next
  }
  
  # Get all cell types for this dataset type
  cell_types <- list.dirs(dataset_type_dir, full.names = FALSE, recursive = FALSE)
  if (length(cell_types) == 0) {
    warning(paste("No cell types found in:", dataset_type_dir))
    next
  }

  # Process each cell type
  for (current_cell_type in cell_types) {
    message(paste("  Processing cell type:", current_cell_type))
    
    cell_type_dir <- file.path(dataset_type_dir, current_cell_type)
    
    # Get all gene types for this cell type
    gene_types <- list.dirs(cell_type_dir, full.names = FALSE, recursive = FALSE)
    if (length(gene_types) == 0) {
      warning(paste("No gene types found in:", cell_type_dir))
      next
    }
    
    # Process each gene type
    for (current_gene_type in gene_types) {
      message(paste("    Processing gene type:", current_gene_type))
      
      gene_type_dir <- file.path(cell_type_dir, current_gene_type)
      meta_scores_file_path <- file.path(gene_type_dir, "meta_scores.csv")
      
      if (!file.exists(meta_scores_file_path)) {
        warning(paste("Meta scores file not found:", meta_scores_file_path))
        next
      }
      
      # Read meta scores
      meta_scores_dt <- fread(meta_scores_file_path)
      
      if (nrow(meta_scores_dt) == 0) {
        warning(paste("Empty meta scores file:", meta_scores_file_path))
        next
      }
      
      # Map genes to Entrez IDs
      genes_for_combination <- unique(meta_scores_dt$gene)
      entrez_df <- map_symbols_to_entrez(genes_for_combination)
      entrez_vector <- unique(na.omit(entrez_df$ENTREZID))
      symbol_entrez_map <- setNames(entrez_df$SYMBOL, entrez_df$ENTREZID)
      
      if (length(entrez_vector) == 0) {
        warning(paste("No valid Entrez IDs found for combination:", current_dataset_type, current_cell_type, current_gene_type))
        next
      }
      
      # Perform KEGG mapping
      kegg_map_obj <- perform_kegg_mapping(entrez_vector)
      
      if (!is.null(kegg_map_obj) && nrow(as.data.frame(kegg_map_obj)) > 0) {
        processed_results <- process_kegg_results(
          as.data.frame(kegg_map_obj),
          symbol_entrez_map,
          meta_scores_dt
        )
        
        # Filter for metabolic pathways
        processed_results <- filter_metabolic_pathways(processed_results, metabolic_pathways)
        
        # Create output directory structure: dataset_type/cell_type/gene_type/
        combination_output_dir <- file.path(
          PATHWAY_MAPPING_OUTPUT_PARENT_DIR, 
          current_dataset_type, 
          current_cell_type, 
          current_gene_type
        )
        dir_create(combination_output_dir, recurse = TRUE)
        
        # Save results
        output_file <- file.path(
          combination_output_dir,
          paste0("kegg_pathway_enrichment_", current_dataset_type, "_", current_cell_type, "_", current_gene_type, ".csv")
        )
        write.csv(processed_results, output_file, row.names = FALSE)
        
        message(paste("      Saved results to:", output_file))
      } else {
        message(paste("      No KEGG enrichment results for combination:", current_dataset_type, current_cell_type, current_gene_type))
      }
    }
  }
}

message("Pathway enrichment analysis complete. Results saved in:", PATHWAY_MAPPING_OUTPUT_PARENT_DIR)
