library(dplyr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fs)
library(stringr)

source("src/0. utils/format_utils.R") 

INTERSECTOR_PARENT_INPUT_DIR <- "output/3. intersector"
PATHWAY_MAPPING_OUTPUT_PARENT_DIR <- "output/5. pathway mapping"
KEGG_ORGANISM_CODE <- "hsa"
METABOLIC_PATHWAYS_FILE <- "output/1. data preprocessing/kegg/hsa01100_rel_pathways.csv"
METABOLIC_PATHWAYS <- read.csv(METABOLIC_PATHWAYS_FILE)$PathwayID


# Parse KEGG output format to get gene symbols
map_kegg_geneid_to_symbols <- function(gene_id_str, entrez_to_symbol_map) {
  entrez_ids <- str_split(gene_id_str, "/")[[1]]
  symbols <- entrez_to_symbol_map[entrez_ids]
  symbols <- ifelse(is.na(symbols), entrez_ids, symbols)
  paste(na.omit(symbols), collapse="/")
}


# Parse KEGG result
process_kegg_results <- function(kegg_result_df, entrez_to_symbol_map, meta_scores_dt = NULL) {
  # These columns are based on what the KEGG API returns
  kegg_result_df$InputGeneSymbolsInPathway <- sapply(
    kegg_result_df$geneID, map_kegg_geneid_to_symbols, entrez_to_symbol_map
  )

  # The number of genes that are part of the pathway
  kegg_result_df$GeneCountFromInput <- sapply(
    str_split(kegg_result_df$geneID, "/"), length
  )

  # Get all unique cell types in the data
  all_cell_types <- unique(meta_scores_dt$sublineage)
  
  # For each pathway count genes per cell type
  for (cell_type in all_cell_types) {
    # Replaces any special characters in cell type names with underscores
    col_name <- paste0("Genes_", gsub("[^A-Za-z0-9]", "_", cell_type))
    
    kegg_result_df[[col_name]] <- sapply(
      kegg_result_df$InputGeneSymbolsInPathway, function(gs_string) {
        current_pathway_gene_symbols <- str_split(gs_string, "/")[[1]]
        # Count genes from this pathway that belong to this cell type
        sum(meta_scores_dt[gene %in% current_pathway_gene_symbols & sublineage == cell_type, .N])
      }
    )
  }
  
  # Calculate total number of contributing cell types
  kegg_result_df$NumberOfContributingCellTypes <- sapply(
    kegg_result_df$InputGeneSymbolsInPathway, function(gs_string) {
      current_pathway_gene_symbols <- str_split(gs_string, "/")[[1]]
      length(unique(meta_scores_dt[gene %in% current_pathway_gene_symbols, sublineage]))
    }
  )
  
  # Calculate total number features contributing to each pathway
  kegg_result_df$TotalGeneCellTypeCombinations <- sapply(
    kegg_result_df$InputGeneSymbolsInPathway, function(gs_string) {
      current_pathway_gene_symbols <- str_split(gs_string, "/")[[1]]
      # Count number of features
      nrow(meta_scores_dt[gene %in% current_pathway_gene_symbols])
    }
  )

  # Calculate average and sum meta score of contributing genes
  pathway_scores_list <- lapply(
    kegg_result_df$InputGeneSymbolsInPathway, function(gs_string) {
      current_pathway_gene_symbols <- str_split(gs_string, "/")[[1]]
      meta_scores_dt[gene %in% current_pathway_gene_symbols, meta_score]
    }
  )

  # vapply is very efficient apparently for mean gene scores
  kegg_result_df$AvgMetaScore <- vapply(
    pathway_scores_list,
    # Strange format that take an anonymous function and then applies it
    function(scores) mean(scores, na.rm = TRUE),
    # This is the return type so a single numeric value
    numeric(1)
  )

  # vapply for the summed gene scores as well
  kegg_result_df$SumMetaScore <- vapply(
    pathway_scores_list,
    function(scores) sum(scores, na.rm = TRUE),
    numeric(1)
  )

  # Create dynamic column selection including all cell type columns
  base_columns <- c("PathwayID", "PathwayName", "GeneCountFromInput", "NumberOfContributingCellTypes", 
                   "TotalGeneCellTypeCombinations", "InputGeneSymbolsInPathway", "AvgMetaScore", "SumMetaScore")

  # If it starts with Genes_ then we want to follow it with cell type names in order to avoid typing everything out
  cell_type_columns <- names(kegg_result_df)[grepl("^Genes_", names(kegg_result_df))]
  all_columns <- c(base_columns, cell_type_columns)
  
  kegg_result_df <- kegg_result_df %>%
    dplyr::select(
      PathwayID = ID,
      PathwayName = Description,
      GeneCountFromInput,
      NumberOfContributingCellTypes,
      TotalGeneCellTypeCombinations,
      InputGeneSymbolsInPathway,
      AvgMetaScore,
      SumMetaScore,
      dplyr::all_of(cell_type_columns)
    ) %>%
    dplyr::arrange(desc(GeneCountFromInput), desc(NumberOfContributingCellTypes))

  kegg_result_df
}


# Retrieve all enriched pathways from KEGG
perform_kegg_mapping <- function(entrez_vector) {
  # Dont do any statistical tests
  enrichKEGG(
    gene         = entrez_vector,
    organism     = KEGG_ORGANISM_CODE,
    pvalueCutoff = 1.0,
    qvalueCutoff = 1.0,
    minGSSize    = 1
  )
}


# Set output dir
dir_create(PATHWAY_MAPPING_OUTPUT_PARENT_DIR, recurse = TRUE)

# Get all normalization strategies from the intersector dir
normalization_methods <- list.dirs(INTERSECTOR_PARENT_INPUT_DIR, full.names = FALSE, recursive = FALSE)

if (length(normalization_methods) == 0) {
  stop("No normalization strategies found in: ", INTERSECTOR_PARENT_INPUT_DIR)
}

# Process each normalization strategy
for (normalization_strategy in normalization_methods) {
  normalization_dir <- file.path(INTERSECTOR_PARENT_INPUT_DIR, normalization_strategy)

  # Get all cell type sets
  cell_type_sets <- list.dirs(normalization_dir, full.names = FALSE, recursive = FALSE)
  
  if (length(cell_type_sets) == 0) {
    warning(paste("No cell type sets found in:", normalization_dir))
    next
  }

  # Process each cell type set like all_clusters or macrophages
  for (cell_type_set in cell_type_sets) {
    cell_type_set_dir <- file.path(normalization_dir, cell_type_set)

    # Only process the metabolic genes
    gene_type <- "metabolic"

    gene_type_dir <- file.path(cell_type_set_dir, gene_type)
    meta_scores_file_path <- file.path(gene_type_dir, "meta_scores.csv")

    if (!file.exists(meta_scores_file_path)) {
      warning(paste("Meta scores file not found:", meta_scores_file_path))
      next
    }

    # Currently uses data tables but switch to tibbles at some point with read_csv
    meta_scores_dt <- fread(meta_scores_file_path)

    # Map genes to Entrez IDs
    genes_for_combination <- unique(meta_scores_dt$gene)
    entrez_df <- map_symbols_to_entrez(genes_for_combination)
    entrez_vector <- unique(na.omit(entrez_df$ENTREZID))

    # Values are symbols and the names are Entrez IDs
    symbol_entrez_map <- setNames(entrez_df$SYMBOL, entrez_df$ENTREZID)

    # Perform KEGG mapping
    kegg_map_obj <- perform_kegg_mapping(entrez_vector)

    if (!is.null(kegg_map_obj) && nrow(as.data.frame(kegg_map_obj)) > 0) {
      # Process KEGG results
      processed_results <- process_kegg_results(
        as.data.frame(kegg_map_obj),
        symbol_entrez_map,
        meta_scores_dt
      )
      
      # Filter for metabolic pathways
      processed_results <- processed_results %>% filter(PathwayID %in% METABOLIC_PATHWAYS)

      # Create output dir structure
      combination_output_dir <- file.path(
        PATHWAY_MAPPING_OUTPUT_PARENT_DIR,
        normalization_strategy,
        cell_type_set,
        gene_type
      )
      dir_create(combination_output_dir, recurse = TRUE)

      # Save results
      output_file <- file.path(
        combination_output_dir,
        paste0("kegg_pathway_enrichment_", normalization_strategy, "_", cell_type_set, "_", gene_type, ".csv")
      )
      write.csv(processed_results, output_file, row.names = FALSE)
    } else {
      message(paste("No KEGG results for combination:", normalization_strategy, cell_type_set, gene_type))
    }
  }
}

message("Metabolic pathway mapping complete.")
