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
DATASET_TYPES_TO_PROCESS <- c("absolute", "relative")
KEGG_ORGANISM_CODE <- "hsa"

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

# Set output dir
dir_create(PATHWAY_MAPPING_OUTPUT_PARENT_DIR, recurse = TRUE)

# Main loop to process each dataset type
for (current_dataset_type in DATASET_TYPES_TO_PROCESS) {
  type_specific_output_dir <- file.path(PATHWAY_MAPPING_OUTPUT_PARENT_DIR, current_dataset_type)
  dir_create(type_specific_output_dir, recurse = TRUE)

  meta_scores_file_path <- file.path(INTERSECTOR_PARENT_INPUT_DIR, current_dataset_type, "meta_scores.csv")
  meta_scores_dt <- fread(meta_scores_file_path)
  
  # Overall KEGG pathway mapping
  all_genes_symbols <- unique(meta_scores_dt$gene)
  entrez_df <- map_symbols_to_entrez(all_genes_symbols)
  entrez_vector <- unique(na.omit(entrez_df$ENTREZID))
  symbol_entrez_map <- setNames(entrez_df$SYMBOL, entrez_df$ENTREZID)
  kegg_map_overall_obj <- perform_kegg_mapping(entrez_vector)

  if (!is.null(kegg_map_overall_obj) && nrow(as.data.frame(kegg_map_overall_obj)) > 0) {
    processed_overall_results <- process_kegg_results(
      as.data.frame(kegg_map_overall_obj),
      symbol_entrez_map,
      meta_scores_dt
    )
    output_file_overall <- file.path(
      type_specific_output_dir,
      paste0(current_dataset_type, "_all_genes_kegg_pathway_map.csv")
    )
    write.csv(processed_overall_results, output_file_overall, row.names = FALSE)
  }

  # Cell type KEGG pathway mapping
  for (current_cell_type in unique(meta_scores_dt$sublineage)) {

    sanitized_cell_type_name <- gsub("[^a-zA-Z0-9_.-]", "_", current_cell_type)
    cell_type_output_dir <- file.path(type_specific_output_dir, sanitized_cell_type_name)
    dir_create(cell_type_output_dir, recurse = TRUE)

    genes_for_cell_type <- meta_scores_dt[sublineage == current_cell_type, unique(gene)]
    entrez_df_cell <- map_symbols_to_entrez(genes_for_cell_type)
    entrez_vector_cell <- unique(na.omit(entrez_df_cell$ENTREZID))
    symbol_entrez_map_cell <- setNames(entrez_df_cell$SYMBOL, entrez_df_cell$ENTREZID)
    kegg_map_cell_type_obj <- perform_kegg_mapping(entrez_vector_cell)

    if (!is.null(kegg_map_cell_type_obj) && nrow(as.data.frame(kegg_map_cell_type_obj)) > 0) {
      processed_cell_type_results <- process_kegg_results(
        as.data.frame(kegg_map_cell_type_obj),
        symbol_entrez_map_cell
      )
      output_file_cell_type <- file.path(
        cell_type_output_dir,
        paste0(sanitized_cell_type_name, "_kegg_pathway_map.csv")
      )
      write.csv(processed_cell_type_results, output_file_cell_type, row.names = FALSE)
    }
  }
}

message("Results are saved in: ", PATHWAY_MAPPING_OUTPUT_PARENT_DIR)
