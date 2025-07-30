library(dplyr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fs)
library(stringr)

source("src/0. utils/format_utils.R") 


# Set all the constants
INTERSECTOR_PARENT_INPUT_DIR <- "output/3. intersector"
PATHWAY_MAPPING_OUTPUT_PARENT_DIR <- "output/5. pathway mapping"
KEGG_ORGANISM_CODE <- "hsa"
METABOLIC_PATHWAYS_FILE <- "output/1. data preprocessing/kegg/hsa01100_rel_pathways.csv"
METABOLIC_PATHWAYS <- read.csv(METABOLIC_PATHWAYS_FILE)$PathwayID


# Parse KEGG output format to get gene symbols
map_kegg_geneid_to_symbols <- function(gene_id_str, entrez_to_symbol_map) {
  # Entrez IDs are returned separated by slashes
  entrez_ids <- str_split(gene_id_str, "/")[[1]]

  # Map them to symbols with our utils
  symbols <- entrez_to_symbol_map[entrez_ids]

  # Keep Entrez IDs if no symbol is found but it shouldnt happen
  symbols <- ifelse(is.na(symbols), entrez_ids, symbols)

  # Join the symbols into a single slash separated string
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

  # If we have meta scores calculate metrics for each cell type
  if (!is.null(meta_scores_dt)) {
    # Get all unique cell types in the data
    all_cell_types <- unique(meta_scores_dt$sublineage)
    
    # For each pathway count genes per cell type
    for (cell_type in all_cell_types) {
      # Use original sublineage format for column names
      col_name <- paste0("Genes_", cell_type)

      # Create a new column for this cell type and count genes
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
  }

  # Create dynamic column selection including all cell type columns
  cell_type_columns <- names(kegg_result_df)[grepl("^Genes_", names(kegg_result_df))]
  
  # Select columns and order by gene count contribution
  kegg_result_df <- kegg_result_df %>%
    dplyr::select(
      PathwayID = ID,
      PathwayName = Description,
      GeneCountFromInput,
      dplyr::any_of(c("NumberOfContributingCellTypes", "TotalGeneCellTypeCombinations")),
      InputGeneSymbolsInPathway,
      dplyr::all_of(cell_type_columns)
    ) %>%
    dplyr::arrange(desc(GeneCountFromInput))

  return(kegg_result_df)
}


# Retrieve all metabolic pathways from KEGG using Entrez IDs
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


# Create output dir
dir_create(PATHWAY_MAPPING_OUTPUT_PARENT_DIR, recurse = TRUE)

# Initialize a variable to hold total gene counts for pathways
pathway_total_gene_counts <- NULL

# Create a baseline where we pass all metabolic genes
all_metabolic_genes_path <- "output/1. data preprocessing/kegg/hsa01100_genes.csv"

# If we have the metabolic genes file then go for it
if (file.exists(all_metabolic_genes_path)) {
  # Read the metabolic genes
  all_metabolic_genes <- read.csv(all_metabolic_genes_path)

  # Map symbols to Entrez IDs
  all_metabolic_entrez_df <- map_symbols_to_entrez(unique(all_metabolic_genes$Symbol))

  # Get unique Entrez IDs list which it should already be
  all_metabolic_entrez_vector <- unique(na.omit(all_metabolic_entrez_df$ENTREZID))
  
  # Perform KEGG mapping to get corresponding pathways
  all_pathways_map_obj <- perform_kegg_mapping(all_metabolic_entrez_vector)
  
  # If we have the KEGG results then process them
  if (!is.null(all_pathways_map_obj) && nrow(as.data.frame(all_pathways_map_obj)) > 0) {
    # Set names to Entrez IDs and keep the symbols as values
    all_pathways_symbol_map <- setNames(all_metabolic_entrez_df$SYMBOL, all_metabolic_entrez_df$ENTREZID)
    
    # Get the baseline count results for all metabolic genes
    processed_all_pathways <- process_kegg_results(as.data.frame(all_pathways_map_obj), all_pathways_symbol_map) %>%
      dplyr::filter(PathwayID %in% METABOLIC_PATHWAYS)
    
    # Write baseline metabolic pathways to file
    output_file <- file.path(PATHWAY_MAPPING_OUTPUT_PARENT_DIR, "all_metabolic_pathways_gene_counts.csv")
    write.csv(processed_all_pathways, output_file, row.names = FALSE)
    message(paste("Saved baseline metabolic pathway file to:", output_file))

    # Load the baseline counts for calculating percentages later
    pathway_total_gene_counts <- processed_all_pathways %>%
      dplyr::select(PathwayID, TotalMetabolicGenesInPathway = GeneCountFromInput)
  }
}

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

  # Process each cell type set
  for (cell_type_set in cell_type_sets) {
    # Set cell type directory
    cell_type_set_dir <- file.path(normalization_dir, cell_type_set)

    # Only process the metabolic genes
    gene_type <- "metabolic"

    # Check if the meta scores file exists
    gene_type_dir <- file.path(cell_type_set_dir, gene_type)
    meta_scores_file_path <- file.path(gene_type_dir, "meta_scores.csv")

    if (!file.exists(meta_scores_file_path)) {
      warning(paste("Meta scores file not found:", meta_scores_file_path))
      next
    }

    # Currently uses data tables but switch to tibbles later
    meta_scores_dt <- fread(meta_scores_file_path)

    # Map genes to Entrez IDs
    genes_for_combination <- unique(meta_scores_dt$gene)
    entrez_df <- map_symbols_to_entrez(genes_for_combination)
    entrez_vector <- unique(na.omit(entrez_df$ENTREZID))

    # Values are symbols and the names are Entrez IDs
    symbol_entrez_map <- setNames(entrez_df$SYMBOL, entrez_df$ENTREZID)

    # Perform KEGG mapping
    kegg_map_obj <- perform_kegg_mapping(entrez_vector)

    # If we have the KEGG results then process them
    if (!is.null(kegg_map_obj) && nrow(as.data.frame(kegg_map_obj)) > 0) {
      # Process KEGG results
      processed_results <- process_kegg_results(
        as.data.frame(kegg_map_obj),
        symbol_entrez_map,
        meta_scores_dt
      )

      # Filter for metabolic pathways
      processed_results <- processed_results %>% dplyr::filter(PathwayID %in% METABOLIC_PATHWAYS)

      # If we have the total gene counts calculate the pathway coverage percentage
      if (!is.null(pathway_total_gene_counts)) {
        processed_results <- processed_results %>%
          dplyr::left_join(pathway_total_gene_counts, by = "PathwayID") %>%
          dplyr::mutate(PercentageOfTotalGenes = (GeneCountFromInput / TotalMetabolicGenesInPathway) * 100) %>%
           # Remove the temporary total count column
          dplyr::select(-TotalMetabolicGenesInPathway)
      }

      # Sort the main results file
      processed_results <- processed_results %>%
        dplyr::arrange(desc(GeneCountFromInput))

      # Create output dir structure
      combination_output_dir <- file.path(
        PATHWAY_MAPPING_OUTPUT_PARENT_DIR,
        normalization_strategy,
        cell_type_set,
        gene_type
      )

      # Create output directory if it doesn't exist
      dir_create(combination_output_dir, recurse = TRUE)

      # Save results for all cell types
      output_file <- file.path(
        combination_output_dir,
        paste0("kegg_pathway_mapping_", normalization_strategy, "_", cell_type_set, "_", gene_type, ".csv")
      )
      write.csv(processed_results, output_file, row.names = FALSE)

      # For each cell type create subdirectory and save pathway stats
      cluster_columns <- names(processed_results)[grepl("^Genes_", names(processed_results))]

      for (cluster_col in cluster_columns) {
        # Only keep pathways where this cell type contributed at least one gene
        idx <- processed_results[[cluster_col]] > 0
        if (!any(idx)) next

        # Get cluster name by removing the genes prefix
        cluster_name <- sub('^Genes_', '', cluster_col)

        # Create subdirectory for this cell type
        cluster_dir <- file.path(combination_output_dir, cluster_name)
        dir_create(cluster_dir, recurse = TRUE)

        # For each pathway extract the contributed gene symbols for this cell type
        contributed_genes <- lapply(seq_len(nrow(processed_results)), function(i) {
          if (!idx[i]) return(NA)

          # Get pathway gene symbols
          pathway_genes <- unlist(strsplit(processed_results$InputGeneSymbolsInPathway[i], "/"))

          # Find which genes in this pathway are present for this cell type
          celltype_genes <- meta_scores_dt$gene[meta_scores_dt$gene %in% pathway_genes & meta_scores_dt$sublineage == cluster_name]
          
          # If no genes found return NA
          if (length(celltype_genes) == 0) return(NA)

          # Join the unique cell type genes into a single string
          paste(unique(celltype_genes), collapse = "/")
        })

        # Create df for this cell type
        celltype_summary <- processed_results[idx, c("PathwayID", "PathwayName", cluster_col, "InputGeneSymbolsInPathway", "PercentageOfTotalGenes"), drop=FALSE]
        names(celltype_summary)[names(celltype_summary) == cluster_col] <- "NumberOfGenesContributed"
        celltype_summary$ContributedGeneSymbols <- unlist(contributed_genes[idx])
        
        # Sort the cell type specific file
        celltype_summary <- celltype_summary %>%
          dplyr::arrange(desc(NumberOfGenesContributed))
        
        # Save file
        cluster_file <- file.path(
          cluster_dir,
          paste0("kegg_pathway_mapping_", normalization_strategy, "_", cell_type_set, "_", gene_type, "_", cluster_name, ".csv")
        )

        # Write the cell type specific file
        write.csv(celltype_summary, cluster_file, row.names = FALSE)
      }
    } else {
      # This should never happen
      message(paste("No KEGG results for combination:", normalization_strategy, cell_type_set, gene_type))
    }
  }
}

message("Metabolic pathway mapping complete.")
