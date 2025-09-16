###############################################################################
# KEGG metabolic gene extraction script (human, hsa01100) — fully annotated
#
# What this script does (high level):
# 1) Fetch the KEGG “overview” metabolic pathway (hsa01100).
# 2) Extract all related (child) metabolic pathway IDs from that overview.
# 3) For each child pathway, retrieve all associated human genes (Entrez IDs).
# 4) Map those Entrez IDs to human gene symbols (via a helper you provide).
# 5) Save both the list of related pathways and the mapped genes to CSV files.
#
# Notes:
# - Requires internet access (KEGG REST API).
# - Uses KEGGREST::keggGet() which returns KEGG records as lists.
# - In KEGGREST, the $GENE field in a pathway is a character vector of PAIRS:
#     c("entrez_id1", "SYMBOL1; Description ...",
#       "entrez_id2", "SYMBOL2; Description ...", ...)
#   We therefore take every 1st, 3rd, 5th, ... element to keep only Entrez IDs.
###############################################################################

# --- Dependencies ------------------------------------------------------------

library(KEGGREST)  # Provides keggGet() to query KEGG entries

# This sourced file defines `map_entrez_to_symbols()`, a helper
# that takes a vector of Entrez Gene IDs (character) and returns a mapped
# representation (vector or data.frame) of gene symbols.
source("src/0. utils/format_utils.R")


# --- Helper function: get genes from a KEGG pathway --------------------------
# Given a KEGG pathway ID (e.g., "hsa00010"), fetch the pathway record and
# extract the Entrez Gene IDs listed under the $GENE field.
get_genes_from_pathway <- function(pathway_id) {
  # Query KEGG for the pathway record; returns a list of length 1
  pathway_info <- keggGet(pathway_id)
  
  # The $GENE field (if present) is a character vector with alternating entries:
  # [1] entrez_id1, [2] "SYMBOL1; Description", [3] entrez_id2, [4] "SYMBOL2; ...", ...
  if (!is.null(pathway_info[[1]]$GENE)) {
    genes <- pathway_info[[1]]$GENE
    
    # Extract only the odd positions: 1, 3, 5, ... (i.e., the Entrez IDs)
    gene_ids <- genes[seq(1, length(genes), by = 2)]
  } else {
    # Some pathways may not list genes (e.g., purely informational maps)
    gene_ids <- character(0)
  }
  
  return(gene_ids)
}


# --- Step 1: Fetch the KEGG overview pathway (human metabolic pathways) ------
# hsa01100 is KEGG’s global “Metabolic pathways” map for Homo sapiens.
pathway_info <- keggGet("hsa01100")


# --- Step 2: Extract child/related pathways ----------------------------------
# The overview map references many specific metabolic pathways (glycolysis,
# TCA, PPP, etc.). KEGGREST exposes these under $REL_PATHWAY as a *named vector*:
#   names($REL_PATHWAY) = pathway IDs (e.g., "hsa00010")
#   values($REL_PATHWAY) = pathway titles (e.g., "Glycolysis / Gluconeogenesis")
related_pathways <- names(pathway_info[[1]]$REL_PATHWAY)


# --- Step 3: Save the list of related pathway IDs for traceability -----------
# This CSV can be used later by a pathway mapper or provenance logs.
pathways_csv_path <- "output/1. data preprocessing/kegg/hsa01100_rel_pathways.csv"

# Ensure the output directory exists (recursive in case parents don’t exist)
dir.create(dirname(pathways_csv_path), recursive = TRUE, showWarnings = FALSE)

# Write the related pathway IDs to disk
write.csv(
  data.frame(PathwayID = related_pathways),
  pathways_csv_path,
  row.names = FALSE
)


# --- Step 4: For each related pathway, collect the genes ---------------------
# lapply() will call get_genes_from_pathway(p) for each p in related_pathways.
# The result is a list of character vectors (Entrez IDs). unlist() flattens
# this into a single character vector (may contain duplicates across pathways).
genes <- unlist(lapply(related_pathways, get_genes_from_pathway))

# OPTIONAL (not changing behavior): If you want unique gene counts downstream,
# consider: unique_genes <- unique(genes)


# --- Step 5: Map Entrez IDs to human gene symbols ----------------------------
# This helper is expected to be defined in your sourced script. It should return
# either a vector of symbols (aligned to `genes`) or a data.frame with columns
# like `entrez_id`, `symbol`, etc. (behavior depends on your implementation).
mapped_metabolic_genes <- map_entrez_to_symbols(genes)


# --- Step 6: Save the mapped genes -------------------------------------------
# Persist the mapped gene set for downstream analysis or reporting.
genes_csv_path <- "output/1. data preprocessing/kegg/hsa01100_genes.csv"

# Ensure the output directory exists
dir.create(dirname(genes_csv_path), recursive = TRUE, showWarnings = FALSE)

# Write out the mapped gene object. If your mapper returns a vector, this will
# produce a one-column CSV; if it returns a data.frame, its columns are saved.
write.csv(mapped_metabolic_genes, file = genes_csv_path, row.names = FALSE)


# --- Step 7: Console summary message -----------------------------------------
# Note: `genes` may include duplicates (same gene in multiple pathways), 
# so print unique genes only.
message("Finished hsa01100 KEGG extracts. Total unique genes: ", length(unique(genes)))