library(KEGGREST)

source("src/0. utils/format_utils.R")

# Function to retrieve gene IDs from a specific pathway
get_genes_from_pathway <- function(pathway_id) {
  pathway_info <- keggGet(pathway_id)
  
  if (!is.null(pathway_info[[1]]$GENE)) {
    genes <- pathway_info[[1]]$GENE
    gene_ids <- genes[seq(1, length(genes), by = 2)]
  } else {
    gene_ids <- character(0)
  }

  return(gene_ids)
}

# hsa01100 (metabolic pathways)
pathway_info <- keggGet("hsa01100")

# hsa01100 is an overview pathway and all the related pathways are the actual gene containing metabolic pathways
related_pathways <- names(pathway_info[[1]]$REL_PATHWAY)

# Set output path and save the related metabolic pathways for the pathway mapper
pathways_csv_path <- "output/1. data preprocessing/kegg/hsa01100_rel_pathways.csv"
dir.create(dirname(pathways_csv_path), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(PathwayID = related_pathways), pathways_csv_path, row.names = FALSE)

# Get the genes from each related pathway using the helper function
genes <- unlist(lapply(related_pathways, get_genes_from_pathway))
mapped_metabolic_genes <- map_entrez_to_symbols(genes)

# Save the results
genes_csv_path <- "output/1. data preprocessing/kegg/hsa01100_genes.csv"
dir.create(dirname(genes_csv_path), recursive = TRUE, showWarnings = FALSE)
write.csv(mapped_metabolic_genes, file = genes_csv_path, row.names = FALSE)

# Not all Entrez IDs are in use or have a corresponding gene symbol so indicate how many were found
message("Finished hsa01100 KEGG extracts. Total unique genes: ", length(genes))
