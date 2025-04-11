# Load the KEGGREST library
library(KEGGREST)

# Function to retrieve gene IDs from a specific pathway
get_genes_from_pathway <- function(pathway_id) {
  # Query the pathway
  pathway_info <- keggGet(pathway_id)
  # Extract genes (if available)
  if (!is.null(pathway_info[[1]]$GENE)) {
    genes <- pathway_info[[1]]$GENE
    gene_ids <- genes[seq(1, length(genes), by = 2)]
  } else {
    gene_ids <- character(0)  # No genes found
  }
  return(gene_ids)
}

# Query hsa01100 (Metabolic Pathways)
pathway_info <- keggGet("hsa01100")

# Extract related pathways
related_pathways <- names(pathway_info[[1]]$REL_PATHWAY)

# Write sub-pathways to file for confirmation
sink("output/kegg/rel_pathways.txt")
print(related_pathways)
sink()

# Extract genes from all related pathways
genes <- unlist(lapply(related_pathways, get_genes_from_pathway))

source("src/utils/format_utils.R")
mapped_metabolic_genes <- map_entrez_to_symbols(genes)

# Write to file
cat("Total unique genes related to hsa01100:", length(genes), "\n")
write.csv(
  mapped_metabolic_genes,
  file = "output/kegg/hsa01100_genes.csv",
  row.names = FALSE
)
cat("Genes saved to 'hsa01100_genes.txt'.\n")
