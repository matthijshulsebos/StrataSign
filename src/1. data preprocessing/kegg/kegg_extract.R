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

# Query hsa01100 (metabolic pathways)
pathway_info <- keggGet("hsa01100")

# Extract related pathways
related_pathways <- names(pathway_info[[1]]$REL_PATHWAY)

# Save related pathways to a csv file
pathways_csv_path <- "output/1. data preprocessing/kegg/hsa01100_rel_pathways.csv"
dir.create(dirname(pathways_csv_path), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(PathwayID = related_pathways), pathways_csv_path, row.names = FALSE)
cat("Pathways saved to:", pathways_csv_path, "\n")

# Extract genes from all related pathways
genes <- unlist(lapply(related_pathways, get_genes_from_pathway))

# Map Entrez IDs to gene symbols
source("src/0. utils/format_utils.R")
mapped_metabolic_genes <- map_entrez_to_symbols(genes)

# Save genes to a csv file
genes_csv_path <- "output/1. data preprocessing/kegg/hsa01100_genes.csv"
dir.create(dirname(genes_csv_path), recursive = TRUE, showWarnings = FALSE)
write.csv(mapped_metabolic_genes, file = genes_csv_path, row.names = FALSE)
cat("Total unique genes related to hsa01100:", length(genes), "\n")
cat("Genes saved to:", genes_csv_path, "\n")
