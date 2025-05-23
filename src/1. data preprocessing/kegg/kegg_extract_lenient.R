# Load the KEGGREST library
library(KEGGREST)

# --- Configuration ---
OUTPUT_BASE_DIR <- "output/1. data preprocessing/kegg"
GENES_FILENAME <- "metabolic_genes_lenient.csv"
PATHWAYS_FILENAME <- "pathways_lenient.txt"

# Ensure output directory exists
if (!dir.exists(OUTPUT_BASE_DIR)) {
  dir.create(OUTPUT_BASE_DIR, recursive = TRUE)
}

# Function to retrieve gene IDs from a specific pathway
get_genes_from_pathway <- function(pathway_id) {
  tryCatch({
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
  }, error = function(e) {
    warning(paste("Failed to retrieve genes for pathway:", pathway_id, "-", e$message))
    return(character(0))
  })
}

# --- Step 1: Retrieve pathways ---
cat("Retrieving pathways...\n")

# Method 1: Pathways related to hsa01100
cat("Retrieving pathways related to hsa01100...\n")
hsa01100_info <- keggGet("hsa01100")
related_pathways <- if (!is.null(hsa01100_info[[1]]$REL_PATHWAY)) {
  names(hsa01100_info[[1]]$REL_PATHWAY)
} else {
  character(0)
}

# Method 2: Pathways matching refined metabolic keywords
cat("Searching for pathways using refined keywords...\n")
all_pathways <- keggList("pathway", "hsa")
keywords <- c(
  "metabolic", "metabolism", "biosynthesis", "degradation", "catabolism", "anabolism",
  "glycolysis", "gluconeogenesis", "TCA cycle", "citrate cycle", "oxidation", "reduction",
  "phosphorylation", "lipid metabolism", "amino acid metabolism", "nucleotide metabolism",
  "carbohydrate metabolism", "secondary metabolism", "polyketide", "terpenoid", "alkaloid",
  "flavonoid", "phenylpropanoid", "energy metabolism", "oxidative phosphorylation",
  "electron transport", "ATP synthesis"
)
pattern <- paste(keywords, collapse = "|")
keyword_pathways <- names(all_pathways)[grepl(pattern, all_pathways, ignore.case = TRUE)]

# Method 3: Pathways under hsa01xxx (Metabolism category)
cat("Retrieving pathways under hsa01xxx...\n")
category_pathways <- names(all_pathways)[grepl("^hsa01", names(all_pathways))]

# Combine all pathways and deduplicate
all_pathways <- unique(c(related_pathways, keyword_pathways, category_pathways))
cat("Total unique pathways identified:", length(all_pathways), "\n")

# Write pathways to file
pathways_file_path <- file.path(OUTPUT_BASE_DIR, PATHWAYS_FILENAME)
writeLines(all_pathways, pathways_file_path)
cat("Pathways written to:", pathways_file_path, "\n")

# --- Step 2: Extract genes from pathways ---
cat("Extracting genes from pathways...\n")
all_genes <- unlist(lapply(all_pathways, get_genes_from_pathway))

# Deduplicate genes
unique_genes <- unique(all_genes)
cat("Total unique genes identified:", length(unique_genes), "\n")

# --- Step 3: Map genes to symbols ---
cat("Mapping genes to symbols...\n")
source("src/0. utils/format_utils.R")
mapped_genes <- map_entrez_to_symbols(unique_genes)

# --- Step 4: Write results to file ---
cat("Writing results to file...\n")
genes_file_path <- file.path(OUTPUT_BASE_DIR, GENES_FILENAME)
write.csv(mapped_genes, file = genes_file_path, row.names = FALSE)
cat("Metabolic genes saved to:", genes_file_path, "\n")