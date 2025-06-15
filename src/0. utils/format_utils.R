library(clusterProfiler)
library(org.Hs.eg.db)

# Function that simulates scRNAseq as bulk
sc_to_bulk <- function(sc_data) {
  # TODO: implement function
  bulk_data <- sc_data
  return(bulk_data)
}

# Function to map gene symbols to Entrez IDs
map_symbols_to_entrez <- function(gene_symbols) {
  # Convert gene symbols to Entrez IDs using org.Hs.eg.db
  entrez_ids <- bitr(gene_symbols,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  # Return the resulting data frame
  return(entrez_ids)
}

# Function to map Entrez IDs to gene symbols
map_entrez_to_symbols <- function(entrez_ids) {
  # Convert Entrez IDs to Gene Symbols using org.Hs.eg.db
  gene_symbols <- bitr(entrez_ids,
                       fromType = "ENTREZID",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)
  column_names <- c("Entrez_ID", "Gene_Symbol")
  colnames(gene_symbols) <- column_names
  
  return(gene_symbols)
}
