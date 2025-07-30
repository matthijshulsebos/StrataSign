library(clusterProfiler)
library(org.Hs.eg.db)


# Gene symbols to Entrez IDs mapping
map_symbols_to_entrez <- function(gene_symbols) {
  entrez_ids <- bitr(gene_symbols,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  
  return(entrez_ids)
}

# Entrez IDs to gene symbols mapping
map_entrez_to_symbols <- function(entrez_ids) {
  gene_symbols <- bitr(entrez_ids,
                       fromType = "ENTREZID",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)

  # Messed up backward compatibility for data preprocessing script
  column_names <- c("Entrez_ID", "Symbol")
  colnames(gene_symbols) <- column_names
  
  return(gene_symbols)
}
