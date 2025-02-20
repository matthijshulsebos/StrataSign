# Load required libraries
library(DESeq2)
library(dplyr)
library(tidyverse)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

# Create output directory
dir.create("src_output/sample_de", recursive = TRUE, showWarnings = FALSE)

# Read the data
counts_data <- read_csv("src_output/Leader/pb_sample_metabolic.csv") %>% as.data.frame()
sample_info <- read_csv("input_tables/table_s1_sample_table.csv") %>% as.data.frame()

# Prepare count matrix
count_matrix <- counts_data %>%
  column_to_rownames(var = names(counts_data)[1]) %>%
  as.matrix()

# Prepare sample information
sample_metadata <- sample_info %>%
  dplyr::select(
    sample_id = 1,  # column containing sample IDs
    condition = 6   # column containing condition (tumor/normal)
  ) %>%
  filter(sample_id %in% colnames(count_matrix)) %>%
  mutate(condition = factor(condition)) %>%
  arrange(match(sample_id, colnames(count_matrix)))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = sample_metadata,
  design = ~ condition
)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Convert results to data frame and add gene names
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Save full results
write_csv(res_df, "src_output/sample_de/differential_expression_results.csv")

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6) +
  scale_color_manual(
    values = c("grey", "red"),
    labels = c("Not significant", "Significant"),
    name = "Differential Expression"
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differential Expression Analysis",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

# Save plot
ggsave("src_output/sample_de/volcano_plot.png", volcano_plot, width = 10, height = 8)

# Save significant genes (padj < 0.05)
significant_genes <- res_df %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange)))

write_csv(significant_genes, "src_output/sample_de/significant_genes.csv")

# Print summary
print("\nAnalysis Summary:")
print(summary(res))
print(paste("\nNumber of significant genes (padj < 0.05):", nrow(significant_genes)))

# Pathway Enrichment Analysis
#' Function to create and save enrichment plots
create_enrichment_plots <- function(enrichment_result, filename_prefix, title_prefix) {
  # Create dot plot
  p1 <- dotplot(enrichment_result, showCategory=20) +
    ggtitle(paste(title_prefix, "- Dot Plot")) +
    theme(axis.text.y = element_text(size=8))
  ggsave(paste0("src_output/sample_de/", filename_prefix, "_dotplot.png"), 
         p1, width=12, height=8)
  
  # Create bar plot
  df <- as.data.frame(enrichment_result)
  if(nrow(df) > 0) {
    p2 <- ggplot(head(df, 20),
                 aes(x=reorder(Description, -p.adjust), y=-log10(p.adjust))) +
      geom_bar(stat="identity", fill="steelblue") +
      coord_flip() +
      theme_minimal() +
      labs(title=paste(title_prefix, "- Top 20 Terms"),
           x="Pathway",
           y="-log10(Adjusted P-value)") +
      theme(axis.text.y=element_text(size=8))
    ggsave(paste0("src_output/sample_de/", filename_prefix, "_barplot.png"), 
           p2, width=12, height=8)
  }
}

# Convert gene IDs
genes <- bitr(significant_genes$gene, 
              fromType = "SYMBOL",
              toType = c("ENTREZID", "ENSEMBL"),
              OrgDb = org.Hs.eg.db)

# Merge with original data
deg_results_with_ids <- significant_genes %>%
  inner_join(genes, by = c("gene" = "SYMBOL"))

# Create ranked gene list
ranked_genes <- deg_results_with_ids$log2FoldChange
names(ranked_genes) <- deg_results_with_ids$ENTREZID
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Perform enrichment analyses
go_bp <- enrichGO(gene = deg_results_with_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

go_mf <- enrichGO(gene = deg_results_with_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

go_cc <- enrichGO(gene = deg_results_with_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

kegg_enrichment <- enrichKEGG(gene = deg_results_with_ids$ENTREZID,
                             organism = 'hsa',
                             pvalueCutoff = 0.05)

reactome_enrichment <- enrichPathway(gene = deg_results_with_ids$ENTREZID,
                                   organism = "human",
                                   pvalueCutoff = 0.05)

# Save enrichment results
write_csv(as.data.frame(go_bp), "src_output/sample_de/GO_biological_process.csv")
write_csv(as.data.frame(go_mf), "src_output/sample_de/GO_molecular_function.csv")
write_csv(as.data.frame(go_cc), "src_output/sample_de/GO_cellular_component.csv")
write_csv(as.data.frame(kegg_enrichment), "src_output/sample_de/KEGG_pathways.csv")
write_csv(as.data.frame(reactome_enrichment), "src_output/sample_de/Reactome_pathways.csv")

# Create enrichment plots
create_enrichment_plots(go_bp, "GO_BP", "GO Biological Process")
create_enrichment_plots(go_mf, "GO_MF", "GO Molecular Function")
create_enrichment_plots(go_cc, "GO_CC", "GO Cellular Component")
create_enrichment_plots(kegg_enrichment, "KEGG", "KEGG Pathways")
create_enrichment_plots(reactome_enrichment, "Reactome", "Reactome Pathways")

# Create network plots
tryCatch({
  kegg_network <- pairwise_termsim(kegg_enrichment)
  p_kegg <- emapplot(kegg_network, showCategory = 30)
  ggsave("src_output/sample_de/KEGG_network.png", p_kegg, width=12, height=12)
  
  go_bp_network <- pairwise_termsim(go_bp)
  p_go <- emapplot(go_bp_network, showCategory = 30)
  ggsave("src_output/sample_de/GO_BP_network.png", p_go, width=12, height=12)
}, error = function(e) {
  message("Note: Could not create some network plots due to insufficient similarities")
  print(e)
})

# Create summary of top pathways
top_pathways <- bind_rows(
  as.data.frame(go_bp)[1:5, c("Description", "p.adjust", "Count")] %>% 
    mutate(Category = "GO Biological Process"),
  as.data.frame(kegg_enrichment)[1:5, c("Description", "p.adjust", "Count")] %>% 
    mutate(Category = "KEGG"),
  as.data.frame(reactome_enrichment)[1:5, c("Description", "p.adjust", "Count")] %>% 
    mutate(Category = "Reactome")
)

write_csv(top_pathways, "src_output/sample_de/top_pathways_summary.csv")

# Print enrichment summary
cat("\nEnrichment Analysis Summary:\n")
cat("GO Biological Process terms:", nrow(as.data.frame(go_bp)), "\n")
cat("GO Molecular Function terms:", nrow(as.data.frame(go_mf)), "\n")
cat("GO Cellular Component terms:", nrow(as.data.frame(go_cc)), "\n")
cat("KEGG pathways:", nrow(as.data.frame(kegg_enrichment)), "\n")
cat("Reactome pathways:", nrow(as.data.frame(reactome_enrichment)), "\n") 
