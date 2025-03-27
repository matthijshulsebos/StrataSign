# Load required libraries
library(DESeq2)
library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)

# Create output directories
dir.create("output/differential_expression", recursive = TRUE, showWarnings = FALSE)

# Load Leader et al. data
if (!exists("lung_ldm")) {
  load("base/data/lung_ldm.rd")
}

# Load metadata
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv")
annots_list <- read_csv("base/input_tables/annots_list.csv")

# If tissue is a list column, fix it
if(is.list(table_s1$tissue)) {
  table_s1$tissue <- unlist(table_s1$tissue)
}

# Filter samples to those used in the clustering model 
table_s1 <- table_s1 %>% 
  filter(Use.in.Clustering.Model. == "Yes") %>%
  filter(!is.na(tissue)) %>%
  mutate(tissue = toupper(tissue))

# Define doublet clusters to exclude
clusters_to_exclude <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Prepare count data
counts <- lung_ldm$dataset$counts
  
# Filter out doublets
all_clusters <- as.numeric(dimnames(counts)[[3]])
clusters_to_keep <- all_clusters[!all_clusters %in% clusters_to_exclude]
counts <- counts[, , as.character(clusters_to_keep), drop = FALSE]

# Filter samples to match table_s1
sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
counts <- counts[dimnames(counts)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]

# Collapse clusters
counts_summed <- apply(counts, c(1, 2), sum)

# Set reference (normal) and comparison (tumor) levels
unique_tissues <- sort(unique(table_s1$tissue))
ref_level <- unique_tissues[1]
comparison_level <- unique_tissues[2]

# Create sample metadata for DESeq2
sample_metadata <- table_s1 %>%
  select(sample_ID, tissue) %>%
  filter(sample_ID %in% rownames(counts_summed)) %>%
  mutate(
    condition = tissue,
    condition = factor(condition, levels = unique_tissues)
  ) %>%
  column_to_rownames("sample_ID")

# Ensure count matrix and metadata are aligned
counts_summed <- counts_summed[rownames(sample_metadata), ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(t(counts_summed)),
  colData = sample_metadata,
  design = ~ condition
)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", comparison_level, ref_level))

# Update the output file names and labels to match actual tissue values
file_str <- paste0(tolower(comparison_level), "_vs_", tolower(ref_level))

# Convert to data frame and add gene name
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Save full results
results_file <- paste0("output/differential_expression/de_results_", file_str, ".csv")
write_csv(res_df, results_file)

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red"), 
                   labels = c("Not significant", "Significant"), 
                   name = "Differential Expression") +
  theme_minimal() +
  labs(title = "Volcano Plot - Tumor vs Normal",
       x = "Log2 Fold Change (tumor vs normal)",
       y = "-Log10 Adjusted P-value")

ggsave("output/differential_expression/volcano_plot.png", volcano_plot, width = 10, height = 8)

# Generate summary statistics
summary_stats <- data.frame(
  Description = c(
    "Total genes analyzed",
    "Significant genes (padj < 0.05)",
    "Upregulated in tumor (padj < 0.05, log2FC > 0)",
    "Downregulated in tumor (padj < 0.05, log2FC < 0)",
    "Highly upregulated (padj < 0.05, log2FC > 2)",
    "Highly downregulated (padj < 0.05, log2FC < -2)"
  ),
  Count = c(
    nrow(res_df),
    sum(res_df$padj < 0.05, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 2, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -2, na.rm = TRUE)
  )
)

# Save summary
write_csv(summary_stats, "output/differential_expression/de_analysis_summary.csv")

# Create gene-level fold changes file 
gene_fold_changes <- res_df %>%
  mutate(Feature = gene) %>%
  select(Feature, Value = log2FoldChange)

write_csv(gene_fold_changes, "output/differential_expression/gene_fold_changes.csv")
message("Created gene_fold_changes.csv - run cluster_differential_expression.r to get complete feature fold changes")
