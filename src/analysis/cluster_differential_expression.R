# Load required libraries
library(DESeq2)
library(dplyr)
library(tidyverse)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

# Load lung_ldm data
load("data/lung_ldm.rd")

# Read data
counts_array <- lung_ldm$dataset$counts
cluster_annot <- read_csv("input_tables/annots_list.csv")

# Print initial dimensions to understand the data structure
print("Dimensions of counts array:")
print(dim(counts_array))
print("Structure of first slice:")
print(dim(counts_array[,,1]))
print("First few gene names:")
print(head(dimnames(counts_array)[[2]]))  # genes should be in second dimension
print("First few sample names:")
print(head(dimnames(counts_array)[[1]]))  # samples should be in first dimension

# Read and process metabolic genes list
metabolic_genes_df <- read_csv("src_output/KEGG/hsa01100_genes.csv")

# Check if we have SYMBOL column, if not convert from ENTREZID
if("SYMBOL" %in% colnames(metabolic_genes_df)) {
  metabolic_genes <- metabolic_genes_df %>% pull(SYMBOL)
} else if("ENTREZID" %in% colnames(metabolic_genes_df)) {
  # Convert ENTREZID to SYMBOL
  gene_conversion <- bitr(
    metabolic_genes_df$ENTREZID, 
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
  )
  metabolic_genes <- gene_conversion$SYMBOL
} else {
  stop("Metabolic genes file must contain either SYMBOL or ENTREZID column")
}

# Print information about metabolic genes
print("Metabolic genes information:")
print(paste("Total number of metabolic genes:", length(metabolic_genes)))
print(paste("First few metabolic genes:", paste(head(metabolic_genes), collapse=", ")))

# Function to perform differential expression analysis for one cluster
perform_de_analysis <- function(count_matrix, sample_info, cluster_name, metabolic_genes, lineage, sub_lineage) {
  # Create output directory for this cluster
  dir.create(file.path("src_output/cluster_de", cluster_name), recursive = TRUE, showWarnings = FALSE)
  
  # Use sub_lineage if available, otherwise use lineage
  display_name <- if(!is.na(sub_lineage) && sub_lineage != "") {
    sub_lineage
  } else {
    lineage
  }
  
  # Print dimensions before processing
  print("Original count matrix dimensions:")
  print(dim(count_matrix))
  
  # Get gene names
  gene_names <- colnames(count_matrix)
  print(paste("Number of genes before filtering:", length(gene_names)))
  
  # Filter for metabolic genes
  metabolic_indices <- gene_names %in% metabolic_genes
  count_matrix <- count_matrix[, metabolic_indices]
  
  print(paste("Number of metabolic genes after filtering:", sum(metabolic_indices)))
  
  # Check if we have any metabolic genes
  if(sum(metabolic_indices) == 0) {
    warning(paste("No metabolic genes found in cluster", cluster_name))
    return(NULL)
  }
  
  # Transpose matrix for DESeq2 (genes should be rows)
  count_matrix <- t(count_matrix)
  
  print("Final count matrix dimensions (genes x samples):")
  print(dim(count_matrix))
  
  # Add a small pseudocount to handle zeros
  count_matrix <- count_matrix + 1
  
  # Pre-filter low count genes before DESeq2
  keep <- rowSums(count_matrix >= 10) >= 5  # At least 5 samples with counts >= 10
  count_matrix <- count_matrix[keep, ]
  
  print(paste("Number of genes after pre-filtering:", nrow(count_matrix)))
  
  # Check if we have enough genes left
  if(nrow(count_matrix) == 0) {
    warning(paste("No genes left after filtering in cluster", cluster_name))
    return(NULL)
  }
  
  # Prepare sample metadata
  sample_metadata <- sample_info %>%
    dplyr::select(
      sample_id = 1,    # adjust column index as needed
      condition = 6     # adjust column index as needed
    ) %>%
    filter(sample_id %in% colnames(count_matrix)) %>%
    mutate(condition = factor(condition)) %>%
    arrange(match(sample_id, colnames(count_matrix)))
  
  # Print condition distribution
  print("Number of samples per condition:")
  print(table(sample_metadata$condition))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = sample_metadata,
    design = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds, fitType = "local")
  res <- results(dds)
  
  # Convert results to data frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)
  
  # Save results
  write_csv(res_df, file.path("src_output/cluster_de", cluster_name, "metabolic_de_results.csv"))
  
  # Create volcano plot with display name
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_point(aes(color = padj < 0.05), alpha = 0.6) +
    scale_color_manual(
      values = c("grey", "red"),
      labels = c("Not significant", "Significant"),
      name = "Differential Expression"
    ) +
    theme_minimal() +
    labs(
      title = paste("Metabolic Genes Volcano Plot -", display_name),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    )
  
  ggsave(file.path("src_output/cluster_de", cluster_name, "metabolic_volcano_plot.png"), 
         volcano_plot, width = 10, height = 8)
  
  # Save significant genes
  significant_genes <- res_df %>%
    filter(padj < 0.05, !is.na(padj)) %>%
    arrange(desc(abs(log2FoldChange)))
  
  write_csv(significant_genes, 
           file.path("src_output/cluster_de", cluster_name, "significant_metabolic_genes.csv"))
  
  return(significant_genes)
}

# Function to perform pathway enrichment for one cluster
perform_pathway_enrichment <- function(significant_genes, cluster_name, lineage, sub_lineage) {
  if(is.null(significant_genes) || nrow(significant_genes) == 0) {
    warning(paste("No significant metabolic genes for", cluster_name))
    return(NULL)
  }
  
  # Use sub_lineage if available, otherwise use lineage
  display_name <- if(!is.na(sub_lineage) && sub_lineage != "") {
    sub_lineage
  } else {
    lineage
  }
  
  # Create output directory
  dir.create(file.path("src_output/cluster_de", cluster_name, "enrichment"), 
             recursive = TRUE, showWarnings = FALSE)
  
  # Convert gene IDs
  genes <- bitr(significant_genes$gene, 
                fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)
  
  # Merge with original data
  deg_results_with_ids <- significant_genes %>%
    dplyr::inner_join(genes, by = c("gene" = "SYMBOL"))
  
  # Create ranked gene list
  ranked_genes <- deg_results_with_ids$log2FoldChange
  names(ranked_genes) <- deg_results_with_ids$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Perform enrichment analyses focusing on metabolic pathways
  kegg_enrichment <- enrichKEGG(
    gene = deg_results_with_ids$ENTREZID,
    organism = 'hsa',
    pvalueCutoff = 0.05,
    universe = unique(bitr(metabolic_genes, 
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)$ENTREZID)
  )
  
  # Save results
  write_csv(as.data.frame(kegg_enrichment), 
           file.path("src_output/cluster_de", cluster_name, "enrichment", "metabolic_KEGG_pathways.csv"))
  
  # Create and save plots with display name
  if(nrow(as.data.frame(kegg_enrichment)) > 0) {
    dotplot(kegg_enrichment, showCategory = 20) +
      ggtitle(paste(display_name, "- Metabolic KEGG Pathways"))
    ggsave(file.path("src_output/cluster_de", cluster_name, "enrichment", "metabolic_KEGG_dotplot.png"), 
           width = 12, height = 8)
  }
  
  return(kegg_enrichment)
}

# Main execution script
main_analysis <- function() {
  # Read sample information
  sample_info <- read_csv("input_tables/table_s1_sample_table.csv")
  
  # Create results summary dataframe
  results_summary <- data.frame(
    cluster = character(),
    cell_type = character(),
    sub_lineage = character(),
    n_metabolic_genes = numeric(),
    n_significant_metabolic_genes = numeric(),
    n_kegg_pathways = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Perform analysis for each cluster
  for(i in 1:dim(counts_array)[3]) {
    cluster_name <- paste0("cluster_", i)
    
    # Get cell type and sub_lineage
    cell_type <- if("lineage" %in% colnames(cluster_annot)) {
      cluster_annot$lineage[i]
    } else {
      "Unknown"
    }
    
    sub_lineage <- if("sub_lineage" %in% colnames(cluster_annot)) {
      cluster_annot$sub_lineage[i]
    } else {
      ""
    }
    
    # Extract count matrix for this cluster
    count_matrix <- counts_array[,,i]
    
    # Count metabolic genes before DE analysis
    n_metabolic <- sum(colnames(count_matrix) %in% metabolic_genes)
    
    # Perform differential expression analysis with lineage and sub_lineage
    sig_genes <- perform_de_analysis(count_matrix, sample_info, cluster_name, metabolic_genes, cell_type, sub_lineage)
    
    # Perform pathway enrichment if there are significant genes
    kegg_results <- NULL
    if(!is.null(sig_genes) && nrow(sig_genes) > 0) {
      kegg_results <- perform_pathway_enrichment(sig_genes, cluster_name, cell_type, sub_lineage)
    }
    
    # Create a single row data frame for this cluster
    cluster_summary <- data.frame(
      cluster = cluster_name,
      cell_type = cell_type,
      sub_lineage = sub_lineage,
      n_metabolic_genes = n_metabolic,
      n_significant_metabolic_genes = if(is.null(sig_genes)) 0 else nrow(sig_genes),
      n_kegg_pathways = if(is.null(kegg_results)) 0 else nrow(as.data.frame(kegg_results)),
      stringsAsFactors = FALSE
    )
    
    # Add to results summary
    results_summary <- rbind(results_summary, cluster_summary)
  }
  
  # Create base directory if it doesn't exist
  dir.create("src_output/cluster_de", showWarnings = FALSE, recursive = TRUE)
  
  # Save results summary
  write_csv(results_summary, "src_output/cluster_de/metabolic_analysis_summary.csv")
  
  # Create summary visualization
  if(nrow(results_summary) > 0) {
    ggplot(results_summary, aes(x = cell_type, y = n_significant_metabolic_genes)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Number of Significant Metabolic DEGs by Cell Type",
           x = "Cell Type",
           y = "Number of Significant Metabolic Genes")
    
    ggsave("src_output/cluster_de/metabolic_deg_summary_by_celltype.png", width = 12, height = 6)
  }
  
  # Print summary
  print("Analysis complete. Summary of results:")
  print(results_summary)
}

# Run the analysis
main_analysis() 
