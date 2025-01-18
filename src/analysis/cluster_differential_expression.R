# Load required libraries
library(DESeq2)
library(dplyr)
library(tidyverse)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(Matrix)

# Load lung_ldm data
load("data/lung_ldm.rd")

# Read data
raw_counts <- lung_ldm$dataset$umitab
cell_to_cluster <- lung_ldm$dataset$cell_to_cluster
cell_to_sample <- lung_ldm$dataset$cell_to_sample
cluster_annot <- read_csv("input_tables/annots_list.csv")

# Function to aggregate counts by cluster and sample
preprocess_counts <- function(raw_counts, cell_to_cluster, cell_to_sample, cluster_id) {
  # Subset cells belonging to this cluster
  cluster_cells <- which(cell_to_cluster == cluster_id)
  cluster_counts <- raw_counts[, cluster_cells]
  cluster_samples <- cell_to_sample[cluster_cells]
  
  # Aggregate counts by sample
  sample_counts <- sapply(unique(cluster_samples), function(sample) {
    sample_cells <- which(cluster_samples == sample)
    if(length(sample_cells) > 0) {
      # Sum counts across cells for each sample
      rowSums(cluster_counts[, sample_cells, drop = FALSE])
    } else {
      rep(0, nrow(cluster_counts))
    }
  })
  
  # Convert to matrix and name columns
  sample_counts <- as.matrix(sample_counts)
  colnames(sample_counts) <- unique(cluster_samples)
  rownames(sample_counts) <- rownames(raw_counts)
  
  return(sample_counts)
}

# Read and process metabolic genes list
metabolic_genes_df <- read_csv("src_output/KEGG/hsa01100_genes.csv")

# Check if we have SYMBOL column, if not convert from ENTREZID
if("SYMBOL" %in% colnames(metabolic_genes_df)) {
  metabolic_genes <- metabolic_genes_df %>% pull(SYMBOL)
} else if("ENTREZID" %in% colnames(metabolic_genes_df)) {
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
  
  # Filter for metabolic genes
  metabolic_indices <- rownames(count_matrix) %in% metabolic_genes
  count_matrix <- count_matrix[metabolic_indices, ]
  
  print(paste("Number of metabolic genes:", sum(metabolic_indices)))
  
  # Check if we have any metabolic genes
  if(sum(metabolic_indices) == 0) {
    warning(paste("No metabolic genes found in cluster", cluster_name))
    return(NULL)
  }
  
  # Prepare sample metadata
  sample_metadata <- sample_info %>%
    dplyr::select(
      sample_id = 1,    # sample IDs from first column
      condition = 6     # tumor/normal condition from sixth column
    ) %>%
    filter(sample_id %in% colnames(count_matrix)) %>%
    mutate(condition = factor(condition)) %>%
    arrange(match(sample_id, colnames(count_matrix)))
  
  # Print condition distribution
  print("Number of samples per condition:")
  print(table(sample_metadata$condition))
  
  # Pre-filter genes with all zeros
  keep <- rowSums(count_matrix > 0) > 0
  count_matrix <- count_matrix[keep, ]
  
  # Try to run DESeq2 analysis
  tryCatch({
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = round(count_matrix),
      colData = sample_metadata,
      design = ~ condition
    )
    
    # Filter low count genes
    keep <- rowSums(counts(dds) >= 10) >= 5
    dds <- dds[keep,]
    
    # Run DESeq2
    dds <- DESeq(dds)
    res <- results(dds)
    
    # Convert results to data frame
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      arrange(padj)
    
    # Save results
    write_csv(res_df, file.path("src_output/cluster_de", cluster_name, "metabolic_de_results.csv"))
    
    # Create volcano plot
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
    
  }, error = function(e) {
    warning(paste("Error in DESeq2 analysis for cluster", cluster_name, ":", e$message))
    return(NULL)
  })
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
  
  # Get number of clusters correctly
  clusters <- sort(unique(cell_to_cluster))
  print(paste("Total number of clusters to process:", length(clusters)))
  
  # Perform analysis for each cluster
  for(i in clusters) {
    print(paste("\n==== Starting analysis for cluster", i, "===="))
    cluster_name <- paste0("cluster_", i)
    
    # Get cell type and sub_lineage
    cell_type <- if("lineage" %in% colnames(cluster_annot)) {
      print(paste("Lineage:", cluster_annot$lineage[i]))
      cluster_annot$lineage[i]
    } else {
      print("No lineage information found")
      "Unknown"
    }
    
    sub_lineage <- if("sub_lineage" %in% colnames(cluster_annot)) {
      print(paste("Sub-lineage:", cluster_annot$sub_lineage[i]))
      cluster_annot$sub_lineage[i]
    } else {
      print("No sub-lineage information found")
      ""
    }
    
    print("Preprocessing counts...")
    # Preprocess counts for this cluster
    count_matrix <- preprocess_counts(raw_counts, cell_to_cluster, cell_to_sample, i)
    print(paste("Count matrix dimensions:", paste(dim(count_matrix), collapse=" x ")))
    
    # Count metabolic genes before DE analysis
    n_metabolic <- sum(rownames(count_matrix) %in% metabolic_genes)
    print(paste("Number of metabolic genes:", n_metabolic))
    
    print("Running differential expression analysis...")
    # Perform differential expression analysis
    sig_genes <- perform_de_analysis(count_matrix, sample_info, cluster_name, 
                                   metabolic_genes, cell_type, sub_lineage)
    
    print("Running pathway enrichment...")
    # Perform pathway enrichment if there are significant genes
    kegg_results <- NULL
    if(!is.null(sig_genes) && nrow(sig_genes) > 0) {
      kegg_results <- perform_pathway_enrichment(sig_genes, cluster_name, 
                                               cell_type, sub_lineage)
    }
    
    print("Creating summary...")
    # Create summary row
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
    
    print(paste("==== Completed analysis for cluster", i, "====\n"))
    
    # Save intermediate results after each cluster
    write_csv(results_summary, "src_output/cluster_de/metabolic_analysis_summary_intermediate.csv")
  }
  
  print("\nAll clusters processed. Creating final outputs...")
  
  # Save final results summary
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
    
    ggsave("src_output/cluster_de/metabolic_deg_summary_by_celltype.png", 
           width = 12, height = 6)
  }
  
  # Print summary
  print("\nAnalysis complete. Summary of results:")
  print(results_summary)
}

# Run the analysis
main_analysis() 
