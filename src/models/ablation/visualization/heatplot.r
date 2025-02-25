# Required packages (tidyverse removed as we only need dplyr and tidyr)
library(dplyr)
library(tidyr)
library(data.table)
library(fs)          
library(ComplexHeatmap)
library(circlize)  # for color scales
library(cluster)  # for silhouette analysis
library(tibble)

#' Calculate distance matrix using scaled Euclidean distance
#' @param x Input matrix
#' @param scale_rows Whether to scale rows before distance calculation
#' @return Distance matrix
calculate_robust_distance <- function(x, scale_rows = FALSE) {  
  x <- as.matrix(x)
  
  if(scale_rows) {
    row_sums <- rowSums(abs(x))
    x <- x / ifelse(row_sums > 0, row_sums, 1)
  }
  
  return(dist(t(x), method = "euclidean"))
}

#' Prepare data for heatmap plotting
#' @param data_path Path to input data file
#' @param n_top_genes Number of top genes to display
#' @return Matrix containing processed data
prepare_heatmap_data <- function(data_path, n_top_genes = 50) {
  # Read data with explicit column types
  data_df <- fread(data_path, colClasses = c("Feature" = "character", "Value" = "numeric"))
  
  # Transform to long format and process features
  data_long <- data_df %>% 
    mutate(
      gene = sapply(strsplit(Feature, "@"), `[`, 1),
      cluster = sapply(strsplit(Feature, "@"), `[`, 2),
      value = as.numeric(Value)
    ) %>% 
    filter(!is.na(cluster)) %>% 
    select(gene, cluster, value)
  
  # Select top n genes by absolute value
  top_genes <- data_long %>% 
    group_by(gene) %>% 
    summarize(total_impact = sum(abs(value))) %>% 
    arrange(desc(total_impact)) %>% 
    head(n_top_genes) %>% 
    pull(gene)
  
  # Create matrix for clustering, ensuring all clusters are represented
  data_matrix <- data_long %>%
    filter(gene %in% top_genes) %>%
    complete(gene = top_genes, cluster = unique(cluster), fill = list(value = 0)) %>%
    pivot_wider(
      names_from = cluster,
      values_from = value,
      values_fill = 0
    ) %>% 
    column_to_rownames("gene")
  
  return(as.matrix(data_matrix))  # Simplified return
}

#' Find optimal number of clusters using silhouette analysis
#' @param dist_matrix Distance matrix
#' @param max_k Maximum number of clusters to try
#' @param min_sil_score Minimum silhouette score to consider clustering valid
#' @return Optimal number of clusters (1 if no clear clusters)
find_optimal_clusters <- function(dist_matrix, max_k = 8, min_sil_score = 0.15) {  # Lowered threshold
  # Convert distance matrix to regular matrix if needed
  if(class(dist_matrix)[1] == "dist") {
    dist_matrix <- as.matrix(dist_matrix)
  }
  
  # Check if we have enough data points for clustering
  n <- nrow(dist_matrix)
  if(n < 3) {
    cat("Warning: Too few items to cluster meaningfully\n")
    return(1)  # No split
  }
  
  # Adjust max_k if necessary
  max_k <- min(max_k, n-1)
  
  # Try different numbers of clusters
  sil_widths <- numeric(max_k - 1)
  for(k in 2:max_k) {
    tryCatch({
      # Perform hierarchical clustering
      hc <- hclust(as.dist(dist_matrix), method = "ward.D2")
      clusters <- cutree(hc, k = k)
      
      # Check if we have enough samples in each cluster
      if(any(table(clusters) < 2)) {
        cat(sprintf("Warning: k=%d produces singleton clusters, skipping\n", k))
        sil_widths[k-1] <- -1
        next
      }
      
      # Calculate silhouette width
      sil <- silhouette(clusters, dist_matrix)
      sil_widths[k-1] <- mean(sil[, "sil_width"])
    }, error = function(e) {
      cat(sprintf("Warning: Failed for k=%d: %s\n", k, e$message))
      sil_widths[k-1] <- -1
    })
  }
  
  # Print silhouette analysis results
  cat(sprintf("Silhouette analysis results:\n"))
  for(k in 2:max_k) {
    if(sil_widths[k-1] >= 0) {
      cat(sprintf("k=%d: %.3f\n", k, sil_widths[k-1]))
    } else {
      cat(sprintf("k=%d: failed\n", k))
    }
  }
  
  # Find best k among valid calculations
  valid_k <- which(sil_widths >= 0)
  if(length(valid_k) > 0) {
    best_score <- max(sil_widths[valid_k])
    if(best_score >= min_sil_score) {
      optimal_k <- valid_k[which.max(sil_widths[valid_k])] + 1
      cat(sprintf("Found meaningful clusters (score: %.3f)\n", best_score))
      return(optimal_k)
    }
  }
  
  # If no good clustering found, return 1 (no split)
  cat("No clear clustering structure found - keeping as single group\n")
  return(1)
}

#' Enhanced heatmap plotting function
#' @param data_path Path to input data file
#' @param n_top_genes Number of top genes to display
#' @param output_path Output path for the plot
#' @return Heatmap object
plot_enhanced_heatmap <- function(data_path, n_top_genes = 50, output_path = NULL) {
  # Prepare data
  data_matrix <- prepare_heatmap_data(data_path, n_top_genes)
  
  # Calculate sparsity
  sparsity <- mean(data_matrix == 0)
  
  # Define color mapping function
  col_fun <- colorRamp2(
    c(-max(abs(data_matrix)), 0, max(abs(data_matrix))),
    c("blue", "white", "red")
  )
  
  # Calculate distances and clustering
  row_dist <- calculate_robust_distance(t(data_matrix), scale_rows = FALSE)
  col_dist <- calculate_robust_distance(data_matrix, scale_rows = TRUE)  # Scale for columns
  
  # Create hierarchical clustering
  row_hc <- hclust(row_dist, method = "ward.D2")
  col_hc <- hclust(col_dist, method = "ward.D2")
  
  # Find optimal clusters using silhouette analysis
  row_k <- find_optimal_clusters(row_dist)
  col_k <- find_optimal_clusters(col_dist)
  
  # Set splits to NULL if k=1 (no clear clusters)
  row_split <- if(row_k > 1) row_k else NULL
  col_split <- if(col_k > 1) col_k else NULL
  
  # Create heatmap with natural clustering
  ht <- Heatmap(data_matrix,
    name = "Value",
    col = col_fun,
    cluster_rows = row_hc,
    cluster_columns = col_hc,
    # Use natural cluster breaks from dendrogram only if k > 1
    row_split = row_split,
    column_split = col_split,
    
    # Visual parameters
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    row_dend_width = unit(3, "cm"),
    column_dend_height = unit(3, "cm"),
    row_dend_gp = gpar(lwd = 0.8),
    column_dend_gp = gpar(lwd = 0.8),
    row_gap = unit(4, "mm"),
    column_gap = unit(4, "mm"),
    border = TRUE,
    rect_gp = gpar(col = "darkgrey", lwd = 0.5),
    heatmap_height = unit(8, "inches"),
    heatmap_width = unit(10, "inches")
  )
  
  # Adjust output dimensions to accommodate larger dendrograms
  if(!is.null(output_path)) {
    png(output_path, width = 14, height = 12, units = "in", res = 300)  # Increased dimensions
    draw(ht)
    dev.off()
  }
  
  return(ht)
}

#' Create heatplots for all input files in heatplot directory
create_all_heatplots <- function(input_dir = "data/ablation/intermediates/heatplot",
                                output_dir = "results/ablation/figures/heatplot") {
  # Get all CSV files recursively
  csv_files <- dir_ls(input_dir, recurse = TRUE, glob = "*.csv")
  
  for (input_file in csv_files) {
    rel_path <- path_rel(input_file, input_dir)
    parts <- path_split(rel_path)[[1]]
    
    # Debug path parsing
    message("\nProcessing file: ", rel_path)
    
    # Get filename (always last part)
    file_name <- parts[length(parts)]
    base_name <- tools::file_path_sans_ext(file_name)
    
    # Get output directory path (preserve full directory structure)
    output_subdir <- file.path(output_dir, dirname(rel_path))
    dir_create(output_subdir, recurse = TRUE)
    
    # Generate unique output path
    output_path <- file.path(
      output_subdir,
      paste0("heatmap_", base_name, ".png")
    )
    
    message("Creating heatplot for: ", rel_path)
    message("Saving to: ", output_path)
    
    plot_enhanced_heatmap(
      data_path = input_file,
      n_top_genes = 50,
      output_path = output_path
    )
  }
}

# Run the function
create_all_heatplots()
