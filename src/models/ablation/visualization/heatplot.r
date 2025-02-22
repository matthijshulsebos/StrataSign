# Required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggdendro)    # Add this
library(grid)
library(gridExtra)   # Add this
library(tidyverse)
library(fs)          # Add this
library(ComplexHeatmap)
library(circlize)  # for color scales

#' Calculate similarity matrix
#' @param data Matrix of values
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @return Similarity matrix
calculate_similarity_matrix <- function(data, method = "pearson") {
  # Check if any rows or columns have zero variance
  row_vars <- apply(data, 1, var, na.rm = TRUE)
  col_vars <- apply(data, 2, var, na.rm = TRUE)
  
  if (any(row_vars == 0) || any(col_vars == 0)) {
    # Add small amount of noise to zero-variance features
    data <- data + matrix(rnorm(nrow(data) * ncol(data), 0, 1e-10), 
                         nrow = nrow(data))
  }
  
  cor_matrix <- cor(data, method = method)
  # Handle NA values in correlation matrix
  cor_matrix[is.na(cor_matrix)] <- 0
  # Convert correlations to distances (1 - correlation)
  dist_matrix <- as.dist(1 - cor_matrix)
  return(dist_matrix)
}

#' Calculate robust distance matrix
#' @param x Input matrix
#' @param method Correlation method
#' @return Distance matrix
calculate_robust_distance <- function(x, method = "pearson") {
  x <- as.matrix(x)
  n <- ncol(x)
  dist_matrix <- matrix(0, n, n)
  
  # Calculate pairwise distances directly
  for(i in 1:n) {
    for(j in i:n) {
      # Get vectors
      vec1 <- x[,i]
      vec2 <- x[,j]
      
      # Calculate cosine similarity
      cos_sim <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
      if(is.na(cos_sim)) cos_sim <- 0
      
      # Calculate magnitude similarity
      mag_ratio <- min(sum(abs(vec1)), sum(abs(vec2))) / max(sum(abs(vec1)), sum(abs(vec1)), sum(abs(vec2)))
      if(is.na(mag_ratio)) mag_ratio <- 0
      
      # Combine both metrics
      similarity <- (cos_sim + 1)/2 * mag_ratio
      dist_matrix[i,j] <- dist_matrix[j,i] <- 1 - similarity
    }
  }
  
  return(as.dist(dist_matrix))
}

#' Perform hierarchical clustering
#' @param dist_matrix Distance matrix
#' @param method Clustering method
#' @return Hierarchical clustering object
perform_clustering <- function(dist_matrix, method = "complete") {
  hclust(dist_matrix, method = method)
}

#' Create dendrogram plot
#' @param hclust_obj Hierarchical clustering object
#' @param direction Direction of the dendrogram ("horizontal" or "vertical")
#' @return ggplot object
create_dendrogram <- function(hclust_obj, direction = "horizontal") {
  dendro_data <- dendro_data(hclust_obj)
  
  if(direction == "horizontal") {
    p <- ggplot() +
      geom_segment(data = segment(dendro_data),
                  aes(x = -y, y = x, xend = -yend, yend = xend)) + # Flip direction
      scale_x_reverse() +  # Reverse scale to point towards heatmap
      coord_cartesian(expand = FALSE) +  # Remove padding
      theme_void()
  } else {
    p <- ggplot() +
      geom_segment(data = segment(dendro_data),
                  aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_cartesian(expand = FALSE) +  # Remove padding
      theme_void()
  }
  return(p)
}

#' Prepare data for heatmap plotting
#' @param data_path Path to input data file
#' @param n_top_genes Number of top genes to display
#' @return List containing processed data and top gene information
#' @export
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
  
  return(list(
    data_matrix = as.matrix(data_matrix)
  ))
}

#' Find natural cluster breaks using dendrogram height
#' @param hc Hierarchical clustering object
#' @return Number of natural clusters
find_natural_clusters <- function(hc, min_height_ratio = 1.25) {  # Lowered from 1.5
  # Get dendrogram heights
  heights <- sort(hc$height, decreasing = TRUE)
  
  # Calculate height ratios between successive levels
  height_ratios <- heights[-length(heights)] / heights[-1]
  
  # Find all significant breaks (where ratio exceeds threshold)
  significant_breaks <- which(height_ratios > min_height_ratio)
  
  if(length(significant_breaks) > 0) {
    # Use number of breaks + 1 for more granular clustering
    n_clusters <- length(significant_breaks) + 1
  } else {
    # Default to 3 clusters if no clear breaks
    n_clusters <- 3
  }
  
  # Cap at higher maximum
  return(min(n_clusters, 8))  # Increased from 5
}

#' Enhanced heatmap plotting function with natural clustering
plot_enhanced_heatmap <- function(data_path, dataset_name, version, analysis_type,
                                n_top_genes = 50, cluster_by_similarity = TRUE,
                                similarity_method = "pearson", output_path = NULL) {
  # Prepare data
  prepared_data <- prepare_heatmap_data(data_path, n_top_genes)
  data_matrix <- as.matrix(prepared_data$data_matrix)
  
  # Calculate sparsity
  sparsity <- mean(data_matrix == 0)
  
  # Define color mapping function
  col_fun <- colorRamp2(
    c(-max(abs(data_matrix)), 0, max(abs(data_matrix))),
    c("blue", "white", "red")
  )
  
  # Calculate distances and clustering
  row_dist <- calculate_robust_distance(t(data_matrix), method = similarity_method)
  col_dist <- calculate_robust_distance(data_matrix, method = similarity_method)
  
  # Create hierarchical clustering
  row_hc <- hclust(row_dist, method = "ward.D2")
  col_hc <- hclust(col_dist, method = "ward.D2")
  
  if(cluster_by_similarity) {
    # Find natural cluster breaks from dendrogram structure
    row_k <- find_natural_clusters(row_hc)
    col_k <- find_natural_clusters(col_hc)
    
    # Debug info
    cat(sprintf("Natural row clusters: %d (based on height differences)\n", row_k))
    cat(sprintf("Natural column clusters: %d (based on height differences)\n", col_k))
  } else {
    row_k <- NULL
    col_k <- NULL
  }
  
  # Create heatmap with natural clustering
  ht <- Heatmap(data_matrix,
    name = "Value",
    col = col_fun,
    cluster_rows = row_hc,
    cluster_columns = col_hc,
    # Use natural cluster breaks from dendrogram
    row_split = if(cluster_by_similarity) row_k else NULL,
    column_split = if(cluster_by_similarity) col_k else NULL,
    
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
  if(is.null(output_path)) {
    output_path <- file.path(output_dir, dataset_name, version,
                            paste0("heatmap_", analysis_type, ".png"))
  }
  
  png(output_path, width = 14, height = 12, units = "in", res = 300)  # Increased dimensions
  draw(ht)
  dev.off()
  
  return(ht)
}

#' Create heatplots for all input files in heatplot directory
#' @param input_dir Base directory containing heatplot input files
#' @param output_dir Base directory for saving heatplot figures
#' @export
create_all_heatplots <- function(input_dir = "data/ablation/intermediates/heatplot",
                                output_dir = "results/ablation/figures/heatplot") {
  # Get all CSV files recursively
  csv_files <- dir_ls(input_dir, recurse = TRUE, glob = "*.csv")
  
  for (input_file in csv_files) {
    # Extract dataset and version from path
    rel_path <- path_rel(input_file, input_dir)
    parts <- path_split(rel_path)[[1]]
    dataset_name <- parts[1]
    version <- parts[2]
    file_name <- parts[3]
    
    # Extract base name without extension
    base_name <- tools::file_path_sans_ext(file_name)
    
    # Create corresponding output directory
    output_subdir <- file.path(output_dir, dataset_name, version)
    dir_create(output_subdir, recurse = TRUE)
    
    # Generate unique output path using the input filename
    output_path <- file.path(
      output_subdir,
      paste0("heatmap_", base_name, ".png")
    )
    
    print(paste("Creating heatplot for:", rel_path))
    print(paste("Saving to:", output_path))
    
    plot_enhanced_heatmap(
      data_path = input_file,
      dataset_name = dataset_name,
      version = version,
      analysis_type = base_name,
      n_top_genes = 50,
      cluster_by_similarity = TRUE,  # Changed to TRUE
      similarity_method = "pearson",
      output_path = output_path
    )
  }
  
  print("Heatplot creation completed successfully!")
}

# Run the function
create_all_heatplots()
