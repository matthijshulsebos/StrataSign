# Required packages
library(dplyr)
library(tidyr)
library(data.table)
library(fs)
library(ComplexHeatmap)
library(circlize)
library(tibble)


#' Cut dendrogram at optimal height to create meaningful clusters
#' @param hc Hierarchical clustering object
#' @param h Height at which to cut (NULL for automatic)
#' @param min_clusters Minimum number of clusters
#' @param max_clusters Maximum number of clusters
#' @param return_factor Whether to return factor or number
#' @return Either cluster factor or number of clusters
cut_dendrogram <- function(hc, h = NULL, min_clusters = 2, max_clusters = 6, return_factor = FALSE) {
  if (is.null(h)) {
    # Calculate height gaps
    heights <- sort(hc$height, decreasing = TRUE)
    
    if (length(heights) >= 4) {
      # Look for significant gaps using more robust method
      height_gaps <- diff(heights) / heights[-length(heights)]
      
      # Apply smoothing to reduce sensitivity to small fluctuations
      if (length(height_gaps) >= 5) {
        smooth_gaps <- stats::filter(height_gaps, rep(1/3, 3), sides = 2)
        smooth_gaps[is.na(smooth_gaps)] <- height_gaps[is.na(smooth_gaps)]
      } else {
        smooth_gaps <- height_gaps
      }
      
      # Look at top 33% of tree for significant gaps
      top_n <- max(3, ceiling(length(heights) * 0.33))
      
      # Find significant gap (at least 30% larger than average gap)
      avg_gap <- mean(smooth_gaps[1:top_n])
      sig_threshold <- 1.3 * avg_gap
      
      sig_gaps <- which(smooth_gaps[1:top_n] > sig_threshold)
      
      if (length(sig_gaps) > 0) {
        # Use the highest significant gap
        max_gap_idx <- sig_gaps[1]
        h <- heights[max_gap_idx + 1]
        cat(sprintf("Significant gap at height %.3f (%.2f times average)\n", 
                  h, smooth_gaps[max_gap_idx]/avg_gap))
      } else {
        # If no significant gap, use fixed percentage height
        h <- heights[ceiling(length(heights) * 0.25)]
        cat(sprintf("No significant gaps found, using height %.3f (25%%)\n", h))
      }
    } else {
      # For small dendrograms use percentage of max height
      h <- 0.7 * max(heights)
      cat(sprintf("Using height: %.3f (70%% of max)\n", h))
    }
  }
  
  # Cut tree at specified height
  clusters <- cutree(hc, h = h)
  n_clusters <- length(unique(clusters))
  
  # Adjust if needed
  if (n_clusters < min_clusters) {
    clusters <- cutree(hc, k = min_clusters)
    cat(sprintf("Adjusted to %d clusters (minimum)\n", min_clusters))
  } else if (n_clusters > max_clusters) {
    clusters <- cutree(hc, k = max_clusters)
    cat(sprintf("Adjusted to %d clusters (maximum)\n", max_clusters))
  }
  
  # Print cluster sizes
  sizes <- table(clusters)
  cat(sprintf("Final clusters (%d): %s\n", 
              length(sizes), paste(sizes, collapse = ", ")))
  
  if (return_factor) {
    return(factor(clusters))
  } else {
    return(length(unique(clusters)))
  }
}

#' Prepare data matrix for heatmap visualization
#' @param data_path Path to feature importance CSV
#' @param n_top_genes Number of top genes to include
#' @return Preprocessed data matrix
prepare_heatmap_data <- function(data_path, n_top_genes = 50) {
  # Read data
  data_df <- fread(data_path, colClasses = c("Feature" = "character", "Value" = "numeric"))
  
  # Extract gene and cluster information
  data_long <- data_df %>% 
    mutate(
      gene = sapply(strsplit(Feature, "@"), `[`, 1),
      cluster = sapply(strsplit(Feature, "@"), `[`, 2),
      value = as.numeric(Value)
    ) %>% 
    filter(!is.na(cluster)) %>% 
    select(gene, cluster, value)
  
  # Print value range information
  cat(sprintf("Original data range: [%.2e, %.2e]\n", min(data_long$value), max(data_long$value)))
  
  # Select top genes by absolute importance
  top_genes <- data_long %>% 
    group_by(gene) %>% 
    summarize(total_impact = sum(abs(value))) %>% 
    arrange(desc(total_impact)) %>% 
    head(n_top_genes) %>% 
    pull(gene)
  
  # Create matrix
  data_matrix <- data_long %>%
    filter(gene %in% top_genes) %>%
    pivot_wider(
      names_from = cluster,
      values_from = value,
      values_fill = 0  # Use zeros for missing combinations
    ) %>% 
    column_to_rownames("gene")
  
  return(as.matrix(data_matrix))
}

#' Generate heatmap with Manhattan distance clustering
#' @param data_path Path to feature importance file
#' @param n_top_genes Number of top genes to include
#' @param output_path Output path for the heatmap image
#' @return Heatmap object
plot_enhanced_heatmap <- function(data_path, n_top_genes = 50, output_path = NULL) {
  # Prepare data
  cat("\nProcessing file:", basename(data_path), "\n")
  data_matrix <- prepare_heatmap_data(data_path, n_top_genes)
  
  # Create color scale
  max_abs <- max(abs(data_matrix))
  col_fun <- colorRamp2(
    c(-max_abs, 0, max_abs),
    c("blue", "white", "red")
  )
  
  # Print matrix properties
  cat(sprintf("Matrix dimensions: %d rows x %d columns\n", 
              nrow(data_matrix), ncol(data_matrix)))
  
  # Use Manhattan distance - less sensitive to outliers
  row_dist <- dist(data_matrix, method = "manhattan")
  col_dist <- dist(t(data_matrix), method = "manhattan")
  
  # Use complete linkage - better boundary separation for feature importance
  # Complete linkage emphasizes separation between positive/negative groups
  row_hc <- hclust(row_dist, method = "complete")
  col_hc <- hclust(col_dist, method = "complete")
  
  # Dynamically determine optimal number of clusters
  cat("Determining optimal clusters based on dendrogram structure:\n")
  row_clusters <- cut_dendrogram(row_hc, min_clusters = 1, max_clusters = 5)
  col_clusters <- cut_dendrogram(col_hc, min_clusters = 1, max_clusters = 5)
  
  cat(sprintf("Optimal clustering: %d row clusters and %d column clusters\n", 
              row_clusters, col_clusters))
  
  # Create heatmap
  ht <- Heatmap(data_matrix,
    name = "Feature\nImportance",
    col = col_fun,
    cluster_rows = row_hc,
    cluster_columns = col_hc,
    row_split = row_clusters,
    column_split = col_clusters,
    
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm")
  )
  
  # Save the plot with transparent background
  if(!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    png(output_path, width = 12, height = 10, units = "in", res = 300, bg = "transparent")
    draw(ht, background = "transparent")
    dev.off()
    cat("Created:", output_path, "\n")
  }
  
  return(ht)
}

#' Create heatplots for all feature importance files
#' @return NULL
create_all_heatplots <- function() {
  # Base directories
  models_base_dir <- "output/models"
  output_base_dir <- "output/figures/heatmap"
  
  # Create output directory
  dir_create(output_base_dir, recurse = TRUE)
  
  # Get model types
  model_types <- dir_ls(models_base_dir, type = "directory") %>%
    path_file()
  
  # Process each model type
  for (model_type in model_types) {
    message(paste0("Processing model type: ", model_type))
    
    # Create output directory
    model_output_dir <- file.path(output_base_dir, model_type)
    dir_create(model_output_dir, recurse = TRUE)
    
    # Find feature importance directories
    feature_dirs <- dir_ls(file.path(models_base_dir, model_type), 
                          type = "directory", 
                          recurse = TRUE, 
                          regexp = "feature_importance$",
                          fail = FALSE)
    
    if (length(feature_dirs) == 0) {
      message(paste0("  No feature importance directories found for ", model_type))
      next
    }
    
    # Process each directory
    for (feat_dir in feature_dirs) {
      # Get CSV files
      csv_files <- dir_ls(feat_dir, regexp = "\\.csv$", fail = FALSE)
      
      if (length(csv_files) == 0) {
        message(paste0("  No CSV files found in ", feat_dir))
        next
      }
      
      # Process each file
      for (csv_file in csv_files) {
        file_basename <- path_file(csv_file) %>% path_ext_remove()
        message(paste0("  Processing file: ", file_basename))
        
        output_file <- file.path(model_output_dir, 
                                paste0("heatmap_", file_basename, ".png"))
        
        # Generate heatmap
        tryCatch({
          plot_enhanced_heatmap(
            data_path = csv_file,
            n_top_genes = 50,
            output_path = output_file
          )
        }, error = function(e) {
          message(paste0("   Error processing ", file_basename, ": ", e$message))
        })
      }
    }
  }
  
  message("Heatplot generation complete.")
}

# Run the function
create_all_heatplots()
