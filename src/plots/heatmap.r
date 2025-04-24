library(dplyr)
library(tidyr)
library(data.table)
library(fs)
library(ComplexHeatmap)
library(circlize)
library(tibble)

cut_dendrogram <- function(hc, h = NULL, min_clusters = 2, max_clusters = 6, return_factor = FALSE) {
  "Cut dendrogram at optimal height to create meaningful clusters
  
  Parameters:
    hc: Hierarchical clustering object from hclust
    h: Height at which to cut (NULL for automatic detection)
    min_clusters: Minimum number of clusters to create
    max_clusters: Maximum number of clusters to create
    return_factor: If TRUE, returns factor of clusters; if FALSE, returns count
  
  Returns:
    Either factor of cluster assignments or number of clusters"
  
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
  
  if (return_factor) {
    return(factor(clusters))
  } else {
    return(length(unique(clusters)))
  }
}

prepare_heatmap_data <- function(data_path, n_top_genes = 50) {
  "Prepare data matrix for heatmap visualization
  
  Parameters:
    data_path: Path to feature importance CSV file
    n_top_genes: Number of top genes to include by absolute importance
  
  Returns:
    Matrix with genes as rows and clusters as columns"
  
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
      message(paste0("No feature importance directories found for ", model_type))
      next
    }
    
    # Process each directory
    for (feat_dir in feature_dirs) {
      # Get CSV files
      csv_files <- dir_ls(feat_dir, regexp = "\\.csv$", fail = FALSE)
      
      if (length(csv_files) == 0) {
        message(paste0("No CSV files found in ", feat_dir))
        next
      }
      
      # Process each file
      for (csv_file in csv_files) {
        file_basename <- path_file(csv_file) %>% path_ext_remove()
        message(paste0("Processing file: ", file_basename))
        
        # Debug print the filename for troubleshooting
        message(paste0("Filename for pattern matching: ", file_basename))
        
        # Default values
        cluster_info <- "all_clusters"  # Default cluster info
        gene_set <- "metabolic"  # Default gene set
        
        if (grepl("lcam_hi", file_basename, ignore.case = TRUE)) {
          cluster_info <- "lcam_hi"
        } else if (grepl("lcam_lo", file_basename, ignore.case = TRUE)) {
          cluster_info <- "lcam_lo"
        } else if (grepl("lcam_both", file_basename, ignore.case = TRUE)) {
          cluster_info <- "lcam_both"
        } else if (grepl("all_clusters", file_basename, ignore.case = TRUE)) {
          cluster_info <- "all_clusters"
        }
        
        # Check for nonmetabolic first (since it contains "metabolic")
        if (grepl("nonmetabolic", file_basename, ignore.case = TRUE)) {
          gene_set <- "nonmetabolic"
          message("Matched nonmetabolic gene set")
        } else if (grepl("metabolic", file_basename, ignore.case = TRUE)) {
          gene_set <- "metabolic"
          message("Matched metabolic gene set")
        } else if (grepl("random", file_basename, ignore.case = TRUE)) {
          gene_set <- "random"
          message("Matched random gene set")
        } else {
          message("WARNING: No gene set matched, using default: metabolic")
        }
        
        message(paste0("Using: ", cluster_info, "/", gene_set))
        
        # Create specific output directories
        specific_output_dir <- file.path(model_output_dir, cluster_info, gene_set)
        dir_create(specific_output_dir, recurse = TRUE)
        
        # Generate comparison heatmaps directly in the main directory
        tryCatch({
          # Clean the basename if it starts with "feature_importance_" to avoid redundancy
          clean_basename <- file_basename
          if (grepl("^feature_importance_", clean_basename)) {
            clean_basename <- sub("^feature_importance_", "", clean_basename)
          }
          
          # Set up the output paths
          fi_heatmap_path <- file.path(specific_output_dir, paste0("fi_", clean_basename, ".png"))
          fc_heatmap_path <- file.path(specific_output_dir, paste0("fc_", clean_basename, ".png"))
          ratio_heatmap_path <- file.path(specific_output_dir, paste0("ratio_", clean_basename, ".png"))
          
          # Load fold change data
          fc_data <- load_fold_changes()
          fold_changes <- fc_data$fold_changes
          gene_fc <- fc_data$gene_fc
          
          # Prepare feature importance data
          message("Loading feature importance data...")
          feature_matrix <- prepare_heatmap_data(csv_file, n_top_genes = 50)
          
          # Create fold change matrix matching feature matrix structure
          fc_matrix <- matrix(0, nrow = nrow(feature_matrix), ncol = ncol(feature_matrix))
          rownames(fc_matrix) <- rownames(feature_matrix)
          colnames(fc_matrix) <- colnames(feature_matrix)
          
          # ALSO CREATE RATIO MATRIX - add these lines
          ratio_matrix <- matrix(0, nrow = nrow(feature_matrix), ncol = ncol(feature_matrix))
          rownames(ratio_matrix) <- rownames(feature_matrix)
          colnames(ratio_matrix) <- colnames(feature_matrix)
          
          # Populate with matching fold changes
          for (i in 1:nrow(feature_matrix)) {
            for (j in 1:ncol(feature_matrix)) {
              gene <- rownames(feature_matrix)[i]
              cluster <- colnames(feature_matrix)[j]
              feature_key <- paste0(gene, "@", cluster)
              
              # Try to find exact feature match
              fc_match <- fold_changes[Feature == feature_key]
              
              if (nrow(fc_match) > 0) {
                fc_matrix[i,j] <- fc_match$Value
              } else {
                # Fall back to gene-level fold change
                gene_match <- gene_fc[gene == gene]
                if (nrow(gene_match) > 0) {
                  fc_matrix[i,j] <- gene_match$value
                } else {
                  # Default to small value if no match found
                  fc_matrix[i,j] <- 0.1
                }
              }
            }
          }
          
          # Print diagnostic info
          cat(sprintf("FC matrix range: [%.4f, %.4f]\n", min(fc_matrix), max(fc_matrix)))
          
          # Replace the ratio calculation with a properly centered fold change distance
          # Calculate ratio matrices with much stronger fold change penalty
          for (i in 1:nrow(feature_matrix)) {
            for (j in 1:ncol(feature_matrix)) {
              fi_value <- feature_matrix[i,j]
              fc_value <- fc_matrix[i,j]  # This is log2 fold change
              
              # Set importance to 0 if absolute log2 fold change is 1 or larger
              if (abs(fc_value) >= 1) {
                ratio_matrix[i,j] <- 0
              } else {
                # Keep original importance if fold change is small
                ratio_matrix[i,j] <- fi_value
              }
            }
          }
          
          # Add diagnostic print statements
          cat(sprintf("Feature matrix range: [%.2e, %.2e]\n", min(feature_matrix), max(feature_matrix)))
          cat(sprintf("Ratio matrix range: [%.2e, %.2e]\n", min(ratio_matrix), max(ratio_matrix)))
          cat(sprintf("Fold changes distribution: min=%.2f, median=%.2f, mean=%.2f, max=%.2f\n", 
                     min(abs(fc_matrix)), median(abs(fc_matrix)), mean(abs(fc_matrix)), max(abs(fc_matrix))))
          
          # Calculate clustering
          row_dist <- dist(feature_matrix, method = "manhattan")
          col_dist <- dist(t(feature_matrix), method = "manhattan")
          row_hc <- hclust(row_dist, method = "complete")
          col_hc <- hclust(col_dist, method = "complete")
          row_clusters <- cut_dendrogram(row_hc, min_clusters = 1, max_clusters = 5, return_factor = FALSE)
          col_clusters <- cut_dendrogram(col_hc, min_clusters = 1, max_clusters = 5, return_factor = FALSE)
          
          # 1. FEATURE IMPORTANCE HEATMAP
          # Get max value to share with ratio heatmap
          fi_max <- max(abs(feature_matrix))
          
          plot_matrix_heatmap(
            data_matrix = feature_matrix,
            title = "Feature\nImportance",
            row_hc = row_hc,
            col_hc = col_hc,
            row_clusters = row_clusters,
            col_clusters = col_clusters,
            output_path = fi_heatmap_path,
            fixed_max = fi_max  # Pass the max value
          )
          
          # 2. FOLD CHANGE HEATMAP - keep its own scale
          plot_matrix_heatmap(
            data_matrix = fc_matrix,
            title = "Log2 Fold\nChange",
            row_hc = row_hc,
            col_hc = col_hc,
            row_clusters = row_clusters,
            col_clusters = col_clusters,
            output_path = fc_heatmap_path
            # No fixed_max, uses its own scale
          )
          
          # Change the title to match the calculation
          plot_matrix_heatmap(
            data_matrix = ratio_matrix,
            title = "Importance to\nFold Change Ratio",
            row_hc = row_hc,
            col_hc = col_hc,
            row_clusters = row_clusters,
            col_clusters = col_clusters,
            output_path = ratio_heatmap_path,
            fixed_max = fi_max  # Use the SAME max as feature importance instead of ratio-specific scaling
          )
          
          # Add Option 3 as a FOURTH heatmap (Key Features highlight)
          highlight_matrix <- matrix(0, nrow = nrow(feature_matrix), ncol = ncol(feature_matrix))
          rownames(highlight_matrix) <- rownames(feature_matrix)
          colnames(highlight_matrix) <- colnames(feature_matrix)
          
          # Set thresholds safely
          if(sum(abs(feature_matrix) > 0) >= 4 && sum(abs(fc_matrix) > 0) >= 4) {
            threshold_importance <- quantile(abs(feature_matrix)[abs(feature_matrix) > 0], 0.75)  # Top 25% of importance
            threshold_fold_change <- quantile(abs(fc_matrix)[abs(fc_matrix) > 0], 0.25)     # Bottom 25% of fold change
          } else {
            threshold_importance <- max(abs(feature_matrix)) * 0.5
            threshold_fold_change <- max(abs(fc_matrix)) * 0.25
          }
          
          # Populate highlight matrix
          for (i in 1:nrow(feature_matrix)) {
            for (j in 1:ncol(feature_matrix)) {
              fi_value <- feature_matrix[i,j]
              fc_value <- fc_matrix[i,j]
              
              # Highlight where importance is high but fold change is low
              if (abs(fi_value) >= threshold_importance && abs(fc_value) <= threshold_fold_change) {
                highlight_matrix[i,j] <- fi_value  # Use importance value for highlighted cells
              } else {
                highlight_matrix[i,j] <- 0  # Others get zero
              }
            }
          }
          
          # Set up the output path for the highlight heatmap
          highlight_heatmap_path <- file.path(specific_output_dir, paste0("highlight_", clean_basename, ".png"))
          
          # 4. HIGHLIGHT HEATMAP - show key features
          plot_matrix_heatmap(
            data_matrix = highlight_matrix,
            title = "Key Features:\nHigh Importance,\nLow Fold Change",
            row_hc = row_hc,
            col_hc = col_hc,
            row_clusters = row_clusters,
            col_clusters = col_clusters,
            output_path = highlight_heatmap_path,
            fixed_max = fi_max  # Use same max as feature importance
          )
          
        }, error = function(e) {
          message(paste0("   Error processing heatmaps: ", e$message))
        })
      }
    }
  }
  
  message("Heatplot generation complete.")
}

# Fast fold change loading with data.table for performance
load_fold_changes <- function() {
  "Load fold change data efficiently using data.table
  
  Loads and processes fold change data from the standard output file.
  Optimizes for performance using data.table indexing.
  
  Returns:
    A list containing:
      - fold_changes: Data table with all fold change entries
      - gene_fc: Data table with gene-level fold changes
      - has_clusters: Boolean indicating if cluster-specific fold changes exist"
  
  fold_changes_path <- "output/differential_expression/feature_fold_changes.csv"
  if (!file.exists(fold_changes_path)) {
    stop("Fold change file not found at: ", fold_changes_path)
  }
  
  # Load fold change data
  message("Loading fold change data...")
  fold_changes <- fread(fold_changes_path, key="Feature")  # Set key for indexing
  message("Loaded ", nrow(fold_changes), " fold change entries")
  
  # Extract gene-level fold changes 
  gene_fc <- fold_changes[!grepl("@", Feature), .(gene = Feature, value = Value)]
  setkey(gene_fc, gene)  # Index for fast lookups
  
  # Return data.tables directly for fastest lookup performance
  message("Fold change mapping ready")
  return(list(
    fold_changes = fold_changes,
    gene_fc = gene_fc,
    has_clusters = any(grepl("@", fold_changes$Feature))
  ))
}

# Helper function for plotting heatmaps with consistent settings
plot_matrix_heatmap <- function(data_matrix, title, row_hc, col_hc, row_clusters, col_clusters, 
                              output_path = NULL, right_annotation = NULL, fixed_max = NULL) {
  "Create a standardized heatmap visualization with consistent formatting
  
  Parameters:
    data_matrix: Numeric matrix to visualize
    title: Title for the heatmap legend
    row_hc: Hierarchical clustering object for rows
    col_hc: Hierarchical clustering object for columns
    row_clusters: Number of row clusters for splitting
    col_clusters: Number of column clusters for splitting
    output_path: Path to save the heatmap image (NULL to not save)
    right_annotation: Optional annotation to display at the right side
    fixed_max: Optional fixed maximum value for color scaling
    
  Returns:
    ComplexHeatmap object"
  
  # Use provided max or calculate from data
  max_abs <- if(!is.null(fixed_max)) fixed_max else max(abs(data_matrix))
  
  col_fun <- colorRamp2(
    c(-max_abs, 0, max_abs),
    c("blue", "white", "red")
  )
  
  # Create heatmap with compatible parameters
  ht <- Heatmap(data_matrix,
    name = title,
    col = col_fun,
    cluster_rows = row_hc,
    cluster_columns = col_hc,
    row_split = row_clusters,
    column_split = col_clusters,
    right_annotation = right_annotation,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    # Use basic border parameters that work with older ComplexHeatmap versions
    gap = unit(1.5, "mm"),
    border = TRUE,
    border_gp = gpar(lwd = 0.8)
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

# Run the function
create_all_heatplots()
