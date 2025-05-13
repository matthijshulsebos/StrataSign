library(dplyr)
library(tidyr)
library(data.table)
library(fs)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# Function to cut a dendrogram to obtain a specified range of clusters
cut_dendrogram <- function(hc, h = NULL, min_clusters = 2, max_clusters = 6, return_factor = FALSE) {
  if (is.null(h)) {
    # Get sorted branch heights from the dendrogram
    heights <- sort(hc$height, decreasing = TRUE) 

    # Calculate relative height gaps
    if (length(heights) > 0) {
        if (length(heights) == 1) {
            height_gaps <- numeric(0)
        } else {
            height_gaps <- diff(heights) / heights[-length(heights)] 
        }
    } else {
        height_gaps <- numeric(0)
    }
      
    # Apply smoothing to height gaps to reduce sensitivity to minor differences
    if (length(height_gaps) >= 5) {
      smooth_gaps <- stats::filter(height_gaps, rep(1/3, 3), sides = 2)
      smooth_gaps[is.na(smooth_gaps)] <- height_gaps[is.na(smooth_gaps)]
    } else {
      # Default to raw gaps if smoothing is not appropriate
      smooth_gaps <- height_gaps
    }
    
    # Examine the top portion of the tree for significant gaps
    effective_length_for_top_n <- length(smooth_gaps)
    top_n <- if (effective_length_for_top_n > 0) max(1, min(effective_length_for_top_n, ceiling(effective_length_for_top_n * 0.33))) else 0


    if (top_n > 0 && effective_length_for_top_n > 0) {
        # Define a significant gap as one that is at least 30% larger than the average
        avg_gap <- mean(smooth_gaps[1:top_n], na.rm = TRUE)
        sig_threshold <- 1.3 * avg_gap
        
        # Find indices of significant gaps
        sig_gaps <- which(smooth_gaps[1:top_n] > sig_threshold)
        
        if (length(sig_gaps) > 0) {
          # Use the height corresponding to the highest significant gap
          max_gap_idx <- sig_gaps[1]
          # We use +1 because the gaps are between the heights
          h <- heights[max_gap_idx + 1]
          cat(sprintf("Significant gap at height %.3f (%.2f times average)\n", 
                    h, smooth_gaps[max_gap_idx]/avg_gap))
        } else {
          # If no significant gap is found use a fixed percentage of the dendrogram height 
          h <- heights[max(1, ceiling(length(heights) * 0.25))]
          cat(sprintf("No significant gaps found, using height %.3f (25th percentile of heights)\n", h))
        }
    } else {
        # Fallback if there are no gaps
        if (length(heights) > 0) {
            # Cut at the highest point
            h <- heights[1] 
            cat(sprintf("Few items/gaps to analyze, using highest merge height: %.3f\n", h))
        } else {
            h <- 0 
            cat("Warning: No heights found in dendrogram\n")
        }
    }
  }
  
  # Cut the tree at the determined height
  clusters <- cutree(hc, h = h)
  n_clusters <- length(unique(clusters))
  
  # Adjust the number of clusters if it falls out of bounds
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

# Function to prepare data
prepare_heatmap_data <- function(data_path, n_top_genes = 50) {
  # Read feature importance csv
  data_df <- fread(data_path, colClasses = c("Feature" = "character", "Value" = "numeric"))
  
  # Transform data from long format
  data_long <- data_df %>% 
    mutate(
      gene = sapply(strsplit(Feature, "@"), `[`, 1),
      cluster = sapply(strsplit(Feature, "@"), `[`, 2),
      value = as.numeric(Value)
    ) %>% 
    filter(!is.na(cluster)) %>%
    select(gene, cluster, value)
  
  # Select top N genes based on their total absolute importance summed across all clusters
  top_genes <- data_long %>% 
    group_by(gene) %>% 
    summarize(total_impact = sum(abs(value), na.rm = TRUE)) %>%
    arrange(desc(total_impact)) %>%
    head(n_top_genes) %>%
    pull(gene)
  
  # Create the final matrix for the heatmap
  data_matrix <- data_long %>%
    filter(gene %in% top_genes) %>%
    pivot_wider(
      names_from = cluster,
      values_from = value,
      values_fill = 0
    ) %>% 
    column_to_rownames("gene")
  
  return(as.matrix(data_matrix))
}

# Main function
create_all_heatplots <- function() {
  # Define base paths
  models_base_dir <- "output/2. models" 
  output_base_dir <- "output/6. plots/heatmap" 
  
  dir_create(output_base_dir, recurse = TRUE)
  
  # Get dataset types
  dataset_types <- dir_ls(models_base_dir, type = "directory") %>%
    path_file()
  
  # If no dataset type directories are found
  if (length(dataset_types) == 0) {
    message(paste0("No dataset type directories (e.g., 'absolute', 'relative') found in ", models_base_dir))
    return()
  }

  # Load fold change data once to be used for all heatmaps
  fc_data_list <- tryCatch({
    load_fold_changes()
  }, error = function(e) {
    # Handle errors during fold change loading
    message(paste0("Critical error loading fold change data: ", e$message))
    message("Heatmap generation will proceed without fold change or ratio plots if feature importance files are found.")
    return(NULL)
  })
  
  # Init fold_changes
  fold_changes <- NULL

  if (!is.null(fc_data_list)) {
    fold_changes <- fc_data_list$fold_changes
  }

  # Iterate over each dataset type
  for (current_dataset_type in dataset_types) {
    message(paste0("\nProcessing dataset type: ", current_dataset_type))
    dataset_type_input_path <- file.path(models_base_dir, current_dataset_type)

    # Find all feature importance CSV files
    all_fi_csv_files <- dir_ls(dataset_type_input_path,
                               recurse = TRUE,
                               regexp = "feature_importance.*\\.csv$",
                               type = "file",
                               fail = FALSE)

    # If no feature importance files are found skip to the next
    if (length(all_fi_csv_files) == 0) {
      message(paste0("  No feature importance CSV files found for dataset type: ", current_dataset_type))
      next
    }

    # Iterate over each found feature importance csv file
    for (csv_file_path in all_fi_csv_files) {
      message(paste0("  Processing file: ", csv_file_path))

      # Parse model name, cluster info, and gene set from the file path
      relative_path_from_dataset_type <- path_rel(csv_file_path, dataset_type_input_path)
      path_parts <- path_split(relative_path_from_dataset_type)[[1]]
      
      model_algorithm <- "unknown_model"
      if (length(path_parts) > 1) {
        model_algorithm <- path_parts[1]
      } else {
         message(paste0("Could not determine model algorithm from path structure for: ", csv_file_path))
      }

      file_basename_no_ext <- path_file(csv_file_path) %>% path_ext_remove()
      
      # String to parse for cluster and gene set info
      parse_string <- sub("^feature_importance_", "", file_basename_no_ext) 
      
      # Default parsed info
      cluster_info_parsed <- "all_clusters"
      gene_set_parsed <- "metabolic" 

      # Parse cluster information from the string
      if (grepl("lcam_hi", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "lcam_hi"
        parse_string <- sub("lcam_hi(?:_)?", "", parse_string, ignore.case = TRUE)
      } else if (grepl("lcam_lo", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "lcam_lo"
        parse_string <- sub("lcam_lo(?:_)?", "", parse_string, ignore.case = TRUE)
      } else if (grepl("lcam_both", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "lcam_both"
        parse_string <- sub("lcam_both(?:_)?", "", parse_string, ignore.case = TRUE)
      } else if (grepl("all_clusters", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "all_clusters"
        parse_string <- sub("all_clusters(?:_)?", "", parse_string, ignore.case = TRUE)
      }

      # Parse gene set information from the remaining string
      if (grepl("nonmetabolic", parse_string, ignore.case = TRUE)) {
        gene_set_parsed <- "nonmetabolic"
      } else if (grepl("metabolic", parse_string, ignore.case = TRUE)) {
        gene_set_parsed <- "metabolic"
      } else if (grepl("random", parse_string, ignore.case = TRUE)) {
        gene_set_parsed <- "random"
      }
      
      message(paste0("Parsed Info: model_algo=", model_algorithm, 
                     ", cluster_info=", cluster_info_parsed, 
                     ", gene_set=", gene_set_parsed))

      # Create specific output directory structure
      specific_output_dir <- file.path(output_base_dir, current_dataset_type, model_algorithm, gene_set_parsed, cluster_info_parsed)
      dir_create(specific_output_dir, recurse = TRUE)
      
      # Used for naming plot output files
      plot_file_identifier <- file_basename_no_ext 

      # Define output paths for different types of heatmaps
      fi_heatmap_path <- file.path(specific_output_dir, paste0("fi_", plot_file_identifier, ".png"))
      fc_heatmap_path <- file.path(specific_output_dir, paste0("fc_", plot_file_identifier, ".png"))
      ratio_heatmap_path <- file.path(specific_output_dir, paste0("ratio_", plot_file_identifier, ".png"))
      
      tryCatch({
        message("Loading feature importance data...")
        # Prepare the feature importance matrix for the heatmap
        feature_matrix <- prepare_heatmap_data(csv_file_path, n_top_genes = 50) 
        
        # Skip if feature matrix is empty or invalid
        if (is.null(feature_matrix) || nrow(feature_matrix) == 0 || ncol(feature_matrix) == 0) {
          message("Skipping heatmap generation due to empty or invalid feature matrix.")
          next
        }

        # Perform hierarchical clustering
        row_dist <- dist(feature_matrix, method = "manhattan")
        col_dist <- dist(t(feature_matrix), method = "manhattan")
        row_hc <- hclust(row_dist, method = "complete")
        col_hc <- hclust(col_dist, method = "complete")
        
        # Determine number of row and column clusters
        max_r_clust <- min(5, nrow(feature_matrix) -1, if(nrow(feature_matrix) < 2) 1 else 5)
        max_c_clust <- min(5, ncol(feature_matrix) -1, if(ncol(feature_matrix) < 2) 1 else 5)

        # Use cut_dendrogram to get cluster counts
        row_clusters <- if (nrow(feature_matrix) > 1) cut_dendrogram(row_hc, min_clusters = 1, max_clusters = max_r_clust, return_factor = FALSE) else 1
        col_clusters <- if (ncol(feature_matrix) > 1) cut_dendrogram(col_hc, min_clusters = 1, max_clusters = max_c_clust, return_factor = FALSE) else 1
        
        # Determine fixed maximum for color scale to ensure consistency across heatmaps
        fi_max <- max(abs(feature_matrix), na.rm = TRUE)
        if (!is.finite(fi_max)) fi_max <- 1

        # Plot Feature Importance (FI) heatmap
        plot_matrix_heatmap(
          data_matrix = feature_matrix, title = "Feature\nImportance",
          row_hc = row_hc, col_hc = col_hc, row_clusters = row_clusters, col_clusters = col_clusters,
          output_path = fi_heatmap_path, fixed_max = fi_max
        )

        # Proceed only if fold change data is available
        if (!is.null(fold_changes)) {
          # Create Fold Change (FC) matrix aligned with the feature_matrix
          fc_matrix <- matrix(0.1, nrow = nrow(feature_matrix), ncol = ncol(feature_matrix))
          rownames(fc_matrix) <- rownames(feature_matrix)
          colnames(fc_matrix) <- colnames(feature_matrix)
          
          # Populate fc_matrix by looking up fold changes
          for (i in 1:nrow(feature_matrix)) {
            for (j in 1:ncol(feature_matrix)) {
              gene_name <- rownames(feature_matrix)[i]
              cluster_name <- colnames(feature_matrix)[j]
              feature_key <- paste0(gene_name, "@", cluster_name)
              
              fc_match <- fold_changes[Feature == feature_key]
              if (nrow(fc_match) > 0) {
                fc_matrix[i,j] <- fc_match$Value[1] 
              }
            }
          }

          # Plot Fold Change (FC) heatmap
          plot_matrix_heatmap(
            data_matrix = fc_matrix, title = "Log2 Fold\nChange",
            row_hc = row_hc, col_hc = col_hc, row_clusters = row_clusters, col_clusters = col_clusters,
            output_path = fc_heatmap_path 
          )

          # Create Ratio matrix (FI values where |log2FC| < 1)
          ratio_matrix <- matrix(0, nrow = nrow(feature_matrix), ncol = ncol(feature_matrix))
          rownames(ratio_matrix) <- rownames(feature_matrix)
          colnames(ratio_matrix) <- colnames(feature_matrix)
          for (i in 1:nrow(feature_matrix)) {
            for (j in 1:ncol(feature_matrix)) {
              if (abs(fc_matrix[i,j]) < 1) {
                ratio_matrix[i,j] <- feature_matrix[i,j]
              }
            }
          }

          # Plot Ratio heatmap
          plot_matrix_heatmap(
            data_matrix = ratio_matrix, title = "FI (where\n|log2FC|<1)", 
            row_hc = row_hc, col_hc = col_hc, row_clusters = row_clusters, col_clusters = col_clusters,
            output_path = ratio_heatmap_path, fixed_max = fi_max
          )
        } else {
          # Fold change data is not available
          message("Fold change data not available. Skipping FC, Ratio, and Highlight heatmaps.")
        }
        
      }, error = function(e) {
        # Error handling for processing a single csv file
        message(paste0("Error processing heatmaps for ", csv_file_path, ": ", e$message))
        message("    Stack trace: ")
        try(print(sys.calls()))
      })
    }
  }
  message("\nHeatplot generation complete.")
}

# Function to load and process fold change data
load_fold_changes <- function() {
  # Define the path to the fold change data file
  fold_changes_path <- "output/4. differential expression/feature_fold_changes.csv"
  
  # Check if the fold change file exists stop if not found
  if (!file.exists(fold_changes_path)) {
    stop("Fold change file not found at: ", fold_changes_path)
  }
  
  # Load fold change data
  message("Loading fold change data...")
  fold_changes <- fread(fold_changes_path, key="Feature")
  
  # Return a list containing the full fold change data
  return(list(
    fold_changes = fold_changes,
    has_clusters = any(grepl("@", fold_changes$Feature))
  ))
}

# Function to create and save a single heatmap using ComplexHeatmap
plot_matrix_heatmap <- function(data_matrix, title, row_hc, col_hc, row_clusters, col_clusters, 
                              output_path = NULL, right_annotation = NULL, fixed_max = NULL) {
  
  # Determine the maximum absolute value for color scaling
  max_abs_val <- if(!is.null(fixed_max) && is.finite(fixed_max)) fixed_max else max(abs(data_matrix), na.rm = TRUE)
  if (!is.finite(max_abs_val) || max_abs_val == 0) max_abs_val <- 1
  
  # Define the color mapping
  col_fun <- colorRamp2(
    c(-max_abs_val, 0, max_abs_val),
    c("blue", "white", "red")
  )
  
  # Create the heatmap object using ComplexHeatmap
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
    gap = unit(1.5, "mm"),
    border = TRUE,
    border_gp = gpar(lwd = 0.8)
  )
  
  # Save the plot to file
  if(!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    png(output_path, width = 12, height = 10, units = "in", res = 300, bg = "transparent")
    draw(ht, background = "transparent")
    dev.off()
  }
  
  return(ht)
}

create_all_heatplots()
