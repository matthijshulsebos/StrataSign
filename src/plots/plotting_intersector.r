# Visualization functions for feature intersection results
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(readr)
library(fs)
library(stringr)

# Create a consistent color mapping for all cell types/sublineages
create_sublineage_color_map <- function(data) {
  # Extract all unique sublineages
  all_sublineages <- data %>%
    mutate(
      cluster_part = sapply(strsplit(as.character(feature_id), "@"), function(x) {
        if(length(x) > 1) x[2] else ""
      }),
      sublineage = ifelse(cluster_part == "", "None", gsub("_[0-9]+$", "", cluster_part))
    ) %>%
    pull(sublineage) %>%
    unique()
  
  # Set1 palette can handle up to 9 distinct colors
  # If we have more than 9, we'll use multiple palettes
  if (length(all_sublineages) <= 9) {
    colors <- brewer.pal(max(3, length(all_sublineages)), "Set1")
  } else {
    # Combine Set1, Set2, and Set3 for more colors
    colors <- c(
      brewer.pal(9, "Set1"),
      brewer.pal(8, "Set2"),
      brewer.pal(12, "Set3")
    )[1:length(all_sublineages)]
  }
  
  # Create named vector for mapping
  color_map <- setNames(colors, all_sublineages)
  return(color_map)
}

#' Plot intersector analysis results from the comprehensive output file
#'
#' @param output_dir Directory containing the intersector output files
#' @param figures_dir Directory where visualizations will be saved
#' @param top_n Number of top features to include in visualizations
#'
plot_intersector_results <- function(output_dir = "output/intersector", 
                                    figures_dir = "output/figures/intersector",
                                    top_n = 50) {
  
  # Create figures directory if it doesn't exist
  dir_create(figures_dir, recursive = TRUE)
  
  # Load the comprehensive analysis file
  feature_analysis_path <- file.path(output_dir, "feature_analysis.csv")
  
  if (!file.exists(feature_analysis_path)) {
    stop("Comprehensive feature analysis file not found: ", feature_analysis_path)
  }
  
  # Load all feature data from the single file
  all_results <- read_csv(feature_analysis_path, show_col_types = FALSE)
  message("Loaded comprehensive analysis data with ", nrow(all_results), " features")
  
  # Generate the color map for all sublineages
  sublineage_colors <- create_sublineage_color_map(all_results)
  
  # Extract model names from column names - handle both prefixes
  importance_cols <- names(all_results)[grepl("^importance_", names(all_results))]
  if (length(importance_cols) == 0) {
    # Try without the prefix if the new format doesn't use it
    model_names <- names(all_results)[!names(all_results) %in% 
                                      c("feature_id", "gene", "n_models", "n_datasets", 
                                        "n_models_nonzero", "meta_score", "sum_rank", 
                                        "sparsity_bonus", "model_breadth", "is_top",
                                        "version_filter", "dataset_filter", "creation_date")]
    importance_cols <- model_names
  } else {
    model_names <- sub("^importance_", "", importance_cols)
  }
  
  message("Found data for ", length(model_names), " models: ", paste(model_names, collapse=", "))
  
  # Get top features
  top_50 <- all_results %>% 
    filter(is_top) %>% 
    head(top_n)
  
  # 2. FEATURE IMPORTANCE HEATMAP
  # Create matrix for heatmap directly from the importance columns
  heatmap_data <- as.matrix(top_50[, importance_cols])
  rownames(heatmap_data) <- top_50$feature_id
  colnames(heatmap_data) <- model_names
  
  # Create annotation for meta-score
  meta_scores <- top_50$meta_score
  names(meta_scores) <- top_50$feature_id
  
  ha_row <- HeatmapAnnotation(
    "Meta-score" = anno_barplot(meta_scores, width = unit(2, "cm")),
    which = "row"
  )
  
  # Create the heatmap with transparent background
  png(file.path(figures_dir, "top_features_heatmap.png"), 
      width = 14, height = 12, units = "in", res = 300, bg = "transparent")
  
  # Fire gradient: yellow-orange-red
  color_palette <- colorRamp2(
    c(0, max(heatmap_data)/3, max(heatmap_data)*2/3, max(heatmap_data)), 
    c("#FFFDE7", "#FFA726", "#F57C00", "#C62828")
  )
  
  hm <- Heatmap(
    heatmap_data,
    name = "Importance",
    col = color_palette,
    row_names_side = "left",
    row_names_max_width = unit(10, "cm"),  # Allow space for long feature IDs
    row_names_gp = gpar(fontsize = 9),     # Slightly smaller font for better fit
    cluster_rows = FALSE,
    cluster_columns = TRUE, 
    show_row_dend = FALSE,
    show_column_dend = TRUE,
    row_title = paste0("Top ", nrow(top_50), " Features"),
    column_title = "Models",
    right_annotation = ha_row,
    rect_gp = gpar(col = "darkgray", lwd = 0.1),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  # Draw with transparent background
  draw(hm, background = NA)
  dev.off()
  message("Created feature importance heatmap")
  
  # 3. LOLLIPOP CHART FOR TOP FEATURES
  top_50_lollipop <- top_50 %>%
    # Extract gene part for all features
    mutate(
      gene_part = sapply(strsplit(as.character(feature_id), "@"), `[`, 1),
      cluster_part = sapply(strsplit(as.character(feature_id), "@"), function(x) {
        if(length(x) > 1) x[2] else ""
      })
    ) %>%
    # Extract sublineage from cluster part
    mutate(
      # Get sublineage part (before the underscore)
      sublineage = ifelse(cluster_part == "", "None", 
                          gsub("_[0-9]+$", "", cluster_part)),
      # Create unique labels
      short_name = paste0(gene_part, ifelse(cluster_part == "", "", paste0("@", cluster_part))),
      # Make sure labels are unique
      short_name = make.unique(short_name, sep = "_")
    ) %>%
    # Create ordered factor for plotting
    mutate(short_name = factor(short_name, levels = short_name[order(meta_score)]))
  
  p3 <- ggplot(top_50_lollipop, aes(x = short_name, y = meta_score)) +
    geom_segment(aes(xend = short_name, yend = 0, color = sublineage), size = 1.2) +
    geom_point(aes(size = n_models, color = sublineage), alpha = 0.8) +
    scale_color_manual(values = sublineage_colors) +  # Use our fixed color map
    scale_size_continuous(range = c(2, 6)) +
    labs(
      title = paste0("Top ", nrow(top_50), " Features by Meta-Score"),
      x = "Feature",
      y = "Meta-Score",
      color = "Sublineage",
      size = "# Models"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )
  
  ggsave(file.path(figures_dir, "top_features_lollipop.png"), p3, width = 10, height = 14)
  message("Created feature lollipop plot")
  
  # Save this color map for the batched lollipop plots to use
  saveRDS(sublineage_colors, file.path(output_dir, "sublineage_colors.rds"))
  
  # 4. MODEL PREVALENCE BARPLOT
  p4 <- ggplot(top_50, aes(x = reorder(feature_id, meta_score), y = n_models, fill = meta_score)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c() +
    labs(
      title = "Number of Models Where Feature Appears",
      x = "Feature",
      y = "Number of Models",
      fill = "Meta-Score"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.title.x = element_blank()
    )
  
  ggsave(file.path(figures_dir, "model_prevalence.png"), p4, width = 14, height = 6)
  message("Created model prevalence plot")
  
  # 5. FEATURE AGREEMENT DOTPLOT WITH FULL FEATURE NAMES
  # Handle different column naming conventions in the new approach
  raw_importance_cols <- names(all_results)[grepl("^raw_importance_", names(all_results))]
  if (length(raw_importance_cols) == 0) {
    # If not using prefixes, fall back to the same columns as the heatmap
    raw_importance_cols <- importance_cols
  }
  
  # Reshape data to long format for the agreement plot
  agreement_data <- top_50 %>%
    select(feature_id, all_of(raw_importance_cols)) %>%
    pivot_longer(
      cols = all_of(raw_importance_cols),
      names_to = "model",
      values_to = "importance"
    ) %>%
    # Clean up model names if they have prefixes
    mutate(model = gsub("^raw_importance_", "", model))
  
  # Get top 20 for a clearer visualization
  top_20_ids <- top_50$feature_id[1:min(20, nrow(top_50))]
  top_20_agreement <- agreement_data %>%
    filter(feature_id %in% top_20_ids) %>%
    mutate(
      gene_part = sapply(strsplit(as.character(feature_id), "@"), `[`, 1),
      cluster_part = sapply(strsplit(as.character(feature_id), "@"), function(x) if(length(x) > 1) x[2] else ""),
      short_name = ifelse(
        cluster_part == "", 
        gene_part, 
        paste0(gene_part, "\n@", cluster_part)
      ),
      short_name = factor(
        short_name, 
        levels = unique(short_name[match(top_20_ids, feature_id)])
      )
    )
  
  p5 <- ggplot(top_20_agreement %>% filter(importance > 0), 
              aes(x = model, y = short_name, size = importance, color = importance)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(range = c(2, 8)) +
    scale_color_viridis_c(option = "plasma") +
    labs(
      title = "Feature Importance Agreement Across Models (Top 20 Features)",
      x = "Model",
      y = "Feature",
      size = "Importance",
      color = "Importance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "gray90")
    )
  
  ggsave(file.path(figures_dir, "feature_agreement.png"), p5, width = 12, height = 12)
  message("Created feature agreement plot")
  
  message("All visualizations created in: ", figures_dir)
}

#' Create multiple lollipop plots for different ranges of top features
#'
#' @param output_dir Directory containing the intersector output files
#' @param figures_dir Directory where visualizations will be saved
#' @param batch_size Number of features in each plot
#' @param max_features Total number of features to plot
#'
create_batched_lollipop_plots <- function(output_dir = "output/intersector", 
                                         figures_dir = "output/figures/intersector/lollipop_batches",
                                         batch_size = 50,
                                         max_features = 300) {
  
  # Create output directory
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load the comprehensive analysis file
  feature_analysis_path <- file.path(output_dir, "feature_analysis.csv")
  
  if (!file.exists(feature_analysis_path)) {
    stop("Comprehensive feature analysis file not found: ", feature_analysis_path)
  }
  
  # Load all feature data
  all_results <- read_csv(feature_analysis_path, show_col_types = FALSE)
  message("Loaded comprehensive analysis data with ", nrow(all_results), " features")
  
  # Try to load saved color map, or create a new one if it doesn't exist
  color_map_path <- file.path(output_dir, "sublineage_colors.rds")
  if (file.exists(color_map_path)) {
    sublineage_colors <- readRDS(color_map_path)
  } else {
    sublineage_colors <- create_sublineage_color_map(all_results)
  }
  
  # Filter for top features only
  top_features <- all_results %>% 
    arrange(desc(meta_score)) %>%  # Sort by meta_score
    head(max_features)           # Take top features regardless of is_top flag
  
  # Get the maximum meta_score for consistent y-axis scaling across all plots
  max_meta_score <- max(top_features$meta_score)
  message(paste("Maximum meta-score:", round(max_meta_score, 3)))
  
  # Add small buffer to y-axis maximum
  y_axis_max <- max_meta_score * 1.05
  
  # Number of batches
  num_batches <- ceiling(max_features / batch_size)
  
  # Create each batch of lollipop plots
  for (batch_num in 1:num_batches) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, nrow(top_features))
    
    if (start_idx > nrow(top_features)) break
    
    # Get features for this batch
    batch_features <- top_features[start_idx:end_idx, ]
    
    # Create the feature labels
    batch_lollipop <- batch_features %>%
      mutate(
        gene_part = sapply(strsplit(as.character(feature_id), "@"), `[`, 1),
        cluster_part = sapply(strsplit(as.character(feature_id), "@"), function(x) {
          if(length(x) > 1) x[2] else ""
        })
      ) %>%
      # Extract sublineage from cluster part
      mutate(
        # Get sublineage part (before the underscore)
        sublineage = ifelse(cluster_part == "", "None", 
                          gsub("_[0-9]+$", "", cluster_part)),
        # Create unique labels
        short_name = paste0(gene_part, ifelse(cluster_part == "", "", paste0("@", cluster_part))),
        # Make sure labels are unique
        short_name = make.unique(short_name, sep = "_")
      ) %>%
      # Create ordered factor for plotting
      mutate(short_name = factor(short_name, levels = short_name[order(meta_score)]))
    
    # Create the lollipop plot with consistent y-axis scale
    p <- ggplot(batch_lollipop, aes(x = short_name, y = meta_score)) +
      geom_segment(aes(xend = short_name, yend = 0, color = sublineage), size = 1.2) +
      geom_point(aes(size = n_models, color = sublineage), alpha = 0.8) +
      scale_color_manual(values = sublineage_colors) +  # Use consistent colors
      scale_size_continuous(range = c(2, 6)) +
      # Use consistent y-axis limits across all plots
      scale_y_continuous(limits = c(0, y_axis_max)) +
      labs(
        title = paste0("Features Ranked ", start_idx, "-", end_idx, " by Meta-Score"),
        x = "Feature",
        y = "Meta-Score",
        color = "Sublineage",
        size = "# Models"
      ) +
      coord_flip() +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        legend.position = "right"
      )
    
    # Save the plot
    filename <- paste0("top_features_lollipop_", start_idx, "_", end_idx, ".png")
    ggsave(file.path(figures_dir, filename), p, width = 10, height = 14)
    message("Created lollipop plot for features ranked ", start_idx, "-", end_idx)
  }
}

plot_cumulative_importance <- function(output_dir = "output/models",
                                      figures_dir = "output/figures/intersector/cumulative",
                                      version_filter = "metabolic",
                                      threshold = 0.7) {
  
  # Create output directory
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Find all feature importance files
  all_fi_files <- dir_ls(output_dir, 
                         recurse = TRUE, 
                         regexp = "feature_importance.*\\.csv$", 
                         fail = FALSE)
  
  # Filter for specific version if specified
  if (!is.null(version_filter)) {
    pattern <- paste0("/", version_filter, "/")
    all_fi_files <- all_fi_files[grep(pattern, all_fi_files, fixed = FALSE)]
  }
  
  # Extract metadata from file paths
  file_metadata <- data.frame(
    file_path = all_fi_files,
    model_type = sapply(all_fi_files, function(path) {
      parts <- strsplit(path_rel(path, output_dir), "/|\\\\")[[1]]
      return(parts[1])
    }),
    dataset = sapply(all_fi_files, function(path) {
      parts <- strsplit(path_rel(path, output_dir), "/|\\\\")[[1]]
      return(parts[2])
    }),
    version = sapply(all_fi_files, function(path) {
      parts <- strsplit(path_rel(path, output_dir), "/|\\\\")[[1]]
      return(parts[3])
    }),
    stringsAsFactors = FALSE
  )
  
  # For each model type, create cumulative importance plot
  for (model in unique(file_metadata$model_type)) {
    # Get files for this model
    model_files <- file_metadata[file_metadata$model_type == model, ]
    
    for (i in 1:nrow(model_files)) {
      # Load feature importance data
      fi_data <- read_csv(model_files$file_path[i], show_col_types = FALSE)
      
      # Skip if empty
      if (nrow(fi_data) == 0) next
      
      # Extract dataset and version info
      current_dataset <- model_files$dataset[i]
      current_version <- model_files$version[i]
      
      # Ensure consistent column names
      if (!("Feature" %in% names(fi_data))) {
        fi_data <- fi_data %>% rename(Feature = 1, Value = 2)
      }
      
      # Take absolute values of importance scores for models like LASSO
      fi_data$Value <- abs(fi_data$Value)
      
      # Sort by importance (descending)
      fi_data <- fi_data %>% 
        arrange(desc(Value))
      
      # Calculate cumulative importance
      total_importance <- sum(fi_data$Value)
      fi_data$cumulative <- cumsum(fi_data$Value)
      fi_data$cumulative_pct <- fi_data$cumulative / total_importance
      fi_data$feature_index <- 1:nrow(fi_data)
      
      # Calculate density information
      total_features <- nrow(fi_data)
      nonzero_features <- sum(fi_data$Value > 0)
      density <- nonzero_features / total_features
      is_dense <- density > 0.1  # Define is_dense before using it
      
      # Add origin point (0,0)
      origin_point <- data.frame(
        Feature = "ORIGIN", 
        Value = 0, 
        cumulative = 0, 
        cumulative_pct = 0, 
        feature_index = 0
      )
      
      # IMPORTANT: Combine with original data and sort only ONCE
      fi_data <- rbind(origin_point, fi_data) %>%
        arrange(feature_index)
      
      # Find threshold index
      threshold_idx <- which(fi_data$cumulative_pct >= threshold)[1]
      if (is.na(threshold_idx)) threshold_idx <- nrow(fi_data)
      
      # Find where cumulative percentage reaches 100% (or very close)
      max_idx <- which(fi_data$cumulative_pct >= 0.999)[1]
      if (is.na(max_idx)) max_idx <- nrow(fi_data)
      
      # Get feature count for title
      feature_count <- ifelse(is_dense, 
                              paste0(threshold_idx, " of ", nonzero_features, " non-zero features (", 
                                    round(threshold_idx / nonzero_features * 100), "%)"),
                              paste0(nonzero_features, " non-zero features"))
      
      # Create plot
      p <- ggplot(fi_data, aes(x = feature_index, y = cumulative_pct)) +
        # Change from blue (#1976D2) to black
        geom_line(color = "black", size = 1) +
        # Red horizontal threshold line (keep as is)
        geom_hline(yintercept = threshold, linetype = "dashed", color = "#D32F2F", size = 1) +
        # Change from green (#388E3C) to red - use a different red shade for distinction
        geom_vline(xintercept = threshold_idx, linetype = "dashed", color = "#E53935", size = 0.8) +
        annotate("text", x = threshold_idx * 0.5, y = threshold + 0.05, 
                label = paste0(round(threshold * 100), "% threshold"), color = "#D32F2F") +
        # Keep text horizontal and black as requested previously
        annotate("text", x = threshold_idx + min(threshold_idx * 0.1, 5), 
                y = threshold/2, 
                label = paste0("Top ", threshold_idx, " features"), 
                color = "black",
                hjust = 0) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        scale_x_continuous(limits = c(0, max_idx * 1.05)) +
        labs(
          title = paste0(model, " Cumulative Feature Importance"),
          subtitle = paste0("Dataset: ", current_dataset, " / ", current_version, "\nModel density: ", 
                           round(density * 100), "% (", feature_count, ")"),
          x = "Feature Rank",
          y = "Cumulative Importance"
        ) +
        theme_minimal() +
        theme(
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 10),
          plot.title = element_text(face = "bold")
        )
      
      # Save plot
      filename <- paste0("cumulative_importance_", model, "_", current_dataset, "_", current_version, ".png")
      ggsave(file.path(figures_dir, filename), p, width = 10, height = 6, dpi = 300)
      message("Created cumulative importance plot for ", model, " / ", current_dataset, " / ", current_version)
    }
  }
}

# Modify the run section at the end of the file to include the new function
plot_intersector_results()  # This creates the original top 50 plot
plot_cumulative_importance()
create_batched_lollipop_plots(batch_size = 50, max_features = 300)
