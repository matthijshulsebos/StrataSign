library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer) 

# Create a consistent color mapping for all cell types
create_sublineage_color_map <- function(data) {
  # Expects to have a sublineage column
  if (!"sublineage" %in% names(data)) {
    message("Error in create_sublineage_color_map: 'sublineage' column not found in input data.")
    # Fallback if 'sublineage' is missing but 'feature_id' is present
    if ("feature_id" %in% names(data)) {
        message("Attempting to derive sublineage from 'feature_id' as a fallback.")
        data <- data %>%
            mutate(
                cluster_part_temp = sapply(strsplit(as.character(feature_id), "@"), function(x) if(length(x) > 1) x[2] else ""),
                sublineage = ifelse(cluster_part_temp == "", "None", gsub("_[0-9]+(_[a-zA-Z]+)?$", "", cluster_part_temp))
            )
    } else {
        message("Cannot derive sublineage information. Color map will be incomplete or NULL.")
        return(NULL)
    }
  }

  all_sublineages <- data %>%
    pull(sublineage) %>%
    unique() %>%
    sort() 

  if (length(all_sublineages) == 0 || (length(all_sublineages) == 1 && all_sublineages[1] == "None" && !any(all_sublineages != "None"))) {
    message("No distinct sublineages found (or only 'None') to create a color map.")
    if (length(all_sublineages) == 1 && all_sublineages[1] == "None") {
        return(setNames("grey50", "None")) # Provide a color for "None" if it's the only one
    }
    return(NULL)
  }
  
  num_colors_needed <- length(all_sublineages)
  available_colors <- character(0)
  
  # Color palette generation logic (same as in intersector.R for consistency)
  if (num_colors_needed <= 9 && num_colors_needed > 0) {
    available_colors <- brewer.pal(max(3, num_colors_needed), "Set1")
  } else if (num_colors_needed <= 17 && num_colors_needed > 0) { 
    available_colors <- c(brewer.pal(9, "Set1"), brewer.pal(max(3, num_colors_needed - 9), "Set2"))
  } else if (num_colors_needed <= 29 && num_colors_needed > 0) { 
    available_colors <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(max(3, num_colors_needed - 17), "Set3"))
  } else if (num_colors_needed > 0) { 
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    available_colors <- rep(col_vector, length.out = num_colors_needed) 
  } else {
     message("No colors needed for sublineage map (num_colors_needed is 0).")
     return(NULL) 
  }
  
  if (length(available_colors) < num_colors_needed && num_colors_needed > 0) {
      available_colors <- rep(available_colors, length.out = num_colors_needed)
  }

  color_map <- setNames(available_colors[1:num_colors_needed], all_sublineages)
  return(color_map)
}


plot_top_features_lollipop <- function(top_n_data, sublineage_colors, figures_dir) {
  # top_n_data comes from meta_scores.csv and should contain:
  # gene, sublineage, meta_score, n_models_occur
  if (nrow(top_n_data) == 0) {
    message("No data for top features lollipop plot.")
    return()
  }
  if (is.null(sublineage_colors) || length(sublineage_colors) == 0) {
    message("Sublineage colors not available for lollipop plot. Using default gray.")
    # Create a fallback color map if none is provided
    unique_sublineages_in_data <- unique(top_n_data$sublineage)
    sublineage_colors <- setNames(rep("gray", length(unique_sublineages_in_data)), unique_sublineages_in_data)
  }
  
  # Prepare data for plotting
  top_n_lollipop_data <- top_n_data %>%
    mutate(
      # Create a display name from gene and sublineage
      short_name = ifelse(sublineage == "None" | sublineage == "" | is.na(sublineage), 
                          as.character(gene), 
                          paste0(as.character(gene), "@", as.character(sublineage))),
      short_name = make.unique(short_name, sep = "_") # Ensure unique names for plot
    ) %>%
    # Order features by meta_score for the plot
    mutate(short_name = factor(short_name, levels = short_name[order(meta_score)])) 

  p <- ggplot(top_n_lollipop_data, aes(x = short_name, y = meta_score)) +
    geom_segment(aes(xend = short_name, yend = 0, color = sublineage), linewidth = 1.2) + 
    geom_point(aes(size = n_models_occur, color = sublineage), alpha = 0.8) + # Uses n_models_occur for size
    scale_color_manual(values = sublineage_colors, name = "Sublineage", na.value = "grey50") + 
    scale_size_continuous(range = c(2, 6), name = "# Models Occur") + # Label reflects n_models_occur
    labs(
      title = paste0("Top ", nrow(top_n_lollipop_data), " Features by Meta-Score"),
      x = "Feature (Gene@Sublineage)",
      y = "Meta-Score"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )
  
  # Dynamic height based on number of features
  plot_height <- max(6, nrow(top_n_lollipop_data) * 0.25)
  ggsave(file.path(figures_dir, "top_features_lollipop.png"), p, width = 10, height = plot_height, limitsize = FALSE) 
  message("Created top features lollipop plot in: ", figures_dir)
}


create_batched_lollipop_plots <- function(all_results_data, 
                                         sublineage_colors,
                                         figures_dir, 
                                         batch_size = 50,
                                         max_features = 200) { 
  # all_results_data comes from meta_scores.csv
  
  fs::dir_create(figures_dir, recurse = TRUE, showWarnings = FALSE)
  
  if (nrow(all_results_data) == 0) {
    message("No features to plot for batched lollipops.")
    return()
  }
   if (is.null(sublineage_colors) || length(sublineage_colors) == 0) {
    message("Sublineage colors not available for batched lollipop plots. Using default gray.")
    unique_sublineages_in_data <- unique(all_results_data$sublineage)
    sublineage_colors <- setNames(rep("gray", length(unique_sublineages_in_data)), unique_sublineages_in_data)
  }
  
  # Select top features based on meta_score
  top_features_data <- all_results_data %>% 
    arrange(desc(meta_score)) %>%
    head(max_features)
  
  if (nrow(top_features_data) == 0) {
    message("No features selected after head(max_features) for batched lollipops.")
    return()
  }
  
  # Consistent y-axis limit for all batches
  y_axis_max <- if (nrow(top_features_data) > 0) max(top_features_data$meta_score, na.rm = TRUE) * 1.05 else 1 
  
  num_batches <- ceiling(nrow(top_features_data) / batch_size)
  
  for (batch_num in 1:num_batches) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, nrow(top_features_data))
    
    if (start_idx > nrow(top_features_data)) break # Should not happen with ceiling
    
    batch_features <- top_features_data[start_idx:end_idx, ]
    
    # Prepare data for the current batch
    batch_lollipop_data <- batch_features %>%
      mutate(
        short_name = ifelse(sublineage == "None" | sublineage == "" | is.na(sublineage), 
                            as.character(gene), 
                            paste0(as.character(gene), "@", as.character(sublineage))),
        short_name = make.unique(short_name, sep = "_")
      ) %>%
      mutate(short_name = factor(short_name, levels = short_name[order(meta_score)]))
    
    p <- ggplot(batch_lollipop_data, aes(x = short_name, y = meta_score)) +
      geom_segment(aes(xend = short_name, yend = 0, color = sublineage), linewidth = 1.2) + 
      geom_point(aes(size = n_models_occur, color = sublineage), alpha = 0.8) + # Uses n_models_occur
      scale_color_manual(values = sublineage_colors, name = "Sublineage", na.value = "grey50") +
      scale_size_continuous(range = c(2, 6), name = "# Models Occur") + # Label reflects n_models_occur
      scale_y_continuous(limits = c(0, y_axis_max)) + # Consistent y-axis
      labs(
        title = paste0("Features Ranked ", start_idx, "-", end_idx, " by Meta-Score"),
        x = "Feature (Gene@Sublineage)",
        y = "Meta-Score"
      ) +
      coord_flip() +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        legend.position = "right"
      )
    
    filename <- paste0("top_features_lollipop_batch_", start_idx, "_to_", end_idx, ".png")
    plot_height_batch <- max(6, nrow(batch_lollipop_data) * 0.25)
    ggsave(file.path(figures_dir, filename), p, width = 10, height = plot_height_batch, limitsize = FALSE) 
    message("Created batched lollipop plot: ", filename, " in: ", figures_dir)
  }
}