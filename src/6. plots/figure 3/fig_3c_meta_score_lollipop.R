library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(fs)
library(readr)

# Source utility functions for feature parsing
source("src/0. utils/feature_name_utils.R") 


# === PLOTTING FUNCTIONS ===

# Create batched lollipop plots
create_batched_lollipop_plots <- function(all_results_data, sublineage_colors, figures_dir, batch_size = 50, max_features = 200) {
  # Create output directory if it doesn't exist
  dir_create(figures_dir, recurse = TRUE)
  
  # Select top features based on meta_score
  top_features_data <- all_results_data %>%
    arrange(desc(meta_score)) %>%
    head(max_features)

  # Create display names using existing feature_id plus fold change indicator
  top_features_data <- top_features_data %>%
    mutate(
      # Create fold change indicator arrows using the existing 'fold_change' column
      fold_change_indicator = ifelse(is.na(fold_change), "",
        ifelse(fold_change > 0, "↑", "↓")
      ),
      
      # Simple display name
      display_name_base = paste0(feature_id, " ", fold_change_indicator),
      
      # Ensure unique display names to avoid ggplot issues
      display_name = make.unique(as.character(display_name_base), sep = "..v"),
      
      # Parse sublineage for coloring
      display_sublineage_for_color = get_simplified_sublineage(feature_id, default_value = NA_character_)
    )

  # Ensure we have color mapping for sublineages
  if (is.null(sublineage_colors) || length(sublineage_colors) == 0) {
    unique_subs <- unique(top_features_data$display_sublineage_for_color[!is.na(top_features_data$display_sublineage_for_color)])
    sublineage_colors <- if (length(unique_subs)) setNames(rep("gray50", length(unique_subs)), unique_subs) else list()
  }

  # Calculate number of batches needed
  num_batches <- ceiling(nrow(top_features_data) / batch_size)
  
  # Generate plots for each batch
  for (batch_num in 1:num_batches) {
    # Calculate batch boundaries
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx   <- min(batch_num * batch_size, nrow(top_features_data))

    # Extract and prepare batch data
    batch_features <- top_features_data[start_idx:end_idx, ] %>%
      mutate(display_name = factor(display_name, levels = display_name[order(meta_score)]))

    # Create lollipop plot for this batch
    p <- ggplot(batch_features, aes(x = display_name, y = meta_score)) +
      # Add lollipop stems (lines from 0 to points)
      geom_segment(aes(xend = display_name, yend = 0, color = display_sublineage_for_color), 
                   linewidth = 0.8, alpha = 0.7) +
      # Add lollipop heads (points sized by model occurrence)
      geom_point(aes(size = n_models_occur, color = display_sublineage_for_color), 
                 alpha = 0.9, stroke = 0.3) +
      # Configure color and size scales with better legends
      scale_color_manual(values = sublineage_colors, 
                        name = "Sublineage", 
                        na.value = "grey60",
                        guide = guide_legend(override.aes = list(size = 3), order = 2)) +
      scale_size_continuous(range = c(1.5, 4), 
                           name = "Models\n(count)",
                           breaks = function(x) pretty(x, n = 4) %>% as.integer() %>% unique(),
                           labels = function(x) as.integer(x),
                           guide = guide_legend(override.aes = list(alpha = 1), order = 1)) +
      # Set y-axis limits with better formatting
      scale_y_continuous(limits = c(0, max(top_features_data$meta_score, na.rm = TRUE) * 1.05),
                        expand = expansion(mult = c(0, 0.02)),
                        breaks = scales::pretty_breaks(n = 6)) +
      # Professional axis labels
      labs(x = "Feature", y = "Meta-Score") +
      # Flip coordinates for horizontal lollipops
      coord_flip() +
      # Apply publication-ready theme
      theme_classic() +
      theme(
        # Text styling
        text = element_text(family = "Arial", color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 7, hjust = 1),
        
        # Axis styling
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.3),
        axis.ticks.length = unit(0.15, "cm"),
        
        # Legend styling
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = "black", linewidth = 0.3),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(5, 5, 5, 5),
        legend.box.spacing = unit(0.3, "cm"),
        
        # Panel styling
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3, linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        
        # Plot margins
        plot.margin = margin(10, 15, 10, 10, "pt"),
        
        # Remove panel border
        panel.border = element_blank()
      )

    # Save plot with higher resolution and better format
    filename <- paste0("top_features_lollipop_batch_", 
                      sprintf("%03d", start_idx), "_to_", sprintf("%03d", end_idx), ".png")
    plot_height <- max(4, nrow(batch_features) * 0.2)  # More compact height
    plot_width <- 8  # Standard width for academic papers
    
    ggsave(
      file.path(figures_dir, filename), 
      plot = p, 
      width = plot_width,
      height = plot_height, 
      dpi = 300, 
      units = "in",
      bg = "transparent", 
      limitsize = FALSE
    )
  }
}


# === MAIN LOOPING FUNCTION ===

# Main unified function to run all intersector lollipop plots
run_all_intersector_lollipop_plots <- function(
    intersector_parent_input_dir = "output/3. intersector",
    figures_parent_output_dir = "output/6. plots/figure 3", 
    lollipop_batch_size = 50,
    lollipop_max_total_features = 200,
    dataset_types_to_process = c("ctnorm_global", "ctnorm_relative", "read_depth"), 
    cell_type_groups_to_process = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
    gene_sets_to_process = c("metabolic", "nonmetabolic", "random")
) {
  
  # Load sublineage color map
  sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
  loaded_sublineage_colors <- NULL
  loaded_sublineage_colors <- readRDS(sublineage_color_map_path)

  # Process all combinations
  for (type in dataset_types_to_process) {
    for (cell_type_group in cell_type_groups_to_process) {
      for (gene_set in gene_sets_to_process) {
        message(paste0("\nProcessing dataset type: ", type, 
                       ", Cell Type Group: ", cell_type_group,
                       ", Gene Set: ", gene_set))

        current_input_dir <- file.path(intersector_parent_input_dir, type, cell_type_group, gene_set)
        current_figures_output_dir <- file.path(figures_parent_output_dir, type, cell_type_group, gene_set)
        dir_create(current_figures_output_dir, recurse = TRUE)

        meta_scores_file <- file.path(current_input_dir, "meta_scores.csv")
        fold_changes_file <- file.path("output", "4. fold changes", type, cell_type_group, gene_set, "fold_changes.csv")

        if (!file.exists(meta_scores_file)) {
          message(paste("Meta scores file not found, skipping:", meta_scores_file))
          next
        }

        all_results <- read_csv(meta_scores_file, show_col_types = FALSE)

        # If fold_change column is missing, try to load and join from fold_changes.csv
        if (!"fold_change" %in% names(all_results)) {
          if (file.exists(fold_changes_file)) {
            fold_changes_df <- read_csv(fold_changes_file, show_col_types = FALSE)
            all_results <- all_results %>%
              left_join(fold_changes_df, by = c("feature_id" = "Feature")) %>%
              rename(fold_change = Value)
            message(paste("Added fold_change column from", fold_changes_file))
          } else {
            all_results <- all_results %>% mutate(fold_change = NA_real_)
            message(paste("Warning: 'fold_change' column not found in", meta_scores_file, "and fold_changes.csv not found - adding as NA."))
          }
        }

        # Create lollipop plots
        if (nrow(all_results) > 0) { 
          lollipop_batches_dir_path <- file.path(current_figures_output_dir, "lollipop_batches")
          create_batched_lollipop_plots(
            all_results_data = all_results,
            sublineage_colors = loaded_sublineage_colors, 
            figures_dir = lollipop_batches_dir_path, 
            batch_size = lollipop_batch_size,
            max_features = lollipop_max_total_features
          )
        }

        message("Lollipop plots for ", type, " normalization", cell_type_group, " cell types", gene_set, "genes.")
      }
    }
  }
  
  message("\nFinished generating all intersector lollipop plots.")
}

# === EXECUTE LOLLIPOP PLOTTING ===

run_all_intersector_lollipop_plots()
  