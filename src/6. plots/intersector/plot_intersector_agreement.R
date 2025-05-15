library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(viridis)


plot_feature_agreement_dotplot <- function(top_n_overall_data, 
                                           scaled_importance_cols_to_use, 
                                           figures_dir, 
                                           top_n_agreement = 20) {

  if (nrow(top_n_overall_data) == 0) {
    message("No data provided for feature agreement dotplot.")
    return()
  }
  
  # Ensure the provided scaled_importance_cols_to_use exist in top_n_overall_data
  valid_scaled_cols <- intersect(scaled_importance_cols_to_use, names(top_n_overall_data))
  if(length(valid_scaled_cols) == 0){
    message("No valid scaled importance columns (e.g., modelname_score) found in the input data for agreement plot.")
    return()
  }

  # Select a subset of features for this specific plot, ordered by meta_score
  # The input top_n_overall_data is already a subset (e.g., top 50)
  agreement_subset_data <- top_n_overall_data %>%
    arrange(desc(meta_score)) %>%
    head(top_n_agreement)

  if (nrow(agreement_subset_data) == 0) {
    message("No features selected for agreement plot after taking top_n_agreement.")
    return()
  }

  # Pivot the data to long format for plotting
  # Include gene and sublineage for creating display names
  agreement_data_long <- agreement_subset_data %>%
    dplyr::select(feature_id, gene, sublineage, meta_score, dplyr::all_of(valid_scaled_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(valid_scaled_cols),
      names_to = "model_col_name", 
      values_to = "scaled_importance"
    ) %>%
    dplyr::mutate(model_display_name = gsub("_score$", "", model_col_name))

  # Create a display name for features using gene and sublineage
  # Order features by their original meta_score for consistent plot appearance
  
  # First, prepare the ordered levels for the factor based on the subset of data
  ordered_feature_levels <- agreement_subset_data %>%
    mutate(
      short_name_display_temp = ifelse(
        sublineage == "None" | sublineage == "" | is.na(sublineage), 
        as.character(gene), 
        paste0(as.character(gene), "\n@", as.character(sublineage)) # Use newline for better y-axis display
      )
    ) %>%
    arrange(desc(meta_score)) %>% # Ensure this order is based on the subset for the plot
    pull(short_name_display_temp) %>%
    unique() # Get unique display names in the correct order

  agreement_plot_data <- agreement_data_long %>%
    mutate(
      short_name_display = ifelse(
        sublineage == "None" | sublineage == "" | is.na(sublineage), 
        as.character(gene), 
        paste0(as.character(gene), "\n@", as.character(sublineage))
      ),
      short_name_display = factor(short_name_display, levels = rev(ordered_feature_levels)) # Reverse for ggplot top-to-bottom
    )

  # Filter for positive scaled_importance to avoid plotting zero-score points if not desired
  # Scaled scores are typically > 0 if the feature was selected by a model.
  plot_data_to_use <- agreement_plot_data %>% filter(scaled_importance > 1e-9) # Use a small epsilon

  if (nrow(plot_data_to_use) == 0) {
    message("No positive scaled importance values to plot for feature agreement.")
    return()
  }
  
  p <- ggplot(plot_data_to_use, 
              aes(x = model_display_name, y = short_name_display, size = scaled_importance, color = scaled_importance)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(range = c(2, 8), name = "Scaled Score") + # Legend label reflects scaled score
    scale_color_viridis_c(option = "plasma", name = "Scaled Score") + # Legend label reflects scaled score
    labs(
      title = paste0("Feature Scaled Score Agreement (Top ", nrow(agreement_subset_data), " Features)"), # Title reflects scaled score
      x = "Model",
      y = "Feature (Gene@Sublineage)" # Y-axis label
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=9),
      axis.text.y = element_text(size = 8), # Adjust size as needed
      panel.grid.major = element_line(color = "gray90"),
      legend.position = "right"
    )
  
  # Dynamic height based on the number of features in this specific plot
  plot_height <- max(6, 0.4 * length(unique(agreement_plot_data$short_name_display))) 
  ggsave(file.path(figures_dir, "feature_agreement_scaled_scores.png"), p, width = 12, height = plot_height, limitsize = FALSE)
  message("Created feature agreement plot (using scaled scores) in: ", figures_dir)
}