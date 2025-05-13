library(ggplot2)
library(dplyr)
library(viridis) # For scale_fill_viridis_c

#' Plot model prevalence barplot for top N features
#'
#' @param top_n_data Data frame of top N features, expected to contain
#'        'gene', 'sublineage', 'meta_score', and 'n_models_occur'.
#' @param figures_dir Directory where the plot will be saved
plot_model_prevalence_barplot <- function(top_n_data, figures_dir) {

  if (nrow(top_n_data) == 0) {
    message("No data for model prevalence barplot.")
    return()
  }

  # Prepare data for plotting: create a display name and ensure correct ordering
  plot_data <- top_n_data %>%
    mutate(
      # Create a display name from gene and sublineage for better readability
      display_name = ifelse(sublineage == "None" | sublineage == "" | is.na(sublineage), 
                            as.character(gene), 
                            paste0(as.character(gene), "@", as.character(sublineage))),
      display_name = make.unique(display_name, sep = "_") # Ensure unique names
    )

  p <- ggplot(plot_data, aes(x = reorder(display_name, meta_score), y = n_models_occur, fill = meta_score)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c() +
    labs(
      title = "Number of Models Where Top Features Appear",
      x = "Feature (Gene@Sublineage)", # Updated x-axis label
      y = "Number of Models Feature Occurs In", # Updated y-axis label to reflect n_models_occur
      fill = "Meta-Score"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), # Keep rotation for long names
      axis.title.x = element_text(size=10) # Add back x-axis title for clarity
    )
  
  # Adjust width and height as needed
  ggsave(file.path(figures_dir, "model_prevalence.png"), p, width = 14, height = 8) 
  message("Created model prevalence plot in: ", figures_dir)
}
