library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(fs)
library(readr)

# Source utility functions for feature parsing
# This is sourced in the main plot_intersector.R, so it might be redundant here
# but doesn't harm if sourced again unless there are side effects in the util script.
source("src/0. utils/feature_name_utils.R") 

# Create batched lollipop plots
# MODIFIED: Accepts processed_fold_changes instead of fold_changes_path
create_batched_lollipop_plots <- function(all_results_data, sublineage_colors, figures_dir, batch_size = 50, max_features = 200, processed_fold_changes) {
  # Removed DEBUG messages for cleanliness
  if (nrow(all_results_data) == 0) {
    message("create_batched_lollipop_plots: No all_results_data provided.")
    return()
  }
  dir_create(figures_dir, recurse = TRUE) 

  if (!"gene" %in% names(all_results_data)) {
    message("ERROR: 'gene' column missing in all_results_data for create_batched_lollipop_plots.")
    return()
  }
  if (!"sublineage" %in% names(all_results_data)) {
    all_results_data$sublineage <- NA_character_
  }
  
  all_results_data <- all_results_data %>%
    mutate(
      gene = str_trim(tolower(gene)),
      sublineage_from_data = str_trim(tolower(get_simplified_sublineage(sublineage, default_value = "")))
    )
  all_results_data$sublineage_from_data[is.na(all_results_data$sublineage_from_data)] <- ""

  top_features_data <- all_results_data %>%
    arrange(desc(meta_score)) %>%
    head(max_features)

  if (nrow(top_features_data) == 0) {
    message("create_batched_lollipop_plots: No features selected after top_n filter.")
    return()
  }
  
  # Use pre-processed fold_changes
  if (is.null(processed_fold_changes) || nrow(processed_fold_changes) == 0) {
    message("create_batched_lollipop_plots: No fold change data provided. Indicators will be missing.")
    top_features_data$fold_change <- NA_real_
  } else {
    top_features_data <- top_features_data %>%
      left_join(processed_fold_changes, by = c("gene" = "gene", "sublineage_from_data" = "sublineage_from_fc"))
  }

  top_features_data <- top_features_data %>%
    mutate(
      fold_change_indicator = ifelse(is.na(fold_change), "", ifelse(fold_change > 0, "↑", ifelse(fold_change < 0, "↓", ""))),
      display_sublineage = ifelse(is.na(sublineage) | str_trim(tolower(sublineage)) %in% c("none", "", "na"), NA_character_, sublineage),
      display_name_base = ifelse(
        is.na(display_sublineage),
        paste0(gene, " ", fold_change_indicator),
        paste0(gene, "@", display_sublineage, " ", fold_change_indicator)
      )
    ) %>%
    mutate(display_name = make.unique(as.character(display_name_base), sep = "_"))

  if (is.null(sublineage_colors) || length(sublineage_colors) == 0) {
    unique_display_sublineages <- unique(top_features_data$display_sublineage[!is.na(top_features_data$display_sublineage)])
    sublineage_colors <- if (length(unique_display_sublineages) > 0) setNames(rep("gray50", length(unique_display_sublineages)), unique_display_sublineages) else list()
  }
  
  num_batches <- ceiling(nrow(top_features_data) / batch_size)

  for (batch_num in 1:num_batches) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, nrow(top_features_data))

    batch_features <- top_features_data[start_idx:end_idx, ] %>%
      mutate(display_name = factor(display_name, levels = display_name[order(meta_score)]))

    p <- ggplot(batch_features, aes(x = display_name, y = meta_score)) +
      geom_segment(aes(xend = display_name, yend = 0, color = display_sublineage), linewidth = 1.2) +
      geom_point(aes(size = n_models_occur, color = display_sublineage), alpha = 0.8) +
      scale_color_manual(values = sublineage_colors, name = "Sublineage", na.value = "grey70") +
      scale_size_continuous(range = c(2, 6), name = "# Models Occur") +
      scale_y_continuous(limits = c(0, max(top_features_data$meta_score, na.rm = TRUE) * 1.05)) + 
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
    save_target_path <- file.path(figures_dir, filename)

    plot_height_batch <- max(6, nrow(batch_features) * 0.25)
    ggsave(save_target_path, p, width = 10, height = plot_height_batch, limitsize = FALSE)
    message("Created batched lollipop plot: ", filename, " in: ", figures_dir)
  }
}
