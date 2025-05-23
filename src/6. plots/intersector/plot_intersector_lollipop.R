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
  if (nrow(all_results_data) == 0) {
    message("create_batched_lollipop_plots: No all_results_data provided.")
    return()
  }
  dir_create(figures_dir, recurse = TRUE)

  # --- Normalize input data to ensure gene, sublineage, cluster_ID columns ---
  if ("feature_id" %in% names(all_results_data)) {
    message("INFO: Parsing from feature_id")
    all_results_data <- all_results_data %>%
      mutate(
        parsed_gene = get_gene_from_feature(feature_id),
        parsed_sublineage = get_simplified_sublineage(feature_id, default_value = NA_character_),
        parsed_cluster_ID = get_cluster_id_from_feature(feature_id)
      ) %>%
      mutate(
        gene = str_trim(tolower(parsed_gene)),
        sublineage = parsed_sublineage,
        cluster_ID = as.character(parsed_cluster_ID)
      ) %>%
      select(-starts_with("parsed_"))
  } else {
    message("INFO: No feature_id; using existing columns")
    if (!"gene" %in% names(all_results_data)) {
      message("ERROR: 'gene' missing. Aborting.")
      return()
    }
    all_results_data <- all_results_data %>% mutate(gene = str_trim(tolower(gene)))
    if (!"sublineage" %in% names(all_results_data)) all_results_data$sublineage <- NA_character_
    if (!"cluster_ID" %in% names(all_results_data)) {
      all_results_data$cluster_ID <- NA_character_
    } else {
      all_results_data$cluster_ID <- as.character(all_results_data$cluster_ID)
    }
  }
  # end normalization

  # prepare for fold‐change join
  all_results_data <- all_results_data %>%
    mutate(sublineage_from_data = str_trim(tolower(
      get_simplified_sublineage(sublineage, default_value = "")
    )))

  top_features_data <- all_results_data %>%
    arrange(desc(meta_score)) %>%
    head(max_features)

  if (nrow(top_features_data) == 0) {
    message("No features after filtering.")
    return()
  }

  if (is.null(processed_fold_changes) || nrow(processed_fold_changes) == 0) {
    message("No fold‐change data; indicators will be missing.")
    top_features_data$fold_change <- NA_real_
  } else {
    top_features_data <- top_features_data %>%
      left_join(processed_fold_changes,
                by = c("gene" = "gene", "sublineage_from_data" = "sublineage_from_fc"))
  }

  top_features_data <- top_features_data %>%
    mutate(
      fold_change_indicator = ifelse(is.na(fold_change), "",
        ifelse(fold_change > 0, "↑", "↓")
      ),

      .clean_sublineage = ifelse(is.na(sublineage) |
                                 str_trim(tolower(sublineage)) %in% c("none","", "na"),
                                 NA_character_, sublineage),
      display_sublineage_for_color = .clean_sublineage,
      temp_sub = .clean_sublineage,
      temp_clus = ifelse(is.na(cluster_ID) |
                         str_trim(tolower(cluster_ID)) %in% c("none","", "na"),
                         "NO_ID", cluster_ID),

      full_context = case_when(
        !is.na(temp_sub) & temp_clus != "NO_ID" ~ paste(temp_sub, temp_clus, sep = "_"),
        !is.na(temp_sub) & temp_clus == "NO_ID" ~ paste0(temp_sub, "_NO_ID"),
        is.na(temp_sub) & temp_clus != "NO_ID" ~ temp_clus,
        TRUE ~ NA_character_
      ),

      display_name_base = if_else(
        is.na(full_context),
        paste0(gene, " ", fold_change_indicator),
        paste0(gene, "@", full_context, " ", fold_change_indicator)
      ),
      display_name = make.unique(as.character(display_name_base), sep = "..v")
    ) %>%
    select(-starts_with("temp"), -.clean_sublineage)

  if (is.null(sublineage_colors) || length(sublineage_colors) == 0) {
    unique_subs <- unique(top_features_data$display_sublineage_for_color[!is.na(top_features_data$display_sublineage_for_color)])
    sublineage_colors <- if (length(unique_subs)) setNames(rep("gray50", length(unique_subs)), unique_subs) else list()
  }

  num_batches <- ceiling(nrow(top_features_data) / batch_size)
  for (batch_num in 1:num_batches) {
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx   <- min(batch_num * batch_size, nrow(top_features_data))

    batch_features <- top_features_data[start_idx:end_idx, ] %>%
      mutate(display_name = factor(display_name, levels = display_name[order(meta_score)]))

    p <- ggplot(batch_features, aes(x = display_name, y = meta_score)) +
      geom_segment(aes(xend = display_name, yend = 0, color = display_sublineage_for_color), linewidth = 1.2) +
      geom_point(aes(size = n_models_occur, color = display_sublineage_for_color), alpha = 0.8) +
      scale_color_manual(values = sublineage_colors, name = "Sublineage", na.value = "grey70") +
      scale_size_continuous(range = c(2, 6), name = "# Models Occur") +
      scale_y_continuous(limits = c(0, max(top_features_data$meta_score, na.rm = TRUE) * 1.05)) +
      labs(title = paste0("Features ", start_idx, "-", end_idx, " by Meta-Score"),
           x = "Feature (Gene@Sublineage_ClusterID)", y = "Meta-Score") +
      coord_flip() +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 9),
            panel.grid.major.y = element_blank(),
            legend.position = "right")

    filename <- paste0("top_features_lollipop_batch_", start_idx, "_to_", end_idx, ".png")
    ggsave(file.path(figures_dir, filename), p, width = 10,
           height = max(6, nrow(batch_features) * 0.25), limitsize = FALSE)
    message("Saved lollipop plot: ", filename)
  }
}
