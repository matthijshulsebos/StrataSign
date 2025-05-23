library(dplyr)
library(readr)
library(fs)
library(stringr)

# Source plotting scripts 
source("src/6. plots/intersector/plot_intersector_lollipop.R")
source("src/6. plots/intersector/plot_intersector_model_prevalence.R")
source("src/6. plots/intersector/plot_intersector_agreement.R")
# Source utility functions as they are used for fold change processing here
source("src/0. utils/feature_name_utils.R") 

run_all_intersector_plots <- function(
    intersector_parent_input_dir,
    figures_parent_output_dir,
    fold_changes_file = "output/4. differential expression/feature_fold_changes.csv",
    top_n_visualizations = 50, # Used for prevalence and agreement plots
    agreement_plot_n_features = 20,
    lollipop_batch_size = 50,
    lollipop_max_total_features = 200, # For batched lollipops
    dataset_types_to_process = c("absolute", "relative"),
    cell_type_groups_to_process = c("all_clusters", "macrophages")
) {
  
  sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
  loaded_sublineage_colors <- NULL
  if (file.exists(sublineage_color_map_path)) {
    loaded_sublineage_colors <- readRDS(sublineage_color_map_path)
    message(paste("Loaded sublineage color map from:", sublineage_color_map_path,
                  "with", length(loaded_sublineage_colors), "entries."))
    if (length(loaded_sublineage_colors) == 0) message("WARNING: Loaded sublineage color map is empty.")
  } else {
    message(paste("WARNING: Sublineage color map not found at:", sublineage_color_map_path, 
                  "Lollipop plots may use default or fallback coloring."))
  }

  # --- Load and process fold changes ONCE ---
  processed_fold_changes <- NULL
  if (file.exists(fold_changes_file)) {
    message(paste("Loading and processing fold changes from:", fold_changes_file, "..."))
    processed_fold_changes <- read_csv(fold_changes_file, show_col_types = FALSE) %>%
      rename(feature_identifier_raw = Feature, fold_change_raw = Value) %>%
      mutate(
        gene_parsed_fc = get_gene_from_feature(feature_identifier_raw),
        sublineage_parsed_fc = get_simplified_sublineage(feature_identifier_raw, default_value = ""),
        gene = str_trim(tolower(gene_parsed_fc)),
        sublineage_from_fc = str_trim(tolower(sublineage_parsed_fc)),
        fold_change = as.numeric(as.character(fold_change_raw))
      ) %>%
      filter(!is.na(gene) & gene != "" & !is.na(fold_change)) %>%
      group_by(gene, sublineage_from_fc) %>%
      summarise(fold_change = first(fold_change), .groups = 'drop')
    message(paste("Finished processing fold changes. Loaded", nrow(processed_fold_changes), "unique gene-sublineage fold changes."))
    if (nrow(processed_fold_changes) == 0) {
        message("WARNING: No fold change data loaded or processed. Lollipop plots may not show fold change indicators.")
        processed_fold_changes <- NULL 
    }
  } else {
    message(paste("WARNING: Fold changes file not found at:", fold_changes_file,
                  "Lollipop plots will not show fold change indicators."))
  }
  # --- End of fold change processing ---

  for (type in dataset_types_to_process) {
    for (cell_type_group in cell_type_groups_to_process) {
      message(paste0("\nProcessing dataset type: ", type, ", Cell Type Group: ", cell_type_group))
      message(paste(rep("=", 50), collapse = ""))

      current_input_dir <- file.path(intersector_parent_input_dir, type, cell_type_group)
      current_figures_output_dir <- file.path(figures_parent_output_dir, type, cell_type_group)
      dir_create(current_figures_output_dir, recurse = TRUE)

      meta_scores_file <- file.path(current_input_dir, "meta_scores.csv")
      if (!file.exists(meta_scores_file)) {
        warning(paste("Meta scores file not found:", meta_scores_file, "- Skipping this group."))
        next
      }
      
      all_results <- read_csv(meta_scores_file, show_col_types = FALSE)
      if (nrow(all_results) == 0) {
        message(paste("No data in meta_scores.csv for", type, cell_type_group, "- Skipping plots for this group."))
        next
      }
      message(paste("Loaded meta scores data for '", type, "' and cell group '", cell_type_group, 
                    "' with", nrow(all_results), "features from:", meta_scores_file))

      model_score_cols <- names(all_results)[str_starts(names(all_results), "scaled_importance_")]
      
      top_n_data_for_plots <- all_results %>% 
        arrange(desc(meta_score)) %>%
        head(top_n_visualizations)

      if (nrow(top_n_data_for_plots) > 0 && "n_models_occur" %in% names(top_n_data_for_plots)) {
        plot_model_prevalence_barplot( 
          top_n_data = top_n_data_for_plots, 
          figures_dir = current_figures_output_dir
        )
      }
      
      if (nrow(top_n_data_for_plots) > 0 && length(model_score_cols) > 0) {
        plot_feature_agreement_dotplot(
          top_n_overall_data = top_n_data_for_plots, 
          scaled_importance_cols_to_use = model_score_cols,
          figures_dir = current_figures_output_dir,
          top_n_agreement = agreement_plot_n_features 
        )
      }
      
      if (nrow(all_results) > 0) { 
        lollipop_batches_dir_path <- file.path(current_figures_output_dir, "lollipop_batches")
        create_batched_lollipop_plots(
          all_results_data = all_results,
          sublineage_colors = loaded_sublineage_colors, 
          figures_dir = lollipop_batches_dir_path, 
          batch_size = lollipop_batch_size,
          max_features = lollipop_max_total_features,
          processed_fold_changes = processed_fold_changes 
        )
      }
      
      message("All intersector visualizations for '", type, "', cell group '", cell_type_group, "' orchestrated. Check: ", current_figures_output_dir)
    }
  }
  
  message("\nFinished processing all dataset types and cell type groups.")
}

# Example call (ensure your paths and parameters are correct for your setup)
run_all_intersector_plots(
  intersector_parent_input_dir = "output/3. intersector", 
  figures_parent_output_dir = "output/6. plots/intersector", 
  fold_changes_file = "output/4. differential expression/feature_fold_changes.csv",
  top_n_visualizations = 50,
  agreement_plot_n_features = 20,
  lollipop_batch_size = 50,
  lollipop_max_total_features = 200,
  dataset_types_to_process = c("absolute", "relative"),
  cell_type_groups_to_process = c("all_clusters", "macrophages")
)
