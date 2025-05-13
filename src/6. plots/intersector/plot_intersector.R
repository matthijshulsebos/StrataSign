library(dplyr)
library(readr)
library(fs)
library(stringr)

# Source plotting scripts 
source("src/6. plots/intersector/plot_intersector_lollipop.R")
source("src/6. plots/intersector/plot_intersector_model_prevalence.R")
source("src/6. plots/intersector/plot_intersector_agreement.R")
source("src/6. plots/intersector/plot_importance_cumulative.R")

run_all_intersector_plots <- function(
    intersector_parent_input_dir,
    model_parent_input_dir_for_cumulative,
    figures_parent_output_dir,
    top_n_visualizations = 50,
    agreement_plot_n_features = 20,
    lollipop_batch_size = 50,
    lollipop_max_total_features = 100,
    cumulative_version_filter = "metabolic",
    cumulative_importance_threshold = 0.7,
    dataset_types_to_process = c("absolute", "relative")
) {
  
  for (type in dataset_types_to_process) {
    message(paste0("\nProcessing dataset type: ", type, "\n", paste(rep("=", 30), collapse = "")))

    current_intersector_input_dir <- file.path(intersector_parent_input_dir, type)
    current_model_input_dir_for_cumulative <- file.path(model_parent_input_dir_for_cumulative, type)
    current_figures_output_dir <- file.path(figures_parent_output_dir, type)
    
    dir_create(current_figures_output_dir, recursive = TRUE)

    # Changed to read meta_scores.csv
    meta_scores_path <- file.path(current_intersector_input_dir, "meta_scores.csv")
    if (!file.exists(meta_scores_path)) {
      warning("Meta scores file not found for type '", type, "': ", meta_scores_path)
      next 
    }
    all_results <- read_csv(meta_scores_path, show_col_types = FALSE)
    message("Loaded meta scores data for '", type, "' with ", nrow(all_results), " features from: ", meta_scores_path)

    if(nrow(all_results) == 0) {
      message("No data in meta_scores.csv for type '", type, "'. Skipping plot generation for this type.")
      next
    }

    # Load sublineage colors
    sublineage_colors_path <- file.path(current_intersector_input_dir, "sublineage_colors.rds")
    sublineage_colors <- NULL
    if (file.exists(sublineage_colors_path)) {
        sublineage_colors <- readRDS(sublineage_colors_path)
        message("Loaded sublineage color map for '", type, "' from: ", sublineage_colors_path)
    } else {
        message("Sublineage color map RDS file not found at: ", sublineage_colors_path)
        message("Attempting to generate sublineage color map from meta_scores.csv data...")
        # Ensure column as expected by the create_sublineage_color_map
        if ("sublineage" %in% names(all_results)) {
            sublineage_colors <- create_sublineage_color_map(all_results)
            if (!is.null(sublineage_colors)) {
                 saveRDS(sublineage_colors, sublineage_colors_path)
                 message("Generated and saved new sublineage color map to: ", sublineage_colors_path)
            } else {
                 message("Failed to generate sublineage color map.")
            }
        } else {
            message("Cannot generate sublineage color map: 'sublineage' column missing in meta_scores.csv.")
        }
    }
    
    if (is.null(sublineage_colors)) {
        message("Warning: Sublineage colors are not available. Plots requiring them may error or look incorrect.")
    }


    # Determine model score columns
    model_score_cols <- names(all_results)[grepl("_score$", names(all_results)) & names(all_results) != "meta_score"]
    
    if (length(model_score_cols) == 0) {
        message("No model-specific score columns (e.g., 'modelname_score') found in meta_scores.csv for type '", type, "'. Some plots may not be generated correctly.")
    } else {
        message("Identified model score columns for '", type, "': ", paste(model_score_cols, collapse=", "))
    }
    
    # Top N data for plots that need a subset
    top_n_data_for_plots <- all_results %>% 
      arrange(desc(meta_score)) %>%
      head(top_n_visualizations)
    
    if(nrow(top_n_data_for_plots) == 0) {
      message("No features available after selecting top_n_visualizations for '", type, "' (or top_n_visualizations is 0). Some plots might be empty or not generated.")
    } else {
      message("Prepared top ", nrow(top_n_data_for_plots), " features for '", type, "' visualization based on meta_score.")
    }

    # Top N Lollipop Chart
    if (nrow(top_n_data_for_plots) > 0 && !is.null(sublineage_colors)) {
      plot_top_features_lollipop(
        top_n_data = top_n_data_for_plots,
        sublineage_colors = sublineage_colors,
        figures_dir = current_figures_output_dir
      )
    } else {
      message("Skipping Top N Lollipop Chart for '", type, "' due to missing data or colors.")
    }
    
    # Model Prevalence Barplot
    if (nrow(top_n_data_for_plots) > 0) {
      plot_model_prevalence_barplot(
        top_n_data = top_n_data_for_plots,
        figures_dir = current_figures_output_dir
      )
    } else {
      message("Skipping Model Prevalence Barplot for '", type, "' due to missing data.")
    }
    
    # Feature Agreement Dotplot
    if (nrow(top_n_data_for_plots) > 0 && length(model_score_cols) > 0) {
      plot_feature_agreement_dotplot(
        top_n_overall_data = top_n_data_for_plots, 
        scaled_importance_cols_to_use = model_score_cols,
        figures_dir = current_figures_output_dir,
        top_n_agreement = agreement_plot_n_features 
      )
    } else {
      message("Skipping Feature Agreement Dotplot for '", type, "' due to missing data or model score columns.")
    }

    # Batched Lollipop Plots
    if (nrow(all_results) > 0 && !is.null(sublineage_colors)) {
      create_batched_lollipop_plots(
        all_results_data = all_results,
        sublineage_colors = sublineage_colors,
        figures_dir = file.path(current_figures_output_dir, "lollipop_batches"),
        batch_size = lollipop_batch_size,
        max_features = lollipop_max_total_features
      )
    } else {
      message("Skipping Batched Lollipop Plots for '", type, "' due to missing data or colors.")
    }
    
    # Cumulative Importance Plots
    plot_cumulative_importance(
      model_output_dir = current_model_input_dir_for_cumulative, 
      figures_dir = file.path(current_figures_output_dir, "cumulative"),
      version_filter = cumulative_version_filter,
      threshold = cumulative_importance_threshold
    )
    
    message("All intersector visualizations for '", type, "' orchestrated. Check: ", current_figures_output_dir)
  } 
  
  message("\nFinished processing all dataset types.")
}

run_all_intersector_plots(
  intersector_parent_input_dir = "output/3. intersector", 
  model_parent_input_dir_for_cumulative = "output/2. models", 
  figures_parent_output_dir = "output/6. plots/intersector", 
  top_n_visualizations = 50,
  agreement_plot_n_features = 20,
  lollipop_batch_size = 50,
  lollipop_max_total_features = 200,
  cumulative_version_filter = "metabolic", 
  cumulative_importance_threshold = 0.7
)
