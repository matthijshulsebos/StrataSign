library(dplyr)
library(tidyr)
library(readr)
library(fs)
library(ggplot2) 
library(RColorBrewer)

# Source the utility functions for feature name parsing
source("src/0. utils/feature_name_utils.R")

# Function to create sublineage color map
create_and_save_sublineage_color_map <- function(data, output_path) {
  # Check for valid input data
  if (is.null(data) || !"feature_id" %in% names(data) || nrow(data) == 0) {
    message("Cannot generate sublineage color map: input data is NULL, empty, or missing 'feature_id'.")
    return(NULL)
  }

  # Extract sublineages from feature identifier using the utils
  all_sublineages <- data %>%
    mutate(
      sublineage = get_simplified_sublineage(as.character(feature_id))
    ) %>%
    pull(sublineage) %>%
    unique() %>%
    sort()

  # Check if sublineages exist for color mapping
  if (length(all_sublineages) == 0 || (length(all_sublineages) == 1 && all_sublineages[1] == "None")) {
    message("No distinct sublineages found, excluding 'None', to create a color map.")
    return(NULL)
  }
  
  num_colors_needed <- length(all_sublineages)
  
  # Color vector generation
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  # Create a single vector of all available unique qualitative colors
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Select the required number of colors
  available_colors <- rep(col_vector, length.out = num_colors_needed)
  
  # Create a named vector mapping sublineages to colors
  color_map <- setNames(available_colors, all_sublineages)
  
  # Save the color map to an RDS file
  tryCatch({
    saveRDS(color_map, file = output_path)
    message("Sublineage color map saved to: ", output_path)
  }, error = function(e) {
    message("Error saving sublineage color map: ", e$message)
  })
  
  return(color_map)
}

# Load feature importance files and extract relevant features
load_and_extract_features <- function(file_path, model_name) {
  # Read csv
  df <- tryCatch(
    readr::read_csv(file_path, show_col_types = FALSE), 
    error = function(e) {
      message(paste("Error reading file:", file_path, "-", e$message)) # Log error
      return(NULL) # Return NULL if reading fails
    }
  )

  # Validate loaded data frame
  if (is.null(df) || nrow(df) == 0 || !"Feature" %in% names(df) || !"Value" %in% names(df)) {
    message(paste("Skipping file due to missing columns, empty content, or read error:", file_path))
    return(NULL)
  }

  # Extract gene, sublineage, absolute importance value
  df_out <- df %>%
    mutate(
      feature_id = as.character(Feature),
      gene = get_gene_from_feature(feature_id),
      sublineage = get_simplified_sublineage(feature_id),
      abs_value = abs(Value),
      model_type = model_name
    ) %>%
    select(feature_id, gene, sublineage, abs_value, model_type) %>%
    filter(!is.na(gene))
  
  return(df_out)
}

# Main function to find feature intersections across models
find_feature_intersection <- function(models_dir = "output/models/absolute",
                                      output_dir = "output/intersector/absolute",
                                      min_models_occurrence = 2) {
  
  # Create output directory
  dir_create(output_dir, recurse = TRUE)

  # Filters are hardcoded for gene and cell type sets
  gene_type_filter_val <- "metabolic" 
  cell_type_filter_val <- "all_clusters" 

  # List all feature importance files recursively within the models_dir
  all_fi_files_raw <- dir_ls(models_dir, recurse = TRUE, regexp = "feature_importance.*\\.csv$", fail = FALSE)
  
  # Construct the base path to filter files based on gene and cell type
  filter_path_segment <- fs::path("", gene_type_filter_val, cell_type_filter_val, "")
  # Filter files based on the constructed path segment using grep
  filtered_fi_files <- grep(filter_path_segment, all_fi_files_raw, value = TRUE, fixed = FALSE)


  # Helper function for initting meta score output df
  create_empty_final_table <- function(processed_model_names = character(0)) {
      cols_def <- list(
          feature_id = character(), 
          gene = character(), 
          sublineage = character(),
          meta_score = numeric(), 
          n_models_occur = integer()
      )
      # Add model score columns
      if (length(processed_model_names) > 0) {
          for (m_name in processed_model_names) {
              cols_def[[paste0(m_name, "_score")]] <- numeric()
          }
      }
      df <- tibble::as_tibble(cols_def)
      return(df)
  }

  # Handle case where no files match the filter criteria
  if (length(filtered_fi_files) == 0) {
    message(paste("No files found for gene_type '", gene_type_filter_val, "' and cell_type '", cell_type_filter_val, "' in ", models_dir))
    empty_df <- create_empty_final_table()
    readr::write_csv(empty_df, file.path(output_dir, "meta_scores.csv"))
    # Ensure path exists before saving empty colors file
    sublineage_color_map_path <- file.path(output_dir, "sublineage_colors.rds")
    if (!dir.exists(dirname(sublineage_color_map_path))) {
        dir.create(dirname(sublineage_color_map_path), recursive = TRUE)
    }
    saveRDS(NULL, file = sublineage_color_map_path)
    message("Empty meta_scores.csv and sublineage_colors.rds (NULL) created as no files were processed.")
    return(empty_df)
  }

  # Create a list of file paths and their corresponding model names
  file_model_list <- lapply(filtered_fi_files, function(f_path) {
    # Extract model name from the file path
    model_name <- path_split(path_rel(f_path, models_dir))[[1]][1]
    list(file_path = f_path, model_name = model_name)
  })
  
  # Sorted list of all model names that were processed
  all_processed_model_names <- sort(unique(sapply(file_model_list, `[[`, "model_name")))

  # Load and extract features from all filtered files
  all_features_raw_list <- lapply(file_model_list, function(fm) {
    load_and_extract_features(fm$file_path, fm$model_name)
  })
  # Combine all loaded feature data frames
  all_features_df <- dplyr::bind_rows(Filter(Negate(is.null), all_features_raw_list))

  # Handle case where no features were loaded
  if (nrow(all_features_df) == 0) {
    message("No features loaded from any files.")
    empty_df <- create_empty_final_table(all_processed_model_names)
    readr::write_csv(empty_df, file.path(output_dir, "meta_scores.csv"))
    return(empty_df)
  }

  # Determine if each model is dense
  model_density_map <- sapply(all_processed_model_names, function(m_name) {
    model_features <- all_features_df %>% filter(model_type == m_name)
    if (nrow(model_features) == 0) return(FALSE)
    num_nonzero <- sum(model_features$abs_value > 1e-10, na.rm = TRUE)
    # Definition of dense is having more than 100 non-zero features
    return(num_nonzero > 100)
  }, USE.NAMES = TRUE)

  # List to store scores for each model
  list_of_model_specific_scores_df <- list() 

  # Loop through each model to calculate and scale its feature scores
  for (current_model_name in all_processed_model_names) {
    # Filter features for the current model
    current_model_data <- all_features_df %>% 
      filter(model_type == current_model_name) %>%
      select(feature_id, gene, sublineage, abs_value)
    
    # Define score column name for this model
    model_score_col_name <- paste0(current_model_name, "_score") 
    # Init score with 0 for all features
    current_model_data <- current_model_data %>%
      mutate(!!rlang::sym(model_score_col_name) := 0.0)

    is_dense <- model_density_map[[current_model_name]]
    if (is.na(is_dense)) is_dense <- FALSE

    # Filter candidate features (those with importance > epsilon) for scaling because of floating point nonsense
    candidate_features <- current_model_data %>% filter(abs_value > 1e-10)
    features_to_scale_scores_for <- tibble()

    if (nrow(candidate_features) > 0) {
      # Select features for scaling based on model density
      if (is_dense) {
        # For dense models select top features contributing to 70% of cumulative importance
        candidate_features_ordered <- candidate_features %>% arrange(desc(abs_value))
        total_abs_importance <- sum(candidate_features_ordered$abs_value, na.rm = TRUE)
        if (total_abs_importance > 0) {
          # Calculate cumulative importance proportion
          candidate_features_ordered <- candidate_features_ordered %>%
            mutate(cumulative_importance_proportion = cumsum(abs_value) / total_abs_importance)
          
          # Select features up to 70% cumulative importance
          selected_indices <- which(candidate_features_ordered$cumulative_importance_proportion <= 0.7)
          # Select at least one feature
          if(length(selected_indices) == 0 && nrow(candidate_features_ordered) > 0) selected_indices <- 1
          features_to_scale_scores_for <- candidate_features_ordered %>% slice(selected_indices)
        }
      } else {
        # For non-dense (sparse) models, scale all candidate features
        features_to_scale_scores_for <- candidate_features
      }

      # Scale scores for selected features to a range of 0.5 to 1.0
      if (nrow(features_to_scale_scores_for) > 0) {
        min_val <- min(features_to_scale_scores_for$abs_value, na.rm = TRUE)
        max_val <- max(features_to_scale_scores_for$abs_value, na.rm = TRUE)
        
        # Perform min-max scaling to 0.5 - 1.0
        scaled_values <- if (max_val == min_val) {
          # If all values are the same which should never happen assign 0.75 to all
          rep(0.75, nrow(features_to_scale_scores_for))
        } else {
          # Standard min-max scaling adjusted to the 0.5-1.0 range
          0.5 + 0.5 * (features_to_scale_scores_for$abs_value - min_val) / (max_val - min_val)
        }
        
        # Create a temporary data frame with feature identifier and their new scaled scores
        scaled_scores_to_update_df <- features_to_scale_scores_for %>%
            select(feature_id) %>%
            mutate(!!rlang::sym(model_score_col_name) := scaled_values)

        # Join scaled scores to the current model data
        current_model_data <- current_model_data %>%
            left_join(scaled_scores_to_update_df, by = "feature_id", suffix = c(".old", ".new"))
        
        new_score_col_suffix_name <- paste0(model_score_col_name, ".new")
        old_score_col_suffix_name <- paste0(model_score_col_name, ".old")

        # If a scaled score exists use it otherwise keep the old 0.0 score 
        current_model_data <- current_model_data %>%
            mutate(
                !!rlang::sym(model_score_col_name) := ifelse(!is.na(!!rlang::sym(new_score_col_suffix_name)), 
                                                 !!rlang::sym(new_score_col_suffix_name), 
                                                 !!rlang::sym(old_score_col_suffix_name))
            ) %>%
            select(-any_of(c(new_score_col_suffix_name, old_score_col_suffix_name)))
      }
    }
    # Store the processed model data frame in the list
    list_of_model_specific_scores_df[[current_model_name]] <- current_model_data %>% 
      select(feature_id, gene, sublineage, !!rlang::sym(model_score_col_name))
  }

  # Edge case if models dont have important features
  if (length(list_of_model_specific_scores_df) == 0) {
    message("No model-specific scores were generated.")
    empty_df <- create_empty_final_table(all_processed_model_names)
    readr::write_csv(empty_df, file.path(output_dir, "meta_scores.csv"))
    return(empty_df)
  }

  # Combine all model score data frames
  final_merged_df <- purrr::reduce(
    list_of_model_specific_scores_df, 
    dplyr::full_join, 
    by = c("feature_id", "gene", "sublineage")
  )

  # Store all model score columns
  model_score_columns <- grep("_score$", names(final_merged_df), value = TRUE)
  
  # Replace NAs in score columns with 0.0
  final_merged_df <- final_merged_df %>%
    mutate(across(all_of(model_score_columns), ~ifelse(is.na(.), 0.0, .)))
  
  # Calculate meta_score sum of model scores
  final_merged_df <- final_merged_df %>%
    mutate(
      meta_score = rowSums(select(., all_of(model_score_columns)), na.rm = TRUE),
      n_models_occur = rowSums(select(., all_of(model_score_columns)) > 1e-10, na.rm = TRUE)
    )

  # Filter features based on the minimum number of model occurrences
  final_results_df <- final_merged_df %>% 
    filter(n_models_occur >= min_models_occurrence)

  # Handle case where no features meet the minimum occurrence threshold
  if (nrow(final_results_df) == 0) {
    message(paste("No features met min_models_occurrence of", min_models_occurrence))
    empty_df <- create_empty_final_table(all_processed_model_names)
    readr::write_csv(empty_df, file.path(output_dir, "meta_scores.csv"))
    if (nrow(all_features_df) > 0) {
        sublineage_color_map_path <- file.path(output_dir, "sublineage_colors.rds")
        create_and_save_sublineage_color_map(all_features_df, sublineage_color_map_path)
    } else {
        sublineage_color_map_path <- file.path(output_dir, "sublineage_colors.rds")
        if (!dir.exists(dirname(sublineage_color_map_path))) {
            dir.create(dirname(sublineage_color_map_path), recursive = TRUE)
        }
        saveRDS(NULL, file = sublineage_color_map_path)
    }
    return(empty_df)
  }

  # Sort final results by meta_score in descending order
  final_results_df <- final_results_df %>% arrange(desc(meta_score))
  
  # Define column groups for ordering the output csv for better readability
  id_cols <- c("feature_id", "gene", "sublineage")
  summary_cols <- c("meta_score", "n_models_occur")

  # Avoid sorting meta_score even though it has the _score suffix
  model_score_columns_for_final_order <- setdiff(model_score_columns, "meta_score")
  sorted_model_score_columns <- sort(model_score_columns_for_final_order)
  
  # Define the column order for the output csv
  final_column_order <- c(id_cols, summary_cols, sorted_model_score_columns)
  # Ensure all columns are in the defined order
  final_column_order <- intersect(final_column_order, names(final_results_df)) 
  
  # Select and reorder columns for the final output
  final_results_df <- final_results_df %>% select(all_of(final_column_order))

  # Write final results to meta_scores.csv
  readr::write_csv(final_results_df, file.path(output_dir, "meta_scores.csv"))
  message("Intersector analysis (meta_scores.csv) saved to ", file.path(output_dir, "meta_scores.csv"))
  
  # Create and save sublineage color map based on the final results
  sublineage_color_map_path <- file.path(output_dir, "sublineage_colors.rds")
  create_and_save_sublineage_color_map(final_results_df, sublineage_color_map_path)
  
  return(final_results_df)
}

# Wrapper function to run intersector for absolute and relative data
run_intersector_typed <- function(
    base_models_dir = "output/2. models",
    base_output_dir = "output/3. intersector",
    min_models_param = 2,
    dataset_types_to_process = c("absolute", "relative") 
) {
  
  all_type_results <- list()
  
  # Loop through each specified dataset type
  for (type in dataset_types_to_process) {
    message(paste0("\nProcessing intersector for dataset type: ", type, "\n", paste(rep("=", 30), collapse = "")))

    # Define input and output directories for the current type
    current_models_dir <- file.path(base_models_dir, type)
    current_output_dir <- file.path(base_output_dir, type)

    # Check if the input directory for the current type exists
    if (!dir.exists(current_models_dir)) {
      warning(paste("Input models directory for type '", type, "' not found: ", current_models_dir, ". Skipping."), immediate. = TRUE)
      next
    }
    
    message(paste("Models_dir:", current_models_dir, " Output_dir:", current_output_dir))

    # Run the main intersector function
    result_df <- find_feature_intersection( 
      models_dir = current_models_dir,
      output_dir = current_output_dir,
      min_models_occurrence = min_models_param
    )
    
    # Store results
    if (!is.null(result_df) && nrow(result_df) > 0) {
      all_type_results[[type]] <- result_df
      message(paste("Finished processing for type ", type))
    } else {
      message(paste("Intersector returned empty for type ", type))
    }
  }
  
  message("\nFinished all dataset types for intersector analysis.")
  return(invisible(all_type_results))
}

run_intersector_typed(
  base_models_dir = "output/2. models", 
  base_output_dir = "output/3. intersector", 
  min_models_param = 2
)
