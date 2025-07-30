library(dplyr)
library(tidyr)
library(readr)
library(fs)
library(RColorBrewer)
library(glue)
library(purrr)
library(data.table)

# Source utility functions
source("src/0. utils/feature_name_utils.R")


# Load an important features file
load_and_extract_features <- function(file_path, model_name) {
  df <- as_tibble(fread(file_path))

  # Data format is not adhered to
  if (is.null(df) || !"Feature"%in%names(df)||!"Value"%in%names(df)||nrow(df)==0) {
    return(NULL)
  }
  
  # Decompose the feature into components using the utils
  df %>%
    transmute(
      feature_id = as.character(Feature),
      gene        = get_gene_from_feature(feature_id),
      sublineage  = get_simplified_sublineage(feature_id),
      abs_value   = abs(Value),
      model_type  = model_name
    ) %>%
    filter(!is.na(gene))
}

# Find and filter files
get_filtered_feature_files <- function(models_dir, target_gene_type, target_cell_type) {
  # Look for cell type directory
  cell_type_dirs <- dir_ls(models_dir, type = "directory")
  
  # No cell type directories found
  if (length(cell_type_dirs) == 0) {
    return(list())
  }

  all_found_files_list <- list()

  # Iterate over each cell type directory
  for (cell_type_dir_path in cell_type_dirs) {
    cell_type_name <- path_file(cell_type_dir_path)
    
    # Skip if this is not the target cell type set
    if (cell_type_name != target_cell_type) {
      next
    }
    
    # Look for gene type directory
    gene_type_dirs <- dir_ls(cell_type_dir_path, type = "directory")
    
    # Iterate over each gene type directory
    for (gene_type_dir_path in gene_type_dirs) {
      gene_type_name <- path_file(gene_type_dir_path)
      
      # Skip if this is not the target gene type
      if (gene_type_name != target_gene_type) {
        next
      }
      
      # Look for model directories
      model_dirs <- dir_ls(gene_type_dir_path, type = "directory")
      
      # Iterate over each model directory
      for (model_dir_path in model_dirs) {
        model_name_from_path <- path_file(model_dir_path)
        
        if (dir_exists(model_dir_path)) {
          # Create regex from glob pattern
          file_name_glob_pattern <- "feature_importance*.csv"
          regex_pattern_for_list_files <- glob2rx(file_name_glob_pattern)
          
          # List all feature importance files in the model directory
          files_in_model_dir <- list.files(path = model_dir_path, 
                                            pattern = regex_pattern_for_list_files, 
                                            full.names = TRUE, 
                                            recursive = FALSE,
                                            ignore.case = TRUE)
          
          # Add feature importance files to list with the model name
          if (length(files_in_model_dir) > 0) {
            for (f_path in files_in_model_dir) {
              # Append file path and model name to list
              all_found_files_list[[length(all_found_files_list) + 1]] <- list(
                file_path = as.character(f_path), 
                model_name = model_name_from_path
              )
            }
          }
        }
      }
    }
  }
  
  # Return list of all feature importance files found
  return(all_found_files_list) 
}


# Scale single model score table based on rank
scale_scores_for_model <- function(df, model_name, is_dense, threshold = 0.7) {
  # Create a score column name based on the model name
  score_col <- paste0(model_name, "_score")

  # Handle empty input dataframe by returning a correctly structured empty tibble
  if (nrow(df) == 0) {
    return(tibble(
      feature_id = character(),
      gene = character(),
      sublineage = character(),
      # This dynamically parses the score column
      !!score_col := numeric()
    ))
  }

  # Base data for processing
  df_out <- df %>%
    select(feature_id, gene, sublineage, abs_value)

  # Only consider features with nonzero importance
  nonzero_features <- df_out %>% filter(abs_value > 0)
  
  # If no features have importance return all scores as 0
  if (nrow(nonzero_features) == 0) {
    return(df_out %>% 
             select(feature_id, gene, sublineage) %>% 
             mutate(!!score_col := 0))
  }

  # For dense models keep only features up to the cumulative threshold which is a percentage of total importance
  if (is_dense) {
    nonzero_features <- nonzero_features %>% arrange(desc(abs_value))
    cum_importance <- cumsum(nonzero_features$abs_value) / sum(nonzero_features$abs_value)
    selected_idx <- which(cum_importance <= threshold)

    # This is sort of arbitrary but should never happen
    if (length(selected_idx) == 0) selected_idx <- 1
    features_to_scale <- nonzero_features[selected_idx, ]
  } else {
    # Default scenario for sparse models
    features_to_scale <- nonzero_features
  }

  # Scale the features to scale between 0.5 and 1 based on rank
  n <- nrow(features_to_scale)
  if (n == 1) {
    scaled_scores <- 1
    ranks <- 1
  } else {
    features_to_scale <- features_to_scale %>% arrange(desc(abs_value))
    ranks <- seq_len(n)
    scaled_scores <- 1.0 - 0.5 * (ranks - 1) / (n - 1)
  }

  # Create a tibble with scores for selected features
  score_update <- tibble(
    feature_id = features_to_scale$feature_id, 
    !!score_col := scaled_scores
  )

  # Left join scores with original base data and replace NA with 0 and select final columns
  df_out %>%
    select(feature_id, gene, sublineage) %>%
    left_join(score_update, by = "feature_id") %>%
    mutate(!!score_col := replace_na(.data[[score_col]], 0))
}


# Merge and compute meta score
merge_and_score <- function(score_tables, all_models, fold_changes_data = NULL) {
  # Get expected score column names from all models that were supposed to be processed
  expected_score_cols <- paste0(all_models, "_score")

  # Perform join on all score tables
  merged <- reduce(score_tables, full_join, by = c("feature_id", "gene", "sublineage"))

  # Identify which score columns are actually present after the merge
  current_score_cols <- grep("_score$", names(merged), value = TRUE)
  
  # Identify and add any missing score columns
  missing_cols <- setdiff(expected_score_cols, current_score_cols)
  if (length(missing_cols) > 0) {
    # Dynamically create a list of new columns to add
    new_cols_to_add <- setNames(rep(list(0), length(missing_cols)), missing_cols)
    merged <- merged %>%
      add_column(!!!new_cols_to_add)
  }

  # Add fold change information
  if (!is.null(fold_changes_data)) {
    merged <- merged %>%
      left_join(fold_changes_data %>% select(Feature, Value), 
                by = c("feature_id" = "Feature")) %>%
      rename(fold_change = Value)
  } else {
    # Fold changes are missing should never happen but isnt critical
    merged <- merged %>%
      mutate(fold_change = NA)
  }
  
  # Replace NA scores with 0 for model scores only
  merged <- merged %>%
    mutate(across(all_of(expected_score_cols), ~replace_na(.x, 0)))

  # Meta score is the sum of all model scores
  merged <- merged %>%
    mutate(meta_score = rowSums(select(., all_of(expected_score_cols))))

  # Model prevalence of feature
  merged <- merged %>%
    mutate(n_models_occur = rowSums(select(., all_of(expected_score_cols)) > 0))

  return(merged)
}


# Filter meta scores and write to file
filter_and_write <- function(df, min_models, all_models, out_dir) {
  # Filter features by minimum model prevalence and sort by meta score
  filtered <- df %>%
    filter(n_models_occur >= min_models) %>%
    arrange(desc(meta_score))

  # Prepare score column names
  score_cols <- paste0(all_models, "_score")
  all_cols <- c("feature_id", "gene", "sublineage", "meta_score", "n_models_occur", "fold_change", sort(score_cols))

  # Return nothing if all features occur in less than minimum model prevalence requirement
  if (nrow(filtered) == 0) {
    message(sprintf("No features met min_models_occurrence of %s.", min_models))
    return(NULL)
  }
  # Write filtered results to file
  write_csv(select(filtered, all_of(all_cols)), file.path(out_dir, "meta_scores.csv"))

  # Return filtered tibble df
  return(filtered[all_cols])
}


# Finds features that are in at least minimum number of models and writes to file
find_feature_intersection <- function(models_dir, output_dir, min_models_occurrence,cell_type_filter_val, gene_type_filter_val, type = NULL) {
  # Create output directory if it doesnt exist
  dir_create(output_dir, recurse = TRUE)

  # Find all relevant feature files for this subset
  file_list <- get_filtered_feature_files(models_dir, gene_type_filter_val, cell_type_filter_val)

  # Throw warning and return nothing if there are no feature importance files found
  if (length(file_list) == 0) {
    message(sprintf("No files for gene type=%s and cell type set=%s", gene_type_filter_val, cell_type_filter_val))
    return(NULL)
  }

  # Get model names from the file list
  models <- file_list %>%
    # Returns a list of model names for file list entries
    map_chr(~ .x$model_name) %>%
    unique() %>%
    sort()

  # Load and combine all features
  all_features <- file_list %>%
    map(~ load_and_extract_features(.x$file_path, .x$model_name)) %>%
    compact() %>%
    bind_rows()

  # Return null if no features were loaded
  if (nrow(all_features) == 0) {
    message("No features found after loading files.")
    return(NULL)
  }
  
  # Load fold change data
  fold_changes_path <- file.path(
    "output", "4. fold changes", type,
    cell_type_filter_val, gene_type_filter_val,
    "fold_changes.csv"
  )

  # Read from path if fold changes file exists
  fold_changes_data <- tryCatch({
    as_tibble(fread(fold_changes_path))
  }, error = function(e) {
    stop(sprintf("Error reading fold change file: %s", e$message))
  })

  # Determine which models are dense
  density_map <- setNames(
    vapply(models, function(m) sum(all_features$abs_value[all_features$model_type == m] > 0, na.rm = TRUE) > 100, logical(1)),
    models
  )

  # Scale feature scores for each model
  score_tables <- map(models, function(m) {
    model_df <- filter(all_features, model_type == m)
    scale_scores_for_model(model_df, m, density_map[[m]])
  })

  # Merge model scores and fold changes
  merged <- merge_and_score(score_tables, models, fold_changes_data)

  # Filter and write results
  result <- filter_and_write(merged, min_models_occurrence, models, output_dir)
  return(result)
}


# Error handling wrapper for a single dataset combination for cell type and gene type
process_dataset_cell_gene_combination <- function(type, cell_type, gene_type, base_models_dir, base_output_dir, min_models_param) {
  # Set up paths
  models_dir <- file.path(base_models_dir, type)
  output_dir <- file.path(base_output_dir, type, cell_type, gene_type)
  
  # Check if models directory exists
  if (!dir.exists(models_dir)) {
    warning(sprintf("Models directory not found: %s.", models_dir))
    return(NULL)
  }
  
  message(sprintf("Processing dataset type: %s, cell type: %s, gene type: %s", type, cell_type, gene_type))

  # Do all the intersector work with this function
  result <- find_feature_intersection(
    models_dir = models_dir,
    output_dir = output_dir,
    min_models_occurrence = min_models_param,
    cell_type_filter_val = cell_type,
    gene_type_filter_val = gene_type,
    type = type
  )
  
  if (is.null(result) || nrow(result) == 0) {
    message(sprintf("No results for dataset type: %s, cell type: %s, gene type: %s.", type, cell_type, gene_type))
    return(NULL)
  }
  
  return(file.path(output_dir, "meta_scores.csv"))
}


generate_and_save_sublineage_color_map <- function(meta_score_files, output_path) {
  # Extract unique sublineages from the meta score files
  sublineages <- meta_score_files %>%
    keep(file.exists) %>%
    map_dfr(~ as_tibble(fread(.x)) %>% select(sublineage)) %>%
    pull(sublineage) %>%
    unique() %>%
    sort() %>%
    discard(is.na)

  # Create directory if it does not exist
  dir_create(dirname(output_path), recurse = TRUE)

  # Number of unique sublineages
  n <- length(sublineages)

  # Use colorblind friendly seed colors
  seed_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")

  # Generate a color palette using polychrome
  col_vector <- Polychrome::createPalette(n, seed_colors)
  color_map <- setNames(col_vector[seq_along(sublineages)], sublineages)

  # Save the color map to an RDS file
  saveRDS(color_map, file = output_path)

  return(color_map)
}


# Main function to run intersector for all combinations
run_intersector <- function(
  base_models_dir = "output/2. models",
  base_output_dir = "output/3. intersector",
  min_models_param = 4,
  dataset_types_to_process = c("ctnorm_global", "ctnorm_relative", "read_depth"),
  cell_types_to_process = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
  gene_types_to_process = c("metabolic", "nonmetabolic", "random"),
  sublineage_color_map_path = "output/3. intersector/sublineage_colors.rds"
) {
  # Create all combinations
  meta_score_files <- expand.grid(
    type = dataset_types_to_process,
    cell_type = cell_types_to_process,
    gene_type = gene_types_to_process,
    stringsAsFactors = FALSE
  ) %>%
    # Process each combination
    pmap(~ process_dataset_cell_gene_combination(..1, ..2, ..3, base_models_dir, base_output_dir, min_models_param)) %>%
    compact()
  
  # For all important sublineages we found create a polychrome color map
  generate_and_save_sublineage_color_map(meta_score_files, sublineage_color_map_path)

  message("Finished intersector analysis.")
}


# Run the intersector
run_intersector(
  base_models_dir = "output/2. models",
  base_output_dir = "output/3. intersector",
  min_models_param = 4,
  cell_types_to_process = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
  gene_types_to_process = c("metabolic", "nonmetabolic", "random"),
  sublineage_color_map_path = "output/3. intersector/sublineage_colors.rds"
)
