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


# Load one file
load_and_extract_features <- function(file_path, model_name) {
  df <- tryCatch({
    as_tibble(fread(file_path))
  }, error = function(e){
    message("Error reading ",file_path,": ",e$message)
    NULL
  })

  if (is.null(df) || !"Feature"%in%names(df)||!"Value"%in%names(df)||nrow(df)==0) {
    return(NULL)
  }
  
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
  
  if (length(cell_type_dirs) == 0) {
    return(list())
  }

  all_found_files_list <- list()

  for (cell_type_dir_path in cell_type_dirs) {
    cell_type_name <- path_file(cell_type_dir_path)
    
    # Skip if this is not the target cell type
    if (cell_type_name != target_cell_type) {
      next
    }
    
    # Look for gene type directory
    gene_type_dirs <- dir_ls(cell_type_dir_path, type = "directory")
    
    for (gene_type_dir_path in gene_type_dirs) {
      gene_type_name <- path_file(gene_type_dir_path) # Get 'metabolic' from '.../all_clusters/metabolic'
      
      # Skip if this is not the target gene type
      if (gene_type_name != target_gene_type) {
        next
      }
      
      # Look for model directories
      model_dirs <- dir_ls(gene_type_dir_path, type = "directory")
      
      for (model_dir_path in model_dirs) {
        model_name_from_path <- path_file(model_dir_path)
        
        if (dir_exists(model_dir_path)) {
          # Create regex from glob pattern
          file_name_glob_pattern <- "feature_importance*.csv"
          regex_pattern_for_list_files <- glob2rx(file_name_glob_pattern)
          
          files_in_model_dir <- list.files(path = model_dir_path, 
                                            pattern = regex_pattern_for_list_files, 
                                            full.names = TRUE, 
                                            recursive = FALSE,
                                            ignore.case = TRUE)
          
          if (length(files_in_model_dir) > 0) {
            for (f_path in files_in_model_dir) {
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
  
  return(all_found_files_list) 
}


# Scale single model score table based on rank
scale_scores_for_model <- function(df, model_name, is_dense, threshold = 0.7) {
  score_col <- paste0(model_name, "_score")

  # Handle empty input dataframe by returning a correctly structured empty tibble
  if (nrow(df) == 0) {
    return(tibble(
      feature_id = character(),
      gene = character(),
      sublineage = character(),
      !!score_col := numeric()
    ))
  }

  # Base data for processing
  df_out <- df %>%
    select(feature_id, gene, sublineage, abs_value)

  # Only consider features with nonzero importance
  nonzero_features <- df_out %>% filter(abs_value > 0)
  
  # If no features have importance, return with all scores as 0
  if (nrow(nonzero_features) == 0) {
    return(df_out %>% 
             select(feature_id, gene, sublineage) %>% 
             mutate(!!score_col := 0))
  }

  # For dense models keep only features up to the cumulative threshold
  if (is_dense) {
    nonzero_features <- nonzero_features %>% arrange(desc(abs_value))
    cum_importance <- cumsum(nonzero_features$abs_value) / sum(nonzero_features$abs_value)
    selected_idx <- which(cum_importance <= threshold)
    if (length(selected_idx) == 0) selected_idx <- 1
    features_to_scale <- nonzero_features[selected_idx, ]
  } else {
    features_to_scale <- nonzero_features
  }

  # Assign scaled scores
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

  # Join scores to the base data, replace NA with 0, and select final columns
  df_out %>%
    select(feature_id, gene, sublineage) %>%
    left_join(score_update, by = "feature_id") %>%
    mutate(!!score_col := replace_na(.data[[score_col]], 0))
}


# Merge and compute meta score
merge_and_score <- function(score_tables, all_models, fold_changes_data = NULL) {
  # Get expected score column names from all models that were supposed to be processed
  expected_score_cols <- paste0(all_models, "_score")

  # Perform the full join on all score tables
  merged <- reduce(score_tables, full_join, by = c("feature_id", "gene", "sublineage"))

  # Identify which score columns are actually present after the merge
  current_score_cols <- grep("_score$", names(merged), value = TRUE)
  
  # Identify and add any missing score columns, initializing them with 0
  missing_cols <- setdiff(expected_score_cols, current_score_cols)
  if (length(missing_cols) > 0) {
    # Dynamically create a list of new columns to add, with 0 as the default value
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
    merged <- merged %>%
      mutate(fold_change = NA)
  }
  
  # Replace NA scores with 0 for all model score columns
  merged <- merged %>%
    mutate(across(all_of(expected_score_cols), ~replace_na(.x, 0)))

  # Calculate meta score as sum of all model scores
  merged <- merged %>%
    mutate(meta_score = rowSums(select(., all_of(expected_score_cols))))

  # Count in how many models the feature was important
  merged <- merged %>%
    mutate(n_models_occur = rowSums(select(., all_of(expected_score_cols)) > 0))

  merged
}


# Filter meta scores and write to file
filter_and_write <- function(df, min_models, all_models, out_dir) {
  # Filter features by minimum model occurrence and sort by meta score
  filtered <- df %>%
    filter(n_models_occur >= min_models) %>%
    arrange(desc(meta_score))

  # Prepare score column names
  score_cols <- paste0(all_models, "_score")
  all_cols <- c("feature_id", "gene", "sublineage", "meta_score", "n_models_occur", "fold_change", sort(score_cols))

  # If no features pass the filter stop the intersector
  if (nrow(filtered) == 0) {
    message(sprintf("No features met min_models_occurrence of %s.", min_models))
    return(NULL)
  }
  # Write filtered results to CSV and return the filtered tibble
  write_csv(select(filtered, all_of(all_cols)), file.path(out_dir, "meta_scores.csv"))
  filtered[all_cols]
}


# Main intersector function
find_feature_intersection <- function(models_dir, output_dir, min_models_occurrence,cell_type_filter_val, gene_type_filter_val, type = NULL) {
  # Ensure output directory exists
  dir_create(output_dir, recurse = TRUE)

  # Find all relevant feature files for this cell/gene type
  file_list <- get_filtered_feature_files(models_dir, gene_type_filter_val, cell_type_filter_val)
  if (length(file_list) == 0) {
    message(sprintf("No files for gene=%s cell=%s", gene_type_filter_val, cell_type_filter_val))
    return(NULL)
  }

  # Get unique model names directly from the file list, ensuring all models are accounted for
  models <- file_list %>%
    map_chr(~ .x$model_name) %>%
    unique() %>%
    sort()

  # Load and combine all features
  all_features <- file_list %>%
    map(~ load_and_extract_features(.x$file_path, .x$model_name)) %>%
    compact() %>%
    bind_rows()
  if (nrow(all_features) == 0) {
    message("No features found after loading files.")
    # Even if no features are found, we may need to create an empty meta_scores file
    # with the correct columns if downstream processing depends on it.
    # For now, we return NULL as before.
    return(NULL)
  }
  # Load fold change data
  fold_changes_path <- file.path(
    "output", "4. fold changes", type,
    cell_type_filter_val, gene_type_filter_val,
    "fold_changes.csv"
  )
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


# Process a single dataset type, cell type, and gene type combination
process_dataset_cell_gene_combination <- function(type, cell_type, gene_type, base_models_dir, base_output_dir, min_models_param) {
  models_dir <- file.path(base_models_dir, type)
  output_dir <- file.path(base_output_dir, type, cell_type, gene_type)
  
  if (!dir.exists(models_dir)) {
    warning(sprintf("Models directory not found: %s.", models_dir))
    return(NULL)
  }
  
  message(sprintf("Processing dataset type: %s, cell type: %s, gene type: %s", type, cell_type, gene_type))
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


# Generate sublineage color map
generate_and_save_sublineage_color_map <- function(meta_score_files, output_path) {
  # Unique sublineages from meta score files
  sublineages <- meta_score_files %>%
    keep(file.exists) %>%
    map_dfr(~ as_tibble(fread(.x)) %>% select(sublineage)) %>%
    pull(sublineage) %>%
    unique() %>%
    sort() %>%
    discard(is.na)

  # Ensure the output directory exists
  dir_create(dirname(output_path), recurse = TRUE)

  # Get qualitative color palettes for categorical data
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rep(col_vector, length.out = length(sublineages))

  # Assign a color to each sublineage
  color_map <- setNames(col_vector[seq_along(sublineages)], sublineages)
  saveRDS(color_map, file = output_path)

  return(color_map)
}


# Main function to run intersector for all combinations
run_intersector <- function(
  base_models_dir = "output/2. models",
  base_output_dir = "output/3. intersector",
  min_models_param = 3,
  dataset_types_to_process = c("ctnorm_global", "ctnorm_relative", "read_depth"),
  cell_types_to_process = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
  gene_types_to_process = c("metabolic", "nonmetabolic", "random"),
  sublineage_color_map_path = "output/3. intersector/sublineage_colors.rds"
) {
  meta_score_files <- expand.grid(
    type = dataset_types_to_process,
    cell_type = cell_types_to_process,
    gene_type = gene_types_to_process,
    stringsAsFactors = FALSE
  ) %>%
    pmap(~ process_dataset_cell_gene_combination(..1, ..2, ..3, base_models_dir, base_output_dir, min_models_param)) %>%
    compact()
  
  generate_and_save_sublineage_color_map(meta_score_files, sublineage_color_map_path)
  message("Finished intersector analysis.")
}


# Run the intersector
run_intersector(
  base_models_dir = "output/2. models",
  base_output_dir = "output/3. intersector",
  min_models_param = 3,
  cell_types_to_process = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
  gene_types_to_process = c("metabolic", "nonmetabolic", "random"),
  sublineage_color_map_path = "output/3. intersector/sublineage_colors.rds"
)
