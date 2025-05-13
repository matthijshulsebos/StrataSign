library(readr)
library(dplyr)
library(fs)

# Source the utility functions for feature name parsing
source("src/0. utils/feature_name_utils.R")

# Define base directory for intersector results
base_intersector_dir <- "output/3. intersector"
output_file_path <- file.path(base_intersector_dir, "common_features_with_ranks.csv")

# Define paths to the meta_scores.csv files
absolute_scores_path <- file.path(base_intersector_dir, "absolute", "meta_scores.csv")
relative_scores_path <- file.path(base_intersector_dir, "relative", "meta_scores.csv")

# Function to load scores, extract gene/cell type, and add rank
load_and_rank_features <- function(file_path, dataset_label) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL) 
  }
  
  scores_df <- readr::read_csv(file_path, show_col_types = FALSE)
  
  if (!"feature_id" %in% names(scores_df) || !"meta_score" %in% names(scores_df)) {
    warning(paste("Required columns ('feature_id', 'meta_score') not found in:", file_path))
    return(NULL)
  }
  
  if (nrow(scores_df) == 0) {
    warning(paste("No data found in:", file_path))
    return(NULL)
  }
  
  # Add rank based on descending meta_score and extract feature components
  ranked_df <- scores_df %>%
    dplyr::select(feature_id, meta_score) %>%
    dplyr::mutate(
      gene = get_gene_from_feature(feature_id),
      simplified_cell_type = get_simplified_sublineage(feature_id)
    ) %>%
    dplyr::arrange(desc(meta_score)) %>%
    dplyr::mutate(rank = min_rank(desc(meta_score))) %>%
    # Rename only score and rank columns with dataset_label
    dplyr::rename_with(~paste0(., "_", dataset_label), c(meta_score, rank))
  
  return(ranked_df)
}

# Load and rank features from the absolute dataset
message(paste("Loading and processing:", absolute_scores_path))
features_absolute <- load_and_rank_features(absolute_scores_path, "absolute")

# Load and rank features from the relative dataset
message(paste("Loading and processing:", relative_scores_path))
features_relative <- load_and_rank_features(relative_scores_path, "relative")

# Proceed only if both dataframes were loaded successfully
if (!is.null(features_absolute) && !is.null(features_relative)) {
  
  # Find common features by joining on feature_id
  common_features_df <- dplyr::inner_join(
    features_absolute,
    features_relative,
    by = "feature_id",
    suffix = c("", "_relative_dup")
  )
  
  if (nrow(common_features_df) > 0) {
    # Select final columns
    comparison_table <- common_features_df %>%
      dplyr::select(
        feature_id, 
        gene,
        simplified_cell_type,
        rank_absolute, 
        meta_score_absolute, 
        rank_relative, 
        meta_score_relative
      ) %>%
      dplyr::arrange(rank_absolute, rank_relative)
      
    message(paste("\nFound", nrow(comparison_table), "common features."))
    
    readr::write_csv(comparison_table, output_file_path)
    message(paste("\nComparison table saved to:", output_file_path))
    
    message("\nFirst few rows of the comparison table:")
    print(head(comparison_table))
    
  } else {
    message("\nNo common features found between 'absolute' and 'relative' intersector results.")
  }
  
} else {
  message("\nCould not perform comparison due to issues loading one or both meta_scores files or parsing feature names.")
}

message("\nComparison script finished.")
