library(dplyr)
library(tidyr)
library(data.table)
library(fs)
library(ggplot2)

find_feature_intersection <- function(models_dir = "output/models",
                                     output_dir = "output/intersector",
                                     figures_dir = "output/figures/intersector",
                                     top_n_features = 100,
                                     min_models = 2,
                                     version_filter = "metabolic",
                                     dataset_filter = "all_clusters") {
  
  # Create output directories
  dir_create(output_dir, recursive = TRUE)
  dir_create(figures_dir, recursive = TRUE)
  
  # Find all feature importance files across model types
  all_fi_files <- dir_ls(models_dir, 
                         recurse = TRUE, 
                         regexp = "feature_importance.*\\.csv$", 
                         fail = FALSE)
  
  # Filter for specific version if specified
  if (!is.null(version_filter)) {
    pattern <- paste0("/", version_filter, "/")
    all_fi_files <- all_fi_files[grep(pattern, all_fi_files, fixed = FALSE)]
  }
  
  # Extract metadata from file paths
  file_metadata <- data.frame(
    file_path = all_fi_files,
    model_type = sapply(all_fi_files, function(path) {
      parts <- strsplit(path_rel(path, models_dir), "/|\\\\")[[1]]
      return(parts[1])
    }),
    dataset = sapply(all_fi_files, function(path) {
      parts <- strsplit(path_rel(path, models_dir), "/|\\\\")[[1]]
      return(parts[2])
    }),
    version = sapply(all_fi_files, function(path) {
      parts <- strsplit(path_rel(path, models_dir), "/|\\\\")[[1]]
      return(parts[3])
    }),
    stringsAsFactors = FALSE
  )
  
  # Filter for specific dataset if specified
  if (!is.null(dataset_filter)) {
    file_metadata <- file_metadata[file_metadata$dataset == dataset_filter, ]
    all_fi_files <- file_metadata$file_path
  }
  
  if (length(all_fi_files) == 0) {
    stop("No feature importance files found matching the filters")
  }
  
  # Initialize combined feature data
  all_features <- data.table()
  
  # First pass: collect all feature importance data and assess model density
  model_density <- list()
  
  for (i in 1:nrow(file_metadata)) {
    file_info <- file_metadata[i, ]
    file_path <- file_info$file_path
    model <- file_info$model_type
    dataset <- file_info$dataset
    version <- file_info$version
    
    # Load feature importance data
    features <- fread(file_path)
    
    # Skip if empty
    if (nrow(features) == 0) {
      next
    }
    
    # Add metadata
    features <- features %>%
      mutate(
        model = model,
        dataset = dataset,
        version = version,
        abs_value = abs(Value),
        gene = sapply(strsplit(Feature, "@"), `[`, 1),
        feature_id = Feature
      )
    
    # Add to combined data
    all_features <- rbind(all_features, features, fill = TRUE)
    
    # Calculate model density (% of features with importance)
    total_features <- nrow(features)
    nonzero_features <- sum(features$abs_value > 1e-10)
    density <- nonzero_features / total_features
    
    # Store density information
    if (!model %in% names(model_density)) {
      model_density[[model]] <- list(
        total_features = total_features,
        nonzero_features = nonzero_features,
        density = density,
        # Change from 50 to 100 non-zero features threshold
        is_dense = nonzero_features > 100  # Now using >100 features instead of >50
      )
    } else {
      # Update if multiple files for same model
      model_density[[model]]$total_features <- model_density[[model]]$total_features + total_features
      model_density[[model]]$nonzero_features <- model_density[[model]]$nonzero_features + nonzero_features
      model_density[[model]]$density <- model_density[[model]]$nonzero_features / model_density[[model]]$total_features
      model_density[[model]]$is_dense <- model_density[[model]]$nonzero_features > 100  # Also updated here
    }
  }
  
  # Process features using density-aware normalization
  normalized_features <- all_features %>%
    group_by(model, dataset, version) %>%
    mutate(
      # Apply the 70% cumulative importance cutoff for dense models
      normalized_score = case_when(
        # For dense models (>100 non-zero features)
        model_density[[first(model)]]$is_dense ~ {
          # Sort by absolute importance (descending)
          sorted_idx <- order(abs_value, decreasing = TRUE)
          sorted_vals <- abs_value[sorted_idx]
          total_importance <- sum(sorted_vals)
          
          # Calculate cumulative importance
          cum_importance <- cumsum(sorted_vals) / total_importance
          
          # Find cutoff at 70% cumulative importance
          cutoff_idx <- which(cum_importance >= 0.7)[1] 
          if (is.na(cutoff_idx)) cutoff_idx <- length(sorted_vals)
          
          # Create mask for important features
          keep_mask <- rep(FALSE, length(sorted_vals))
          keep_mask[1:cutoff_idx] <- TRUE
          keep <- rep(FALSE, n())
          keep[sorted_idx] <- keep_mask
          
          # Normalize scores for kept features, zero out others
          result <- rep(0, n())
          if (sum(keep) > 0) {
            min_val <- min(abs_value[keep])
            max_val <- max(abs_value[keep])
            range_val <- max_val - min_val
            if (range_val > 0) {
              result[keep] <- 0.5 + 0.5 * (abs_value[keep] - min_val) / range_val
            } else {
              result[keep] <- 0.75 # Default if all values are the same
            }
          }
          result
        },
        # For sparse models (≤100 non-zero features)
        abs_value > 0 ~ {
          min_val <- min(abs_value[abs_value > 0])
          max_val <- max(abs_value[abs_value > 0])
          range_val <- max_val - min_val
          if (range_val > 0) {
            0.5 + 0.5 * (abs_value - min_val) / range_val
          } else {
            ifelse(abs_value > 0, 0.75, 0)
          }
        },
        TRUE ~ 0
      )
    ) %>%
    ungroup()
  
  # Calculate meta-scores based on summing normalized importance
  feature_scores <- normalized_features %>%
    group_by(feature_id) %>%
    summarise(
      gene = first(gene),
      # Count models and non-zero appearances
      n_models = length(unique(model)),
      n_datasets = length(unique(dataset)),
      n_models_nonzero = sum(normalized_score > 0),
      
      # Sum normalized importance across all models
      meta_score = sum(normalized_score),
      
      # For compatibility with original output
      sum_rank = sum(normalized_score),
      sparsity_bonus = sum(normalized_score),
      model_breadth = log1p(n_models_nonzero),
      
      .groups = "drop"
    ) %>%
    filter(n_models >= min_models) %>%
    arrange(desc(meta_score))
  
  # Get the top N features
  top_features <- feature_scores %>%
    head(top_n_features)
  
  # Create a single comprehensive CSV
  unified_results <- feature_scores %>%
    # Add top feature flag
    mutate(is_top = row_number() <= top_n_features) %>%
    # Join with summary statistics for each model
    left_join(
      normalized_features %>%
        group_by(feature_id, model) %>%
        summarise(
          importance = mean(normalized_score),
          raw_importance = mean(abs_value),
          direction = sign(sum(Value)),
          .groups = "drop"
        ) %>%
        # Create wide format for easier plotting
        pivot_wider(
          id_cols = feature_id,
          names_from = model,
          values_from = c(importance, raw_importance, direction),
          names_sep = "_",
          values_fill = list(importance = 0, raw_importance = 0, direction = 0)
        ),
      by = "feature_id"
    ) %>%
    mutate(
      version_filter = version_filter,
      dataset_filter = dataset_filter,
      creation_date = as.character(Sys.Date())
    )
  
  # Save as a single CSV file
  write.csv(unified_results, file.path(output_dir, "feature_analysis.csv"), row.names = FALSE)
  
  # Return the unified results for further processing if needed
  return(unified_results)
}

# Simplified model sparsity calculation - just for backward compatibility
calculate_model_sparsity <- function(models_dir = "output/models", version_filter = NULL) {
  # Find all feature importance files
  all_fi_files <- dir_ls(models_dir, 
                       recurse = TRUE, 
                       regexp = "feature_importance.*\\.csv$", 
                       fail = FALSE)
  
  # Filter for specific version if provided
  if (!is.null(version_filter)) {
    pattern <- paste0("/", version_filter, "/")
    all_fi_files <- all_fi_files[grep(pattern, all_fi_files, fixed = FALSE)]
  }
  
  # Extract model types from file paths
  model_types <- sapply(all_fi_files, function(path) {
    parts <- strsplit(path_rel(path, models_dir), "/|\\\\")[[1]]
    return(parts[1])
  })
  
  # Initialize sparsity tracking
  model_stats <- data.frame(
    model = unique(model_types),
    total_features = 0,
    nonzero_features = 0,
    sparsity_ratio = 0,
    sparsity_score = 0.5,
    stringsAsFactors = FALSE
  )
  
  return(model_stats)
}

run_intersector <- function(top_n = 100, version_filter = "metabolic") {
  # Run the feature intersection algorithm with version filter
  result <- find_feature_intersection(
    top_n_features = top_n,
    min_models = 2,
    version_filter = version_filter
  )
  
  return(result)
}

# Run the analysis when script is executed
run_intersector()
