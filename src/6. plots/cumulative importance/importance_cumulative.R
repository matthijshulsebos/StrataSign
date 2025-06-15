library(ggplot2)
library(dplyr)
library(readr)
library(fs)
library(stringr)
library(scales)

# Fallback: use ‘a’ if it’s a non-empty string, otherwise ‘b’
`%||%` <- function(a, b) {
  if (!is.null(a) && !is.na(a) && nzchar(trimws(a))) a else b
}

# Read and perform initial validation of feature importance data
read_fi_data <- function(file_path) {
  fi_data_raw <- tryCatch({
    readr::read_csv(file_path, show_col_types = FALSE)
  }, error = function(e_read) {
    warning(paste0("Error reading CSV '", file_path, "': ", e_read$message))
    return(NULL)
  })
  if (is.null(fi_data_raw) || nrow(fi_data_raw) == 0) {
    message("Skipping empty or unreadable file: ", file_path)
    return(NULL)
  }
  return(fi_data_raw)
}

# Extract metadata from file path
extract_path_metadata <- function(file_path, model_output_dir) {
  relative_path_str <- as.character(fs::path_rel(file_path, start = model_output_dir))
  path_parts <- strsplit(relative_path_str, "/|\\\\")[[1]]
  
  # New path format: {normalization}/{cell_type_group}/{gene_set}/{model_algo}/feature_importance_*.csv
  metadata <- list(
    data_type = if (length(path_parts) >= 1) path_parts[1] else NA_character_,  # normalization type
    cell_type_group = if (length(path_parts) >= 2) path_parts[2] else NA_character_,
    gene_set_ext = if (length(path_parts) >= 3) path_parts[3] else NA_character_,
    model_algo = if (length(path_parts) >= 4) path_parts[4] else NA_character_  # model algorithm
  )
  
  essential_components <- c(metadata$data_type, metadata$cell_type_group, metadata$gene_set_ext, metadata$model_algo)
  if(any(is.na(essential_components)) || any(sapply(essential_components, function(x) trimws(x) == ""))) {
      warning(paste0("Path component extraction failed for: '", file_path, "'. ",
                     "Extracted: dataType='", metadata$data_type, 
                     "', cellGroup='", metadata$cell_type_group, 
                     "', geneSet='", metadata$gene_set_ext,
                     "', modelAlgo='", metadata$model_algo, 
                     "'. Skipping this file."))
      return(NULL)
  }
  return(metadata)
}

# Data preparation for plotting
prepare_fi_plot_data <- function(fi_data_raw, threshold_val) {
  fi_data <- fi_data_raw
  if (!("Feature" %in% names(fi_data)) && ncol(fi_data) >= 1) {
    if (ncol(fi_data) >= 2 && !("Value" %in% names(fi_data))) {
      fi_data <- fi_data %>% dplyr::rename(Feature = 1, Value = 2)
    } else if (!("Value" %in% names(fi_data))) {
      fi_data <- fi_data %>% dplyr::rename(Feature = 1)
    }
  }
  if (!("Feature" %in% names(fi_data)) || !("Value" %in% names(fi_data))) {
    message("Skipping file due to missing 'Feature' or 'Value' column.")
    return(NULL)
  }
  
  fi_data$Value <- abs(as.numeric(fi_data$Value))
  fi_data <- fi_data %>% dplyr::filter(!is.na(Value)) %>% dplyr::arrange(dplyr::desc(Value))
  
  if (nrow(fi_data) == 0) {
    message("No valid data after processing NA values.")
    return(NULL)
  }
  
  total_importance <- sum(fi_data$Value, na.rm = TRUE)
  
  if (total_importance == 0) {
    message("Total importance is zero, preparing minimal plot data.")
    fi_data_plot <- tibble::tibble(Feature = "ORIGIN", Value = 0, cumulative = 0, cumulative_pct = 0, feature_index = 0)
    x_max_rank <- 0
    threshold_idx_in_plot <- NA_integer_
  } else {
    fi_data$cumulative <- cumsum(fi_data$Value)
    fi_data$cumulative_pct <- fi_data$cumulative / total_importance
    fi_data$feature_index <- 1:nrow(fi_data)
    origin_point <- tibble::tibble(Feature = "ORIGIN", Value = 0, cumulative = 0, cumulative_pct = 0, feature_index = 0)
    fi_data_plot <- dplyr::bind_rows(origin_point, fi_data) %>% dplyr::arrange(feature_index)
    
    threshold_idx_in_plot <- which(fi_data_plot$cumulative_pct >= threshold_val)[1]
    
    idx_100_pct_rows <- which(fi_data_plot$cumulative_pct >= 0.99999)
    x_max_rank <- if (length(idx_100_pct_rows) > 0) {
      fi_data_plot$feature_index[min(idx_100_pct_rows)]
    } else {
      max(fi_data_plot$feature_index, na.rm = TRUE)
    }
  }
  
  if (nrow(fi_data_plot) > 1 && (is.na(x_max_rank) || x_max_rank < 1)) {
    x_max_rank <- max(fi_data_plot$feature_index, na.rm = TRUE)
  }

  if (is.na(x_max_rank) || (!is.finite(x_max_rank)) || (x_max_rank < 0) || (x_max_rank == 0 && nrow(fi_data_plot) > 1) ) {
    message("WARNING: x_max_rank ('", x_max_rank, "') is NA, non-finite, negative, or zero when data exists beyond origin. Attempting to reset.")
    if (x_max_rank == 0 && nrow(fi_data_plot) > 1) {
        x_max_rank <- max(fi_data_plot$feature_index, na.rm = TRUE)
        message("INFO: x_max_rank reset to max feature_index: ", x_max_rank)
    }
    if (is.na(x_max_rank) || !is.finite(x_max_rank) || x_max_rank < 0 ) {
         message("ERROR: x_max_rank is still invalid after fallback. Skipping plot.")
         return(NULL)
    }
  }

  return(list(
    fi_data_plot = fi_data_plot, 
    total_importance = total_importance,
    threshold_idx_in_plot = threshold_idx_in_plot,
    x_max_rank = x_max_rank
  ))
}

# Generate plots for all matching files
plot_cumulative_importance <- function(model_output_dir,
                                       figures_dir,
                                       threshold) {
  # Ensure figures_dir is a valid path
  if (!is.character(figures_dir) || length(figures_dir) != 1 ||
      !nzchar(trimws(figures_dir)) || identical(trimws(figures_dir), "FALSE")) {
    warning("Invalid figures_dir.")
    return()
  }
  # Create base dir
  tryCatch(dir_create(figures_dir, recurse = TRUE),
           error = function(e) {
             warning("Cannot create figures_dir: ", e$message)
             return()
           })
  message("Plots will be saved under: ", figures_dir)

  # Find all feature_importance files
  all_fi_files <- dir_ls(model_output_dir,
                         recurse = TRUE,
                         regexp = "feature_importance.*\\.csv$",
                         fail = FALSE)

  message("Processing all ", length(all_fi_files), " found files.")

  if (length(all_fi_files) == 0) {
    message("No files found.")
    return()
  }

  # Process each file
  for (file_path in all_fi_files) {
    message("\nProcessing: ", file_path)
    
    fi_data_raw <- read_fi_data(file_path)
    if (is.null(fi_data_raw)) next
    
    metadata <- extract_path_metadata(file_path, model_output_dir)
    if (is.null(metadata)) next
    
    effective_gene_set_label <- metadata$gene_set_ext # Directly use gene_set_ext from metadata
    effective_gene_set_label <- effective_gene_set_label %||% "unknown_geneset"

    fi_data_prepared <- prepare_fi_plot_data(fi_data_raw, threshold)
    if (is.null(fi_data_prepared)) next
    
    fi_data_plot <- fi_data_prepared$fi_data_plot
    total_importance <- fi_data_prepared$total_importance
    threshold_idx_in_plot <- fi_data_prepared$threshold_idx_in_plot
    x_max_rank <- fi_data_prepared$x_max_rank
    
    # Determine coordinates for 70% line
    vline_x_coord <- NA_real_
    if (!is.na(threshold_idx_in_plot) && threshold_idx_in_plot > 0 && threshold_idx_in_plot <= nrow(fi_data_plot)) {
        vline_x_coord <- fi_data_plot$feature_index[threshold_idx_in_plot]
    } else if (total_importance > 0 && !is.na(threshold)) { # If threshold not met but data exists, place at end
        vline_x_coord <- max(fi_data_plot$feature_index, na.rm = TRUE)
    }

    # Subtitle construction
    total_features_in_original_fi_data <- nrow(fi_data_raw %>% filter(!is.na(Value) & Value !=0))
    density_in_fi_data <- ifelse(nrow(fi_data_raw) > 0, total_features_in_original_fi_data / nrow(fi_data_raw), 0)
    features_at_threshold_label_count_val <- if(!is.na(vline_x_coord)) vline_x_coord else "N/A"
    
    feature_count_label_for_subtitle <- ifelse(density_in_fi_data > 0.1 && total_features_in_original_fi_data > 0, 
                                   paste0(features_at_threshold_label_count_val, " of ", total_features_in_original_fi_data, " non-zero (", 
                                          if(!is.na(vline_x_coord) && total_features_in_original_fi_data > 0 && vline_x_coord > 0) round(vline_x_coord / total_features_in_original_fi_data * 100) else "N/A", "%)"),
                                   paste0(total_features_in_original_fi_data, " non-zero features"))
    
    plot_subtitle <- paste0("Cell Group: ", metadata$cell_type_group, " / Gene Set: ", effective_gene_set_label,
                           "\nModel density: ", round(density_in_fi_data * 100), "% (", feature_count_label_for_subtitle, ")")

    if (is.na(vline_x_coord) || !is.finite(vline_x_coord)) {
        message("INFO: Vertical line for threshold ", threshold*100, "% (at x=", vline_x_coord, ") might not be drawn.")
    }

    # Create the plot
    p <- ggplot(fi_data_plot, aes(x = feature_index, y = cumulative_pct)) +
      geom_line(color = "black", size = 1) +
      {
        if (!is.na(vline_x_coord) && is.finite(vline_x_coord) && vline_x_coord >= 0 && 
            (x_max_rank == 0 || vline_x_coord <= x_max_rank) ) {
          geom_vline(xintercept = vline_x_coord, linetype = "dashed", color = "firebrick", size = 1)
        }
      } +
      scale_y_continuous(
        labels = scales::percent, 
        limits = c(0, 1.01), 
        expand = expansion(mult = c(0, 0.01)) 
      ) +
      scale_x_continuous(
        limits = c(0, if(x_max_rank > 0) x_max_rank else if(nrow(fi_data_plot) > 1) 1 else 0.1),
        expand = expansion(mult = c(0.01, 0.01)) 
      ) +
      labs(
        title = paste0("Cumulative Feature Importance: ", metadata$data_type, " - ", metadata$model_algo),
        subtitle = plot_subtitle,
        x = "Feature Rank (by importance)",
        y = "Cumulative Importance (%)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10)
      )

    # Write file to output directory - include model algorithm to prevent overwriting
    sane_data_type <- gsub("[^a-zA-Z0-9_.-]", "_", metadata$data_type)
    sane_model_algo <- gsub("[^a-zA-Z0-9_.-]", "_", metadata$model_algo)
    sane_cell_group <- gsub("[^a-zA-Z0-9_.-]", "_", metadata$cell_type_group)
    sane_gs_label <- gsub("[^a-zA-Z0-9_.-]", "_", effective_gene_set_label)
    
    # Include model algorithm in the path to prevent overwriting
    specific_output_subdir <- file.path(figures_dir, sane_data_type, sane_cell_group, 
                                        sane_gs_label, sane_model_algo)
    
    dir_creation_success <- FALSE
    tryCatch({
      if (!dir_exists(specific_output_subdir)) dir_create(specific_output_subdir, recurse = TRUE)
      dir_creation_success <- TRUE 
    }, error = function(e_dir) {
      warning(paste0("Error creating output directory '", specific_output_subdir, "': ", e_dir$message))
    })

    if (!dir_creation_success) next 

    base_input_filename <- gsub("\\.csv$", "", basename(file_path)) 
    
    # Include model algorithm in filename as additional safeguard
    plot_filename <- paste0("cumulative_", sane_model_algo, "_", base_input_filename, ".png")
    output_filepath <- file.path(specific_output_subdir, plot_filename)
    
    message("Attempting to save plot to: ", output_filepath)
    tryCatch({
      ggsave(filename = output_filepath, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
      message("Successfully saved plot: ", output_filepath)
    }, error = function(e_save) { 
      warning(paste0("Error during ggsave for '", output_filepath, "': ", e_save$message))
    })
  }
  message("\nFinished processing all files.")
}

# Run plotting function
plot_cumulative_importance(
  model_output_dir = "output/2. models", 
  figures_dir = "output/6. plots/cumulative importance",
  threshold = 0.7 
)
