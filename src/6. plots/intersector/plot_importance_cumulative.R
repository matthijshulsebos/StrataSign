library(ggplot2)
library(dplyr)
library(readr)
library(fs)
library(stringr)


plot_cumulative_importance <- function(model_output_dir = "output/models",
                                      figures_dir = "output/figures/intersector/cumulative",
                                      version_filter = "metabolic",
                                      threshold = 0.7) {
  
  dir_create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_fi_files <- dir_ls(model_output_dir, 
                         recurse = TRUE, 
                         regexp = "feature_importance.*\\.csv$", 
                         fail = FALSE)
  
  if (!is.null(version_filter) && version_filter != "") {
    pattern <- paste0("/", version_filter, "/")
    all_fi_files <- all_fi_files[grepl(pattern, all_fi_files, fixed = FALSE)]
  }
  
  if (length(all_fi_files) == 0) {
    message("No feature importance files found matching criteria.")
    return()
  }
  
  file_metadata <- tibble(
    file_path = all_fi_files
  ) %>%
  mutate(
    relative_path = path_rel(file_path, model_output_dir),
    path_parts = strsplit(as.character(relative_path), "/|\\\\"),
    model_type = sapply(path_parts, function(x) if(length(x) >= 1) x[1] else NA_character_),
    dataset = sapply(path_parts, function(x) if(length(x) >= 2) x[2] else NA_character_),
    version = sapply(path_parts, function(x) if(length(x) >= 3) x[3] else NA_character_)
  ) %>% select(-path_parts, -relative_path)


  for (model_t in unique(file_metadata$model_type)) {
    model_files_subset <- file_metadata %>% filter(model_type == model_t)
    
    for (i in 1:nrow(model_files_subset)) {
      fi_data <- read_csv(model_files_subset$file_path[i], show_col_types = FALSE)
      
      if (nrow(fi_data) == 0) next
      
      current_dataset <- model_files_subset$dataset[i]
      current_version <- model_files_subset$version[i]
      
      if (!("Feature" %in% names(fi_data)) && ncol(fi_data) >= 2) {
        fi_data <- fi_data %>% rename(Feature = 1, Value = 2)
      } else if (!("Feature" %in% names(fi_data))) {
        message("Skipping file due to missing 'Feature' column: ", model_files_subset$file_path[i])
        next
      }
      
      fi_data$Value <- abs(fi_data$Value)
      fi_data <- fi_data %>% arrange(desc(Value))
      
      total_importance <- sum(fi_data$Value, na.rm = TRUE)
      if (total_importance == 0) next

      fi_data$cumulative <- cumsum(fi_data$Value)
      fi_data$cumulative_pct <- fi_data$cumulative / total_importance
      fi_data$feature_index <- 1:nrow(fi_data)
      
      total_features <- nrow(fi_data)
      nonzero_features <- sum(fi_data$Value > 0, na.rm = TRUE)
      density <- ifelse(total_features > 0, nonzero_features / total_features, 0)
      
      origin_point <- tibble(
        Feature = "ORIGIN", Value = 0, cumulative = 0, 
        cumulative_pct = 0, feature_index = 0
      )
      
      fi_data_plot <- bind_rows(origin_point, fi_data) %>%
        arrange(feature_index)
      
      threshold_idx <- which(fi_data_plot$cumulative_pct >= threshold)[1]
      if (is.na(threshold_idx)) threshold_idx <- nrow(fi_data_plot)
      
      max_idx <- which(fi_data_plot$cumulative_pct >= 0.999)[1]
      if (is.na(max_idx)) max_idx <- nrow(fi_data_plot)
      
      feature_count_label <- ifelse(density > 0.1 && nonzero_features > 0, 
                                   paste0(threshold_idx -1, " of ", nonzero_features, " non-zero (", 
                                          round((threshold_idx-1) / nonzero_features * 100), "%)"),
                                   paste0(nonzero_features, " non-zero features"))
      
      p <- ggplot(fi_data_plot, aes(x = feature_index, y = cumulative_pct)) +
        geom_line(color = "black", size = 1) +
        geom_hline(yintercept = threshold, linetype = "dashed", color = "#D32F2F", size = 1) +
        geom_vline(xintercept = threshold_idx -1, linetype = "dashed", color = "#E53935", size = 0.8) + 
        annotate("text", x = (threshold_idx-1) * 0.5, y = threshold + 0.05, 
                label = paste0(round(threshold * 100), "% threshold"), color = "#D32F2F") +
        annotate("text", x = (threshold_idx-1) + min((threshold_idx-1) * 0.1, 5), 
                y = threshold/2, 
                label = paste0("Top ", threshold_idx-1, " features"), 
                color = "black", hjust = 0) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        scale_x_continuous(limits = c(0, max_idx * 1.05)) +
        labs(
          title = paste0(model_t, " Cumulative Feature Importance"),
          subtitle = paste0("Dataset: ", current_dataset, " / ", current_version, 
                           "\nModel density: ", round(density * 100), "% (", feature_count_label, ")"),
          x = "Feature Rank",
          y = "Cumulative Importance"
        ) +
        theme_minimal() +
        theme(
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 10),
          plot.title = element_text(face = "bold")
        )
      
      filename <- paste0("cumulative_importance_", model_t, "_", current_dataset, "_", current_version, ".png")
      ggsave(file.path(figures_dir, filename), p, width = 10, height = 6, dpi = 300)
      message("Created cumulative plot: ", filename, " in: ", figures_dir)
    }
  }
}
