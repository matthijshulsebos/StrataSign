library(ComplexHeatmap)
library(circlize)
library(grid)
library(stats)

# --- Constants ---
COLOR_DARK_BLUE <- "#2166AC"
COLOR_MEDIUM_BLUE <- "#67A9CF"
COLOR_WHITE <- "#FFFFFF"
COLOR_MEDIUM_RED <- "#EF8A62"
COLOR_DARK_RED <- "#B2182B"
SMALLEST_POSITIVE_STEP <- 1e-9
MEDIUM_SATURATION_BAND_PROP <- 0.15

# --- Helper Functions for Color Scale ---

# Default color scale for all-zero or invalid data
default_color_scale <- function() {
  list(breaks = c(0, 0), colors = c(COLOR_WHITE, COLOR_WHITE))
}

# Positive data color scale
positive_color_scale <- function(scale_max) {
  if (round(scale_max, 9) <= 0) return(default_color_scale())
  medium_band_upper <- scale_max * MEDIUM_SATURATION_BAND_PROP
  if (SMALLEST_POSITIVE_STEP >= medium_band_upper) {
    breaks <- c(0, SMALLEST_POSITIVE_STEP, scale_max)
    colors <- c(COLOR_WHITE, COLOR_MEDIUM_RED, COLOR_DARK_RED)
  } else {
    breaks <- c(0, SMALLEST_POSITIVE_STEP, medium_band_upper, scale_max)
    colors <- c(COLOR_WHITE, COLOR_MEDIUM_RED, COLOR_MEDIUM_RED, COLOR_DARK_RED)
  }
  list(breaks = breaks, colors = colors)
}

# Negative data color scale
negative_color_scale <- function(scale_min) {
  if (round(scale_min, 9) >= 0) return(default_color_scale())
  medium_band_lower <- scale_min * (1 - MEDIUM_SATURATION_BAND_PROP)
  if (-SMALLEST_POSITIVE_STEP <= medium_band_lower) {
    breaks <- c(scale_min, -SMALLEST_POSITIVE_STEP, 0)
    colors <- c(COLOR_DARK_BLUE, COLOR_MEDIUM_BLUE, COLOR_WHITE)
  } else {
    breaks <- c(scale_min, medium_band_lower, -SMALLEST_POSITIVE_STEP, 0)
    colors <- c(COLOR_DARK_BLUE, COLOR_MEDIUM_BLUE, COLOR_MEDIUM_BLUE, COLOR_WHITE)
  }
  list(breaks = breaks, colors = colors)
}

# Mixed data color scale
mixed_color_scale <- function(abs_limit) {
  if (!is.finite(abs_limit) || round(abs_limit, 9) == 0) abs_limit <- 1
  pos_medium_band_upper <- abs_limit * MEDIUM_SATURATION_BAND_PROP
  neg_medium_band_value <- -abs_limit * MEDIUM_SATURATION_BAND_PROP
  if (-SMALLEST_POSITIVE_STEP <= neg_medium_band_value) {
    breaks_neg_part <- c(-abs_limit, -SMALLEST_POSITIVE_STEP)
    colors_neg_part <- c(COLOR_DARK_BLUE, COLOR_MEDIUM_BLUE)
  } else {
    breaks_neg_part <- c(-abs_limit, neg_medium_band_value, -SMALLEST_POSITIVE_STEP)
    colors_neg_part <- c(COLOR_DARK_BLUE, COLOR_MEDIUM_BLUE, COLOR_MEDIUM_BLUE)
  }
  if (SMALLEST_POSITIVE_STEP >= pos_medium_band_upper) {
    breaks_pos_part <- c(SMALLEST_POSITIVE_STEP, abs_limit)
    colors_pos_part <- c(COLOR_MEDIUM_RED, COLOR_DARK_RED)
  } else {
    breaks_pos_part <- c(SMALLEST_POSITIVE_STEP, pos_medium_band_upper, abs_limit)
    colors_pos_part <- c(COLOR_MEDIUM_RED, COLOR_MEDIUM_RED, COLOR_DARK_RED)
  }
  list(breaks = c(breaks_neg_part, 0, breaks_pos_part),
       colors = c(colors_neg_part, COLOR_WHITE, colors_pos_part))
}

# Generate color scale based on data
create_color_scale <- function(data_matrix, fixed_max = NULL, title_for_warning = "heatmap") {
  min_data_val <- min(data_matrix, na.rm = TRUE)
  max_data_val <- max(data_matrix, na.rm = TRUE)
  raw_scale_list <- if (!is.finite(min_data_val) && !is.finite(max_data_val)) {
    default_color_scale()
  } else {
    if (!is.finite(min_data_val)) min_data_val <- if (max_data_val < 0) max_data_val else 0
    if (!is.finite(max_data_val)) max_data_val <- if (min_data_val > 0) min_data_val else 0
    if (round(min_data_val, 9) == 0 && round(max_data_val, 9) == 0) {
      default_color_scale()
    } else if (round(min_data_val, 9) >= 0) {
      scale_max <- if (!is.null(fixed_max) && is.finite(fixed_max)) fixed_max else max_data_val
      positive_color_scale(scale_max)
    } else if (round(max_data_val, 9) <= 0) {
      scale_min <- if (!is.null(fixed_max) && is.finite(fixed_max)) -abs(fixed_max) else min_data_val
      negative_color_scale(scale_min)
    } else {
      abs_limit <- if (!is.null(fixed_max) && is.finite(fixed_max)) abs(fixed_max) else max(abs(c(min_data_val, max_data_val)), na.rm = TRUE)
      mixed_color_scale(abs_limit)
    }
  }
  breaks <- raw_scale_list$breaks
  colors <- raw_scale_list$colors
  if (length(breaks) == 0 || length(breaks) != length(colors)) {
    warning(paste("Initial breaks/colors issue for", title_for_warning, ". Defaulting."))
    return(default_color_scale())
  }
  df <- data.frame(original_breaks = breaks, rounded_breaks = round(breaks, 9), colors = colors)
  df <- df[order(df$rounded_breaks, df$original_breaks), ]
  unique_indices <- !duplicated(df$rounded_breaks)
  df_unique <- df[unique_indices, ]
  final_breaks <- df_unique$original_breaks
  final_colors <- df_unique$colors
  if (length(final_breaks) < 2) {
    warning(paste("Insufficient unique breaks (<2) after sanitization for:", title_for_warning, ". Defaulting."))
    return(default_color_scale())
  }
  if (length(final_breaks) != length(final_colors)) {
    warning(paste("CRITICAL: Mismatch in final_breaks/final_colors length for:", title_for_warning, ". Defaulting."))
    return(default_color_scale())
  }
  list(breaks = final_breaks, colors = final_colors)
}

# Cut dendrogram into clusters
cut_dendrogram <- function(hc, h = NULL, min_clusters = 2, max_clusters = 6, return_factor = FALSE) {
  if (is.null(h)) {
    heights <- sort(hc$height, decreasing = TRUE)
    height_gaps <- if (length(heights) > 1) diff(heights) / heights[-length(heights)] else numeric(0)
    smooth_gaps <- if (length(height_gaps) >= 5) {
      stats::filter(height_gaps, rep(1/3, 3), sides = 2)
    } else {
      height_gaps
    }
    smooth_gaps[is.na(smooth_gaps)] <- height_gaps[is.na(smooth_gaps)]
    avg_gap <- mean(smooth_gaps[1:min(5, length(smooth_gaps))], na.rm = TRUE)
    sig_threshold <- 1.3 * avg_gap
    sig_gaps <- which(smooth_gaps > sig_threshold)
    h <- if (length(sig_gaps) > 0) heights[sig_gaps[1] + 1] else heights[max(1, ceiling(length(heights) * 0.25))]
  }
  clusters <- cutree(hc, h = h)
  n_clusters <- length(unique(clusters))
  if (n_clusters < min_clusters) clusters <- cutree(hc, k = min_clusters)
  if (n_clusters > max_clusters) clusters <- cutree(hc, k = max_clusters)
  if (return_factor) factor(clusters) else length(unique(clusters))
}

# Create and save heatmap
plot_matrix_heatmap <- function(data_matrix, title, row_hc, col_hc, row_clusters, col_clusters,
                                output_path = NULL, right_annotation = NULL, fixed_max = NULL) {
  if (!is.matrix(data_matrix)) data_matrix <- as.matrix(data_matrix)
  if (!is.numeric(data_matrix)) {
    warning(paste("Data matrix for '", title, "' is not numeric. Attempting to convert.", sep=""))
    mode(data_matrix) <- "numeric"
  }
  color_scale <- create_color_scale(data_matrix, fixed_max, title_for_warning = title)
  if (length(color_scale$breaks) != length(color_scale$colors) || length(color_scale$breaks) < 2) {
    warning(paste("FATAL: Invalid breaks/colors for heatmap '", title, "'. Using default.", sep=""))
    color_scale <- default_color_scale()
  }
  col_fun <- circlize::colorRamp2(color_scale$breaks, color_scale$colors)
  ht <- Heatmap(data_matrix, name = title, col = col_fun,
                cluster_rows = row_hc, cluster_columns = col_hc,
                row_split = row_clusters, column_split = col_clusters,
                right_annotation = right_annotation,
                show_row_names = TRUE, show_column_names = TRUE,
                column_names_rot = 45,
                row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                row_dend_width = unit(2, "cm"), column_dend_height = unit(2, "cm"),
                gap = unit(1.5, "mm"), border = TRUE, border_gp = gpar(lwd = 0.8))
  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    png(output_path, width = 12, height = 10, units = "in", res = 300, bg = "white")
    draw(ht, background = "white")
    dev.off()
  }
  ht
}
