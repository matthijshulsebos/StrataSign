# Probability histogram for all models (all_clusters, metabolic) for each normalization method
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

glob_norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")

glob_pred_files <- function(norm) {
  list.files(
    path = file.path("output/2. models", norm, "all_clusters", "metabolic"),
    pattern = "predictions_all_clusters_metabolic.csv$", 
    recursive = TRUE, full.names = TRUE
  )
}

get_model_name <- function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  idx <- which(parts == "metabolic")
  if (length(idx) > 0 && length(parts) > idx) {
    return(parts[idx+1])
  } else {
    return(basename(dirname(path)))
  }
}

theme_unified <- theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"),
    legend.box.background = element_rect(color = "black", size = 1)
  )

plot_probability_histograms_for_norm <- function(norm) {
  pred_files <- glob_pred_files(norm)
  if (length(pred_files) == 0) {
    message(sprintf("No prediction files found for all_clusters & metabolic in %s.", norm))
    return(NULL)
  }
  # Read all predictions into one data frame
  prob_df <- map_dfr(pred_files, function(f) {
    dat <- read_csv(f, show_col_types = FALSE)
    model <- get_model_name(f)
    if ("y_pred_prob" %in% names(dat)) {
      data.frame(model = model, y_pred_prob = dat$y_pred_prob)
    } else {
      warning(sprintf("File %s has no y_pred_prob column, skipping.", f))
      return(NULL)
    }
  })
  if (nrow(prob_df) == 0) {
    message(sprintf("No valid probability data found for %s.", norm))
    return(NULL)
  }
  # Assign colors (primary palette, repeat as needed)
  primary_colors <- c(
    "#E41A1C", # red
    "#377EB8", # blue
    "#4DAF4A", # green
    "#000000", # black
    "#FF7F00", # orange
    "#984EA3", # purple
    "#00CED1", # cyan
    "#FFD700"  # yellow
  )
  model_levels <- unique(prob_df$model)
  color_vals <- setNames(rep(primary_colors, length.out = length(model_levels)), model_levels)
  prob_df$model <- factor(prob_df$model, levels = model_levels)
  plt <- ggplot(prob_df, aes(x = model, y = y_pred_prob, fill = model, color = model)) +
    geom_violin(trim = FALSE, alpha = 0.5, width = 0.8) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = color_vals) +
    scale_color_manual(values = color_vals) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Model", y = "Predicted Probability", color = "Model", fill = "Model"
    ) +
    theme_unified
  outdir <- file.path("output/6. plots/figure 2", norm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, "probability_violin_all_clusters_metabolic.png"), plt, width=8, height=6)
  return(plt)
}

if (interactive()) {
  for (norm in glob_norms) plot_probability_histograms_for_norm(norm)
}
