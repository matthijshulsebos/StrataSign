# Accuracy heatmap comparing all models and normalization methods (all_clusters, metabolic)
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

# Aggregate accuracy for all models and normalization methods
gather_accuracy <- function() {
  bind_rows(lapply(glob_norms, function(norm) {
    pred_files <- glob_pred_files(norm)
    if (length(pred_files) == 0) return(NULL)
    map_dfr(pred_files, function(f) {
      dat <- read_csv(f, show_col_types = FALSE)
      if (!all(c("y_test", "y_pred") %in% names(dat))) return(NULL)
      model <- get_model_name(f)
      acc <- mean((dat$y_pred >= 0.5) == dat$y_test)
      data.frame(model = model, norm = norm, accuracy = acc)
    })
  }))
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
    legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"),
    legend.box.background = element_rect(color = "black", size = 1)
  )

plot_accuracy_heatmap_all_norms <- function() {
  acc_df <- gather_accuracy()
  if (nrow(acc_df) == 0) {
    message("No valid accuracy data found for any normalization method.")
    return(NULL)
  }
  acc_df$norm <- factor(acc_df$norm, levels = glob_norms)
  plt <- ggplot(acc_df, aes(x = model, y = norm, fill = accuracy)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", accuracy)), color = "black", size = 6, fontface = "bold") +
    scale_fill_viridis_c(option = "C", name = "Accuracy") +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Model", y = "Normalization"
    ) +
    theme_unified
  outdir <- file.path("output/6. plots/figure 2")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, "accuracy_heatmap_all_clusters_metabolic.png"), plt, width=8, height=4)
  return(plt)
}

if (interactive()) {
  plot_accuracy_heatmap_all_norms()
}
