# ROC curve plotting for all models (all_clusters, metabolic) for each normalization method
library(ggplot2)
library(dplyr)
library(readr)
library(pROC)
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

plot_roc_curves_for_norm <- function(norm) {
  pred_files <- glob_pred_files(norm)
  if (length(pred_files) == 0) {
    message(sprintf("No prediction files found for all_clusters & metabolic in %s.", norm))
    return(NULL)
  }
  roc_list <- list()
  for (f in pred_files) {
    dat <- read_csv(f, show_col_types = FALSE)
    if (!all(c("y_test", "y_pred_prob") %in% names(dat))) {
      warning(sprintf("File %s does not have y_pred_prob, skipping.", f))
      next
    }
    model <- get_model_name(f)
    roc_obj <- tryCatch({
      pROC::roc(dat$y_test, dat$y_pred_prob, quiet=TRUE)
    }, error=function(e) NULL)
    if (!is.null(roc_obj)) {
      roc_list[[model]] <- roc_obj
    }
  }
  if (length(roc_list) == 0) {
    message(sprintf("No valid ROC curves could be computed for %s.", norm))
    return(NULL)
  }
  # Assign the same color to models with identical ROC curves
  roc_data_list <- lapply(roc_list, function(roc) paste0(round(roc$sensitivities, 6), collapse = ",")
    |> paste0("|", paste0(round(roc$specificities, 6), collapse = ",")))
  group_ids <- match(roc_data_list, unique(roc_data_list))
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
  group_colors <- rep(primary_colors, length.out = length(unique(group_ids)))
  model_colors <- setNames(group_colors[group_ids], names(roc_list))
  auc_labels <- sapply(names(roc_list), function(name) {
    auc_val <- pROC::auc(roc_list[[name]])
    sprintf("%s (AUC=%.3f)", name, auc_val)
  })
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
  plt <- ggroc(roc_list, legacy.axes=TRUE) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "1 - Specificity", y = "Sensitivity", color = "Model"
    ) +
    scale_color_manual(values = model_colors, labels = auc_labels) +
    theme_unified
  plt <- plt + coord_cartesian(xlim = c(0, 1.05))
  outdir <- file.path("output/6. plots/figure 2", norm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, "roc_auc_all_clusters_metabolic.png"), plt, width=8, height=6)
  return(plt)
}

if (interactive()) {
  for (norm in glob_norms) plot_roc_curves_for_norm(norm)
}
