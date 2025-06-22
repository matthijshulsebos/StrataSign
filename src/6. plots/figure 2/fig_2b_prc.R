# PRC curve plotting for all models (all_clusters, metabolic) for each normalization method
library(ggplot2)
library(dplyr)
library(readr)
library(PRROC)
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

plot_prc_curves_for_norm <- function(norm) {
  pred_files <- glob_pred_files(norm)
  if (length(pred_files) == 0) {
    message(sprintf("No prediction files found for all_clusters & metabolic in %s.", norm))
    return(NULL)
  }
  prc_list <- list()
  for (f in pred_files) {
    dat <- read_csv(f, show_col_types = FALSE)
    if (!all(c("y_test", "y_pred_prob") %in% names(dat))) {
      warning(sprintf("File %s does not have y_pred_prob, skipping.", f))
      next
    }
    model <- get_model_name(f)
    prc_obj <- tryCatch({
      PRROC::pr.curve(scores.class0 = dat$y_pred_prob[dat$y_test == 1],
                     scores.class1 = dat$y_pred_prob[dat$y_test == 0],
                     curve = TRUE)
    }, error=function(e) NULL)
    if (!is.null(prc_obj)) {
      prc_list[[model]] <- prc_obj
    }
  }
  if (length(prc_list) == 0) {
    message(sprintf("No valid PRC curves could be computed for %s.", norm))
    return(NULL)
  }
  # Assign the same color to models with identical PRC curves
  prc_data_list <- lapply(prc_list, function(prc) paste0(round(prc$curve[,1], 6), collapse = ",")
    |> paste0("|", paste0(round(prc$curve[,2], 6), collapse = ",")))
  group_ids <- match(prc_data_list, unique(prc_data_list))
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
  model_colors <- setNames(group_colors[group_ids], names(prc_list))
  auc_labels <- sapply(names(prc_list), function(name) {
    auc_val <- prc_list[[name]]$auc.integral
    sprintf("%s (AUC=%.3f)", name, auc_val)
  })
  # Build PRC data frame for plotting (no forced endpoints)
  prc_df <- do.call(rbind, lapply(names(prc_list), function(model) {
    data.frame(
      recall = prc_list[[model]]$curve[,1],
      precision = prc_list[[model]]$curve[,2],
      model = model
    )
  }))
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
  plt <- ggplot(prc_df, aes(x = recall, y = precision, color = model)) +
    geom_line(size = 0.7) +
    labs(
      x = "Recall", y = "Precision", color = "Model"
    ) +
    scale_color_manual(values = model_colors, labels = auc_labels) +
    theme_unified
  plt <- plt + coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.05))
  outdir <- file.path("output/6. plots/figure 2", norm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, "prc_auc_all_clusters_metabolic.png"), plt, width=8, height=6)
  return(plt)
}

if (interactive()) {
  for (norm in glob_norms) plot_prc_curves_for_norm(norm)
}
