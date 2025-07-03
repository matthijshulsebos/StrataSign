library(ggplot2)
library(dplyr)
library(readr)
library(pROC)
library(stringr)
library(purrr)


# Define all combinations
glob_norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")
celltype_sets <- c("all_clusters", "macrophages")
gene_types <- c("metabolic", "nonmetabolic", "random")


# Generalized prediction file finder
glob_pred_files <- function(norm, celltype, genetype) {
  list.files(
    path = file.path("output/2. models", norm, celltype, genetype),
    pattern = paste0("predictions_", celltype, "_", genetype, ".csv$"),
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


# Clean model names for plotting (shared with other scripts)
clean_model_name <- function(x) {
  x_low <- tolower(x)
  if (x_low == "lasso") return("LASSO")
  if (x_low == "lightgbm") return("LightGBM")
  if (x_low == "xgboost") return("XGBoost")
  if (x_low == "randomforest") return("Random Forest")
  if (x_low == "elasticnet") return("Elastic Net")
  if (x_low == "spls") return("sPLS")
  if (x_low == "svmrbf") return("SVM (RBF)")
  if (x_low == "svmlinear") return("SVM (Linear)")
  return(x)
}


# ROC plotting func
plot_roc_curves_for_combo <- function(norm, celltype, genetype) {
  pred_files <- glob_pred_files(norm, celltype, genetype)
  if (length(pred_files) == 0) {
    message(sprintf("No prediction files found for %s, %s, %s.", norm, celltype, genetype))
    return(NULL)
  }

  roc_list <- list()
  for (f in pred_files) {
    dat <- read_csv(f, show_col_types = FALSE)
    model_raw <- get_model_name(f)
    model <- clean_model_name(model_raw)
    roc_list[[model]] <- pROC::roc(dat$y_test, dat$y_pred_prob, quiet=TRUE)
  }

  # Generate unique identifier for ROC curves
  auc_vals <- sapply(roc_list, function(roc) as.numeric(pROC::auc(roc)))
  # If auc is 1 then all curves are the same
  is_perfect <- abs(auc_vals - 1) < 1e-6
  # Create unique id based on sensitivities and specificities
  roc_data_list <- lapply(seq_along(roc_list), function(i) {
    if (is_perfect[i]) {
      "PERFECT"
    } else {
      roc <- roc_list[[i]]
      paste0(
        paste0(round(roc$sensitivities, 6), collapse = ","),
        "|",
        paste0(round(roc$specificities, 6), collapse = ",")
      )
    }
  })

  # Create a mapping of unique identifiers to colors
  group_ids <- match(roc_data_list, unique(roc_data_list))
  primary_colors <- c(
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#000000",
    "#FF7F00",
    "#984EA3",
    "#00CED1",
    "#FFD700"
  )
  group_colors <- rep(primary_colors, length.out = length(unique(group_ids)))
  model_colors <- setNames(group_colors[group_ids], names(roc_list))
  auc_labels <- sapply(names(roc_list), function(name) {
    auc_val <- pROC::auc(roc_list[[name]])
    sprintf("%s (AUC=%.2f)", name, auc_val)
  })

  theme_unified <- theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "plain", hjust = 0.5, size = 18),
      plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15),
      axis.title = element_text(face = "plain", size = 16),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title = element_text(face = "plain", size = 15),
      legend.text = element_text(size = 13),
      panel.grid = element_blank(),
      legend.background = element_rect(color = NA, fill = NA),
      legend.box.background = element_blank()
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

  # Save all to figure s1/roc
  outdir <- file.path("output/6. plots/figure s1", "roc")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0("roc_", norm, "_", celltype, "_", genetype, ".png")
  ggsave(file.path(outdir, filename), plt, width=8, height=6)

  # Copy the figure if this is the global normalization, all_clusters, metabolic
  if (norm == "ctnorm_global" && celltype == "all_clusters" && genetype == "metabolic") {
    outdir_main <- "output/6. plots/figure 4"
    dir.create(outdir_main, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(outdir_main, "fig_4a_roc_all_metabolic.png"), plt, width=8, height=6)

    message("Completed writing figure 4a to file.")
  }
  return(plt)
}


# Loop over all combinations
for (norm in glob_norms) {
  for (celltype in celltype_sets) {
    for (genetype in gene_types) {
      plot_roc_curves_for_combo(norm, celltype, genetype)
    }
  }
}
