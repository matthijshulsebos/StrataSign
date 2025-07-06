library(ggplot2)
library(dplyr)
library(readr)
library(pROC)
library(stringr)
library(purrr)


# Define all combinations
glob_norms <- c("ctnorm_global", "ctnorm_global_zscaled", "ctnorm_relative", "read_depth")
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

  # Use RColorBrewer Set1 palette for all model colors
  set1_colors <- RColorBrewer::brewer.pal(9, "Set1")
  n_models <- length(roc_list)
  model_colors <- setNames(rep(set1_colors, length.out = n_models), names(roc_list))
  auc_labels <- sapply(names(roc_list), function(name) {
    auc_val <- pROC::auc(roc_list[[name]])
    sprintf("%s (AUC=%.2f)", name, auc_val)
  })

  theme_unified <- theme_bw(base_size = 16, base_family = "sans") +
    theme(
      plot.title = element_text(face = "plain", hjust = 0.5, size = 18, color = "black"),
      plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15, color = "black"),
      axis.title = element_text(face = "plain", size = 16, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      legend.title = element_text(face = "plain", size = 15, color = "black"),
      legend.text = element_text(size = 13, color = "black"),
      # Enable grid lines for publication clarity
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      legend.background = element_rect(color = NA, fill = NA),
      legend.box.background = element_blank()
    )
  plt <- ggroc(roc_list, legacy.axes=TRUE, size = 1) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "1 - Specificity", y = "Sensitivity", color = "Model"
    ) +
    scale_color_manual(values = model_colors, labels = auc_labels) +
    theme_unified
  plt <- plt + coord_fixed(ratio = 1, xlim = c(0, 1.05), ylim = c(0, 1.05))
  # Move legend to the right outside the plot area and use a wide output to maximize the square plot panel
  plt <- plt + theme(
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.box.margin = margin(0, 0, 0, 20),
    plot.margin = margin(5, 5, 5, 5)
  )

  # Save all to figure s5/roc/<norm>/<celltype>/<genetype>/roc.png
  s5_roc_root <- file.path("output/6. plots/figure s5/roc")
  safe_norm <- gsub("[^a-zA-Z0-9_]+", "_", norm)
  safe_celltype <- gsub("[^a-zA-Z0-9_]+", "_", celltype)
  safe_genetype <- gsub("[^a-zA-Z0-9_]+", "_", genetype)
  combo_dir <- file.path(s5_roc_root, safe_norm, safe_celltype, safe_genetype)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
  s5_roc_file <- file.path(combo_dir, "roc.png")
  ggsave(s5_roc_file, plt, width=10, height=6)

  # Do not save a copy to figure 4 or print a message for figure 4
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
