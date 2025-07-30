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


# Retrieves model name from the file path
get_model_name <- function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  idx <- which(parts == "metabolic")
  if (length(idx) > 0 && length(parts) > idx) {
    return(parts[idx+1])
  } else {
    return(basename(dirname(path)))
  }
}


# Clean model names for plotting
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


# ROC for each model
plot_roc_curves_for_combo <- function(norm, celltype, genetype) {
  pred_files <- glob_pred_files(norm, celltype, genetype)
  if (length(pred_files) == 0) {
    message(sprintf("No prediction files found for %s, %s, %s.", norm, celltype, genetype))
    return(NULL)
  }

  for (f in pred_files) {
    dat <- read_csv(f, show_col_types = FALSE)
    model_raw <- get_model_name(f)
    model <- clean_model_name(model_raw)
    roc_obj <- pROC::roc(dat$y_test, dat$y_pred_prob, quiet=TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    # Prepare output path
    roc_root <- file.path("output/6. plots/figure 6/roc")
    safe_norm <- gsub("[^a-zA-Z0-9_]+", "_", norm)
    safe_celltype <- gsub("[^a-zA-Z0-9_]+", "_", celltype)
    safe_genetype <- gsub("[^a-zA-Z0-9_]+", "_", genetype)
    safe_model <- gsub("[^a-zA-Z0-9_]+", "_", model)
    combo_dir <- file.path(roc_root, safe_norm, safe_celltype, safe_genetype)
    dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
    roc_file <- file.path(combo_dir, sprintf("roc_%s.png", safe_model))

    # Create the plot with ggroc
    p <- ggroc(roc_obj, color = "#1B9E77", size = 1.5) +
      ggtitle(sprintf("Model: %s (AUC: %.2f)", model, auc_val)) +
      xlab("False positive rate") +
      ylab("True positive rate") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 40, hjust = 0, face = "plain"),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        axis.text.x = element_text(size = 34),
        axis.text.y = element_text(size = 34),
        text = element_text(family = "sans"),

        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),

        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )

    # Save the plot
    ggsave(
      roc_file,
      plot = p,
      width = 12,
      height = 12,
      units = "in",
      dpi = 300,
      bg = "transparent"
    )
  }
  return(NULL)
}


# Loop over all combinations
for (norm in glob_norms) {
  for (celltype in celltype_sets) {
    for (genetype in gene_types) {
      plot_roc_curves_for_combo(norm, celltype, genetype)
    }
  }
}
