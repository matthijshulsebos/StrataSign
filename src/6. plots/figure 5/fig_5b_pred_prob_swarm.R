library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggbeeswarm)


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

# Main plotting function
plot_probability_swarm <- function(norm, cell_type, gene_type) {
  pred_files <- list.files(
    path = file.path("output/2. models", norm, cell_type, gene_type),
    pattern = "predictions_.*.csv$",
    recursive = TRUE, full.names = TRUE
  )

  # Check if there are no files
  if (length(pred_files) == 0) {
    message(sprintf("No prediction files found for %s/%s/%s.", norm, cell_type, gene_type))
    return(NULL)
  }

  # Read all prediction files and extract y_pred_prob, color by y_test (0=Normal, 1=Tumor)
  prob_df <- map_dfr(pred_files, function(f) {
    model <- basename(dirname(f))
    dat <- read_csv(f, show_col_types = FALSE)
    # y_test: 0=Normal, 1=Tumor
    if ("y_test" %in% names(dat)) {
      tissue <- ifelse(dat$y_test == 0, "Normal",
                       ifelse(dat$y_test == 1, "Tumor", as.character(dat$y_test)))
    } else {
      tissue <- rep(NA_character_, nrow(dat))
    }
    data.frame(model = model, y_pred_prob = dat$y_pred_prob, tissue = tissue)
  })


  # Map raw model names to clean names
  prob_df$clean_model <- vapply(prob_df$model, clean_model_name, character(1))
  model_levels <- sort(unique(prob_df$clean_model))
  prob_df$clean_model <- factor(prob_df$clean_model, levels = model_levels)

  # Set tissue color mapping
  tissue_colors <- c("Normal" = "#377EB8", "Tumor" = "#E41A1C")
  prob_df$tissue <- ifelse(is.na(prob_df$tissue), "Unknown", prob_df$tissue)

  plt <- ggplot(prob_df, aes(x = clean_model, y = y_pred_prob, color = tissue, fill = tissue)) +
    geom_beeswarm(size = 2.2, cex = 1.1, priority = "density", dodge.width = 0.3, alpha = 0.7) +
    scale_fill_manual(values = tissue_colors, breaks = names(tissue_colors), labels = names(tissue_colors)) +
    scale_color_manual(values = tissue_colors, breaks = names(tissue_colors), labels = names(tissue_colors)) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL, y = "Predicted probability",
      color = "Tissue",
      fill = "Tissue"
    ) +
    swarm_plot_theme +
    theme(legend.position = "right") +
    scale_x_discrete(expand = expansion(add = 0.5))

  # Save all to s5 pred_prob_swarm directory in subdirectories by norm/cell_type/gene_type
  s5_swarm_root <- "output/6. plots/figure s5/pred_prob_swarm"
  safe_norm <- gsub("[^a-zA-Z0-9_]+", "_", norm)
  safe_celltype <- gsub("[^a-zA-Z0-9_]+", "_", cell_type)
  safe_genetype <- gsub("[^a-zA-Z0-9_]+", "_", gene_type)
  combo_dir <- file.path(s5_swarm_root, safe_norm, safe_celltype, safe_genetype)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
  s5_swarm_file <- file.path(combo_dir, "probability_swarm.png")
  ggsave(s5_swarm_file, plt, width=8, height=6)

  return(plt)
}


# Theme for swarm plots
swarm_plot_theme <- theme_bw(base_size = 16, base_family = "sans") +
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5, size = 18, color = "black"),
    plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15, color = "black"),
    axis.title = element_text(face = "plain", size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    legend.title = element_text(face = "plain", size = 15, color = "black"),
    legend.text = element_text(size = 13, color = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.background = element_rect(color = NA, fill = NA),
    legend.box.background = element_blank(),
    strip.text = element_text(color = "black", face = "bold")
  )

# Define combination components
norms <- c("ctnorm_global", "ctnorm_global_zscaled", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters", "lcam_both", "lcam_lo", "lcam_hi", "macrophages")
gene_types <- c("metabolic", "nonmetabolic", "random")

# Loop over combinations using nested for loops
for (norm in norms) {
  for (cell_type in cell_types) {
    for (gene_type in gene_types) {
      pred_files <- list.files(
        path = file.path("output/2. models", norm, cell_type, gene_type),
        pattern = "predictions_.*.csv$",
        recursive = TRUE, full.names = TRUE
      )
      if (length(pred_files) > 0) {
        plot_probability_swarm(norm, cell_type, gene_type)
      }
    }
  }
}
