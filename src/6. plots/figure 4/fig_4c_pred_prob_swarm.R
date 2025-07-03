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

  # Read all prediction files and extract y_pred_prob
  prob_df <- map_dfr(pred_files, function(f) {
    model <- basename(dirname(f))
    dat <- read_csv(f, show_col_types = FALSE)
    data.frame(model = model, y_pred_prob = dat$y_pred_prob)
  })

  # Map raw model names to clean names
  prob_df$clean_model <- vapply(prob_df$model, clean_model_name, character(1))

  primary_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#000000", "#FF7F00", "#984EA3", "#00CED1", "#FFD700"
  )

  # Get vars from the data
  model_levels <- sort(unique(prob_df$clean_model))
  color_vals <- setNames(rep(primary_colors, length.out = length(model_levels)), model_levels)
  prob_df$clean_model <- factor(prob_df$clean_model, levels = model_levels)

  plt <- ggplot(prob_df, aes(x = clean_model, y = y_pred_prob, color = clean_model, fill = clean_model)) +
    geom_beeswarm(size = 2.2, cex = 1.1, priority = "density", dodge.width = 0.7, alpha = 0.7) +
    scale_fill_manual(values = color_vals, breaks = levels(prob_df$clean_model), labels = levels(prob_df$clean_model)) +
    scale_color_manual(values = color_vals, breaks = levels(prob_df$clean_model), labels = levels(prob_df$clean_model)) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL, y = "Predicted probability"
    ) +
    swarm_plot_theme +
    theme(legend.position = "none")

  # Save all to s1 swarm directory
  s1_swarm_dir <- "output/6. plots/figure s1/swarm"
  dir.create(s1_swarm_dir, recursive = TRUE, showWarnings = FALSE)
  s1_swarm_file <- file.path(s1_swarm_dir, paste0("probability_swarm_", norm, "_", cell_type, "_", gene_type, ".png"))
  ggsave(s1_swarm_file, plt, width=8, height=6)

  # Save global normalization, all_clusters, metabolic to figure 4c
  if (norm == "ctnorm_global" && cell_type == "all_clusters" && gene_type == "metabolic") {
    fig4_dir <- "output/6. plots/figure 4"
    dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(fig4_dir, "fig_4c_probability_swarm.png"), plt, width=8, height=6)
  }
  return(plt)
}


# Theme for swarm plots
swarm_plot_theme <- theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5, size = 18),
    plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15),
    axis.title = element_text(face = "plain", size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_text(face = "plain", size = 15),
    legend.text = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.background = element_rect(color = NA, fill = NA),
    legend.box.background = element_blank()
  )

# Define combination components
norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters", "lcam_both", "lcam_lo", "lcam_hi", "macrophages")
gene_types <- c("metabolic", "nonmetabolic", "random")

# Loop over combinations
combos <- expand.grid(norm = norms, cell_type = cell_types, gene_type = gene_types, stringsAsFactors = FALSE)
for (i in seq_len(nrow(combos))) {
  norm <- combos$norm[i]
  cell_type <- combos$cell_type[i]
  gene_type <- combos$gene_type[i]
  pred_files <- list.files(
    path = file.path("output/2. models", norm, cell_type, gene_type),
    pattern = "predictions_.*.csv$",
    recursive = TRUE, full.names = TRUE
  )
  if (length(pred_files) > 0) {
    plot_probability_swarm(norm, cell_type, gene_type)
  }
}
