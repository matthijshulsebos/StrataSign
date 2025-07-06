library(ggplot2)
library(dplyr)
library(readr)


# File paths for normalization methods
global_path <- "output/3. intersector/ctnorm_global/all_clusters/metabolic/meta_scores.csv"
relative_path <- "output/3. intersector/ctnorm_relative/all_clusters/metabolic/meta_scores.csv"
read_depth_path <- "output/3. intersector/read_depth/all_clusters/metabolic/meta_scores.csv"

# Load data
meta_global <- read_csv(global_path, show_col_types = FALSE) %>% select(meta_score) %>% mutate(method = "ctnorm_global")
meta_relative <- read_csv(relative_path, show_col_types = FALSE) %>% select(meta_score) %>% mutate(method = "ctnorm_relative")
meta_read_depth <- read_csv(read_depth_path, show_col_types = FALSE) %>% select(meta_score) %>% mutate(method = "read_depth")

# Combine for plotting
meta_all <- bind_rows(meta_global, meta_relative, meta_read_depth) %>%
  mutate(method = recode(method,
                        ctnorm_global = "global",
                        ctnorm_relative = "relative",
                        read_depth = "read depth"))


# Plot and save each density separately
plot_and_save_density <- function(df, method_label, outdir, fname) {
  p <- ggplot(df, aes(x = meta_score)) +
    geom_density(fill = "#377EB8", alpha = 1, color = "black", size = 0.6, adjust = 1.1) +
    theme_bw(base_size = 16, base_family = "sans") +
    theme(
      plot.title = element_text(face = "plain", hjust = 0.5, size = 18, color = "black"),
      plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15, color = "black"),
      axis.title = element_text(face = "plain", size = 16, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, color = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      strip.text = element_text(color = "black", face = "bold")
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Meta score",
      y = "Density"
    )

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, fname), p, width = 7, height = 5)
}

# Define all combinations
norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters")
gene_types <- c("metabolic")

# Output root for figure s6
s6_root <- file.path("output/6. plots/figure s6/meta_score_dist")
dir.create(s6_root, recursive = TRUE, showWarnings = FALSE)

# Helper to load meta_scores for a normalization method
load_meta_scores <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  read_csv(file_path, show_col_types = FALSE) %>% select(meta_score)
}

# Iterate over all combinations and save density plots
for (norm in norms) {
  for (cell_type in cell_types) {
    for (gene_type in gene_types) {
      file_path <- file.path("output/3. intersector", norm, cell_type, gene_type, "meta_scores.csv")
      meta_df <- load_meta_scores(file_path)
      if (is.null(meta_df)) next
      method_label <- norm
      # Save to figure s6/meta_score_dist/<norm>/<cell_type>/<gene_type>/meta_score_distribution.png
      safe_norm <- gsub("[^a-zA-Z0-9_]+", "_", norm)
      safe_celltype <- gsub("[^a-zA-Z0-9_]+", "_", cell_type)
      safe_genetype <- gsub("[^a-zA-Z0-9_]+", "_", gene_type)
      combo_dir <- file.path(s6_root, safe_norm, safe_celltype, safe_genetype)
      dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
      fname <- "meta_score_distribution.png"
      plot_and_save_density(meta_df, method_label, combo_dir, fname)
    }
  }
}
