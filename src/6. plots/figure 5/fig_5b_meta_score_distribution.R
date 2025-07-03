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
    geom_density(fill = "#377EB8", alpha = 0.7, color = NA, size = 0.7, adjust = 1.1) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "plain", hjust = 0.5, size = 18),
      plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15),
      axis.title = element_text(face = "plain", size = 16),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
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

# Run plotting function
plot_and_save_density(meta_all %>% filter(method == "global"),
                     "global",
                     "output/6. plots/figure 5",
                     "fig_5b_meta_score_distribution.png")
plot_and_save_density(meta_all %>% filter(method == "relative"),
                     "relative",
                     "output/6. plots/figure s2",
                     "meta_score_distribution_relative.png")
plot_and_save_density(meta_all %>% filter(method == "read depth"),
                     "read depth",
                     "output/6. plots/figure s2",
                     "meta_score_distribution_read_depth.png")
