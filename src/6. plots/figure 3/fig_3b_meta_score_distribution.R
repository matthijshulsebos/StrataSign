library(ggplot2)
library(dplyr)
library(readr)

# Helper to load meta_score for a normalization method
get_meta_scores <- function(file_path, method_name) {
  df <- read_csv(file_path, show_col_types = FALSE)
  df %>% select(meta_score) %>% mutate(method = method_name)
}

# File paths for all three normalization methods (all_clusters/metabolic)
global_path <- "output/3. intersector/ctnorm_global/all_clusters/metabolic/meta_scores.csv"
relative_path <- "output/3. intersector/ctnorm_relative/all_clusters/metabolic/meta_scores.csv"
read_depth_path <- "output/3. intersector/read_depth/all_clusters/metabolic/meta_scores.csv"

# Gather data
meta_global <- get_meta_scores(global_path, "ctnorm_global")
meta_relative <- get_meta_scores(relative_path, "ctnorm_relative")
meta_read_depth <- get_meta_scores(read_depth_path, "read_depth")

# Combine for plotting
meta_all <- bind_rows(meta_global, meta_relative, meta_read_depth) %>%
  mutate(method = recode(method,
                        ctnorm_global = "global",
                        ctnorm_relative = "relative",
                        read_depth = "read depth"))

# Density plot for meta_score by normalization method
p <- ggplot(meta_all, aes(x = meta_score, fill = method)) +
  geom_density(alpha = 0.4, color = NA, size = 0.7, adjust = 1.1) +
  scale_fill_manual(
    name = "Normalization",
    values = c(
      "global" = "#56B4E9",      # blue
      "relative" = "#E69F00",    # orange
      "read depth" = "#009E73"   # green
    )
  ) +
  guides(fill = guide_legend(override.aes = list(color = "black"))) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.background = element_rect(color = "black", fill = "white", size = 0.8, linetype = "solid")
  ) +
  labs(x = "Meta score", y = "Density")

ggsave("output/6. plots/figure 3/fig_3b_meta_score_distribution.png", p, width = 7, height = 5)
