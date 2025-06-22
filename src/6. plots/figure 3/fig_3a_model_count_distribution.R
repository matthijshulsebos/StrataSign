library(ggplot2)
library(dplyr)
library(readr)

# Helper to load n_models_occur for a normalization method
gather_n_models_occur <- function(file_path, method_name) {
  df <- read_csv(file_path, show_col_types = FALSE)
  df %>% select(n_models_occur) %>% mutate(method = method_name)
}

# File paths for all three normalization methods (all_clusters/metabolic)
global_path <- "output/3. intersector/ctnorm_global/all_clusters/metabolic/meta_scores.csv"
relative_path <- "output/3. intersector/ctnorm_relative/all_clusters/metabolic/meta_scores.csv"
read_depth_path <- "output/3. intersector/read_depth/all_clusters/metabolic/meta_scores.csv"

# Gather data
nmodels_global <- gather_n_models_occur(global_path, "ctnorm_global")
nmodels_relative <- gather_n_models_occur(relative_path, "ctnorm_relative")
nmodels_read_depth <- gather_n_models_occur(read_depth_path, "read_depth")

# Combine for plotting
nmodels_all <- bind_rows(nmodels_global, nmodels_relative, nmodels_read_depth) %>%
  mutate(method = recode(method,
                        ctnorm_global = "global",
                        ctnorm_relative = "relative",
                        read_depth = "read depth"))

# Summarize counts for each n_models_occur and method
nmodels_summary <- nmodels_all %>%
  group_by(method, n_models_occur) %>%
  summarise(count = n(), .groups = 'drop')

# Bar plot (dodged) for discrete n_models_occur values
p <- ggplot(nmodels_summary, aes(x = factor(n_models_occur), y = count, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  scale_fill_manual(
    name = "Normalization",
    values = c(
      "global" = "#56B4E9",
      "relative" = "#E69F00",
      "read depth" = "#009E73"
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.background = element_rect(color = "black", fill = "white", size = 0.8, linetype = "solid")
  ) +
  labs(x = "Number of models feature occurs in", y = "Number of features")

ggsave("output/6. plots/figure 3/fig_3a_model_count_distribution.png", p, width = 7, height = 5)
