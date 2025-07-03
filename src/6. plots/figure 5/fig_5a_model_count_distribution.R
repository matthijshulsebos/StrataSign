library(ggplot2)
library(dplyr)
library(readr)

# Helper to load n_models_occur for a normalization method
load_n_models_occur <- function(file_path, method_name) {
  df <- read_csv(file_path, show_col_types = FALSE)
  df %>% select(n_models_occur) %>% mutate(method = method_name)
}

# File paths for all normalization methods 
global_path <- "output/3. intersector/ctnorm_global/all_clusters/metabolic/meta_scores.csv"
relative_path <- "output/3. intersector/ctnorm_relative/all_clusters/metabolic/meta_scores.csv"
read_depth_path <- "output/3. intersector/read_depth/all_clusters/metabolic/meta_scores.csv"

# Load data
nmodels_global <- load_n_models_occur(global_path, "ctnorm_global")
nmodels_relative <- load_n_models_occur(relative_path, "ctnorm_relative")
nmodels_read_depth <- load_n_models_occur(read_depth_path, "read_depth")

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

# Bar plot for model occurence count values
p <- ggplot(nmodels_summary, aes(x = factor(n_models_occur), y = count, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  scale_fill_brewer(
    name = "Normalization",
    palette = "Set1",
    labels = c("global", "relative", "read depth")
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(face = "plain", size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.title = element_text(face = "plain", size = 15),
    legend.text = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.background = element_rect(color = NA, fill = NA),
    legend.box.background = element_blank()
  ) +
  labs(x = "Model prevalence", y = "Number of features")

ggsave("output/6. plots/figure 5/fig_5a_model_count_distribution.png", p, width = 7, height = 5)
