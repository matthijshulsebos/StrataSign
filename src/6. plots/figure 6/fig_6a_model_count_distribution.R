library(ggplot2)
library(dplyr)
library(readr)


# Define all combinations
norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters")
gene_types <- c("metabolic")

# Helper to load n_models_occur for a normalization method
load_n_models_occur <- function(file_path, method_name, cell_type, gene_type) {
  if (!file.exists(file_path)) return(NULL)
  df <- read_csv(file_path, show_col_types = FALSE)
  df %>% select(n_models_occur) %>% mutate(method = method_name, cell_type = cell_type, gene_type = gene_type)
}

# Output root for figure s6
s6_root <- file.path("output/6. plots/figure s6/meta_prevalence")
dir.create(s6_root, recursive = TRUE, showWarnings = FALSE)

# Loop over all combinations
for (norm in norms) {
  for (cell_type in cell_types) {
    for (gene_type in gene_types) {
      # File path for this combination
      file_path <- file.path("output/3. intersector", norm, cell_type, gene_type, "meta_scores.csv")
      nmodels_df <- load_n_models_occur(file_path, norm, cell_type, gene_type)
      if (is.null(nmodels_df)) next

      # Summarize counts for each n_models_occur (single normalization method)
      nmodels_summary <- nmodels_df %>%
        group_by(n_models_occur) %>%
        summarise(count = n(), .groups = 'drop')

      # Bar plot for model occurrence count values (single normalization method)
      p <- ggplot(nmodels_summary, aes(x = factor(n_models_occur), y = count)) +
        geom_bar(stat = "identity", fill = "#377EB8", color = "black", width = 0.7) +
        theme_bw(base_size = 16, base_family = "sans") +
        theme(
          axis.title = element_text(face = "plain", size = 16, color = "black"),
          axis.text = element_text(size = 14, color = "black"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          legend.position = "none",
          panel.grid = element_blank(),
          legend.background = element_rect(color = NA, fill = NA),
          legend.box.background = element_blank()
        ) +
        labs(x = "Model prevalence", y = "Number of features")

      # Save to figure s6/meta_prevalence/<norm>/<cell_type>/<gene_type>/model_count_distribution.png
      safe_norm <- gsub("[^a-zA-Z0-9_]+", "_", norm)
      safe_celltype <- gsub("[^a-zA-Z0-9_]+", "_", cell_type)
      safe_genetype <- gsub("[^a-zA-Z0-9_]+", "_", gene_type)
      combo_dir <- file.path(s6_root, safe_norm, safe_celltype, safe_genetype)
      dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
      s6_file <- file.path(combo_dir, "model_count_distribution.png")
      ggsave(s6_file, p, width = 7, height = 5)
    }
  }
}
