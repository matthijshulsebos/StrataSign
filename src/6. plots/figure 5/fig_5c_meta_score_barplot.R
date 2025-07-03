library(ggplot2)
library(dplyr)
library(readr)

# Define combinations
norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both")
gene_sets <- c("metabolic", "nonmetabolic", "random")


# Load sublineage color mapping
sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
sublineage_colors <- readRDS(sublineage_color_map_path)

for (norm in norms) {
  for (cell_type in cell_types) {
    for (gene_set in gene_sets) {
      meta_scores_file <- file.path("output/3. intersector", norm, cell_type, gene_set, "meta_scores.csv")
      if (!file.exists(meta_scores_file)) next
      all_results <- read_csv(meta_scores_file, show_col_types = FALSE)

      # Select top features by meta_score
      top_features_data <- all_results %>%
        arrange(desc(meta_score)) %>%
        head(15)

      # Use sublineage color mapping
      sublineage_colors_plot <- sublineage_colors
      top_features_data <- top_features_data %>%
        mutate(
          display_name = make.unique(as.character(feature_id), sep = "..v"),
          sublineage = as.character(sublineage),
          n_models_occur_f = as.factor(n_models_occur)
        )

      p <- ggplot(top_features_data, aes(
          x = reorder(display_name, meta_score), 
          y = meta_score, 
          fill = sublineage
        )) +
        geom_col(color = "black", width = 0.7) +
        geom_text(aes(label = n_models_occur, y = meta_score + max(meta_score) * 0.02), 
                  color = "black", fontface = "plain", size = 3.5, hjust = 0) +
        scale_fill_manual(values = sublineage_colors_plot, name = "Sublineage", na.value = "grey60") +
        guides(
          fill = guide_legend(title = "Sublineage")
        ) +
        labs(x = "Feature", y = "Meta score") +
        coord_flip() +
        theme_bw(base_size = 16) +
        theme(
          plot.title = element_text(face = "plain", hjust = 0.5, size = 18),
          plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15),
          axis.title = element_text(face = "plain", size = 16),
          axis.text = element_text(size = 13),
          axis.text.y = element_text(size = 11, hjust = 1),
          legend.position = "right",
          legend.background = element_rect(fill = NA, color = NA),
          legend.title = element_text(size = 13, face = "plain"),
          legend.text = element_text(size = 12),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = margin(10, 15, 10, 10, "pt")
        )
        
      # Save main figure 5c for global/all_clusters/metabolic
      if (norm == "ctnorm_global" && cell_type == "all_clusters" && gene_set == "metabolic") {
        outdir <- "output/6. plots/figure 5"
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        ggsave(file.path(outdir, "fig_5c_meta_score_barplot.png"), plot = p, width = 8, height = 6, dpi = 300, units = "in", bg = "white", limitsize = FALSE)
      } else {
        outdir <- "output/6. plots/figure s2/barplot"
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        fname <- sprintf("barplot_%s_%s_%s.png", norm, cell_type, gene_set)
        ggsave(file.path(outdir, fname), plot = p, width = 8, height = 6, dpi = 300, units = "in", bg = "white", limitsize = FALSE)
      }
    }
  }
}
