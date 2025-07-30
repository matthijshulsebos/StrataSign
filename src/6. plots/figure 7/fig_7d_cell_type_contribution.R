library(ggplot2)
library(dplyr)
library(readr)
library(treemapify)

# Load sublineage color mapping as in the barplot script
sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
sublineage_colors <- readRDS(sublineage_color_map_path)

# Define all combinations to plot
combinations <- expand.grid(
  normalization_method = c("ctnorm_global", "ctnorm_relative", "read_depth"),
  meta_celltype_set = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
  gene_type = c("metabolic", "nonmetabolic", "random"),
  stringsAsFactors = FALSE
)

# Add a column for output type
combinations$output_type <- ifelse(
  combinations$normalization_method == "ctnorm_global" &
    combinations$meta_celltype_set == "all_clusters" &
    combinations$gene_type == "metabolic",
  "figure5", "s2"
)

for (norm in c("ctnorm_global", "ctnorm_relative", "read_depth")) {
  for (celltype in c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both")) {
    for (gene_set in c("metabolic", "nonmetabolic", "random")) {
      meta_scores_file <- file.path(
        "output/3. intersector",
        norm,
        celltype,
        gene_set,
        "meta_scores.csv"
      )
      if (!file.exists(meta_scores_file)) {
        message(sprintf("meta_scores.csv not found at %s", meta_scores_file))
        next
      }
      meta_scores <- suppressMessages(read_csv(meta_scores_file, show_col_types = FALSE))
      celltype_col <- if ("sublineage" %in% names(meta_scores)) "sublineage" else if ("cell_type" %in% names(meta_scores)) "cell_type" else NULL
      if (is.null(celltype_col)) {
        message(sprintf("No cell type column found in %s", meta_scores_file))
        next
      }
      contrib <- meta_scores %>%
        count(celltype = !!sym(celltype_col), name = "n_features") %>%
        mutate(proportion = n_features / sum(n_features))

      # Use sublineage color mapping if available
      if (!is.null(sublineage_colors)) {
        fill_scale <- scale_fill_manual(values = sublineage_colors, name = NULL, na.value = "grey80")
      } else {
        fill_scale <- scale_fill_brewer(palette = "Set1", name = NULL)
      }

      p <- ggplot(contrib, aes(area = n_features, fill = celltype, label = celltype)) +
        geom_treemap(color = "white") +
        geom_treemap_text(colour = "black", place = "centre", grow = TRUE, reflow = TRUE, min.size = 3, family = "sans") +
        fill_scale +
        theme_bw(base_family = "sans") +
        theme(
          legend.position = "none", 
          plot.title = element_blank(), 
          plot.subtitle = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA)
        )

      # Save to figure 7
      type_treemap_root <- file.path("output/6. plots/figure 7/type_treemap")
      combo_dir <- file.path(type_treemap_root, norm, celltype, gene_set)
      if (!dir.exists(combo_dir)) dir.create(combo_dir, recursive = TRUE)
      out_name <- file.path(combo_dir, sprintf("celltype_contribution_%s_%s_%s.png", norm, celltype, gene_set))
      ggsave(out_name, p, width = 7, height = 5, dpi = 300, bg = "transparent")
    }
  }
}
