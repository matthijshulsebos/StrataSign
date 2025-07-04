library(ggplot2)
library(dplyr)
library(readr)
library(treemapify)

# Load sublineage color mapping as in the barplot script
sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
sublineage_colors <- readRDS(sublineage_color_map_path)

# Define all combinations to plot
combinations <- expand.grid(
  normalization_method = c("ctnorm_global"),
  meta_celltype_set = c("all_clusters"),
  gene_type = c("metabolic"),
  stringsAsFactors = FALSE
)

# Add a column for output type (figure5 or s2)
combinations$output_type <- ifelse(
  combinations$normalization_method == "ctnorm_global" &
    combinations$meta_celltype_set == "all_clusters" &
    combinations$gene_type == "metabolic",
  "figure5", "s2"
)

for (i in seq_len(nrow(combinations))) {
  combo <- combinations[i, ]
  meta_scores_file <- file.path(
    "output/3. intersector",
    combo$normalization_method,
    combo$meta_celltype_set,
    combo$gene_type,
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
    geom_treemap_text(colour = "black", place = "centre", grow = TRUE, reflow = TRUE, min.size = 3) +
    fill_scale +
    theme(legend.position = "none", plot.title = element_blank(), plot.subtitle = element_blank())

  if (combo$output_type == "figure5") {
    out_base <- "output/6. plots/figure 5/"
    if (!dir.exists(out_base)) dir.create(out_base, recursive = TRUE)
    out_name <- paste0(out_base, "fig_5d_cell_type_contribution_treemap.png")
    ggsave(out_name, p, width = 7, height = 5, dpi = 300)
    message(sprintf("Wrote fig 5d to %s", out_name))
  } else {
    out_base <- "output/6. plots/figure s2/contribution"
    if (!dir.exists(out_base)) dir.create(out_base, recursive = TRUE)
    out_name <- paste0(out_base, "/celltype_contribution_",
      combo$normalization_method, "_",
      combo$meta_celltype_set, "_",
      combo$gene_type, ".png")
    ggsave(out_name, p, width = 7, height = 5, dpi = 300)
  }
}
