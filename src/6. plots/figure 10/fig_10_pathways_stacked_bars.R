library(tidyverse)
library(fs)
library(RColorBrewer)

# Load sublineage color mapping
sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
sublineage_colors <- if (file.exists(sublineage_color_map_path)) readRDS(sublineage_color_map_path) else NULL

# Directory containing overview pathway mapping files
overview_dir <- "output/5. pathway mapping"
plot_output_dir <- "output/6. plots/figure 10"
dir_create(plot_output_dir, recurse = TRUE)

# Find all overview pathway mapping files
overview_files <- dir_ls(overview_dir, recurse = TRUE, regexp = "kegg_pathway_mapping_.*\\.csv$")

for (csv_file in overview_files) {
  df <- read_csv(csv_file, show_col_types = FALSE)
  gene_cols <- grep("^Genes_", names(df), value = TRUE)
  if (length(gene_cols) == 0) next
  # Reshape to long format
  df_long <- df %>%
    select(PathwayID, PathwayName, all_of(gene_cols)) %>%
    pivot_longer(
      cols = all_of(gene_cols),
      names_to = "CellType",
      values_to = "GeneCount"
    ) %>%
    mutate(
      sublineage = sub("^Genes_", "", CellType)
    ) %>%
    filter(GeneCount > 0)

  # Order pathways by total gene count and keep only top 10
  pathway_order <- df_long %>% group_by(PathwayID, PathwayName) %>% summarise(Total = sum(GeneCount), .groups = "drop") %>% arrange(desc(Total))
  top_pathways <- head(pathway_order$PathwayName, 10)
  pathway_order <- pathway_order %>% filter(PathwayName %in% top_pathways)
  df_long <- df_long %>% filter(PathwayName %in% top_pathways)
  df_long$PathwayName <- factor(df_long$PathwayName, levels = pathway_order$PathwayName)

  p <- ggplot(df_long, aes(x = PathwayName, y = GeneCount, fill = sublineage)) +
    geom_bar(stat = "identity", color = NA, width = 0.6) +
    geom_text(
      data = pathway_order,
      aes(x = PathwayName, y = Total, label = Total),
      inherit.aes = FALSE,
      vjust = -0.5,
      size = 8,
      family = "sans",
      fontface = "plain"
    ) +
    (
      if (!is.null(sublineage_colors)) {
        scale_fill_manual(values = sublineage_colors, name = "Sublineage", na.value = "grey80")
      } else {
        scale_fill_brewer(palette = "Set1", name = "Sublineage")
      }
    ) +
    scale_x_discrete(expand = c(0, 0), position = "bottom") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_bw(base_family = "sans", base_size = 16) +
    labs(
      x = "Pathway",
      y = "Number of Genes Contributed"
    ) +
    theme(
      axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 26, family = "sans", color = "black"),
      axis.text.y = element_text(size = 26, family = "sans", color = "black"),
      axis.title.x = element_text(size = 30, family = "sans", color = "black", face = "plain"),
      axis.title.y = element_text(size = 30, family = "sans", color = "black", face = "plain"),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 26, family = "sans", color = "black"),
      legend.text = element_text(size = 24, family = "sans", color = "black"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      panel.grid = element_blank(),
      strip.text = element_text(size = 26, face = "plain", family = "sans", color = "black"),
      plot.title = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 80, l = 60),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA)
    )

  # Check for missing sublineage colors
  if (!is.null(sublineage_colors)) {
    missing_sublineages <- setdiff(unique(df_long$sublineage), names(sublineage_colors))
    if (length(missing_sublineages) > 0) {
      warning(sprintf(
        "The following sublineages are present in the data but missing from the color map: %s",
        paste(missing_sublineages, collapse = ", ")
      ))
    }
  }
  plot_base <- path_ext_remove(path_file(csv_file))
  ggsave(file.path(plot_output_dir, paste0(plot_base, ".png")), plot = p, width = 18, height = 18, dpi = 400, bg = "transparent")
}
