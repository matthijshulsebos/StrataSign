library(dplyr)
library(ggplot2)
library(readr)
library(Matrix)
library(ggpubr)


# Configuration
output_figure_dir <- "output/6. plots/figure 3"
output_plot_name <- "fig_3a.png"

# Load cell metadata from cp10k
cell_metadata_final <- read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>% mutate(sample_ID = as.character(sample_ID))
annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)

# Create output directory
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
output_plot_path <- file.path(output_figure_dir, output_plot_name)

# Calculate cell counts per sample per cell type
cell_counts <- cell_metadata_final %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
  mutate(
    celltype = ifelse(!is.na(sub_lineage) & sub_lineage != "", sub_lineage, lineage),
    celltype = gsub("/", "-", celltype)
  ) %>%
  group_by(sample_ID, celltype) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(sample_ID = as.character(sample_ID))

# Add tissue information for plotting
cell_counts <- cell_counts %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")

# Create the boxplot (per sample)
p <- ggplot(cell_counts, aes(x = celltype, y = cell_count)) +
  geom_boxplot(aes(fill = tissue), alpha = 0.7, width = 0.6) + 
  scale_y_log10(labels = scales::comma_format()) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02), add = c(0.5, 0.5))) +
  labs(
    x = NULL,
    y = "Cell count (log10 scale)",
    fill = "Tissue type"
  ) +
  theme_bw(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, color = "black"),
    legend.position = "none",
    strip.text = element_text(color = "black", face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  ) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00"))

# Save the plot
ggsave(output_plot_path, plot = p, width = 12, height = 6, dpi = 300)

message("Completed writing figure 2A to file.")
