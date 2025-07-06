library(ggplot2)
library(Matrix)
library(dplyr)
library(readr)
library(ggpubr)
library(cowplot) # Add at the top if not already loaded
library(grid)

source("src/1. data preprocessing/training datasets/data_loader.R")

datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
raw_umitab <- lung_ldm$dataset$umitab
table_s1 <- select_representative_samples(datasets$table_s1)


# Define doublets for filtering
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Prepare cell metadata and filter umitab
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, doublets = DOUBLETS)
umitab_filtered <- filter_umitab(raw_umitab, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# Map each cell to its sample
cell_sample_map <- cell_metadata_final %>% select(cell_ID, sample_ID)

# Calculate total UMIs per cell using filtered umitab
cell_umis <- data.frame(
  cell_ID = colnames(umitab_filtered),
  total_umis = Matrix::colSums(umitab_filtered)
)

# Map each cell to its sample and sum per sample
sample_counts <- cell_umis %>%
  left_join(cell_sample_map, by = "cell_ID") %>%
  filter(!is.na(sample_ID)) %>%
  group_by(sample_ID) %>%
  summarise(total_umis = sum(total_umis), .groups = 'drop')

# Add tissue type
sample_counts <- sample_counts %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
  filter(!is.na(tissue))

# Scale total UMIs to 10K units for readability
sample_counts <- sample_counts %>%
  mutate(total_umis_10k = total_umis / 10000)

# Density plot of total UMIs per sample by tissue type
p <- ggplot(sample_counts, aes(x = total_umis_10k, fill = tissue, color = tissue)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  scale_color_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(x = "Total UMIs per sample (x10³)", y = "Density", fill = "Tissue", color = "Tissue")

ggsave("output/6. plots/figure 1/fig_1a_count_distribution_density.png", p, width = 7, height = 5, dpi = 300)

# Create and save a separate legend for tissue type (shared across all figure 1 plots)
legend_plot <- ggplot(sample_counts, aes(x = total_umis_10k, fill = tissue, color = tissue)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  scale_color_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.8, linetype = "solid"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.margin = margin(6, 6, 6, 6, unit = "pt") # Add spacing inside legend border
  ) +
  labs(fill = "Tissue", color = "Tissue")

# Extract just the legend
ggpubr_legend <- ggpubr::get_legend(legend_plot)
# Calculate legend size in inches
legend_width_in <- convertWidth(grobWidth(ggpubr_legend), "in", valueOnly = TRUE)
legend_height_in <- convertHeight(grobHeight(ggpubr_legend), "in", valueOnly = TRUE)
# Save the legend as a separate file, tightly cropped
ggsave(
  filename = "output/6. plots/figure 1/shared_tissue_legend.png",
  plot = grid::grid.grabExpr(grid::grid.draw(ggpubr_legend)),
  width = legend_width_in,
  height = legend_height_in,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

message("Completed writing figure 1A to file.")
message("Completed writing shared tissue legend to file.")
