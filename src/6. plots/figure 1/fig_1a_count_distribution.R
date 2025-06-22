# Figure 1A: Distribution of counts (UMIs) per sample, split by tissue type (Normal vs Tumor)
library(ggplot2)
library(Matrix)
library(dplyr)
library(readr)

# Load data (assumes working directory is project root)
source("src/1. data preprocessing/training datasets/data_loader.R")
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
ds_matrix <- lung_ldm$dataset$ds[[1]]
table_s1 <- datasets$table_s1

# Define DOUBLETS constant
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Prepare cell metadata (already QC-filtered for downsampled ds_matrix)
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, doublets = DOUBLETS)
cell_metadata <- cell_metadata %>% filter(cell_ID %in% colnames(ds_matrix))

# Map each cell to its sample
cell_sample_map <- cell_metadata %>% select(cell_ID, sample_ID)

# Calculate total UMIs per cell once
cell_umis <- data.frame(
  cell_ID = colnames(ds_matrix),
  total_umis = Matrix::colSums(ds_matrix)
)

# Map each cell to its sample and sum per sample
sample_counts <- cell_umis %>%
  left_join(cell_sample_map, by = "cell_ID") %>%
  filter(!is.na(sample_ID)) %>%
  group_by(sample_ID) %>%
  summarise(total_umis = sum(total_umis), .groups = 'drop')

# Add tissue type (Normal/Tumor) from table_s1
sample_counts <- sample_counts %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
  filter(!is.na(tissue))

# Scale total UMIs to 10K units for better readability
sample_counts <- sample_counts %>%
  mutate(total_umis_10k = total_umis / 10000)

# Density plot of total UMIs per sample (in 10K), colored by tissue type
ggplot(sample_counts, aes(x = total_umis_10k, fill = tissue, color = tissue)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  scale_color_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.background = element_rect(color = "black", fill = "white", size = 0.8, linetype = "solid")
  ) +
  labs(x = "Total UMIs per sample (x10³)", y = "Density", fill = "Tissue", color = "Tissue")

ggsave("output/6. plots/figure 1/fig_1a_count_distribution_density.png", width = 7, height = 5)
