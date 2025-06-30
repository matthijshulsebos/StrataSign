library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)
library(Matrix)

# Load required datasets
source("src/1. data preprocessing/training datasets/data_loader.R")

# Configuration
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
output_figure_dir <- "output/6. plots/figure 1"
output_plot_name <- "fig_1d_metabolic_proportion_boxplot.png"

# Load datasets
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
raw_umitab <- lung_ldm$dataset$umitab
table_s1 <- select_representative_samples(datasets$table_s1)
hsa01100_genes <- datasets$hsa01100_genes

# Create output directory
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
output_plot_path <- file.path(output_figure_dir, output_plot_name)

# Prepare cell metadata and filter umitab
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, doublets = DOUBLETS)
umitab_filtered <- filter_umitab(raw_umitab, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# Get metabolic genes list
metabolic_genes_list <- hsa01100_genes$Symbol

# Calculate metabolic proportions per sample
# Aggregate counts by sample and gene
sample_gene_counts <- umitab_filtered %>%
  summary() %>%
  as.data.frame() %>%
  setNames(c("gene_idx", "cell_idx", "count")) %>%
  mutate(
    gene = rownames(umitab_filtered)[gene_idx],
    cell_ID = colnames(umitab_filtered)[cell_idx]
  ) %>%
  inner_join(cell_metadata_final %>% select(cell_ID, sample_ID), by = "cell_ID") %>%
  group_by(sample_ID, gene) %>%
  summarise(total_count = sum(count), .groups = 'drop')

# Calculate total and metabolic counts per sample
total_counts_per_sample <- sample_gene_counts %>%
  group_by(sample_ID) %>%
  summarise(total_sample_counts = sum(total_count), .groups = 'drop')

metabolic_counts_per_sample <- sample_gene_counts %>%
  filter(gene %in% metabolic_genes_list) %>%
  group_by(sample_ID) %>%
  summarise(total_metabolic_counts = sum(total_count), .groups = 'drop')

# Calculate proportions
proportion_data <- total_counts_per_sample %>%
  left_join(metabolic_counts_per_sample, by = "sample_ID") %>%
  mutate(
    total_metabolic_counts = ifelse(is.na(total_metabolic_counts), 0, total_metabolic_counts),
    proportion_metabolic = ifelse(total_sample_counts > 0, total_metabolic_counts / total_sample_counts, 0)
  ) %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
  filter(!is.na(tissue), tissue %in% c("Normal", "Tumor"))

# Create the boxplot
p <- ggplot(proportion_data, aes(x = tissue, y = proportion_metabolic, fill = tissue)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  labs(
    x = "Tissue Type",
    y = "Proportion of Metabolic Gene Counts",
    fill = "Tissue Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

# Save the plot
ggsave(output_plot_path, plot = p, width = 7, height = 6, dpi = 300)

message("Completed writing figures 1DEF to file.")
