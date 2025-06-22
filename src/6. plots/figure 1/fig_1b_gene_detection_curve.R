# Figure 1B: Gene detection ridgeline plot by lineage and tissue (Normal vs Tumor)
library(ggplot2)
library(Matrix)
library(dplyr)
library(readr)
library(ggridges)

# Load data (assumes working directory is project root)
source("src/1. data preprocessing/training datasets/data_loader.R")
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
ds_matrix <- lung_ldm$dataset$ds[[1]]
table_s1 <- datasets$table_s1
annots_list <- datasets$annots_list

# Define DOUBLETS constant
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Prepare cell metadata (already QC-filtered for downsampled ds_matrix)
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, doublets = DOUBLETS)
cell_metadata <- cell_metadata %>% filter(cell_ID %in% colnames(ds_matrix))

# Add lineage annotation to each cell
cell_metadata <- cell_metadata %>%
  left_join(annots_list %>% select(cluster = cluster, lineage), by = c("cluster_ID" = "cluster"))

# Add tissue annotation to each cell
cell_metadata <- cell_metadata %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")

# Calculate number of detected genes per cell (genes with >=1 UMI)
cell_gene_counts <- data.frame(
  cell_ID = colnames(ds_matrix),
  n_genes = Matrix::colSums(ds_matrix > 0)
)

# Merge with metadata
cell_gene_counts <- cell_gene_counts %>%
  left_join(cell_metadata, by = "cell_ID") %>%
  filter(!is.na(lineage), !is.na(tissue))

# Ridgeline plot: detected genes per cell, by lineage, colored by tissue
p <- ggplot(cell_gene_counts, aes(x = n_genes, y = lineage, fill = tissue, color = tissue)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 1.1, rel_min_height = 0.01, size = 0.4) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  scale_color_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.position = "right",
    legend.background = element_rect(color = "black", fill = "white", size = 0.8, linetype = "solid")
  ) +
  labs(x = "Detected genes per cell", fill = "Tissue", color = "Tissue")

ggsave("output/6. plots/figure 1/fig_1b_gene_detection_ridgeline.png", p, width = 8, height = 6)
