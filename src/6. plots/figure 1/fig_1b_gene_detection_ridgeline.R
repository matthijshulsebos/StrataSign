library(ggplot2)
library(Matrix)
library(dplyr)
library(readr)
library(ggridges)

source("src/1. data preprocessing/training datasets/data_loader.R")

datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
raw_umitab <- lung_ldm$dataset$umitab
table_s1 <- select_representative_samples(datasets$table_s1)
annots_list <- datasets$annots_list


# Define doublets for filtering
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Prepare cell metadata and filter umitab
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, doublets = DOUBLETS)
umitab_filtered <- filter_umitab(raw_umitab, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# Add (sub)lineage annotations to each cell
cell_metadata_final <- cell_metadata_final %>%
  left_join(annots_list %>% select(cluster = cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage))

# Add tissue annotation to each cell
cell_metadata_final <- cell_metadata_final %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")

# Calculate number of detected genes per cell
cell_gene_counts <- data.frame(
  cell_ID = colnames(umitab_filtered),
  n_genes = Matrix::colSums(umitab_filtered > 0)
)

# Merge with metadata
cell_gene_counts <- cell_gene_counts %>%
  left_join(cell_metadata_final, by = "cell_ID") %>%
  filter(!is.na(cell_type), !is.na(tissue))

# Set upper limit percentile for x-axis
x_upper_limit <- quantile(cell_gene_counts$n_genes, 0.98)

# Ridgeline plot of detected genes per cell, by cell type, colored by tissue
p <- ggplot(cell_gene_counts, aes(x = n_genes, y = cell_type, fill = tissue, color = tissue)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 0.95, rel_min_height = 0.01, size = 0.4) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  scale_color_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  scale_x_continuous(limits = c(0, x_upper_limit)) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.position = "right",
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.8, linetype = "solid")
  ) +
  labs(x = "Detected genes per cell", fill = "Tissue", color = "Tissue")

ggsave("output/6. plots/figure 1/fig_1b_gene_detection_ridgeline.png", p, width = 10, height = 10)

message("Completed writing figure 1B to file.")
