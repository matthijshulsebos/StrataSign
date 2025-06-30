library(dplyr)
library(ggplot2)
library(readr)
library(Matrix)

# Load required datasets
source("src/1. data preprocessing/training datasets/data_loader.R")

# Configuration
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
output_figure_dir <- "output/6. plots/figure 1"
output_plot_name <- "fig_1c_cell_type_counts_per_patient_boxplot.png"

# Load datasets
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
raw_umitab <- lung_ldm$dataset$umitab
table_s1 <- select_representative_samples(datasets$table_s1)
annots_list <- datasets$annots_list

# Create output directory
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
output_plot_path <- file.path(output_figure_dir, output_plot_name)

# Prepare cell metadata and filter umitab
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, doublets = DOUBLETS)
umitab_filtered <- filter_umitab(raw_umitab, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# Calculate cell counts per sample per cell type (simpler approach)
cell_counts <- cell_metadata_final %>%
  group_by(sample_ID, cluster_ID) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(
    sample_ID = as.character(sample_ID),
    cluster_ID = as.numeric(cluster_ID)
  )

# Calculate cell counts per patient per cell type
patient_celltype_counts <- cell_counts %>%
  # Add patient information
  left_join(table_s1 %>% select(sample_ID, patient_ID, tissue), by = "sample_ID") %>%
  filter(!is.na(patient_ID)) %>%
  # Add cell type annotations
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    cluster_label = paste0(cluster_name, " (", cluster_ID, ")")
  )

# Create the boxplot
p <- ggplot(patient_celltype_counts, aes(x = reorder(cluster_label, cluster_ID), y = cell_count)) +
  geom_boxplot(aes(fill = tissue), alpha = 0.7) +
  scale_y_log10(labels = scales::comma_format()) +
  labs(
    x = "Cell Type",
    y = "Cell Count (log10 scale)",
    fill = "Tissue Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00"))

# Save the plot
ggsave(output_plot_path, plot = p, width = 16, height = 8, dpi = 300)

message("Completed writing figure 1C to file.")
