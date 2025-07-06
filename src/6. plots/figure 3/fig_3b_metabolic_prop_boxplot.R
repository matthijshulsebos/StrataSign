library(dplyr)
library(ggplot2)
library(readr)
library(Matrix)
library(ggpubr)


# Configuration
output_figure_dir <- "output/6. plots/figure 3"
output_plot_name <- "fig_3b.png"

# Load cp10k normalized umitab and metadata
umitab_filtered <- readRDS("output/6. plots/data/cp10k/cp10k_normalized_umitab.rds")
cell_metadata_final <- read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>% mutate(sample_ID = as.character(sample_ID))
annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)

# Load metabolic genes
metabolic_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", show_col_types = FALSE)

# Create output directory
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
output_plot_path <- file.path(output_figure_dir, output_plot_name)

# Filter on metabolic genes
genes_in_data <- rownames(umitab_filtered)
metabolic_genes_filtered <- metabolic_genes %>%
  filter(Symbol %in% genes_in_data) %>%
  pull(Symbol)

# Calculate metabolic gene proportion per sample per cell type
umitab_metabolic <- umitab_filtered[metabolic_genes_filtered, ]

# Add sample and patient information to cell metadata
cell_metadata_with_sample <- cell_metadata_final %>%
  mutate(sample_ID = as.character(sample_ID)) %>%
  left_join(table_s1 %>% select(sample_ID, patient_ID, tissue), by = "sample_ID")

# Calculate metabolic gene proportion per sample per cell type
cell_metadata_grouped <- cell_metadata_with_sample %>%
  mutate(sample_celltype = paste(sample_ID, cluster_ID, sep = "_")) %>%
  group_by(sample_celltype, sample_ID, patient_ID, cluster_ID, tissue) %>%
  summarise(
    n_cells = n(),
    cell_indices = list(match(cell_ID, colnames(umitab_filtered))),
    .groups = 'drop'
  )

# Create sample cell type identifier
cell_metadata_with_sample <- cell_metadata_with_sample %>%
  mutate(sample_celltype = paste(sample_ID, cluster_ID, sep = "_"))

# Create a mapping for cell to cell type
cell_to_sample_celltype <- setNames(cell_metadata_with_sample$sample_celltype, cell_metadata_with_sample$cell_ID)

# Get cells that are in our filtered matrix
cells_in_matrix <- intersect(names(cell_to_sample_celltype), colnames(umitab_filtered))
sample_celltype_vector <- cell_to_sample_celltype[cells_in_matrix]

# Calculate sums using tapply
total_umi_sums <- tapply(Matrix::colSums(umitab_filtered[, cells_in_matrix]), 
                        sample_celltype_vector, sum)
metabolic_umi_sums <- tapply(Matrix::colSums(umitab_metabolic[, cells_in_matrix]), 
                            sample_celltype_vector, sum)

# Create results dataframe
sample_celltype_counts <- data.frame(
  sample_celltype = names(total_umi_sums),
  total_umi_celltype = as.numeric(total_umi_sums),
  metabolic_umi_celltype = as.numeric(metabolic_umi_sums)
) %>%
  mutate(
    metabolic_proportion = metabolic_umi_celltype / total_umi_celltype,
    # Extract sample and cluster id from sample_celltype
    sample_ID = sub("_[^_]*$", "", sample_celltype),
    cluster_ID = as.numeric(sub(".*_", "", sample_celltype))
  ) %>%
  # Add back other metadata
  left_join(cell_metadata_with_sample %>% 
           distinct(sample_ID, cluster_ID, patient_ID, tissue), 
           by = c("sample_ID", "cluster_ID")) %>%
  select(-sample_celltype) %>%
  filter(!is.na(metabolic_proportion))


# Aggregate by cell type and patient
patient_celltype_metabolic <- sample_celltype_counts %>%
  left_join(annots_list %>% select(cluster = cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
         cell_type = gsub("/", "-", cell_type))

# Create boxplot
p <- ggplot(patient_celltype_metabolic, aes(x = cell_type, y = metabolic_proportion)) +
  geom_boxplot(aes(fill = tissue), alpha = 0.7, width = 0.6) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02), add = c(0.5, 0.5))) +
  labs(
    x = NULL,
    y = "Metabolic gene expression proportion",
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
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))

# Save plot
ggsave(output_plot_path, plot = p, width = 12, height = 6, dpi = 300)

# Create and save a separate legend for tissue type
legend_plot <- ggplot(patient_celltype_metabolic, aes(x = cell_type, y = metabolic_proportion)) +
  geom_boxplot(aes(fill = tissue), alpha = 0.7) +
  scale_fill_manual(values = c("Normal" = "#56B4E9", "Tumor" = "#E69F00")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  labs(fill = "Tissue type")

# Extract just the legend
legend_grob <- get_legend(legend_plot)

# Save the legend as a separate file
legend_output_path <- file.path(output_figure_dir, "fig_3_shared_tissue_legend.png")
ggsave(legend_output_path, legend_grob, width = 1, height = 0.8, dpi = 300)

message("Completed writing figure 3B to file.")
