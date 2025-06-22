library(dplyr)
library(ggplot2)
library(readr)
library(Matrix)

# Load required datasets
source("src/1. data preprocessing/training datasets/data_loader.R")

# Configuration
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
output_figure_dir <- "output/6. plots/cell type abundance"
output_plot_name <- "cell_type_counts_per_patient_boxplot.png"

# Load datasets
message("Loading datasets...")
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
table_s1 <- datasets$table_s1
annots_list <- datasets$annots_list

# Create output directory
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
output_plot_path <- file.path(output_figure_dir, output_plot_name)

# Prepare data for analysis
message("Preparing data...")

# Extract downsampled UMI matrix and prepare cell metadata
ds_matrix <- extract_downsampled_umitab(lung_ldm)
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)
umitab_filtered <- filter_umitab(ds_matrix, cell_metadata)
cell_metadata_final <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# Convert UMI matrix to long format and aggregate by cell type
message("Converting UMI matrix to long format...")

# Ensure input is sparse
if (!inherits(umitab_filtered, "Matrix")) {
  umitab_filtered <- as(umitab_filtered, "sparseMatrix")
}

# Convert sparse matrix to triplet format
umitab_triplet <- summary(umitab_filtered)
colnames(umitab_triplet) <- c("gene_idx", "cell_idx", "umi_count")

# Map indices to names
umitab_triplet$gene <- rownames(umitab_filtered)[umitab_triplet$gene_idx]
umitab_triplet$cell_ID <- colnames(umitab_filtered)[umitab_triplet$cell_idx]

# Join with metadata and aggregate by sample and cluster
counts_long <- umitab_triplet %>%
  select(gene, cell_ID, umi_count) %>%
  inner_join(cell_metadata_final %>% select(cell_ID, sample_ID, cluster_ID), by = "cell_ID") %>%
  group_by(sample_ID, cluster_ID) %>%
  summarise(total_counts = sum(umi_count), .groups = 'drop') %>%
  mutate(
    sample_ID = as.character(sample_ID),
    cluster_ID = as.numeric(cluster_ID)
  )

# Check data completeness before any filtering
message("Checking data completeness...")
all_samples <- unique(counts_long$sample_ID)
all_clusters <- unique(counts_long$cluster_ID)
expected_combinations <- length(all_samples) * length(all_clusters)
actual_combinations <- nrow(counts_long %>% distinct(sample_ID, cluster_ID))

message(sprintf("Expected sample-cluster combinations: %d", expected_combinations))
message(sprintf("Actual sample-cluster combinations: %d", actual_combinations))
message(sprintf("Data completeness: %.1f%%", (actual_combinations / expected_combinations) * 100))

# Check how many combinations have zero counts
zero_count_combinations <- counts_long %>%
  filter(total_counts == 0) %>%
  nrow()

message(sprintf("Sample-cluster combinations with zero total counts: %d (%.1f%%)", 
                zero_count_combinations, 
                (zero_count_combinations / actual_combinations) * 100))

# Calculate total counts per patient per cell type
patient_celltype_counts <- counts_long %>%
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

# Show distribution of zero vs non-zero counts
message("\nDistribution of counts:")
count_distribution <- patient_celltype_counts %>%
  mutate(count_category = ifelse(total_counts == 0, "Zero", "Non-zero")) %>%
  count(count_category) %>%
  mutate(percentage = n / sum(n) * 100)
print(count_distribution)

# Filter out zero counts for visualization (but keep the analysis above)
patient_celltype_counts_nonzero <- patient_celltype_counts %>%
  filter(total_counts > 0)

# Create the boxplot
message("Creating boxplot...")
p <- ggplot(patient_celltype_counts_nonzero, aes(x = reorder(cluster_label, cluster_ID), y = total_counts)) +
  geom_boxplot(aes(fill = tissue), alpha = 0.7) +
  scale_y_log10(labels = scales::comma_format()) +
  labs(
    title = "Raw Counts per Patient by Cell Type",
    subtitle = "Distribution of total raw counts across patients for each cell type (excluding zero counts)",
    x = "Cell Type",
    y = "Total Raw Counts (log10 scale)",
    fill = "Tissue Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("Normal" = "#4CAF50", "Tumor" = "#F44336"))

# Save the plot
message(paste("Saving plot to:", output_plot_path))
ggsave(output_plot_path, plot = p, width = 16, height = 8, dpi = 300)

# Print summary statistics
message("Summary statistics:")
summary_stats <- patient_celltype_counts %>%
  group_by(cluster_ID, cluster_label, tissue) %>%
  summarise(
    n_patients = n(),
    median_counts = median(total_counts),
    mean_counts = mean(total_counts),
    min_counts = min(total_counts),
    max_counts = max(total_counts),
    .groups = 'drop'
  ) %>%
  arrange(cluster_ID)

print(summary_stats)

message("Cell type prevalence analysis complete!")
