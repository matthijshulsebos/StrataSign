# Violin plot: number of cells per cell type per tissue type (across selected samples)
library(dplyr)
library(ggplot2)
library(readr)

# Load data loader functions
data_loader_path <- "src/1. data preprocessing/training datasets/data_loader.R"
source(data_loader_path)

# Load main dataset and sample metadata (using the same selection as preprocessing pipeline)
datasets <- load_all_datasets()
table_s1 <- select_representative_samples(datasets$table_s1)
lung_ldm <- datasets$lung_ldm

ds_matrix <- lung_ldm$dataset$ds[[1]]

# Prepare cell metadata (same as preprocessing pipeline)
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)

# Filter cell_metadata to only those cells present in ds_matrix (as in preprocessing pipeline)
cell_metadata <- cell_metadata %>% filter(cell_ID %in% colnames(ds_matrix))

# Join with sample metadata to get tissue type
cell_metadata <- cell_metadata %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")

# Remove any cells without tissue annotation (shouldn't happen, but safe)
cell_metadata <- cell_metadata %>% filter(!is.na(tissue))

# Count number of cells per (sample, cell type, tissue)
cell_counts <- cell_metadata %>%
  group_by(sample_ID, cluster_ID, tissue) %>%
  summarise(n_cells = n(), .groups = 'drop')

# Optionally, join with annotation to get cell type names
annots_list <- datasets$annots_list
cell_counts <- cell_counts %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
  mutate(cell_type = ifelse(!is.na(sub_lineage) & sub_lineage != "", sub_lineage, lineage))

# For plotting, create a combined label for cell type and tissue
group_label <- function(cell_type, tissue) paste(cell_type, tissue, sep = " | ")
cell_counts <- cell_counts %>% mutate(group = group_label(cell_type, tissue))

# Save cell counts per cell type per sample with tissue and sample ID to CSV
output_csv_path <- "output/6. plots/figure s1/cell_counts_per_sample.csv"
dir.create(dirname(output_csv_path), recursive = TRUE, showWarnings = FALSE)
write.csv(cell_counts, output_csv_path, row.names = FALSE)

# Calculate and save total counts per sample, including clustering model flag
total_counts_per_sample <- cell_counts %>%
  group_by(sample_ID, tissue) %>%
  summarise(total_cells = sum(n_cells), .groups = 'drop') %>%
  left_join(table_s1 %>% select(sample_ID, Use.in.Clustering.Model.), by = "sample_ID")
write.csv(total_counts_per_sample, "output/6. plots/figure s1/total_counts_per_sample.csv", row.names = FALSE)

# Calculate and save total counts per cell type (lineage) per sample
cell_counts_lineage <- cell_counts %>%
  group_by(sample_ID, tissue, lineage) %>%
  summarise(total_cells = sum(n_cells), .groups = 'drop')
write.csv(cell_counts_lineage, "output/6. plots/figure s1/total_counts_per_lineage_per_sample.csv", row.names = FALSE)

# Calculate and print total counts for the entire dataset
total_cells_all <- sum(cell_counts$n_cells)
cat("Total number of cells in the dataset:", total_cells_all, "\n")
write.csv(data.frame(total_cells_all = total_cells_all), "output/6. plots/figure s1/total_cells_all.csv", row.names = FALSE)

# Plot: violin plot of n_cells per (cell type, tissue) across samples, remove outliers
ggplot(cell_counts, aes(x = group, y = n_cells, fill = tissue)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(x = "Cell type | Tissue", y = "Number of cells per sample", title = "Cell counts per cell type and tissue (selected samples)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to file
ggsave(
  filename = "output/6. plots/figure s1/cell_counts_violin.png",
  plot = last_plot(),
  width = 14, height = 7, dpi = 300
)

# For each sublineage, check if it has >=10 cells in >=50% of samples
sublineage_counts <- cell_counts %>%
  filter(!is.na(sub_lineage)) %>%
  group_by(sub_lineage, sample_ID) %>%
  summarise(n_cells = sum(n_cells), .groups = 'drop')

sublineage_coverage <- sublineage_counts %>%
  group_by(sub_lineage) %>%
  summarise(
    n_samples_with_10 = sum(n_cells >= 10),
    n_samples_present = n(),
    prop_samples_with_10 = n_samples_with_10 / n_samples_present,
    .groups = 'drop'
  ) %>%
  mutate(
    meets_criteria = prop_samples_with_10 >= 0.5
  )

write.csv(sublineage_coverage, "output/6. plots/figure s1/sublineage_coverage_10in50.csv", row.names = FALSE)

# For each lineage, check if it has >=10 cells in >=50% of samples
lineage_counts <- cell_counts %>%
  filter(!is.na(lineage)) %>%
  group_by(lineage, sample_ID) %>%
  summarise(n_cells = sum(n_cells), .groups = 'drop')

lineage_coverage <- lineage_counts %>%
  group_by(lineage) %>%
  summarise(
    n_samples_with_10 = sum(n_cells >= 10),
    n_samples_present = n(),
    prop_samples_with_10 = n_samples_with_10 / n_samples_present,
    .groups = 'drop'
  ) %>%
  mutate(
    meets_criteria = prop_samples_with_10 >= 0.5
  )

write.csv(lineage_coverage, "output/6. plots/figure s1/lineage_coverage_10in50.csv", row.names = FALSE)

# Print summary
cat("Sublineages with >=10 cells in >=50% of samples:\n")
print(sublineage_coverage %>% filter(meets_criteria))
cat("Lineages with >=10 cells in >=50% of samples:\n")
print(lineage_coverage %>% filter(meets_criteria))
