# Calculate average number of cells per sublineage per sample and coverage criteria
library(dplyr)
library(readr)

# Load data loader and main dataset
source("src/1. data preprocessing/training datasets/data_loader.R")
datasets <- load_all_datasets()
table_s1 <- select_representative_samples(datasets$table_s1)
lung_ldm <- datasets$lung_ldm

ds_matrix <- lung_ldm$dataset$ds[[1]]
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)
cell_metadata <- cell_metadata %>% filter(cell_ID %in% colnames(ds_matrix))

# Join with annotation to get sublineage
annots_list <- datasets$annots_list
cell_metadata <- cell_metadata %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster"))

# Remove NA sublineages
cell_metadata <- cell_metadata %>% filter(!is.na(sub_lineage))

# Count number of cells per sublineage per sample
sublineage_counts <- cell_metadata %>%
  group_by(sub_lineage, sample_ID) %>%
  summarise(n_cells = n(), .groups = 'drop')

# Calculate average number of cells per sublineage per sample
sublineage_avg <- sublineage_counts %>%
  group_by(sub_lineage) %>%
  summarise(
    avg_cells_per_sample = mean(n_cells),
    n_samples_with_10 = sum(n_cells >= 10),
    n_samples_present = n(),
    prop_samples_with_10 = n_samples_with_10 / n_samples_present,
    .groups = 'drop'
  ) %>%
  mutate(
    has_10_in_50pct = prop_samples_with_10 >= 0.5,
    has_10_in_70pct = prop_samples_with_10 >= 0.7
  )

# Map sublineage to cluster IDs
sublineage_to_clusters <- cell_metadata %>%
  group_by(sub_lineage) %>%
  summarise(cluster_IDs = paste(sort(unique(cluster_ID)), collapse=","), .groups = 'drop')

# Add cluster_IDs to sublineage_avg
sublineage_avg <- sublineage_avg %>%
  left_join(sublineage_to_clusters, by = "sub_lineage")

# Ensure output directory exists before writing CSV
output_dir <- "output/1. data preprocessing/sublineage_mask"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(sublineage_avg, file.path(output_dir, "cells_per_sublineage_per_sample.csv"), row.names = FALSE)

# Write arrays of cluster IDs for 10_in_50pct and 10_in_70pct sublineages
clusters_10_in_50pct <- sublineage_to_clusters$cluster_IDs[sublineage_avg$has_10_in_50pct]
clusters_10_in_50pct_vec <- sort(unique(unlist(strsplit(clusters_10_in_50pct, ","))))
clusters_10_in_50pct_vec <- as.integer(clusters_10_in_50pct_vec)
clusters_10_in_70pct <- sublineage_to_clusters$cluster_IDs[sublineage_avg$has_10_in_70pct]
clusters_10_in_70pct_vec <- sort(unique(unlist(strsplit(clusters_10_in_70pct, ","))))
clusters_10_in_70pct_vec <- as.integer(clusters_10_in_70pct_vec)
writeLines(paste(clusters_10_in_50pct_vec, collapse=", "), file.path(output_dir, "clusters_10_in_50pct.txt"))
writeLines(paste(clusters_10_in_70pct_vec, collapse=", "), file.path(output_dir, "clusters_10_in_70pct.txt"))
