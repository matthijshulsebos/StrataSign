library(Matrix)
library(readr)
library(dplyr)
library(Matrix)
library(readr)

source("src/1. data preprocessing/training datasets/data_loader.R")

DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

annots_list <- load_annotations()

# Load filtered cell metadata
cell_metadata <- read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)

# Filter out doublets
cell_metadata_filtered <- cell_metadata %>% filter(!cluster_ID %in% DOUBLETS)

# Count number of cells per cluster per sample
cells_per_cluster_sample <- cell_metadata_filtered %>%
  group_by(cluster_ID, sample_ID) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  arrange(cluster_ID, sample_ID)

# Coerce cluster_ID and cluster to character immediately after filtering and loading
cells_per_cluster_sample <- cells_per_cluster_sample %>% mutate(cluster_ID = as.character(cluster_ID))
annots_list_full <- annots_list %>% mutate(cluster = as.character(cluster), sub_lineage = as.character(sub_lineage))

# Write summary to output/7. validation
output_dir <- "output/7. validation"
output_file <- file.path(output_dir, "cells_per_cluster_per_sample.csv")
write_csv(cells_per_cluster_sample, output_file)

# Calculate minimum number of cells per cluster
min_cells_per_cluster <- cells_per_cluster_sample %>%
  group_by(cluster_ID) %>%
  summarise(min_cells = min(n_cells), .groups = 'drop') %>%
  arrange(cluster_ID)

# Write min cells summary
min_cells_file <- file.path(output_dir, "min_cells_per_cluster.csv")
write_csv(min_cells_per_cluster, min_cells_file)

# Load cluster annotations and coerce to character immediately
annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
annots_list_full <- annots_list %>% mutate(cluster = as.character(cluster), sub_lineage = as.character(sub_lineage))

matching_clusters <- intersect(cells_per_cluster_sample$cluster_ID, annots_list_full$cluster)

# Join to get lineage for each cluster
cells_with_lineage <- cells_per_cluster_sample %>%
  left_join(annots_list_full %>% select(cluster, lineage), by = c("cluster_ID" = "cluster"))

# Aggregate per lineage per sample
cells_per_lineage_sample <- cells_with_lineage %>%
  group_by(lineage, sample_ID) %>%
  summarise(n_cells = sum(n_cells), .groups = 'drop') %>%
  arrange(lineage, sample_ID)

# Write lineage summary
lineage_file <- file.path(output_dir, "cells_per_lineage_per_sample.csv")
write_csv(cells_per_lineage_sample, lineage_file)

# Load Table S1 sample metadata
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE)

# Join cluster annotations and Table S1 tissue
cells_with_annots <- cells_per_cluster_sample %>%
  left_join(annots_list_full %>% select(cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")

# Strip whitespace from sublineage before checking
cells_with_annots <- cells_with_annots %>%
  mutate(sublineage_clean = trimws(sub_lineage)) %>%
  mutate(sublineage_or_lineage = ifelse(!is.na(sublineage_clean) & sublineage_clean != "" & sublineage_clean != "NA", sublineage_clean, lineage))

cells_per_tissue_sublineage <- cells_with_annots %>%
  group_by(tissue, sublineage_or_lineage) %>%
  summarise(n_cells = sum(n_cells), .groups = 'drop') %>%
  arrange(tissue, sublineage_or_lineage)

# Pivot to wide format and calculate increase
tissue_comparison <- cells_per_tissue_sublineage %>%
  tidyr::pivot_wider(names_from = tissue, values_from = n_cells, values_fill = 0) %>%
  rename_with(tolower) %>%
  mutate(
    count_increase = tumor - normal,
    percentage_increase = case_when(
      normal == 0 & tumor > 0 ~ Inf,
      normal == 0 & tumor == 0 ~ 0,
      TRUE ~ ((tumor - normal) / normal) * 100
    )
  ) %>%
  arrange(desc(abs(percentage_increase)))

# Write summary
tissue_sublineage_file <- file.path(output_dir, "cells_per_tissue_per_sublineage_or_lineage.csv")
write_csv(tissue_comparison, tissue_sublineage_file)
