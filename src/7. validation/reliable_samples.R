library(dplyr)
library(Matrix)
library(readr)

# Load lung_ldm and umitab (raw, filtered)
load("base/data/lung_ldm.rd")
lung_ldm <- get("lung_ldm")
raw_umitab <- lung_ldm$dataset$umitab

# Load sample metadata and annotations

table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>% mutate(sample_ID = as.character(sample_ID))
annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)

# Doublet clusters (from preprocessing)
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Prepare cell metadata (as in preprocessing)
cell_to_sample <- lung_ldm$dataset$cell_to_sample
cell_to_cluster <- lung_ldm$dataset$cell_to_cluster
cell_metadata <- data.frame(
  cell_ID = names(cell_to_sample),
  sample_ID = unname(cell_to_sample),
  cluster_ID = as.numeric(unname(cell_to_cluster)),
  stringsAsFactors = FALSE
) %>%
  left_join(table_s1 %>% dplyr::select(sample_ID, tissue, patient_ID), by = "sample_ID") %>%
  left_join(annots_list %>% dplyr::select(cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    cluster_name = ifelse(is.na(cluster_name), paste0("Cluster_", cluster_ID), cluster_name),
    tissue = ifelse(is.na(tissue), "Unknown", tissue)
  ) %>%
  filter(!cluster_ID %in% DOUBLETS)

# Filter umitab to cells in cell_metadata
umitab_filtered <- raw_umitab[, cell_metadata$cell_ID]

# Only keep cells present in both
cell_metadata <- cell_metadata %>% filter(cell_ID %in% colnames(umitab_filtered))

# For each sample and cell type, count number of cells
cell_counts <- cell_metadata %>%
  group_by(sample_ID, cluster_name) %>%
  summarise(n_cells = n(), .groups = 'drop')

# For each cell type, get the minimum number of cells per sample
min_cells_per_celltype <- cell_counts %>%
  group_by(cluster_name) %>%
  summarise(min_cells = min(n_cells), n_samples = n(), .groups = 'drop')

# Cell types with at least 50 cells in every sample
reliable_celltypes <- min_cells_per_celltype %>% filter(min_cells >= 50)

print(reliable_celltypes)
