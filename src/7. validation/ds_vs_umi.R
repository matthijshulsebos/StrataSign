library(dplyr)
library(Matrix)
library(readr)

# Load lung_ldm and raw umitab
load("base/data/lung_ldm.rd")
lung_ldm <- get("lung_ldm")
raw_umitab <- lung_ldm$dataset$umitab

# Load sample metadata and annotations
cell_to_sample <- lung_ldm$dataset$cell_to_sample
cell_to_cluster <- lung_ldm$dataset$cell_to_cluster
annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
sample_metadata <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>% mutate(sample_ID = as.character(sample_ID))

# Prepare cell metadata (add tissue)
cell_metadata <- data.frame(
  cell_ID = names(cell_to_sample),
  sample_ID = unname(cell_to_sample),
  cluster_ID = as.numeric(unname(cell_to_cluster)),
  stringsAsFactors = FALSE
) %>%
  left_join(sample_metadata %>% dplyr::select(sample_ID, tissue), by = "sample_ID") %>%
  left_join(annots_list %>% dplyr::select(cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    cluster_name = ifelse(is.na(cluster_name), paste0("Cluster_", cluster_ID), cluster_name)
  )

# Calculate total UMIs per cell
umis_per_cell <- Matrix::colSums(raw_umitab)
umis_df <- data.frame(cell_ID = names(umis_per_cell), nUMI = as.numeric(umis_per_cell))

# Join with cell metadata to get cell type and tissue
umis_df <- umis_df %>% left_join(cell_metadata, by = "cell_ID")

# Calculate statistics per cell type and tissue
stats_per_celltype_tissue <- umis_df %>%
  group_by(cluster_name, tissue) %>%
  summarise(
    n_cells = n(),
    mean_UMI = mean(nUMI),
    median_UMI = median(nUMI),
    sd_UMI = sd(nUMI),
    min_UMI = min(nUMI),
    max_UMI = max(nUMI),
    q25_UMI = quantile(nUMI, 0.25),
    q75_UMI = quantile(nUMI, 0.75),
    .groups = 'drop'
  ) %>%
  arrange(cluster_name, desc(mean_UMI))

print(stats_per_celltype_tissue, n = Inf)

# Number of cells per cell type per sample
cells_per_type_sample <- umis_df %>%
  group_by(cluster_name, sample_ID) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  arrange(cluster_name, sample_ID)

#print(cells_per_type_sample, n = Inf)

# Number of cells per cell type per sample and tissue
cells_per_type_sample_tissue <- umis_df %>%
  group_by(cluster_name, tissue, sample_ID) %>%
  summarise(n_cells = n(), .groups = 'drop')

# For each cell type and tissue, get the minimum number of cells per sample
min_cells_per_type_tissue <- cells_per_type_sample_tissue %>%
  group_by(cluster_name, tissue) %>%
  summarise(min_cells = min(n_cells), .groups = 'drop') %>%
  arrange(cluster_name, tissue)

#print(min_cells_per_type_tissue, n = Inf)
