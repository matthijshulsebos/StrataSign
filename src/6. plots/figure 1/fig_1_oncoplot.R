library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(viridis)
library(tibble)

output_figure_dir <- "output/6. plots/figure 1"
output_plot_name <- "fig_1a_patient_celltype_proportion_oncoplot.png"


# Load the sample filtered CP10K cell metadata
cell_metadata_final <- read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)

# Load annotation files
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>%
  mutate(sample_ID = as.character(sample_ID))
annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)

# Create output directory
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
output_plot_path <- file.path(output_figure_dir, output_plot_name)

# Add sample and cell type to cell metadata
cell_metadata_annotated <- cell_metadata_final %>%
  # Ensure sample_ID is character
  mutate(sample_ID = as.character(sample_ID)) %>%
  left_join(table_s1 %>% select(sample_ID, patient_ID, tissue), by = "sample_ID") %>%
  left_join(annots_list %>% select(cluster = cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cell_type = gsub("/", "-", cell_type),
    cluster_label = paste0(cell_type, " (", cluster_ID, ")")
  ) %>%
  filter(!is.na(patient_ID), !is.na(tissue), !is.na(cell_type))

# Get all unique cell types ordered by cluster_ID for the x axis labels
all_cell_types <- cell_metadata_annotated %>%
  select(cluster_ID, cell_type, cluster_label) %>%
  distinct() %>%
  arrange(cluster_ID) %>%
  # Create a dense rank for cell types to ensure no gaps in plotting
  mutate(cell_type_order = dense_rank(cluster_ID))

# Get all unique patients
all_patients <- cell_metadata_annotated %>%
  select(patient_ID) %>%
  distinct() %>%
  pull(patient_ID)

# Calculate cell type proportions per patient per tissue
patient_celltype_counts <- cell_metadata_annotated %>%
  group_by(patient_ID, tissue, cluster_ID, cell_type, cluster_label) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(patient_ID, tissue) %>%
  mutate(
    total_cells_patient_tissue = sum(cell_count),
    proportion = cell_count / total_cells_patient_tissue
  ) %>%
  ungroup() %>%
  select(patient_ID, tissue, cluster_ID, proportion)

# Get patients who have each tissue type
patients_with_normal <- cell_metadata_annotated %>%
  filter(tissue == "Normal") %>%
  distinct(patient_ID) %>%
  pull(patient_ID)

patients_with_tumor <- cell_metadata_annotated %>%
  filter(tissue == "Tumor") %>%
  distinct(patient_ID) %>%
  pull(patient_ID)

# Create complete matrix
complete_matrix <- expand_grid(
  patient_ID = all_patients,
  tissue = c("Normal", "Tumor"),
  # Use the ordered cell_types
  cell_type_order = all_cell_types$cell_type_order
) %>%
  # Join cell type on the cluster ids
  left_join(all_cell_types, by = "cell_type_order") %>%
  # Join proportion data
  left_join(patient_celltype_counts, by = c("patient_ID", "tissue", "cluster_ID")) %>%
  mutate(
    proportion = ifelse(is.na(proportion), 0, proportion),
    # Mark missing biopt
    is_missing_tissue = case_when(
      tissue == "Normal" & !patient_ID %in% patients_with_normal ~ TRUE,
      tissue == "Tumor" & !patient_ID %in% patients_with_tumor ~ TRUE,
      # Sets rows to true or false based on previous conditions
      TRUE ~ FALSE
    )
  )

# Get patient order based on subtype
patient_order <- table_s1 %>%
  select(patient_ID, subtype = disease) %>%
  distinct() %>%
  # Order by cancer subtype and patient id
  arrange(subtype, patient_ID) %>%
  # Create simple sequential numbering
  mutate(patient_order_val = row_number())

# Order cell types by average proportion in normal sample
celltype_avg_normal <- complete_matrix %>%
  filter(tissue == "Normal") %>%
  group_by(cell_type) %>%
  summarise(avg_prop = mean(proportion, na.rm = TRUE), .groups = 'drop') %>%
  arrange(avg_prop)
celltype_levels <- celltype_avg_normal$cell_type

# Prepare data for complete heatmap
plot_data_facet <- complete_matrix %>%
  left_join(patient_order, by = "patient_ID") %>%
  mutate(
    patient_label = paste0("P", patient_ID),
    cell_type = factor(cell_type, levels = celltype_levels),
    tissue = factor(tissue, levels = c("Normal", "Tumor")),
    subtype = factor(subtype, levels = sort(unique(patient_order$subtype))),
    patient_ID = factor(patient_ID, levels = patient_order$patient_ID)
  )

# Complete heatmap
p <- ggplot(plot_data_facet, aes(x = cell_type, y = patient_ID)) +
  geom_tile(aes(fill = ifelse(is_missing_tissue, NA, sqrt(proportion))), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "gray",
    name = "cell type count\n / total cell count",
    guide = guide_colorbar(
      barwidth = 1.5,
      barheight = 8,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  facet_grid(subtype ~ tissue, scales = "free_y", space = "free_y") +
  scale_x_discrete(
    expand = c(0, 0),
    position = "bottom"
  ) +
  scale_y_discrete(
    labels = plot_data_facet %>% distinct(patient_ID, patient_label) %>% deframe(),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, family = "sans"),
    axis.text.y = element_text(size = 10, family = "sans"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 10, family = "sans"),
    legend.text = element_text(size = 8, family = "sans"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 11, face = "plain", family = "sans"),
    plot.margin = margin(t = 10, r = 5, b = 10, l = 30)
  )

# Save the plot
ggsave(output_plot_path, plot = p, width = 12, height = 6, dpi = 300, bg = "white", device = "png")

message("Completed writing figure 1 oncoplot to file.")
