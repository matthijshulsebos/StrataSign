# Figure 2C: Metabolic gene detection by tissue (Normal vs Tumor, only metabolic genes)
library(ggplot2)
library(Matrix)
library(dplyr)
library(readr)
library(ggridges)

# Load CP10K normalized counts and metadata
umitab_filtered <- readRDS("output/6. plots/data/cp10k/cp10k_normalized_umitab.rds")
cell_metadata_final <- readr::read_csv("output/6. plots/data/cp10k/cell_metadata_final.csv", show_col_types = FALSE)
annots_list <- readr::read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
table_s1 <- readr::read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE)

# Add (sub)lineage and tissue annotations to each cell
cell_metadata_final <- cell_metadata_final %>%
  left_join(annots_list %>% select(cluster = cluster, lineage, sub_lineage), by = c("cluster_ID" = "cluster")) %>%
  mutate(cell_type = ifelse(!is.na(sub_lineage), sub_lineage, lineage)) %>%
  left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")

# Load metabolic genes
metabolic_genes <- readr::read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", show_col_types = FALSE)
all_genes <- rownames(umitab_filtered)
metabolic_genes_filtered <- metabolic_genes %>%
  dplyr::filter(Symbol %in% all_genes) %>%
  dplyr::pull(Symbol)

n_metabolic_total <- length(metabolic_genes_filtered)
metabolic_gene_idx <- which(all_genes %in% metabolic_genes_filtered)


# Helper: add cell_type_group for 2-column facet (kept for 2c faceting)
add_cell_type_group <- function(df, cell_type_col = "cell_type") {
  unique_cell_types <- sort(unique(df[[cell_type_col]]))
  n_types <- length(unique_cell_types)
  df %>% mutate(
    cell_type_group = ifelse(
      .data[[cell_type_col]] %in% unique_cell_types[1:ceiling(n_types/2)],
      "Column 1", "Column 2"
    ),
    cell_type_group = factor(cell_type_group, levels = c("Column 1", "Column 2"))
  )
}

cell_metabolic_tissue <- data.frame(
  cell_ID = colnames(umitab_filtered),
  n_metabolic_genes = Matrix::colSums(umitab_filtered[metabolic_gene_idx, , drop = FALSE] > 0) / n_metabolic_total
) %>%
  left_join(cell_metadata_final, by = "cell_ID") %>%
  filter(!is.na(cell_type), !is.na(tissue)) %>%
  add_cell_type_group()

x_upper_limit_2c <- quantile(cell_metabolic_tissue$n_metabolic_genes, 0.98)
color_scheme_tissue <- c("Normal" = "#009E73", "Tumor" = "#D55E00")
p_2c <- ggplot(cell_metabolic_tissue, aes(x = n_metabolic_genes, y = cell_type, fill = tissue, color = tissue)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 0.95, rel_min_height = 0.01, size = 0.4) +
  scale_fill_manual(values = color_scheme_tissue) +
  scale_color_manual(values = color_scheme_tissue) +
  scale_x_continuous(limits = c(0, x_upper_limit_2c)) +
  facet_wrap(~ cell_type_group, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    legend.position = "right",
    legend.background = element_rect(color = NA, fill = "white", linewidth = 0, linetype = "blank"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.margin = margin(16, 16, 16, 16, unit = "pt"),
    strip.text = element_blank()
  ) +
  labs(x = "Detected metabolic genes per cell", fill = "Tissue", color = "Tissue")
ggsave("output/6. plots/figure 2/fig_2c_metabolic_gene_detection_by_tissue_ridgeline.png", p_2c, width = 12, height = 6, dpi = 300)

message("Completed writing metabolic gene detection ridgeline plot (2C) to file.")
