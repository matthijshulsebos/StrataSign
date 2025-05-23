library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales) # For percent_format in ggplot

# --- Configuration ---
# Path to the RData file containing lung_ldm
lung_ldm_rdata_path <- "base/data/lung_ldm.rd" 
# Path for metabolic genes
metabolic_genes_csv_path <- "output/1. data preprocessing/kegg/hsa01100_genes.csv" 
# Path for sample metadata
table_s1_path <- "base/input_tables/table_s1_sample_table.csv"
# Output directory and plot name
output_figure_dir <- "output/6. plots/exploratory" 
output_plot_name <- "metabolic_gene_proportions_by_tissue_from_lung_ldm.png"

# Doublet clusters to be removed (from preprocess_relative.R)
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Create output directory if it doesn't exist
dir.create(output_figure_dir, showWarnings = FALSE, recursive = TRUE)
output_plot_path <- file.path(output_figure_dir, output_plot_name)

# --- Load Data ---
message(paste("Loading Leader et al. data (lung_ldm) from:", lung_ldm_rdata_path))
if (!file.exists(lung_ldm_rdata_path)) {
  stop(paste("lung_ldm.rd file not found at:", lung_ldm_rdata_path))
}
if (!exists("lung_ldm")) { 
  load(lung_ldm_rdata_path)
}

if (!exists("lung_ldm") || !("dataset" %in% names(lung_ldm)) || !("counts" %in% names(lung_ldm$dataset))) {
  stop("lung_ldm object or lung_ldm$dataset$counts not found after loading RData.")
}
counts_3d_raw <- lung_ldm$dataset$counts  
message("Original 3D counts dimensions: ", paste(dim(counts_3d_raw), collapse = " x "))

message(paste("Loading sample metadata (table_s1) from:", table_s1_path))
if (!file.exists(table_s1_path)) {
  stop(paste("table_s1 file not found at:", table_s1_path))
}
table_s1 <- read_csv(table_s1_path) %>%
  mutate(sample_ID = as.character(sample_ID)) # Ensure sample_ID is character

# --- Filter Count Data (similar to preprocess_relative.R) ---
# 1. Filter out doublet clusters
message(paste("Filtering out doublet clusters:", paste(DOUBLETS, collapse=", ")))
all_clusters_raw <- dimnames(counts_3d_raw)[[3]] 
clusters_to_keep <- all_clusters_raw[!all_clusters_raw %in% as.character(DOUBLETS)]
counts_3d_no_doublets <- counts_3d_raw[, , clusters_to_keep, drop = FALSE]
message("3D counts dimensions after doublet filtering: ", paste(dim(counts_3d_no_doublets), collapse = " x "))

# 2. Filter samples to match those present in table_s1
message("Filtering samples to match table_s1...")
sample_ids_from_counts <- dimnames(counts_3d_no_doublets)[[1]]
sample_ids_in_table_s1 <- table_s1 %>% pull(sample_ID) %>% unique()
samples_to_keep_in_counts <- intersect(sample_ids_from_counts, sample_ids_in_table_s1)

if(length(samples_to_keep_in_counts) == 0) {
  stop("No common samples found between counts data and table_s1 after initial filtering.")
}
counts_3d <- counts_3d_no_doublets[samples_to_keep_in_counts, , , drop = FALSE]
message("Final 3D counts dimensions for analysis: ", paste(dim(counts_3d), collapse = " x "))

# --- Process Counts for Plotting ---
# Aggregate counts across clusters to get a 2D matrix (sample x gene)
message("Aggregating filtered 3D counts to 2D (sample x gene)...")
counts_2d_sample_x_gene <- apply(counts_3d, MARGIN = c(1, 2), FUN = sum, na.rm = TRUE)
message("2D (sample x gene) counts dimensions: ", paste(dim(counts_2d_sample_x_gene), collapse = " x "))

# Transpose to get genes as rows and samples as columns, then convert to data frame
counts_data_df <- as.data.frame(t(counts_2d_sample_x_gene))
counts_data_df <- tibble::rownames_to_column(counts_data_df, var = "gene_symbol") 

gene_col_name <- "gene_symbol" 
counts_data <- counts_data_df %>% mutate(across(all_of(gene_col_name), as.character))

message(paste("Loading metabolic genes list from:", metabolic_genes_csv_path))
if (!file.exists(metabolic_genes_csv_path)) {
  stop(paste("Metabolic genes CSV file not found:", metabolic_genes_csv_path))
}
metabolic_genes_df <- read_csv(metabolic_genes_csv_path)

if (!"SYMBOL" %in% names(metabolic_genes_df)) {
  stop(paste("The file", metabolic_genes_csv_path, "does not contain a 'SYMBOL' column."))
}
metabolic_genes_list <- metabolic_genes_df %>%
  pull(SYMBOL) %>%
  unique() %>% 
  na.omit() 

if (length(metabolic_genes_list) == 0) {
  stop("Metabolic genes list is empty after processing the CSV.")
}
message(paste("Loaded", length(metabolic_genes_list), "unique metabolic gene symbols."))

# --- Calculate Proportions ---
message("Calculating metabolic gene proportions...")
counts_long <- counts_data %>%
  pivot_longer(cols = -all_of(gene_col_name), names_to = "sample", values_to = "count") %>%
  mutate(count = as.numeric(count))

total_counts_per_sample <- counts_long %>%
  group_by(sample) %>%
  summarise(total_sample_counts = sum(count, na.rm = TRUE), .groups = 'drop')

metabolic_counts_per_sample <- counts_long %>%
  filter(.data[[gene_col_name]] %in% metabolic_genes_list) %>%
  group_by(sample) %>%
  summarise(total_metabolic_counts = sum(count, na.rm = TRUE), .groups = 'drop')

proportion_data <- total_counts_per_sample %>%
  left_join(metabolic_counts_per_sample, by = "sample") %>%
  mutate(
    total_metabolic_counts = ifelse(is.na(total_metabolic_counts), 0, total_metabolic_counts), 
    proportion_metabolic = ifelse(total_sample_counts > 0, total_metabolic_counts / total_sample_counts, 0)
  )

# --- Join with Tissue Information ---
message("Joining proportion data with tissue information from table_s1...")
proportion_data_with_tissue <- proportion_data %>%
  left_join(
    table_s1 %>% 
      dplyr::select(sample_ID, tissue) %>% 
      dplyr::distinct(),
    by = c("sample" = "sample_ID")
  )

if(any(is.na(proportion_data_with_tissue$tissue))) {
  warning("Some samples in proportion_data did not have matching tissue information in table_s1. This should not happen if sample filtering was effective.")
  # print(proportion_data_with_tissue %>% filter(is.na(tissue)))
}
# Convert tissue to factor for ordered legend
proportion_data_with_tissue$tissue <- factor(proportion_data_with_tissue$tissue, levels = c("Normal", "Tumor", "Metastasis"))

message("Proportion data with tissue calculated:")
print(head(proportion_data_with_tissue))

# --- Plot Proportions ---
message("Generating plot...")
plot_title <- "Proportion of Metabolic Gene Counts per Sample by Tissue"

p <- ggplot(proportion_data_with_tissue, aes(x = reorder(sample, proportion_metabolic), y = proportion_metabolic, fill = tissue)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(name = "Tissue Type", 
                    values = c("Normal" = "skyblue", "Tumor" = "salmon", "Metastasis" = "lightgreen"),
                    na.value = "grey50") + # Should not have NAs if filtering worked
  labs(
    title = plot_title,
    x = "Sample",
    y = "Proportion of Metabolic Gene Counts"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7), # Reduced size for many samples
    plot.title = element_text(hjust = 0.5, size=14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "top"
  )

# Save the plot
ggsave(output_plot_path, plot = p, width = 14, height = 8, dpi = 300) # Increased width for more samples
message(paste("Plot saved to:", output_plot_path))

message("Script finished.")