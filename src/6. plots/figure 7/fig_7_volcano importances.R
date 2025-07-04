library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(ggrepel) # Add ggrepel library
# library(viridis) # No longer strictly needed
# library(RColorBrewer) # No longer strictly needed

# Source utility functions
source("src/0. utils/feature_name_utils.R")

# Create output directory
dir.create("output/6. plots/volcano", recursive = TRUE, showWarnings = FALSE)


# Read meta importance scores and dataset-specific fold changes
message("Loading meta importance scores.")
meta_scores_files <- list.files("output/3. intersector", 
                                pattern = "meta_scores.csv", 
                                recursive = TRUE, 
                                full.names = TRUE)

# Filter to only all_clusters/metabolic meta scores
all_clusters_metabolic_files <- meta_scores_files[grepl("all_clusters.*metabolic", meta_scores_files)]

# Load only the all_clusters/metabolic meta scores
all_meta_scores <- map_dfr(all_clusters_metabolic_files, function(file) {
  # Split path for both / and \\ (Windows)
  path_parts <- strsplit(file, "/|\\\\")[[1]]
  # Find the index of "3. intersector" to be robust
  idx <- which(path_parts == "3. intersector")
  dataset_type <- path_parts[idx + 1]
  cell_type_group <- path_parts[idx + 2]
  gene_set <- path_parts[idx + 3]
  read_csv(file, show_col_types = FALSE) %>%
    mutate(source_file = basename(dirname(file)),
           dataset_type = dataset_type,
           cell_type_group = cell_type_group,
           gene_set = gene_set)
})

# For each meta_scores file, find the corresponding fold_changes file
get_fc_path <- function(dataset_type, cell_type_group, gene_set) {
  file.path("output/4. fold changes", dataset_type, cell_type_group, gene_set, "fold_changes.csv")
}

# Load and join fold changes for each dataset
message("Loading and joining fold changes for each dataset.")
all_meta_scores_with_fc <- all_meta_scores %>%
  group_by(dataset_type, cell_type_group, gene_set) %>%
  group_modify(~ {
    fc_path <- get_fc_path(.y$dataset_type, .y$cell_type_group, .y$gene_set)
    if (file.exists(fc_path)) {
      fc <- read_csv(fc_path, show_col_types = FALSE)
      left_join(.x, fc %>% dplyr::rename(log2_fold_change = Value), by = c("feature_id" = "Feature"))
    } else {
      message("Warning: Fold change file not found for ", fc_path, ". Filling with NA.")
      .x$log2_fold_change <- NA_real_
      .x
    }
  }) %>%
  ungroup()

# Add gene and sublineage columns
combined_data <- all_meta_scores_with_fc %>%
  mutate(
    gene = get_gene_from_feature(feature_id),
    sublineage = get_simplified_sublineage(feature_id)
  )

# Load metabolic genes for classification
hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", 
                          show_col_types = FALSE)

# Filter combined data to only include metabolic genes
message("Rows in combined_data before metabolic gene filter: ", nrow(combined_data))
message("Sample combined_data before metabolic gene filter:")
print(head(combined_data, 5))
combined_data <- combined_data %>%
  filter(gene %in% hsa01100_genes$Symbol) %>%
  mutate(gene_type = "Metabolic")  # All genes are metabolic now
message("Rows in combined_data after metabolic gene filter: ", nrow(combined_data))
message("Sample combined_data after metabolic gene filter:")
print(head(combined_data, 5))

# Filter out extreme fold changes
message("Rows before fold change range filter: ", nrow(combined_data))
combined_data <- combined_data %>%
  filter(abs(log2_fold_change) < 10)
message("Rows after fold change range filter: ", nrow(combined_data))
message("Sample combined_data after fold change range filter:")
print(head(combined_data, 5))
# Remove rows with NAs in columns critical for plotting
message("Rows before NA filter: ", nrow(combined_data))
combined_data <- combined_data %>%
  filter(!is.na(log2_fold_change) & !is.na(meta_score) & !is.na(sublineage))
message("Rows after NA filter: ", nrow(combined_data))
message("Sample combined_data after NA filter:")
print(head(combined_data, 5))

# Calculate the 100th highest meta score as cutoff
meta_score_cutoff <- combined_data %>%
  arrange(desc(meta_score)) %>%
  slice(100) %>%
  pull(meta_score)

# Create significance categories
combined_data <- combined_data %>%
  mutate(
    significance = case_when(
      meta_score >= meta_score_cutoff & abs(log2_fold_change) >= 1 ~ "High importance & High FC", # Changed 0.5 to 1
      meta_score >= meta_score_cutoff & abs(log2_fold_change) < 1 ~ "High importance & Low FC",   # Changed 0.5 to 1
      meta_score < meta_score_cutoff & abs(log2_fold_change) >= 1 ~ "Low importance & High FC",    # Changed 0.5 to 1
      TRUE ~ "Low importance & Low FC"
    ),
    significance = factor(significance, levels = c(
      "High importance & High FC",
      "High importance & Low FC", 
      "Low importance & High FC",
      "Low importance & Low FC"
    ))
  )

# Create volcano plot
message("Creating volcano plot.")
p <- ggplot(combined_data, aes(x = log2_fold_change, y = meta_score)) +
  # Points colored by significance category
  geom_point(aes(color = significance), alpha = 0.7, size = 1.5) + 
  
  # Add cutoff lines
  geom_hline(yintercept = meta_score_cutoff, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate(
    "text",
    x = -Inf, y = meta_score_cutoff, hjust = -0.1, vjust = -0.5, size = 4, color = "red",
    label = sprintf("100th highest meta-score = %.2f", meta_score_cutoff)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") + # Changed xintercept to -1 and 1
  
  # Add labels for top features based on highest importance and lowest absolute fold change (0 < abs(FC) < 1) globally
  geom_text_repel(
    data = combined_data %>% 
             filter(abs(log2_fold_change) > 0 & abs(log2_fold_change) < 1) %>% # Filter for 0 < abs(FC) < 1
             arrange(desc(meta_score), abs(log2_fold_change)) %>% # Sort filtered data
             slice_head(n = 5), # Label top 5 from this filtered & sorted set
    aes(label = paste(gene, sublineage, sep = ", ")), 
    size = 3, 
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'black', 
    segment.size = 0.6,      
    min.segment.length = 0,  
    max.overlaps = Inf 
  ) +
  
  # Styling
  # Revert to scale_color_manual for significance categories
  scale_color_manual(
    name = "Significance",
    values = c(
      "High importance & High FC" = "#e74c3c",
      "High importance & Low FC" = "#f39c12", 
      "Low importance & High FC" = "#3498db",
      "Low importance & Low FC" = "#95a5a6"
    )
  ) +
  
  labs(
    x = "Log2 Fold Change (Tumor vs Normal)",
    y = "Meta Importance Score"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank()
  ) +
  
  guides(
    # Update guide for the significance color scale
    color = guide_legend(title = "Significance", override.aes = list(size = 3, alpha = 1))
  )

# Save the plot
output_path <- "output/6. plots/figure 7/volcano_importance_vs_foldchange.png"
dir.create("output/6. plots/figure 7", recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, p, width = 12, height = 8, dpi = 300, bg = "white")

# Create summary statistics
message("Creating summary statistics.")
summary_stats <- combined_data %>%
  group_by(significance) %>%
  summarise(
    count = n(),
    mean_meta_score = mean(meta_score),
    mean_fold_change = mean(abs(log2_fold_change)),
    .groups = 'drop'
  ) %>%
  arrange(desc(count))

# Top features in each quadrant
top_features_high_importance_high_fc <- combined_data %>%
  filter(significance == "High importance & High FC") %>%
  arrange(desc(meta_score)) %>%
  slice_head(n = 10) %>%
  select(feature_id, gene, sublineage, meta_score, log2_fold_change, gene_type)

top_features_high_importance_low_fc <- combined_data %>%
  filter(significance == "High importance & Low FC") %>%
  arrange(desc(meta_score)) %>%
  slice_head(n = 10) %>%
  select(feature_id, gene, sublineage, meta_score, log2_fold_change, gene_type)

# Save summary tables
write_csv(summary_stats, "output/6. plots/figure 7/volcano_summary_statistics.csv")
write_csv(top_features_high_importance_high_fc, "output/6. plots/figure 7/top_high_importance_high_fc.csv")
write_csv(top_features_high_importance_low_fc, "output/6. plots/figure 7/top_high_importance_low_fc.csv")

# Print key findings
high_imp_high_fc_count <- sum(combined_data$significance == "High importance & High FC")
total_features <- nrow(combined_data)

message("Volcano plot analysis complete.")
