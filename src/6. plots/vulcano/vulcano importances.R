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

# Read meta importance scores and fold changes
message("Loading meta importance scores.")
meta_scores_files <- list.files("output/3. intersector", 
                                pattern = "meta_scores.csv", 
                                recursive = TRUE, 
                                full.names = TRUE)

# Filter to only all_clusters/metabolic meta scores
all_clusters_metabolic_files <- meta_scores_files[grepl("all_clusters.*metabolic", meta_scores_files)]

# Load only the all_clusters/metabolic meta scores
all_meta_scores <- map_dfr(all_clusters_metabolic_files, function(file) {
  read_csv(file, show_col_types = FALSE) %>%
    mutate(source_file = basename(dirname(file)))
})

message("Loading fold changes.")
fold_changes <- read_csv("output/4. fold changes/feature_fold_changes.csv", 
                        show_col_types = FALSE)

# Prepare fold_changes data by renaming 'Value' column before the join
fold_changes_processed <- fold_changes %>%
  dplyr::rename(log2_fold_change = "Value") # Explicit dplyr::rename and quote "Value"

# Combine meta scores with fold changes
combined_data <- all_meta_scores %>%
  inner_join(fold_changes_processed, by = c("feature_id" = "Feature")) %>%
  # 'Value' column is now 'log2_fold_change' from the join, so no further rename needed here for that column.
  # Add gene classification
  mutate(
    gene = get_gene_from_feature(feature_id),
    sublineage = get_simplified_sublineage(feature_id)
  )

# Load metabolic genes for classification
hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", 
                          show_col_types = FALSE)

# Filter combined data to only include metabolic genes
combined_data <- combined_data %>%
  filter(gene %in% hsa01100_genes$SYMBOL) %>%
  mutate(gene_type = "Metabolic")  # All genes are metabolic now

# Remove rows with NAs in columns critical for plotting
combined_data <- combined_data %>%
  filter(!is.na(log2_fold_change) & !is.na(meta_score) & !is.na(sublineage))

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
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank()
  ) +
  
  guides(
    # Update guide for the significance color scale
    color = guide_legend(title = "Significance", override.aes = list(size = 3, alpha = 1))
  )

# Save the plot
output_path <- "output/6. plots/volcano/volcano_importance_vs_foldchange.png"
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
write_csv(summary_stats, "output/6. plots/volcano/volcano_summary_statistics.csv")
write_csv(top_features_high_importance_high_fc, "output/6. plots/volcano/top_high_importance_high_fc.csv")
write_csv(top_features_high_importance_low_fc, "output/6. plots/volcano/top_high_importance_low_fc.csv")

# Print key findings
high_imp_high_fc_count <- sum(combined_data$significance == "High importance & High FC")
total_features <- nrow(combined_data)

message("Volcano plot analysis complete.")
