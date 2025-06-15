library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(ggrepel) # Add ggrepel for labeling
# library(RColorBrewer) # No longer strictly needed if using scale_color_hue

# Source utility functions
source("src/0. utils/feature_name_utils.R")

# Create output directory
dir.create("output/6. plots/quadrant", recursive = TRUE, showWarnings = FALSE)

# Read meta importance scores and fold changes
message("Loading meta importance scores.")
meta_scores_files <- list.files("output/3. intersector", 
                                pattern = "meta_scores.csv", 
                                recursive = TRUE, 
                                full.names = TRUE)

if (length(meta_scores_files) == 0) {
  stop("No meta_scores.csv files found in output/3. intersector")
}

# Filter to only all_clusters and metabolic meta scores
message("Filtering to all_clusters/metabolic meta scores.")
all_clusters_metabolic_files <- meta_scores_files[grepl("all_clusters.*metabolic", meta_scores_files)]

# Load only the all_clusters/metabolic meta scores
all_meta_scores <- map_dfr(all_clusters_metabolic_files, function(file) {
  read_csv(file, show_col_types = FALSE) %>%
    mutate(source_file = basename(dirname(file)))
})

message("Loading fold changes.")
fold_changes <- read_csv("output/4. fold changes/feature_fold_changes.csv", 
                        show_col_types = FALSE)

# Combine meta scores with fold changes
message("Combining meta scores with fold changes.")
combined_data <- all_meta_scores %>%
  inner_join(fold_changes, by = c("feature_id" = "Feature")) %>%
  rename(log2_fold_change = Value) %>%
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

# Create quadrant categories (requires abs_log2_fold_change)
combined_data <- combined_data %>%
  mutate(
    abs_log2_fold_change = abs(log2_fold_change)
  )

# Remove rows with NAs in columns critical for plotting before quadrant assignment and plotting
combined_data <- combined_data %>%
  filter(!is.na(abs_log2_fold_change) & !is.na(meta_score) & !is.na(sublineage) & !is.na(log2_fold_change))



meta_score_cutoff <- combined_data %>%
  arrange(desc(meta_score)) %>%
  slice(100) %>%
  pull(meta_score)

# Create log2 fold change cutoff
fold_change_cutoff <- 1.0

# Create quadrant categories
combined_data <- combined_data %>%
  mutate(
    abs_log2_fold_change = abs(log2_fold_change),
    quadrant = case_when(
      meta_score >= meta_score_cutoff & abs_log2_fold_change >= fold_change_cutoff ~ "High Importance & High FC",
      meta_score >= meta_score_cutoff & abs_log2_fold_change < fold_change_cutoff ~ "High Importance & Low FC",
      meta_score < meta_score_cutoff & abs_log2_fold_change >= fold_change_cutoff ~ "Low Importance & High FC",
      TRUE ~ "Low Importance & Low FC"
    ),
    quadrant = factor(quadrant, levels = c(
      "High Importance & High FC",
      "High Importance & Low FC", 
      "Low Importance & High FC",
      "Low Importance & Low FC"
    ))
  )

# Calculate quadrant boundaries for visualization
x_min <- 0
x_max <- max(combined_data$abs_log2_fold_change) + 0.1
y_min <- min(combined_data$meta_score) - 0.1
y_max <- max(combined_data$meta_score) + 0.1

# Create true quadrant rectangles
quadrant_rects <- data.frame(
  xmin = c(0, fold_change_cutoff, 0, fold_change_cutoff),
  xmax = c(fold_change_cutoff, x_max, fold_change_cutoff, x_max),
  ymin = c(meta_score_cutoff, meta_score_cutoff, y_min, y_min),
  ymax = c(y_max, y_max, meta_score_cutoff, meta_score_cutoff),
  quadrant = c("High Importance & Low FC", "High Importance & High FC", 
               "Low Importance & Low FC", "Low Importance & High FC"),
  stringsAsFactors = FALSE
)

# Create the quadrant plot
message("Creating quadrant plot.")
p <- ggplot() +
  # Add quadrant rectangles
  geom_rect(data = quadrant_rects, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = quadrant),
            alpha = 0.1, color = "black", size = 0.5) +
  
  # Add points for features, colored by sublineage (cell type)
  geom_point(data = combined_data, 
             aes(x = abs_log2_fold_change, y = meta_score, color = sublineage), # Changed color to sublineage
             alpha = 0.7, size = 1.5) +
  
  # Add labels for top features based on highest importance and lowest absolute fold change (0 < abs(FC) < 1) globally
  ggrepel::geom_text_repel(
    data = combined_data %>%
             filter(abs(log2_fold_change) > 0 & abs(log2_fold_change) < 1) %>% # Filter for 0 < abs(FC) < 1
             arrange(desc(meta_score), abs(log2_fold_change)) %>% # Sort filtered data
             slice_head(n = 5), # Label top 5 from this filtered & sorted set
    aes(x = abs_log2_fold_change, y = meta_score, label = paste(gene, sublineage, sep = ", ")), 
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = 'grey50',
    min.segment.length = 0 
  ) +
  
  # Add cutoff lines
  geom_hline(yintercept = meta_score_cutoff, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = fold_change_cutoff, linetype = "dashed", color = "black", size = 1) +
  
  # Color scheme
  scale_fill_manual(
    name = "Quadrant",
    values = c(
      "High Importance & High FC" = "#e74c3c",
      "High Importance & Low FC" = "#f39c12", 
      "Low Importance & High FC" = "#3498db",
      "Low Importance & Low FC" = "#95a5a6"
    ),
    guide = "legend" 
  ) +
  
  # Add a color scale for sublineage (cell type) using ggplot2's default hue scale
  scale_color_hue(name = "Cell Type") + 
  
  labs(
    x = "Absolute Log2 Fold Change",
    y = "Meta Importance Score"
  ) +
  
  theme_bw() + # Changed to theme_bw()
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "lightgray", size = 0.2)
  ) +
  
  guides(
    fill = guide_legend(title = "Quadrant", override.aes = list(alpha = 0.5)), 
    color = guide_legend(title = "Cell Type", override.aes = list(alpha = 1, size = 3)) # Guide for cell type color
  )

# Save the plot
output_path <- "output/6. plots/quadrant/quadrant_importance_vs_foldchange.png"
ggsave(output_path, p, width = 12, height = 8, dpi = 300, bg = "white") # Changed bg to white

# Create summary statistics
message("Creating summary statistics.")
summary_stats <- combined_data %>%
  group_by(quadrant) %>%
  summarise(
    count = n(),
    mean_meta_score = mean(meta_score),
    mean_abs_fold_change = mean(abs_log2_fold_change),
    .groups = 'drop'
  ) %>%
  arrange(desc(count))

# Top features in each quadrant
top_features_per_quadrant <- combined_data %>%
  group_by(quadrant) %>%
  arrange(desc(meta_score)) %>%
  slice_head(n = 5) %>%
  select(quadrant, feature_id, gene, sublineage, meta_score, log2_fold_change) %>%
  ungroup()

# Save summary tables
write_csv(summary_stats, "output/6. plots/quadrant/quadrant_summary_statistics.csv")
write_csv(top_features_per_quadrant, "output/6. plots/quadrant/top_features_per_quadrant.csv")

message("Quadrant plot analysis complete.")
