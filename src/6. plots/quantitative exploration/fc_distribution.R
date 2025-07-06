library(ggplot2)
library(dplyr)
library(readr)

# Create output directory
dir.create("output/6. plots/fold changes", recursive = TRUE, showWarnings = FALSE)

# Read data
message("Reading fold change data...")
fold_changes <- read_csv("output/4. fold changes/feature_fold_changes.csv", 
                        show_col_types = FALSE)

message("Loading metabolic genes...")
hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv", 
                          show_col_types = FALSE)

# Extract gene names and classify as metabolic or non-metabolic
fold_changes <- fold_changes %>%
  mutate(
    gene = gsub("@.*", "", Feature),
    gene_type = ifelse(gene %in% hsa01100_genes$SYMBOL, "Metabolic", "Non-metabolic")
  )

message(sprintf("Classified %d metabolic and %d non-metabolic features", 
                sum(fold_changes$gene_type == "Metabolic"),
                sum(fold_changes$gene_type == "Non-metabolic")))

# Calculate statistics by gene type
fc_stats_by_type <- fold_changes %>%
  group_by(gene_type) %>%
  summarise(
    count = n(),
    mean_fc = mean(Value, na.rm = TRUE),
    median_fc = median(Value, na.rm = TRUE),
    sd_fc = sd(Value, na.rm = TRUE),
    .groups = 'drop'
  )

message("Fold change statistics by gene type:")
print(fc_stats_by_type)

# Create histogram
p_hist <- ggplot(fold_changes, aes(x = Value, fill = gene_type)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "white") +
  facet_wrap(~gene_type, scales = "free_y") +
  scale_fill_manual(values = c("Metabolic" = "#e74c3c", "Non-metabolic" = "#3498db"),
                    name = "Gene type") +
  labs(
    x = "Log2 Fold Change (Tumor vs Normal)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 11)
  )

# Save results
ggsave("output/6. plots/fold changes/fold_change_histogram.png", 
       p_hist, width = 10, height = 6, dpi = 300, bg = "white")

write_csv(fc_stats_by_type, "output/6. plots/fold changes/fold_change_statistics_by_gene_type.csv")

message("Files saved: fold_change_histogram.png, fold_change_statistics_by_gene_type.csv")
