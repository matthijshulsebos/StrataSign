library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(scales)
library(ggplot2)

DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)
GENE_SET_TYPES <- c("metabolic", "nonmetabolic", "random")
BASE_OUTPUT_DIR <- file.path("output", "1. data preprocessing", "training datasets", "absolute") 
CLUSTER_DEFINITIONS <- list(
  lcam_hi = c(44, 9, 17, 28, 46, 11, 42),
  lcam_lo = c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16),
  lcam_both = c(c(44, 9, 17, 28, 46, 11, 42), c(56, 34, 53, 10, 25, 54, 55, 57, 45, 14, 16)),
  macrophages = c(5, 8, 10, 11, 25, 32, 33, 35, 38, 42, 47, 54, 55, 57),
  all_clusters = NULL
)

# Load the Leader dataset
message("Loading Leader et al. data...")
if (!exists("lung_ldm")) {
  load("base/data/lung_ldm.rd")
}

# Load additional datasets
table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv") %>% 
  mutate(sample_ID = as.character(sample_ID))
annots_list <- read_csv("base/input_tables/annots_list.csv")
hsa01100_genes <- read_csv("output/1. data preprocessing/kegg/hsa01100_genes.csv")

# 3D matrix: sample × gene × cluster
counts <- lung_ldm$dataset$counts

# Filter out doublets
all_clusters_initial <- as.numeric(dimnames(counts)[[3]])
clusters_to_keep_initial <- all_clusters_initial[!all_clusters_initial %in% DOUBLETS]
counts_filtered <- counts[, , as.character(clusters_to_keep_initial), drop = FALSE]

# Filter counts to match samples present in table_s1
sample_ids_to_keep <- table_s1 %>% pull(sample_ID)
counts_filtered <- counts_filtered[dimnames(counts_filtered)[[1]] %in% sample_ids_to_keep, , , drop = FALSE]

# Convert filtered 3D array to long format
counts_long <- as.data.frame(as.table(counts_filtered))
colnames(counts_long) <- c("sample_ID", "gene", "cluster_ID", "count")

counts_long <- counts_long %>%
  mutate(
    sample_ID = as.character(sample_ID),
    cluster_ID = as.numeric(as.character(cluster_ID)),
    gene = as.character(gene)
  ) %>%
  filter(count > 0) # Filter zero counts for efficiency

# Calculate total counts per cell type per sample
total_counts_celltype_sample <- counts_long %>%
  group_by(sample_ID, cluster_ID) %>%
  summarise(total_counts_in_celltype = sum(count), .groups = 'drop')

# Calculate the global average of cell type totals
global_avg_celltype_total <- total_counts_celltype_sample %>%
  summarise(mean_total = mean(total_counts_in_celltype)) %>%
  pull(mean_total)

# Calculate gene fraction and normalized counts
counts_normalized <- counts_long %>%
  left_join(total_counts_celltype_sample, by = c("sample_ID", "cluster_ID")) %>%
  mutate(
    gene_fraction_in_celltype = ifelse(total_counts_in_celltype > 0, count / total_counts_in_celltype, 0),
    normalized_count = gene_fraction_in_celltype * global_avg_celltype_total
  ) %>%
  select(sample_ID, gene, cluster_ID, normalized_count)

# Log transform
counts_log_transformed <- counts_normalized %>%
  mutate(log_normalized_count = log1p(normalized_count))

# Create feature names
counts_processed_long <- counts_log_transformed %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>% 
  mutate(
    cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
    cluster_name = gsub("/", "-", cluster_name),
    cluster_identifier = paste(cluster_name, cluster_ID, sep = "_"),
    gene_cluster = paste(gene, cluster_identifier, sep = "@") 
  ) %>%
  select(sample_ID, gene, cluster_ID, gene_cluster, log_normalized_count)

# Process a specific subset defined by cluster and gene type
process_and_save_subset <- function(globally_processed_data_long, 
                                    cluster_subset_ids,
                                    subset_name,
                                    gene_set_type,
                                    hsa01100_genes_df,
                                    table_s1_df,
                                    base_output_path) {

  message(paste("Processing Subset:", subset_name, "| Gene Set:", gene_set_type))

  # Filter by cluster subset
  if (!is.null(cluster_subset_ids)) {
    subset_data_long <- globally_processed_data_long %>% filter(cluster_ID %in% cluster_subset_ids)
  } else {
    subset_data_long <- globally_processed_data_long
  }

  # Filter by gene set type
  all_genes_in_subset <- unique(subset_data_long$gene)
  metabolic_genes <- hsa01100_genes_df$SYMBOL[hsa01100_genes_df$SYMBOL %in% all_genes_in_subset]

  # Define suffix based on gene type
  if (gene_set_type == "metabolic") {
    genes_to_keep <- metabolic_genes
    suffix <- "_metabolic"
  } else if (gene_set_type == "nonmetabolic") {
    nonmetabolic_genes <- setdiff(all_genes_in_subset, metabolic_genes) # Filter from genes *in the subset*
     if (length(nonmetabolic_genes) > length(metabolic_genes)) {
       set.seed(42)
       genes_to_keep <- sample(nonmetabolic_genes, length(metabolic_genes))
     } else {
        genes_to_keep <- nonmetabolic_genes
     }
    suffix <- "_nonmetabolic"
  } else if (gene_set_type == "random") {
    set.seed(43)
    sample_size <- min(length(metabolic_genes), length(all_genes_in_subset))
    genes_to_keep <- sample(all_genes_in_subset, sample_size)
    suffix <- "_random"
  } else {
    stop("Invalid gene set type.")
  }

  # Apply gene filter
  subset_data_long <- subset_data_long %>% filter(gene %in% genes_to_keep)

  # Check if filtering resulted in empty data
  if (nrow(subset_data_long) == 0) {
    warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after filtering."))
    return()
  }

  # Pivot to wide format
  message("Pivoting data to wide format...")
  counts_wide <- subset_data_long %>%
    # Use the log_normalized_count column created globally
    select(sample_ID, gene_cluster, log_normalized_count) %>%
    pivot_wider(names_from = gene_cluster,
                values_from = log_normalized_count,
                values_fill = 0)

  # Merge metadata and split train/test
  counts_wide_meta <- counts_wide %>%
    left_join(
      table_s1_df %>% select(sample_ID, patient_ID, tissue),
      by = "sample_ID"
    ) %>%
    filter(!is.na(tissue))

  # Check if merge resulted in empty data
  if (nrow(counts_wide_meta) == 0) {
    warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in zero rows after merging metadata/filtering NAs."))
    return()
  }

  set.seed(42)
  patients <- unique(counts_wide_meta$patient_ID)
  train_patients <- sample(patients, size = round(length(patients) * 0.7))

  train_df <- counts_wide_meta %>% filter(patient_ID %in% train_patients)
  test_df <- counts_wide_meta %>% filter(!patient_ID %in% train_patients)

  if (nrow(train_df) == 0 || nrow(test_df) == 0) {
     warning(paste("Subset", subset_name, "|", gene_set_type, "resulted in empty train/test set. Skipping."))
     return()
  }

  y_train <- train_df$tissue
  y_test <- test_df$tissue

  # Prepare metadata for saving
  metadata_train <- table_s1_df %>%
    filter(sample_ID %in% train_df$sample_ID) %>%
    arrange(match(sample_ID, train_df$sample_ID))

  # Prepare feature matrices (remove metadata)
  X_train <- train_df %>% select(-sample_ID, -patient_ID, -tissue)
  X_test <- test_df %>% select(-sample_ID, -patient_ID, -tissue)

  # Ensure numeric and handle NAs (should be 0 from pivot)
  X_train <- X_train %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))
  X_test <- X_test %>% mutate(across(everything(), ~replace_na(as.numeric(.), 0)))

  # Save outputs
  output_dir <- file.path(base_output_path, subset_name, gene_set_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Use the suffix defined earlier in the function
  write_csv(metadata_train, file.path(output_dir, paste0("metadata_", subset_name, suffix, ".csv")))
  write_csv(X_train, file.path(output_dir, paste0("X_train_", subset_name, suffix, ".csv")))
  write_csv(X_test, file.path(output_dir, paste0("X_test_", subset_name, suffix, ".csv")))
  write_csv(data.frame(x = y_train), file.path(output_dir, paste0("y_train_", subset_name, suffix, ".csv")))
  write_csv(data.frame(x = y_test), file.path(output_dir, paste0("y_test_", subset_name, suffix, ".csv")))

  message(paste("Finished:", subset_name, "|", gene_set_type))
}

# Main Loop
message("\n Preprocessing all defined subsets...")
for (subset_name in names(CLUSTER_DEFINITIONS)) {
  cluster_ids <- CLUSTER_DEFINITIONS[[subset_name]]
  
  for (gene_type in GENE_SET_TYPES) {
    
    # Call the processing function for the current combination
    process_and_save_subset(
      globally_processed_data_long = counts_processed_long, 
      cluster_subset_ids = cluster_ids,
      subset_name = subset_name,
      gene_set_type = gene_type,
      hsa01100_genes_df = hsa01100_genes,
      table_s1_df = table_s1,
      base_output_path = BASE_OUTPUT_DIR
    )
    
  }
}

message("\n All preprocessing complete!")


# --- Configuration ---
counts_file_path <- "path/to/your/counts_dataset.csv" # Replace with your actual counts file path
# Example: metabolic_genes_list <- c("PDHA1", "ACO2", "IDH1", "SUCLG1", "HK1", "PFKFB3") 
# Or load from a file:
metabolic_genes_file_path <- "path/to/your/metabolic_genes_list.txt" # Replace if loading from file, one gene per line
output_plot_path <- "output/metabolic_gene_proportions.png" # Define where to save the plot
figure_dir <- "output/1. preprocessing_qc" # Define a directory for plots

# Create directory if it doesn't exist
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
output_plot_path <- file.path(figure_dir, "metabolic_gene_proportions.png")


# --- Load Data ---
message(paste("Loading counts data from:", counts_file_path))
counts_data <- read_csv(counts_file_path) # Or read_tsv if tab-separated

# Assuming the first column is gene names, and subsequent columns are samples
# Let's rename the first column to "gene" for clarity if it's not already named appropriately
# If your gene column has a different name, adjust "GENE_COLUMN_NAME"
# counts_data <- counts_data %>% rename(gene = GENE_COLUMN_NAME) 

# Ensure gene names are character
counts_data <- counts_data %>% mutate(across(1, as.character)) 
gene_col_name <- names(counts_data)[1]


message("Loading metabolic genes list...")
# Option 1: Define directly
# metabolic_genes_list <- c("PDHA1", "ACO2", "IDH1", "SUCLG1", "HK1", "PFKFB3") # Make sure these match gene names in counts_data

# Option 2: Load from a file (one gene per line)
if (file.exists(metabolic_genes_file_path)) {
  metabolic_genes_list <- read_lines(metabolic_genes_file_path)
  metabolic_genes_list <- metabolic_genes_list[metabolic_genes_list != ""] # Remove empty lines
} else {
  stop(paste("Metabolic genes file not found:", metabolic_genes_file_path))
  # Or define a default list if the file is optional:
  # metabolic_genes_list <- c("GENE1", "GENE2") 
}

message(paste("Loaded", length(metabolic_genes_list), "metabolic genes."))

# --- Preprocess and Calculate Proportions ---

# Convert to long format for easier processing: gene, sample, count
counts_long <- counts_data %>%
  pivot_longer(cols = -all_of(gene_col_name), names_to = "sample", values_to = "count")

# Calculate total counts per sample
total_counts_per_sample <- counts_long %>%
  group_by(sample) %>%
  summarise(total_sample_counts = sum(count, na.rm = TRUE), .groups = 'drop')

# Filter for metabolic genes and calculate their total counts per sample
metabolic_counts_per_sample <- counts_long %>%
  filter(.data[[gene_col_name]] %in% metabolic_genes_list) %>%
  group_by(sample) %>%
  summarise(total_metabolic_counts = sum(count, na.rm = TRUE), .groups = 'drop')

# Join total and metabolic counts, then calculate proportion
proportion_data <- total_counts_per_sample %>%
  left_join(metabolic_counts_per_sample, by = "sample") %>%
  mutate(
    total_metabolic_counts = ifelse(is.na(total_metabolic_counts), 0, total_metabolic_counts), # Handle samples with no metabolic genes found
    proportion_metabolic = ifelse(total_sample_counts > 0, total_metabolic_counts / total_sample_counts, 0)
  )

message("Proportion data calculated:")
print(head(proportion_data))

# --- Plot Proportions ---
plot_title <- "Proportion of Metabolic Gene Counts per Sample"

p <- ggplot(proportion_data, aes(x = reorder(sample, -proportion_metabolic), y = proportion_metabolic)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = plot_title,
    x = "Sample",
    y = "Proportion of Metabolic Gene Counts"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5)
  )

# Save the plot
ggsave(output_plot_path, plot = p, width = 10, height = 6, dpi = 300)
message(paste("Plot saved to:", output_plot_path))

# Display the plot (optional, if running interactively)
# print(p)
