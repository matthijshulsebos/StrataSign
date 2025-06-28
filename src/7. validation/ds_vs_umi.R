# Compare downsampled dataset (ds) vs UMI data
# Show cell counts per cell type per sample and totals

library(dplyr)
library(readr)
library(Matrix)
library(knitr)
library(ggplot2)
library(tidyr)

# Load data loader functions
source("src/1. data preprocessing/training datasets/data_loader.R")

# Configuration
DOUBLETS <- c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)

# Load datasets
message("Loading lung_ldm dataset...")
datasets <- load_all_datasets()
lung_ldm <- datasets$lung_ldm
table_s1 <- select_representative_samples(datasets$table_s1)
annots_list <- datasets$annots_list

# Extract the two matrices for comparison
message("Extracting matrices for comparison...")

# DS (downsampled) matrix - this is what extract_downsampled_umitab() extracts
ds_matrix <- lung_ldm$dataset$ds[[1]]
message(sprintf("DS matrix dimensions: %d genes x %d cells", nrow(ds_matrix), ncol(ds_matrix)))

# UMI matrix (original counts) - this is the full UMI count matrix before downsampling
umitab_matrix <- lung_ldm$dataset$umitab
message(sprintf("UMI matrix dimensions: %d genes x %d cells", nrow(umitab_matrix), ncol(umitab_matrix)))

# Prepare cell metadata (same filtering as preprocessing pipeline)
message("Preparing cell metadata...")
cell_metadata <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)

# Filter cell metadata to only cells present in ds_matrix
cell_metadata_ds <- cell_metadata %>% filter(cell_ID %in% colnames(ds_matrix))
message(sprintf("Cell metadata after filtering: %d cells", nrow(cell_metadata_ds)))

# Join with sample metadata and annotations
cell_metadata_complete <- cell_metadata_ds %>%
  left_join(table_s1 %>% select(sample_ID, tissue, patient_ID), by = "sample_ID") %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cell_type = ifelse(!is.na(sub_lineage) & sub_lineage != "", sub_lineage, lineage),
    cell_type = ifelse(is.na(cell_type), paste0("Cluster_", cluster_ID), cell_type)
  )

# ===== ANALYSIS 1: DS DATASET CELL COUNTS =====
message("\n=== ANALYZING DS DATASET ===")

# Count cells per cell type per sample in DS dataset
ds_cell_counts <- cell_metadata_complete %>%
  group_by(sample_ID, tissue, patient_ID, cluster_ID, cell_type) %>%
  summarise(n_cells_ds = n(), .groups = 'drop') %>%
  arrange(sample_ID, cluster_ID)

# Summary per cell type (total across all samples)
ds_celltype_totals <- ds_cell_counts %>%
  group_by(cluster_ID, cell_type) %>%
  summarise(
    total_cells_ds = sum(n_cells_ds),
    n_samples_present = n(),
    avg_cells_per_sample = mean(n_cells_ds),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_cells_ds))

# Summary per sample (total across all cell types)  
ds_sample_totals <- ds_cell_counts %>%
  group_by(sample_ID, tissue, patient_ID) %>%
  summarise(
    total_cells_ds = sum(n_cells_ds),
    n_cell_types = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_cells_ds))

# ===== ANALYSIS 2: UMI DATASET CELL COUNTS =====
message("\n=== ANALYZING UMI DATASET (FULL UMITAB) ===")

# The UMI dataset contains all cells before downsampling
# Get cell metadata for all cells in umitab matrix
message("Preparing UMI cell metadata...")
cell_metadata_umi_full <- prepare_cell_metadata(lung_ldm, table_s1, DOUBLETS)

# Filter to cells actually present in umitab matrix
cell_metadata_umi <- cell_metadata_umi_full %>% filter(cell_ID %in% colnames(umitab_matrix))
message(sprintf("UMI matrix cells: %d", ncol(umitab_matrix)))
message(sprintf("UMI cell metadata: %d cells", nrow(cell_metadata_umi)))
message(sprintf("DS matrix cells: %d", ncol(ds_matrix)))
message(sprintf("Cells removed in downsampling: %d", ncol(umitab_matrix) - ncol(ds_matrix)))

# Join with sample metadata and annotations for UMI dataset
cell_metadata_umi_complete <- cell_metadata_umi %>%
  left_join(table_s1 %>% select(sample_ID, tissue, patient_ID), by = "sample_ID") %>%
  left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
  mutate(
    cell_type = ifelse(!is.na(sub_lineage) & sub_lineage != "", sub_lineage, lineage),
    cell_type = ifelse(is.na(cell_type), paste0("Cluster_", cluster_ID), cell_type)
  )

# Count cells per cell type per sample in UMI dataset
umi_cell_counts <- cell_metadata_umi_complete %>%
  group_by(sample_ID, tissue, patient_ID, cluster_ID, cell_type) %>%
  summarise(n_cells_umi = n(), .groups = 'drop') %>%
  arrange(sample_ID, cluster_ID)

# Summary per cell type (total across all samples) for UMI
umi_celltype_totals <- umi_cell_counts %>%
  group_by(cluster_ID, cell_type) %>%
  summarise(
    total_cells_umi = sum(n_cells_umi),
    n_samples_present = n(),
    avg_cells_per_sample = mean(n_cells_umi),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_cells_umi))

# Summary per sample (total across all cell types) for UMI
umi_sample_totals <- umi_cell_counts %>%
  group_by(sample_ID, tissue, patient_ID) %>%
  summarise(
    total_cells_umi = sum(n_cells_umi),
    n_cell_types = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_cells_umi))

# ===== COMBINED COMPARISON =====
message("\n=== CREATING COMPARISON TABLES ===")

# Compare cell type totals
celltype_comparison <- ds_celltype_totals %>%
  full_join(umi_celltype_totals, by = c("cluster_ID", "cell_type"), suffix = c("_ds", "_umi")) %>%
  mutate(
    total_cells_ds = ifelse(is.na(total_cells_ds), 0, total_cells_ds),
    total_cells_umi = ifelse(is.na(total_cells_umi), 0, total_cells_umi),
    cells_removed_in_ds = total_cells_umi - total_cells_ds,
    pct_retained_in_ds = ifelse(total_cells_umi > 0, (total_cells_ds / total_cells_umi) * 100, 0)
  ) %>%
  arrange(desc(total_cells_umi))

# Compare sample totals  
sample_comparison <- ds_sample_totals %>%
  full_join(umi_sample_totals, by = c("sample_ID", "tissue", "patient_ID"), suffix = c("_ds", "_umi")) %>%
  mutate(
    total_cells_ds = ifelse(is.na(total_cells_ds), 0, total_cells_ds),
    total_cells_umi = ifelse(is.na(total_cells_umi), 0, total_cells_umi),
    cells_removed_in_ds = total_cells_umi - total_cells_ds,
    pct_retained_in_ds = ifelse(total_cells_umi > 0, (total_cells_ds / total_cells_umi) * 100, 0)
  ) %>%
  arrange(desc(total_cells_umi))

# ===== CELL TYPE COMPARISON ANALYSIS =====
message("\n=== ANALYZING CELL TYPE OVERLAP BETWEEN DATASETS ===")

# Get unique cell types from each dataset
ds_cell_types <- unique(ds_cell_counts$cell_type)
umi_cell_types <- unique(umi_cell_counts$cell_type)

# Find overlaps and differences
cell_types_in_both <- intersect(ds_cell_types, umi_cell_types)
cell_types_only_in_ds <- setdiff(ds_cell_types, umi_cell_types)
cell_types_only_in_umi <- setdiff(umi_cell_types, ds_cell_types)

# Create summary table
cell_type_overlap_summary <- data.frame(
  metric = c("Total cell types in DS", "Total cell types in UMI", 
             "Cell types in both datasets", "Cell types only in DS", 
             "Cell types only in UMI", "Overlap percentage"),
  count = c(length(ds_cell_types), length(umi_cell_types), 
            length(cell_types_in_both), length(cell_types_only_in_ds),
            length(cell_types_only_in_umi), 
            round((length(cell_types_in_both) / length(umi_cell_types)) * 100, 1))
)

cat("CELL TYPE OVERLAP ANALYSIS:\n")
print(knitr::kable(cell_type_overlap_summary, format = "pipe"))

if (length(cell_types_only_in_umi) > 0) {
  cat("\nCell types only in UMI dataset:\n")
  cat(paste(cell_types_only_in_umi, collapse = ", "))
  cat("\n")
}

if (length(cell_types_only_in_ds) > 0) {
  cat("\nCell types only in DS dataset:\n")
  cat(paste(cell_types_only_in_ds, collapse = ", "))
  cat("\n")
}

# ===== NORMAL VS TUMOR ANALYSIS =====
message("\n=== ANALYZING NORMAL VS TUMOR CELL COUNTS ===")

# First, check what tissue types we actually have
tissue_types <- unique(c(umi_cell_counts$tissue, ds_cell_counts$tissue))
message(sprintf("Tissue types found: %s", paste(tissue_types, collapse = ", ")))

# UMI dataset: Normal vs Tumor counts per cell type
umi_normal_tumor <- umi_cell_counts %>%
  group_by(cell_type, tissue) %>%
  summarise(total_cells = sum(n_cells_umi), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = total_cells, values_fill = 0)

# Get the actual column names after pivot_wider
tissue_cols <- setdiff(names(umi_normal_tumor), "cell_type")
message(sprintf("Column names after pivot_wider: %s", paste(tissue_cols, collapse = ", ")))

# Add calculations using the actual column names
if (length(tissue_cols) >= 2) {
  col1 <- tissue_cols[1]
  col2 <- tissue_cols[2]
  
  umi_normal_tumor <- umi_normal_tumor %>%
    mutate(
      total_across_tissues = .data[[col1]] + .data[[col2]],
      pct_tissue1 = ifelse(total_across_tissues > 0, (.data[[col1]] / total_across_tissues) * 100, 0),
      pct_tissue2 = ifelse(total_across_tissues > 0, (.data[[col2]] / total_across_tissues) * 100, 0),
      fold_change = ifelse(.data[[col1]] > 0, .data[[col2]] / .data[[col1]], ifelse(.data[[col2]] > 0, Inf, 0))
    ) %>%
    arrange(desc(total_across_tissues))
  
  # Rename columns to be more descriptive
  names(umi_normal_tumor)[names(umi_normal_tumor) == "pct_tissue1"] <- paste0("pct_", col1)
  names(umi_normal_tumor)[names(umi_normal_tumor) == "pct_tissue2"] <- paste0("pct_", col2)
  names(umi_normal_tumor)[names(umi_normal_tumor) == "fold_change"] <- paste0("fold_change_", col2, "_vs_", col1)
}

# DS dataset: Normal vs Tumor counts per cell type  
ds_normal_tumor <- ds_cell_counts %>%
  group_by(cell_type, tissue) %>%
  summarise(total_cells = sum(n_cells_ds), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = total_cells, values_fill = 0)

# Add calculations for DS dataset using the same approach
if (length(tissue_cols) >= 2) {
  col1 <- tissue_cols[1]
  col2 <- tissue_cols[2]
  
  ds_normal_tumor <- ds_normal_tumor %>%
    mutate(
      total_across_tissues = .data[[col1]] + .data[[col2]],
      pct_tissue1 = ifelse(total_across_tissues > 0, (.data[[col1]] / total_across_tissues) * 100, 0),
      pct_tissue2 = ifelse(total_across_tissues > 0, (.data[[col2]] / total_across_tissues) * 100, 0),
      fold_change = ifelse(.data[[col1]] > 0, .data[[col2]] / .data[[col1]], ifelse(.data[[col2]] > 0, Inf, 0))
    ) %>%
    arrange(desc(total_across_tissues))
  
  # Rename columns to be more descriptive
  names(ds_normal_tumor)[names(ds_normal_tumor) == "pct_tissue1"] <- paste0("pct_", col1)
  names(ds_normal_tumor)[names(ds_normal_tumor) == "pct_tissue2"] <- paste0("pct_", col2)
  names(ds_normal_tumor)[names(ds_normal_tumor) == "fold_change"] <- paste0("fold_change_", col2, "_vs_", col1)
}

cat("NORMAL VS TUMOR ANALYSIS - UMI DATASET (Top 15 cell types):\n")
if (length(tissue_cols) >= 2) {
  col1 <- tissue_cols[1]
  col2 <- tissue_cols[2]
  display_cols <- c("cell_type", col1, col2, "total_across_tissues", paste0("pct_", col1), paste0("pct_", col2))
  print(knitr::kable(head(umi_normal_tumor %>% select(all_of(display_cols)), 15),
                     format = "pipe", digits = 1))
} else {
  print("Not enough tissue types for comparison")
}

cat("\n\nNORMAL VS TUMOR ANALYSIS - DS DATASET (Top 15 cell types):\n")
if (length(tissue_cols) >= 2) {
  col1 <- tissue_cols[1]
  col2 <- tissue_cols[2]
  display_cols <- c("cell_type", col1, col2, "total_across_tissues", paste0("pct_", col1), paste0("pct_", col2))
  print(knitr::kable(head(ds_normal_tumor %>% select(all_of(display_cols)), 15),
                     format = "pipe", digits = 1))
} else {
  print("Not enough tissue types for comparison")
}

# ===== COMBINED CELLS PER CELLTYPE PER SAMPLE =====
message("\n=== CREATING COMBINED CELLS PER CELLTYPE PER SAMPLE TABLE ===")

# Create combined table with DS cell counts and UMI cell counts per sample per cell type
combined_cells_per_sample_celltype <- ds_cell_counts %>%
  # Start with DS data (sample_ID, tissue, patient_ID, cluster_ID, cell_type, n_cells_ds)
  full_join(
    umi_cell_counts %>% 
      select(sample_ID, tissue, patient_ID, cluster_ID, cell_type, n_cells_umi),
    by = c("sample_ID", "tissue", "patient_ID", "cluster_ID", "cell_type")
  ) %>%
  # Fill missing values with 0
  mutate(
    n_cells_ds = ifelse(is.na(n_cells_ds), 0, n_cells_ds),
    n_cells_umi = ifelse(is.na(n_cells_umi), 0, n_cells_umi),
    cells_removed_in_ds = n_cells_umi - n_cells_ds,
    pct_retained_in_ds = ifelse(n_cells_umi > 0, (n_cells_ds / n_cells_umi) * 100, 0)
  ) %>%
  # Remove rows where both counts are 0
  filter(n_cells_ds > 0 | n_cells_umi > 0) %>%
  arrange(sample_ID, cluster_ID)

# ===== DISPLAY RESULTS =====
message("\n=== RESULTS SUMMARY ===")

cat("OVERALL COMPARISON:\n")
cat(sprintf("DS dataset: %d total cells across %d samples\n", 
            sum(ds_sample_totals$total_cells_ds), nrow(ds_sample_totals)))
cat(sprintf("UMI dataset: %d total cells across %d samples\n",
            sum(umi_sample_totals$total_cells_umi), nrow(umi_sample_totals)))
cat(sprintf("Cells removed in downsampling: %d (%.1f%%)\n\n",
            sum(umi_sample_totals$total_cells_umi) - sum(ds_sample_totals$total_cells_ds),
            ((sum(umi_sample_totals$total_cells_umi) - sum(ds_sample_totals$total_cells_ds)) / sum(umi_sample_totals$total_cells_umi)) * 100))

# Top 20 cell types by cell count
cat("TOP 20 CELL TYPES BY CELL COUNT (UMI dataset):\n")
print(knitr::kable(head(celltype_comparison %>% 
                        select(cell_type, cluster_ID, total_cells_umi, total_cells_ds, cells_removed_in_ds, pct_retained_in_ds), 20),
                   format = "pipe", digits = 1))

cat("\n\nTOP 20 SAMPLES BY CELL COUNT:\n")
print(knitr::kable(head(sample_comparison %>%
                        select(sample_ID, tissue, patient_ID, total_cells_umi, total_cells_ds, cells_removed_in_ds, pct_retained_in_ds), 20),
                   format = "pipe", digits = 1))

# Summary by tissue type
tissue_summary <- sample_comparison %>%
  group_by(tissue) %>%
  summarise(
    n_samples = n(),
    total_cells_umi = sum(total_cells_umi, na.rm = TRUE),
    total_cells_ds = sum(total_cells_ds, na.rm = TRUE),
    cells_removed = sum(cells_removed_in_ds, na.rm = TRUE),
    avg_cells_per_sample_umi = mean(total_cells_umi, na.rm = TRUE),
    avg_cells_per_sample_ds = mean(total_cells_ds, na.rm = TRUE),
    pct_retained_overall = ifelse(total_cells_umi > 0, (total_cells_ds / total_cells_umi) * 100, 0),
    .groups = 'drop'
  )

cat("\n\nSUMMARY BY TISSUE TYPE:\n")
print(knitr::kable(tissue_summary, format = "pipe", digits = 1))

# Display cell type overlap analysis
cat("\n\nCELL TYPE OVERLAP ANALYSIS:\n")
print(knitr::kable(cell_type_overlap_summary, format = "pipe", digits = 1))

if (length(cell_types_only_in_umi) > 0) {
  cat("\n\nCell types ONLY in UMI dataset:\n")
  cat(paste(cell_types_only_in_umi, collapse = ", "))
  cat("\n")
}

if (length(cell_types_only_in_ds) > 0) {
  cat("\nCell types ONLY in DS dataset:\n") 
  cat(paste(cell_types_only_in_ds, collapse = ", "))
  cat("\n")
}

# Display Normal vs Tumor analysis
cat("\n\nTOP 15 CELL TYPES: NORMAL VS TUMOR (UMI dataset):\n")
if (length(tissue_cols) >= 2) {
  col1 <- tissue_cols[1]
  col2 <- tissue_cols[2]
  display_cols <- c("cell_type", col1, col2, "total_across_tissues", paste0("pct_", col1), paste0("pct_", col2))
  print(knitr::kable(head(umi_normal_tumor %>% select(all_of(display_cols)), 15),
                     format = "pipe", digits = 1))
} else {
  print("Not enough tissue types for comparison")
}

cat("\n\nTOP 15 CELL TYPES: NORMAL VS TUMOR (DS dataset):\n")
if (length(tissue_cols) >= 2) {
  col1 <- tissue_cols[1]
  col2 <- tissue_cols[2]
  display_cols <- c("cell_type", col1, col2, "total_across_tissues", paste0("pct_", col1), paste0("pct_", col2))
  print(knitr::kable(head(ds_normal_tumor %>% select(all_of(display_cols)), 15),
                     format = "pipe", digits = 1))
} else {
  print("Not enough tissue types for comparison")
}

# ===== UMI THRESHOLD ANALYSIS =====
message("\n=== ANALYZING UMI THRESHOLD EFFECTS ===")

# Calculate UMI counts per cell for the UMI matrix
message("Calculating UMI counts per cell...")
umi_counts_per_cell <- Matrix::colSums(umitab_matrix)

# Get cells that are in the UMI dataset but not in DS dataset
cells_in_umi_not_ds <- setdiff(colnames(umitab_matrix), colnames(ds_matrix))
cells_in_both <- intersect(colnames(umitab_matrix), colnames(ds_matrix))

message(sprintf("Cells in UMI but not in DS: %d", length(cells_in_umi_not_ds)))
message(sprintf("Cells in both datasets: %d", length(cells_in_both)))

# Analyze UMI counts for cells removed in downsampling
if (length(cells_in_umi_not_ds) > 0) {
  removed_cells_umi_counts <- umi_counts_per_cell[cells_in_umi_not_ds]
  
  # Count how many removed cells have < 2000 UMIs
  cells_below_2000 <- sum(removed_cells_umi_counts < 2000)
  cells_below_1000 <- sum(removed_cells_umi_counts < 1000)
  cells_below_500 <- sum(removed_cells_umi_counts < 500)
  
  # Summary statistics for removed cells
  umi_threshold_summary <- data.frame(
    metric = c("Total cells removed in DS", 
               "Removed cells with < 500 UMIs", 
               "Removed cells with < 1000 UMIs",
               "Removed cells with < 2000 UMIs",
               "% removed cells with < 500 UMIs",
               "% removed cells with < 1000 UMIs", 
               "% removed cells with < 2000 UMIs",
               "Median UMIs in removed cells",
               "Mean UMIs in removed cells"),
    count = c(length(cells_in_umi_not_ds),
              cells_below_500,
              cells_below_1000, 
              cells_below_2000,
              round((cells_below_500 / length(cells_in_umi_not_ds)) * 100, 1),
              round((cells_below_1000 / length(cells_in_umi_not_ds)) * 100, 1),
              round((cells_below_2000 / length(cells_in_umi_not_ds)) * 100, 1),
              round(median(removed_cells_umi_counts), 0),
              round(mean(removed_cells_umi_counts), 0))
  )
  
  cat("UMI THRESHOLD ANALYSIS:\n")
  print(knitr::kable(umi_threshold_summary, format = "pipe"))
  
  # Compare with cells retained in DS
  retained_cells_umi_counts <- umi_counts_per_cell[cells_in_both]
  
  cat("\n\nUMI COUNT DISTRIBUTIONS:\n")
  
  # Summary statistics comparison
  umi_distribution_comparison <- data.frame(
    metric = c("Median UMIs", "Mean UMIs", "Min UMIs", "Max UMIs", 
               "Cells with < 500 UMIs", "Cells with < 1000 UMIs", "Cells with < 2000 UMIs"),
    removed_cells = c(median(removed_cells_umi_counts),
                      mean(removed_cells_umi_counts),
                      min(removed_cells_umi_counts),
                      max(removed_cells_umi_counts),
                      sum(removed_cells_umi_counts < 500),
                      sum(removed_cells_umi_counts < 1000),
                      sum(removed_cells_umi_counts < 2000)),
    retained_cells = c(median(retained_cells_umi_counts),
                       mean(retained_cells_umi_counts), 
                       min(retained_cells_umi_counts),
                       max(retained_cells_umi_counts),
                       sum(retained_cells_umi_counts < 500),
                       sum(retained_cells_umi_counts < 1000),
                       sum(retained_cells_umi_counts < 2000))
  )
  
  print(knitr::kable(umi_distribution_comparison, format = "pipe", digits = 1))
  
  # Check if the difference matches the < 2000 UMI threshold exactly
  total_cells_removed <- ncol(umitab_matrix) - ncol(ds_matrix)
  cells_below_2000_in_removed <- sum(removed_cells_umi_counts < 2000)
  
  cat(sprintf("\n\nTHRESHOLD MATCH ANALYSIS:\n"))
  cat(sprintf("Total cells removed in downsampling: %d\n", total_cells_removed))
  cat(sprintf("Removed cells with < 2000 UMIs: %d\n", cells_below_2000_in_removed))
  cat(sprintf("Match (difference = 0): %s\n", ifelse(total_cells_removed == cells_below_2000_in_removed, "YES", "NO")))
  cat(sprintf("Difference: %d cells\n", abs(total_cells_removed - cells_below_2000_in_removed)))
  
  # Create UMI count bins for removed cells
  umi_bins <- cut(removed_cells_umi_counts, 
                  breaks = c(0, 500, 1000, 1500, 2000, 3000, 5000, Inf),
                  labels = c("0-500", "500-1000", "1000-1500", "1500-2000", "2000-3000", "3000-5000", "5000+"),
                  include.lowest = TRUE)
  
  umi_bin_summary <- table(umi_bins)
  
  cat("\n\nUMI COUNT BINS FOR REMOVED CELLS:\n")
  print(knitr::kable(data.frame(UMI_Range = names(umi_bin_summary), 
                                Cell_Count = as.numeric(umi_bin_summary)), 
                     format = "pipe"))
  
} else {
  cat("No cells were removed in downsampling - DS and UMI datasets have the same cells.\n")
}

# Create output directory and save results
output_dir <- "output/7. validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save detailed tables
write_csv(celltype_comparison, file.path(output_dir, "celltype_comparison_ds_vs_umi.csv"))
write_csv(sample_comparison, file.path(output_dir, "sample_comparison_ds_vs_umi.csv"))
write_csv(ds_cell_counts, file.path(output_dir, "ds_cells_per_celltype_per_sample.csv"))
write_csv(umi_cell_counts, file.path(output_dir, "umi_cells_per_celltype_per_sample.csv"))
write_csv(combined_cells_per_sample_celltype, file.path(output_dir, "combined_cells_per_celltype_per_sample_ds_vs_umi.csv"))

# Save new analyses
write_csv(cell_type_overlap_summary, file.path(output_dir, "cell_type_overlap_summary.csv"))
write_csv(umi_normal_tumor, file.path(output_dir, "umi_celltype_counts_normal_vs_tumor.csv"))
write_csv(ds_normal_tumor, file.path(output_dir, "ds_celltype_counts_normal_vs_tumor.csv"))

# Save UMI threshold analysis if cells were removed
if (exists("umi_threshold_summary")) {
  write_csv(umi_threshold_summary, file.path(output_dir, "umi_threshold_analysis.csv"))
  write_csv(umi_distribution_comparison, file.path(output_dir, "umi_distribution_comparison.csv"))
  
  # Save detailed UMI counts for removed and retained cells
  if (length(cells_in_umi_not_ds) > 0) {
    removed_cells_data <- data.frame(
      cell_ID = cells_in_umi_not_ds,
      umi_count = umi_counts_per_cell[cells_in_umi_not_ds],
      stringsAsFactors = FALSE
    )
    write_csv(removed_cells_data, file.path(output_dir, "removed_cells_umi_counts.csv"))
    
    retained_cells_data <- data.frame(
      cell_ID = cells_in_both,
      umi_count = umi_counts_per_cell[cells_in_both],
      stringsAsFactors = FALSE
    )
    write_csv(retained_cells_data, file.path(output_dir, "retained_cells_umi_counts.csv"))
  }
}

# Save cell type lists
if (length(cell_types_only_in_umi) > 0) {
  writeLines(cell_types_only_in_umi, file.path(output_dir, "cell_types_only_in_umi.txt"))
}
if (length(cell_types_only_in_ds) > 0) {
  writeLines(cell_types_only_in_ds, file.path(output_dir, "cell_types_only_in_ds.txt"))
}
writeLines(cell_types_in_both, file.path(output_dir, "cell_types_in_both_datasets.txt"))

message(sprintf("\nDetailed results saved to: %s", output_dir))
message("Files created:")
message("- celltype_comparison_ds_vs_umi.csv")
message("- sample_comparison_ds_vs_umi.csv") 
message("- ds_cells_per_celltype_per_sample.csv")
message("- umi_cells_per_celltype_per_sample.csv")
message("- combined_cells_per_celltype_per_sample_ds_vs_umi.csv")
message("- cell_type_overlap_summary.csv")
message("- umi_celltype_counts_normal_vs_tumor.csv")
message("- ds_celltype_counts_normal_vs_tumor.csv")
message("- cell_types_only_in_umi.txt")
message("- cell_types_only_in_ds.txt") 
message("- cell_types_in_both_datasets.txt")
if (exists("umi_threshold_summary")) {
  message("- umi_threshold_analysis.csv")
  message("- umi_distribution_comparison.csv")
  message("- removed_cells_umi_counts.csv")
  message("- retained_cells_umi_counts.csv")
}