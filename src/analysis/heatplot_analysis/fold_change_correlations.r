library(dplyr)
library(tidyr)
library(data.table)
library(Matrix)
library(fs)

# Load data and annotations
if (!exists("lung_ldm")) {
  message("Loading Mt. Sinai scRNAseq data into R")
  load("base/data/lung_ldm.rd")
}

table_s1 <- fread("base/input_tables/table_s1_sample_table.csv")
annots_list <- fread("base/input_tables/annots_list.csv")
hsa01100_genes <- fread("data/kegg/hsa01100_genes.csv")

# Filter samples
table_s1 <- table_s1 %>% 
  filter(`Use.in.Clustering.Model.` == "Yes") %>%
  mutate(sample_ID = as.character(sample_ID))

calculate_fold_changes <- function() {
  # Load and filter counts - remove metabolic gene filtering
  counts <- lung_ldm$dataset$counts
  sample_ids <- table_s1$sample_ID
  counts <- counts[sample_ids, , , drop = FALSE]  # Only filter samples
  
  # Debug dimensions
  message("Dimensions after filtering:")
  message("Samples: ", nrow(counts))
  message("Genes: ", ncol(counts))
  message("Clusters: ", dim(counts)[3])
  
  # Convert to long format and add tissue information
  counts_long <- as.data.frame(as.table(counts)) %>%
    setnames(c("sample_ID", "gene", "cluster", "count")) %>%
    mutate(count = as.numeric(count)) %>%
    left_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID")
  
  # Add cluster names
  counts_long <- counts_long %>%
    mutate(cluster_num = as.numeric(cluster)) %>%
    left_join(
      annots_list %>% 
        mutate(
          cluster_num = as.numeric(cluster),
          cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
          cluster_name = paste(cluster_name, cluster, sep = "_")
        ) %>%
        select(cluster_num, cluster_name),
      by = "cluster_num"
    ) %>%
    filter(!is.na(cluster_name))
  
  # Calculate fold changes between Normal and Tumor
  fold_changes <- counts_long %>%
    # Calculate mean expression per tissue condition
    group_by(gene, cluster_name, tissue) %>%
    summarize(mean_expr = mean(count), .groups = 'drop') %>%
    pivot_wider(
      names_from = tissue,
      values_from = mean_expr
    ) %>%
    # Calculate log2 fold change
    mutate(
      Value = log2((Tumor + 1) / (Normal + 1)),  # Add pseudocount to avoid division by zero
      Feature = paste0(gene, "@", cluster_name)
    ) %>%
    select(Feature, Value)
  
  return(fold_changes)
}

calculate_similarities <- function(fold_changes) {
  # Create similarities directory
  dir.create("data/ablation/intermediates/heatplot/similarities", recursive = TRUE, showWarnings = FALSE)
  
  # Recursively find all CSV files in heatplot directory (excluding similarities folder)
  csv_files <- dir_ls("data/ablation/intermediates/heatplot", recurse = TRUE, glob = "*.csv") %>%
    .[!grepl("similarities|fold_changes\\.csv|variances\\.csv", .)]
  
  # Initialize results dataframe
  similarities <- data.frame()
  
  # Process each file
  for(file_path in csv_files) {
    # Read feature importance file
    importance <- fread(file_path)
    message("\nProcessing: ", basename(file_path))
    message("Features in importance file: ", nrow(importance))
    
    # Filter fold changes to only include features present in importance file
    filtered_fold_changes <- fold_changes %>%
      filter(Feature %in% importance$Feature)
    
    message("Matching features found: ", nrow(filtered_fold_changes))
    
    if(nrow(filtered_fold_changes) == 0) {
      warning("No matching features found for ", basename(file_path))
      next
    }
    
    # Merge fold changes with importance scores
    merged <- filtered_fold_changes %>%
      inner_join(importance, by = "Feature", suffix = c("_fold", "_imp"))
    
    message("Features after merge: ", nrow(merged))
    
    # Calculate similarity score for each feature
    feature_similarities <- merged %>%
      mutate(
        # First check sign agreement
        sign_agreement = sign(Value_fold) == sign(Value_imp),
        
        # Calculate normalized magnitude similarity (scale-invariant)
        max_magnitude = pmax(abs(Value_fold), abs(Value_imp)),
        min_magnitude = pmin(abs(Value_fold), abs(Value_imp)),
        magnitude_similarity = min_magnitude / max_magnitude,
        
        # Combine sign and magnitude (ranges from -1 to 1)
        Value = ifelse(sign_agreement, magnitude_similarity, -magnitude_similarity)
      ) %>%
      select(Feature, Value)
    
    # Create corresponding output path in similarities directory
    rel_path <- path_rel(file_path, "data/ablation/intermediates/heatplot")
    output_dir <- file.path("data/ablation/intermediates/heatplot/similarities", dirname(rel_path))
    dir_create(output_dir, recurse = TRUE)
    
    # Save with same structure as input
    output_path <- file.path(output_dir, basename(file_path))
    fwrite(feature_similarities, output_path)
    
    # Store summary statistics
    similarities <- bind_rows(similarities, 
                            data.frame(File = rel_path,
                                     MeanSimilarity = mean(feature_similarities$Value)))
  }
  
  return(similarities)
}

# Create output directories
dir.create("data/ablation/intermediates/analysis", recursive = TRUE, showWarnings = FALSE)

# Calculate fold changes and save to analysis directory
fold_changes <- calculate_fold_changes()
fwrite(fold_changes, "data/ablation/intermediates/analysis/fold_changes.csv")

# Calculate similarities
similarities <- calculate_similarities(fold_changes)

# Print summary statistics
message("\nMean similarities with feature importance files:")
print(similarities)
