# Load necessary libraries
library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(foreach)

# Load and validate data
load("data/lung_ldm.rd")
umitab <- lung_ldm$dataset$umitab

# Create mappings
cell_to_sample <- data.frame(
    cell = names(lung_ldm$dataset$cell_to_sample),
    sample = unname(lung_ldm$dataset$cell_to_sample),
    stringsAsFactors = FALSE
)

cell_to_cluster <- data.frame(
    cell = names(lung_ldm$dataset$cell_to_cluster),
    cluster = unname(lung_ldm$dataset$cell_to_cluster),
    stringsAsFactors = FALSE
)

# Get dimensions
samples <- unique(cell_to_sample$sample)
clusters <- unique(cell_to_cluster$cluster)
genes <- rownames(umitab)

# Print diagnostic information
print("Dataset dimensions:")
print(sprintf("Number of cells: %d", ncol(lung_ldm$dataset$umitab)))
print(sprintf("Number of genes: %d", nrow(lung_ldm$dataset$umitab)))
print(sprintf("Number of samples: %d", length(samples)))
print(sprintf("Number of clusters: %d", length(clusters)))

# Print dimensions of umitab for debugging
print("UMI table dimensions:")
print(dim(umitab))
print("First few cell names in umitab:")
print(head(colnames(umitab)))
print("First few cell names in cell_to_sample:")
print(head(cell_to_sample$cell))

# Update memory calculation
n_samples <- length(samples)
n_genes <- length(genes)
n_clusters <- length(clusters)

print(sprintf("Memory needed for full array: %.2f GB", 
              n_samples * n_genes * n_clusters * 8 / 1e9))

# Define paths
output_dir <- "src_output/reconstructed_counts"
counts_file <- file.path(output_dir, "counts_array.rds")

# Check if counts exist
if (!file.exists(counts_file)) {
    message("Reconstructing counts array...")
    
    # Initialize array
    counts_array <- array(
        0, 
        dim = c(length(samples), length(genes), length(clusters)),
        dimnames = list(samples, genes, clusters)
    )
    
    # Process clusters sequentially
    total_clusters <- length(clusters)
    for (i in seq_along(clusters)) {
        cluster <- clusters[i]
        message(sprintf("Processing cluster %d/%d: %s", i, total_clusters, cluster))
        
        cluster_cells <- cell_to_cluster$cell[cell_to_cluster$cluster == cluster]
        
        for (sample in samples) {
            sample_cells <- intersect(
                cell_to_sample$cell[cell_to_sample$sample == sample],
                cluster_cells
            )
            
            if (length(sample_cells) > 0) {
                counts_array[sample, , cluster] <- rowSums(umitab[, sample_cells, drop=FALSE])
            }
        }
        
        gc()
    }
    
    # Save results
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(counts_array, counts_file)
} else {
    message("Loading existing counts array...")
    counts_array <- readRDS(counts_file)
    message("Counts array loaded successfully")
    message(sprintf("Array dimensions: %d samples x %d genes x %d clusters", 
                   dim(counts_array)[1], 
                   dim(counts_array)[2], 
                   dim(counts_array)[3]))
}

# Validate counts with progress tracking
message("\nStarting validation...")
validation_dir <- "src_output/validation/counts"
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

# Validate reconstructed counts against original counts
validate_counts <- function(reconstructed_counts, true_counts, clusters) {
    # Check dimensions match
    if (!identical(dim(reconstructed_counts), dim(true_counts))) {
        stop("Arrays have different dimensions")
    }
    
    # Check and align dimensions
    true_samples <- dimnames(true_counts)[[1]]
    true_genes <- dimnames(true_counts)[[2]]
    true_clusters <- dimnames(true_counts)[[3]]
    
    recon_samples <- dimnames(reconstructed_counts)[[1]]
    recon_genes <- dimnames(reconstructed_counts)[[2]]
    recon_clusters <- dimnames(reconstructed_counts)[[3]]
    
    # Verify all names exist in both arrays
    if (!identical(sort(true_samples), sort(recon_samples)) ||
        !identical(sort(true_genes), sort(recon_genes)) ||
        !identical(sort(true_clusters), sort(recon_clusters))) {
        stop("Arrays have different dimension names")
    }
    
    # Reorder reconstructed array to match true counts
    aligned_counts <- reconstructed_counts[true_samples, true_genes, true_clusters]
    
    # Compare aligned arrays
    is_identical <- identical(aligned_counts, true_counts)
    total_elements <- length(true_counts)
    matching_elements <- sum(aligned_counts == true_counts)
    
    # Print results
    message("\nValidation Results:")
    message(sprintf("Arrays completely identical: %s", if(is_identical) "YES" else "NO"))
    message(sprintf("Matching elements: %d/%d (%.2f%%)", 
                   matching_elements, 
                   total_elements,
                   100 * matching_elements/total_elements))
    
    if (!is_identical) {
        mismatch_indices <- which(aligned_counts != true_counts, arr.ind = TRUE)
        message(sprintf("Number of mismatches: %d", nrow(mismatch_indices)))
    }
    
    return(is_identical)
}

# Run validation using counts_array instead of final_array
validation_results <- validate_counts(counts_array, 
                                    lung_ldm$dataset$counts, 
                                    clusters)
message("Validation complete")
print(validation_results)
