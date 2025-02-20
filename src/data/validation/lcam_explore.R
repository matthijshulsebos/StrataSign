# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(GGally)

# Load and prepare data
lcam_data <- read_csv("src_output/Leader/lcam_sample.csv")
sample_annots <- read.csv("input_tables/table_s1_sample_table.csv", row.names=1)

# Merge with tissue and patient_ID
merged_data <- merge(
    lcam_data,
    sample_annots[, c("tissue", "patient_ID"), drop=FALSE],
    by.x="...1",
    by.y="row.names"
) %>%
select(-"...1") %>%
select(LCAMhi, LCAMlo, difference, tissue, patient_ID)

# Filter for patients with both tissue types
filtered_patients <- merged_data %>%
    group_by(patient_ID) %>%
    summarize(
        has_both = all(c("Tumor", "Normal") %in% tissue)
    ) %>%
    filter(has_both) %>%
    pull(patient_ID)

# Update merged data
merged_data <- merged_data %>%
    filter(patient_ID %in% filtered_patients)

message(sprintf("Retained %d patients with both tissue types", length(filtered_patients)))

# Create output directory
output_dir <- "src_output/exploration_out"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Plot distributions by tissue type
plot_distributions <- function(data, output_dir) {
    for (var in c("LCAMhi", "LCAMlo", "difference")) {
        p <- ggplot(data, aes(x = .data[[var]], fill = tissue)) +
            geom_density(alpha = 0.5) +
            labs(title = paste("Distribution of", var),
                 x = var,
                 y = "Density") +
            theme_minimal()
        
        ggsave(file.path(output_dir, paste0("dist_", var, ".pdf")), p)
    }
}

# Plot correlations separately for each tissue type
plot_correlations <- function(data, output_dir) {
    for (type in unique(data$tissue)) {
        subset_data <- data %>%
            filter(tissue == type) %>%
            select(LCAMhi, LCAMlo, difference)
        
        p <- ggpairs(subset_data) +
            theme_minimal() +
            labs(title = paste("Correlations in", type, "tissue"))
        
        ggsave(file.path(output_dir, paste0("correlations_", type, ".pdf")), p)
    }
}

# Add correlation plots with tissue type
plot_tissue_correlations <- function(data, output_dir) {
    var_pairs <- list(
        c("LCAMhi", "LCAMlo"),
        c("LCAMhi", "difference"),
        c("LCAMlo", "difference")
    )
    
    for (pair in var_pairs) {
        p <- ggplot(data, aes_string(x = pair[1], y = pair[2], color = "tissue")) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "lm", se = TRUE) +
            theme_minimal() +
            labs(title = paste("Correlation:", pair[1], "vs", pair[2]))
        
        ggsave(file.path(output_dir, paste0("tissue_corr_", pair[1], "_", pair[2], ".pdf")), p)
    }
}

# Add PCA visualization with validation
plot_pca <- function(data, output_dir) {
    # Prepare data for PCA
    lcam_vars <- scale(data[, c("LCAMhi", "LCAMlo", "difference")])
    pca_result <- prcomp(lcam_vars)
    
    pca_data <- data.frame(
        PC1 = pca_result$x[,1],
        PC2 = pca_result$x[,2],
        tissue = data$tissue,
        patient_ID = data$patient_ID
    )
    
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                             color = patient_ID,
                             shape = tissue)) +
        geom_point(size = 4, alpha = 0.7) +
        scale_shape_manual(values = c("Tumor" = 4, "Normal" = 16)) +
        theme_minimal() +
        labs(title = "PCA of LCAM measurements",
             x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
             y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")) +
        guides(color = guide_legend(title = "Patient ID"),
               shape = guide_legend(title = "Tissue Type")) +
        theme(legend.position = "right")
    
    ggsave(file.path(output_dir, "pca_plot.pdf"), p, width = 12, height = 8)
}

# Generate visualizations
plot_distributions(merged_data, output_dir)
plot_correlations(merged_data, output_dir)
plot_tissue_correlations(merged_data, output_dir)
plot_pca(merged_data, output_dir)
