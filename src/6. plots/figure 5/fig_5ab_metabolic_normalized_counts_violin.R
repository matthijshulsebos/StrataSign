library(dplyr)
library(ggplot2)
library(readr)
library(ggpubr)

# Set basic variables
cell_type_set <- "all_clusters"
gene_type <- "metabolic"
output_figure_dir <- "output/6. plots/figure 5"

# Counts data path
base_dir <- file.path("output", "6. plots", "data", "cell_type_norm_counts", gene_type, cell_type_set)

# Create output directory
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read and prepare data for a normalization method
prepare_violin_data <- function(normalization_method) {
  file_path <- file.path(base_dir, paste0(normalization_method, "_total_counts.csv"))
  
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      # Remove the cluster id from the cell type label
      cell_type_clean = paste0(gsub("_[0-9]+$", "", cell_type_label), " (", cluster_ID, ")"),
      # Log transform the total counts for visualization
      log_total_count = log10(total_count + 1)
    ) %>%
    filter(!is.na(total_count), total_count > 0)
  
  return(data)
}

# Create violin plot function
create_violin_plot <- function(data, title, output_name, reference_line_value = NULL) {
  # Determine cell types with at least 3 samples in each tissue type
  cell_type_summary <- data %>%
    group_by(cell_type_clean, tissue, cluster_ID) %>%
    summarise(n_samples = n(), .groups = 'drop') %>%
    group_by(cell_type_clean, cluster_ID) %>%
    summarise(
      total_samples = sum(n_samples),
      tissue_types = n_distinct(tissue),
      .groups = 'drop'
    ) %>%
    filter(total_samples >= 6, tissue_types == 2) %>%
    arrange(cluster_ID)

  # Match the data to the summary
  plot_data <- data %>%
    filter(cell_type_clean %in% cell_type_summary$cell_type_clean) %>%
    arrange(cluster_ID) %>%
    mutate(cell_type_clean = factor(cell_type_clean, levels = unique(cell_type_clean)))

  # Create violin plot
  p <- ggplot(plot_data, aes(x = cell_type_clean, y = log_total_count)) +
    geom_violin(aes(fill = cell_type_clean), alpha = 0.9, color = "gray20", linewidth = 0.5) +
    geom_jitter(color = "gray50", alpha = 0.4, width = 0.15, size = 0.8) +
    scale_fill_manual(
      values = c(
        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
        "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
        "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
        "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
        "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
        "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
        "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
        "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
        "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"
      )
    ) +
    labs(
      x = "",
      y = "Log10(total metabolic gene counts + 1)"
    ) +
    theme_bw(base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, color = "black"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none",
      strip.text = element_text(color = "black", face = "bold"),
      strip.background = element_rect(fill = "transparent", color = "black")
    )

  # If we pass the relative normalization value then plot the horizontal line
  if (!is.null(reference_line_value)) {
    p <- p + geom_hline(yintercept = reference_line_value, 
                        color = "#D73027", linewidth = 1.2, alpha = 0.8, linetype = "dashed")
  }

  # Save the plot
  output_path <- file.path(output_figure_dir, paste0(output_name, ".png"))
  ggsave(output_path, plot = p, width = 12, height = 6, dpi = 300, bg = "transparent")

  return(p)
}

# Read data for each normalization method
read_depth_data <- prepare_violin_data("read_depth")
global_norm_data <- prepare_violin_data("ctnorm_global")
relative_norm_data <- prepare_violin_data("ctnorm_relative")

# Get relative normalization value
relative_norm_value <- unique(relative_norm_data$log_total_count)[1]

# Create figures 4a and 4b
fig_4a <- create_violin_plot(
  read_depth_data, 
  "Metabolic gene expression",
  "fig_5a_read_depth_violin"
)

fig_4b <- create_violin_plot(
  global_norm_data, 
  "Metabolic gene expression",
  "fig_5b_global_norm_violin",
  relative_norm_value
)

message(paste("Completed writing 5AB to file."))
