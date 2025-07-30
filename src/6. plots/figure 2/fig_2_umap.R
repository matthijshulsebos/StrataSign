library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(fs)
library(matrixStats) # Added for rowMedians

# Load utils
source("src/0. utils/format_utils.R")

# Load cp10k normalized data and metadata
load_cp10k_datasets <- function() {
  message("Loading cp10k-normalized data and metadata for UMAP analysis...")

  # Load cp10k normalized data
  cp10k_path <- "output/6. plots/data/cp10k/cp10k_normalized_umitab.rds"
  if (!file.exists(cp10k_path)) {
    stop("cp10k-normalized data file not found: ", cp10k_path)
  }
  cp10k_matrix <- readRDS(cp10k_path)
  if (!inherits(cp10k_matrix, "Matrix")) {
    stop("Loaded object from ", cp10k_path, " is not a Matrix.")
  }
  message(sprintf("Loaded cp10k-normalized matrix with %d genes and %d cells.", nrow(cp10k_matrix), ncol(cp10k_matrix)))

  # Load final cell metadata
  metadata_path <- "output/6. plots/data/cp10k/cell_metadata_final.csv"
  if (!file.exists(metadata_path)) {
    stop("Final cell metadata file not found: ", metadata_path)
  }
  cell_metadata <- read_csv(metadata_path, show_col_types = FALSE)
  cell_metadata <- cell_metadata %>% mutate(sample_ID = as.character(sample_ID))
  message(sprintf("Loaded metadata for %d cells.", nrow(cell_metadata)))

  # Load supplementary tables for annotations
  annots_list <- read_csv("base/input_tables/annots_list.csv", show_col_types = FALSE)
  table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", show_col_types = FALSE) %>% mutate(sample_ID = as.character(sample_ID))

  # Join metadata with annotations and sample info
  cell_metadata <- cell_metadata %>%
    left_join(annots_list, by = c("cluster_ID" = "cluster")) %>%
    left_join(table_s1, by = "sample_ID")

  # Ensure cluster_name exists which it initially doesnt
  if (!"cluster_name" %in% colnames(cell_metadata) || all(is.na(cell_metadata$cluster_name))) {
    message("Column cluster_name not found or is empty in metadata, creating it from lineage/sub_lineage.")
    cell_metadata <- cell_metadata %>%
      mutate(
        cluster_name = ifelse(!is.na(sub_lineage) & sub_lineage != "", sub_lineage, lineage),
        cluster_name = gsub("/", "-", cluster_name),
        cluster_name = ifelse(is.na(cluster_name), paste0("Cluster_", cluster_ID), cluster_name)
      )
  }

  return(list(
    cp10k = cp10k_matrix,
    cell_metadata = cell_metadata
  ))
}


# Run Seurat UMAP on cp10k-normalized data
run_cp10k_umap <- function(cp10k_matrix, cell_metadata, max_cells = 250000) {
  message("Running Seurat UMAP analysis using cp10k-normalized matrix.")
  cells_to_use <- intersect(cell_metadata$cell_ID, colnames(cp10k_matrix))
  if (length(cells_to_use) == 0) {
    warning("No cells to use after intersecting metadata with cp10k colnames.")
    return(data.frame(UMAP1=numeric(), UMAP2=numeric(), cell_ID=character()) %>% 
             left_join(cell_metadata, by="cell_ID"))
  }
  expression_matrix <- cp10k_matrix[, cells_to_use, drop = FALSE]
  cell_metadata <- cell_metadata %>% filter(cell_ID %in% cells_to_use)
  if (nrow(cell_metadata) == 0 || ncol(expression_matrix) == 0) {
    warning("No cells or features remaining for Seurat object creation.")
    return(data.frame(UMAP1=numeric(), UMAP2=numeric(), cell_ID=character()) %>% 
             left_join(cell_metadata, by="cell_ID"))
  }
  seurat_obj <- CreateSeuratObject(
    counts = expression_matrix,
    meta.data = cell_metadata %>% tibble::column_to_rownames("cell_ID") %>% dplyr::select(cluster_ID, cluster_name, tissue),
    min.cells = 3,
    min.features = 200
  ) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>%
    ScaleData(features = VariableFeatures(.), verbose = FALSE) %>%
    RunPCA(features = VariableFeatures(.), npcs = 50, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$cell_ID <- rownames(umap_coords)
  plot_data <- umap_coords %>% left_join(cell_metadata, by = "cell_ID")
  message("Seurat UMAP analysis complete!")
  return(plot_data)
}


# Create normal vs tumor cell type UMAP plot only
create_cp10k_umap_plots <- function(plot_data, output_dir) {
  message("Creating UMAP plots...")
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning("Plot data is empty. Skipping plot generation.")
    return(NULL)
  }
  unique_clusters <- sort(unique(plot_data$cluster_name))
  n_clusters <- length(unique_clusters)
  
  # Load sublineage color mapping
  sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
  colors <- NULL
  if (file.exists(sublineage_color_map_path)) {
    sublineage_colors <- readRDS(sublineage_color_map_path)
    colors <- sublineage_colors[unique_clusters]
  } else {
    warning("Sublineage colors file not found, no color mapping will be applied.")
  }
  cluster_centroids <- plot_data %>%
    group_by(cluster_name, tissue) %>%
    summarise(
      UMAP1_center = median(UMAP1, na.rm = TRUE),
      UMAP2_center = median(UMAP2, na.rm = TRUE),
      n_cells = n(),
      .groups = 'drop'
    ) %>%
    filter(n_cells >= 10)

  # Base plot
  base_plot <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = cluster_name)) +
    scale_color_manual(values = colors) +
    facet_wrap(~ tissue, ncol = 2) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      strip.text = element_text(size = 14, face = "plain", color = "black", hjust = 0)
    ) +
    labs(title = NULL, subtitle = NULL, color = "Transcriptional State", x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1))

  # Scatter plot with text labels (original)
  p_scatter <- base_plot +
    geom_point(size = 0.4, alpha = 0.7) +
    geom_text(data = cluster_centroids, aes(x = UMAP1_center, y = UMAP2_center, label = cluster_name),
              color = "black", size = 2.5, fontface = "bold", inherit.aes = FALSE)
  ggsave(file.path(output_dir, "figure_2_scatter.png"), p_scatter, width = 16, height = 8, dpi = 300)

  # Scatter plot with contours and labels
  p_scatter_contour <- base_plot +
    geom_point(size = 0.4, alpha = 0.5) +
    geom_density_2d(aes(color = cluster_name), contour_var = "ndensity", alpha = 0.8) +
    geom_text(data = cluster_centroids, aes(x = UMAP1_center, y = UMAP2_center, label = cluster_name),
              color = "black", size = 2.5, fontface = "bold", inherit.aes = FALSE)
  ggsave(file.path(output_dir, "figure_2_scatter_contour.png"), p_scatter_contour, width = 16, height = 8, dpi = 300)

  # Contours only plot with labels
  p_contour_only <- base_plot +
    geom_density_2d(aes(color = cluster_name), contour_var = "ndensity") +
    geom_text(data = cluster_centroids, aes(x = UMAP1_center, y = UMAP2_center, label = cluster_name),
              color = "black", size = 2.5, fontface = "bold", inherit.aes = FALSE)
  ggsave(file.path(output_dir, "figure_2_contour.png"), p_contour_only, width = 16, height = 8, dpi = 300)

  return(list(scatter = p_scatter, scatter_contour = p_scatter_contour, contour_only = p_contour_only))
}

# Main function
generate_cp10k_umap_plots <- function(output_dir = "output/6. plots/figure 2", max_cells = 250000, use_cache = TRUE) {
  dir_create(output_dir, recurse = TRUE)
  cached_umap_path <- file.path(output_dir, "umap_plot_data.csv")
  if (use_cache && file.exists(cached_umap_path)) {
    message("Loading cached UMAP data from: ", cached_umap_path)
    plot_data <- read_csv(cached_umap_path, show_col_types = FALSE)
  } else {
    message("Generating new UMAP data.")
    datasets <- load_cp10k_datasets()
    cell_metadata <- datasets$cell_metadata
    plot_data <- run_cp10k_umap(datasets$cp10k, cell_metadata, max_cells)
    write_csv(plot_data, cached_umap_path)
    message("Saved UMAP data to cache: ", cached_umap_path)
  }
  plot_list <- create_cp10k_umap_plots(plot_data, output_dir)
  message("UMAP analysis complete!")
  message(sprintf("Plots saved to: %s", output_dir))
  return(list(plot_data = plot_data, plots = plot_list))
}

umap_results <- generate_cp10k_umap_plots(output_dir = "output/6. plots/figure 2", max_cells = 50000, use_cache = TRUE)
