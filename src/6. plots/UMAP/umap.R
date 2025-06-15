library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(fs)
library(matrixStats) # Added for rowMedians

# Load utils
source("src/0. utils/format_utils.R")

# Load datasets
load_umap_datasets <- function() {
  message("Loading datasets for UMAP analysis...")
  
  env_lung_ldm <- new.env()
  load("base/data/lung_ldm.rd", envir = env_lung_ldm)
  
  if (!"lung_ldm" %in% ls(env_lung_ldm)) {
    stop("Failed to load 'lung_ldm' object from base/data/lung_ldm.rd.")
  }
  lung_ldm_from_file <- env_lung_ldm$lung_ldm
  
  if (is.null(lung_ldm_from_file)) {
    stop("'lung_ldm' object is NULL after loading.")
  }
  if (!is.environment(lung_ldm_from_file) && !is.list(lung_ldm_from_file)) {
     stop(sprintf("'lung_ldm' is of class: %s. Expected a list or environment.", class(lung_ldm_from_file)[1]))
  }

  dataset_obj <- NULL
  if (is.environment(lung_ldm_from_file)) {
    if (exists("dataset", envir = lung_ldm_from_file, inherits = FALSE)) {
      dataset_obj <- lung_ldm_from_file$dataset
    }
  } else { # is.list(lung_ldm_from_file)
    if ("dataset" %in% names(lung_ldm_from_file)) {
      dataset_obj <- lung_ldm_from_file$dataset
    }
  }

  if (is.null(dataset_obj)) {
    stop("'dataset' object/element not found in 'lung_ldm'.")
  }
  if (!is.environment(dataset_obj) && !is.list(dataset_obj)) {
    stop(sprintf("'lung_ldm$dataset' is of class: %s. Expected a list or environment.", class(dataset_obj)[1]))
  }
  
  potential_ds <- NULL
  if (is.environment(dataset_obj)) {
    if (exists("ds", envir = dataset_obj, inherits = FALSE)) {
      potential_ds <- dataset_obj$ds
    }
  } else { # is.list(dataset_obj)
    if ("ds" %in% names(dataset_obj)) {
      potential_ds <- dataset_obj$ds
    }
  }

  if (is.null(potential_ds)) {
    stop("'ds' object/element not found in 'lung_ldm$dataset'.")
  }

  actual_ds_matrix <- NULL
  if (inherits(potential_ds, "Matrix")) {
    actual_ds_matrix <- potential_ds
  } else if (is.list(potential_ds)) {
    message("'lung_ldm$dataset$ds' is a list. Attempting to extract Matrix object...")
    # Try common scenarios for a matrix within a list
    if (length(potential_ds) == 1 && inherits(potential_ds[[1]], "Matrix")) {
      actual_ds_matrix <- potential_ds[[1]]
    } else if ("counts" %in% names(potential_ds) && inherits(potential_ds$counts, "Matrix")) {
      actual_ds_matrix <- potential_ds$counts
    } else if ("matrix" %in% names(potential_ds) && inherits(potential_ds$matrix, "Matrix")) {
      actual_ds_matrix <- potential_ds$matrix
    } else if ("data" %in% names(potential_ds) && inherits(potential_ds$data, "Matrix")) {
      actual_ds_matrix <- potential_ds$data
    }
    
    if (!is.null(actual_ds_matrix)) {
        message("Successfully extracted Matrix from 'lung_ldm$dataset$ds' list.")
        # Update the 'ds' in the original structure to be the matrix itself
        if (is.environment(dataset_obj)) {
            dataset_obj$ds <- actual_ds_matrix
        } else { # dataset_obj is a list
            dataset_obj[['ds']] <- actual_ds_matrix
            # If lung_ldm_from_file is a list, dataset_obj was a copy of lung_ldm_from_file$dataset (a list)
            # So, we need to update lung_ldm_from_file$dataset$ds
            if(is.list(lung_ldm_from_file) && "dataset" %in% names(lung_ldm_from_file) && is.list(lung_ldm_from_file$dataset)) {
                 lung_ldm_from_file$dataset[['ds']] <- actual_ds_matrix
            }
        }
    }
  }

  if (is.null(actual_ds_matrix) || !inherits(actual_ds_matrix, "Matrix")) {
    stop(sprintf("'lung_ldm$dataset$ds' is of class: %s and could not be resolved to a Matrix object. Expected a Matrix or a list containing one.", 
                 class(potential_ds)[1]))
  }
  
  message(sprintf("Successfully loaded and verified 'ds' Matrix with %d rows and %d columns.", 
                  nrow(actual_ds_matrix), ncol(actual_ds_matrix)))
  
  # Load sample metadata
  table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv", 
                       show_col_types = FALSE) %>% 
    mutate(sample_ID = as.character(sample_ID))
  
  # Load annotations
  annots_list <- read_csv("base/input_tables/annots_list.csv", 
                         show_col_types = FALSE)
  
  return(list(
    lung_ldm = lung_ldm_from_file, 
    table_s1 = table_s1,
    annots_list = annots_list
  ))
}

# Load intersector meta scores (ctnorm_global/all_clusters/metabolic only)
load_meta_scores <- function(intersector_path = "output/3. intersector") {
  message("Loading meta scores from intersector results...")
  
  meta_score_files <- list.files(intersector_path, 
                                pattern = "meta_scores.csv", 
                                recursive = TRUE, 
                                full.names = TRUE)
  
  # Look for ctnorm_global/all_clusters/metabolic
  ctnorm_global_metabolic_files <- meta_score_files[
    grepl("ctnorm_global/all_clusters/metabolic", meta_score_files)
  ]
  
  if (length(ctnorm_global_metabolic_files) == 0) {
    warning("No meta score files found for ctnorm_global/all_clusters/metabolic")
    return(NULL)
  }
  
  all_meta_scores <- map_dfr(ctnorm_global_metabolic_files, ~ {
    tryCatch({
      read_csv(.x, show_col_types = FALSE) %>%
        mutate(source_file = basename(dirname(.x)))
    }, error = function(e) {
      message("Error reading ", .x, ": ", e$message)
      NULL
    })
  })
  
  if (nrow(all_meta_scores) == 0) {
    warning("No meta scores loaded from ctnorm_global/all_clusters/metabolic")
    return(NULL)
  }
  
  message(sprintf("Loaded %d meta scores from %d ctnorm_global/all_clusters/metabolic files", nrow(all_meta_scores), length(ctnorm_global_metabolic_files)))
  return(all_meta_scores)
}

# Prepare single cell metadata
prepare_cell_metadata <- function(lung_ldm, table_s1, annots_list, 
                                 doublets = c(3, 4, 6, 12, 15, 21, 22, 24, 26, 27, 60)) {
  
  current_dataset <- lung_ldm$dataset # Assumes load_umap_datasets validated lung_ldm and lung_ldm$dataset
  if (is.null(current_dataset)) { # Basic check
      stop("lung_ldm$dataset is NULL in prepare_cell_metadata.")
  }

  message("Preparing single cell metadata...")
  
  # Ensure cell_to_sample and cell_to_cluster exist and are named vectors/lists
  if (is.null(current_dataset$cell_to_sample) || is.null(names(current_dataset$cell_to_sample))) {
      stop("'cell_to_sample' is missing or not a named vector/list in lung_ldm$dataset.")
  }
  if (is.null(current_dataset$cell_to_cluster) || is.null(names(current_dataset$cell_to_cluster))) {
      stop("'cell_to_cluster' is missing or not a named vector/list in lung_ldm$dataset.")
  }

  cell_to_sample_df <- data.frame(
    cell_ID = names(current_dataset$cell_to_sample), # Use current_dataset
    sample_ID = unname(current_dataset$cell_to_sample), # Use current_dataset
    stringsAsFactors = FALSE
  )
  
  cell_to_cluster_df <- data.frame(
    cell_ID = names(current_dataset$cell_to_cluster), # Use current_dataset
    cluster_ID = as.numeric(unname(current_dataset$cell_to_cluster)), # Use current_dataset
    stringsAsFactors = FALSE
  )
  
  # Combine metadata
  cell_metadata <- inner_join(cell_to_sample_df, cell_to_cluster_df, by = "cell_ID") %>%
    filter(!cluster_ID %in% doublets) %>%
    left_join(table_s1 %>% dplyr::select(sample_ID, tissue, patient_ID), by = "sample_ID") %>%
    left_join(annots_list %>% dplyr::select(cluster, lineage, sub_lineage), 
              by = c("cluster_ID" = "cluster")) %>%
    mutate(
      cluster_name = ifelse(!is.na(sub_lineage), sub_lineage, lineage),
      cluster_name = gsub("/", "-", cluster_name),
      cluster_name = ifelse(is.na(cluster_name), paste0("Cluster_", cluster_ID), cluster_name),
      tissue = ifelse(is.na(tissue), "Unknown", tissue)
    ) %>%
    filter(tissue %in% c("Normal", "Tumor"))  # Only keep Normal and Tumor
  

  message(sprintf("Prepared metadata for %d single cells", nrow(cell_metadata)))
  return(cell_metadata)
}

# Apply meta scores to cells
apply_meta_scores <- function(lung_ldm, cell_metadata, meta_scores_data) {
  message("Applying meta scores to single cells using ds (downsampled matrix)...")
  
  if (is.null(meta_scores_data) || nrow(meta_scores_data) == 0) {
    warning("No meta scores data available")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  
  # Parse feature_id: gene@sublineage_clusterid format
  parsed_meta_scores <- meta_scores_data %>%
    tidyr::separate(feature_id, into = c("gene", "sublineage_cluster"), sep = "@", remove = FALSE) %>%
    mutate(sublineage_clean = gsub("_[0-9]+$", "", sublineage_cluster)) %>%
    filter(!is.na(gene), !is.na(sublineage_clean), !is.na(meta_score)) %>%
    dplyr::select(gene, sublineage_clean, meta_score) %>%
    distinct()
  
  if (nrow(parsed_meta_scores) == 0) {
    warning("No valid gene-sublineage meta scores found")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  
  # Create mapping between cluster names and sublineages
  unique_sublineages <- unique(parsed_meta_scores$sublineage_clean)
  unique_cluster_names <- unique(cell_metadata$cluster_name)
  
  cluster_to_sublineage_map <- data.frame(
    cluster_name = unique_cluster_names,
    sublineage_match = NA_character_
  )
  
  # Find exact and partial matches
  for (i in 1:nrow(cluster_to_sublineage_map)) {
    cluster_name <- cluster_to_sublineage_map$cluster_name[i]
    if (cluster_name %in% unique_sublineages) {
      cluster_to_sublineage_map$sublineage_match[i] <- cluster_name
    } else {
      matches <- unique_sublineages[grepl(cluster_name, unique_sublineages, ignore.case = TRUE)]
      if (length(matches) == 0) {
        matches <- unique_sublineages[grepl(paste0("\\b", gsub("[^A-Za-z0-9]", "", cluster_name), "\\b"), 
                                           gsub("[^A-Za-z0-9]", "", unique_sublineages), 
                                           ignore.case = TRUE)]
      }
      if (length(matches) > 0) {
        cluster_to_sublineage_map$sublineage_match[i] <- matches[1]
      }
    }
  }
  
  genes_with_meta_scores <- unique(parsed_meta_scores$gene)
  ds_matrix_obj <- lung_ldm$dataset$ds 
  
  if (is.null(ds_matrix_obj) || !inherits(ds_matrix_obj, "Matrix")) {
    warning("ds matrix is not available or not a Matrix object in lung_ldm$dataset. Cannot apply meta scores.")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  if (is.null(rownames(ds_matrix_obj))) {
    warning("ds (downsampled matrix) rownames are not available. Cannot apply meta scores.")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  genes_in_expression_data <- intersect(genes_with_meta_scores, rownames(ds_matrix_obj))
  if (length(genes_in_expression_data) == 0) {
    warning("No genes with meta scores found in ds (downsampled matrix) expression data")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  if (is.null(colnames(ds_matrix_obj))) {
    warning("ds (downsampled matrix) colnames are not available. Cannot apply meta scores.")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  cells_to_use <- intersect(cell_metadata$cell_ID, colnames(ds_matrix_obj))
  expression_matrix <- ds_matrix_obj[genes_in_expression_data, cells_to_use, drop = FALSE]
  
  # Add sublineage mapping to cell metadata
  cell_metadata_with_mapping <- cell_metadata %>%
    filter(cell_ID %in% cells_to_use) %>%
    left_join(cluster_to_sublineage_map, by = "cluster_name")
  
  # Prepare meta score lookup table (gene, sublineage) -> meta_score
  meta_score_lookup <- parsed_meta_scores %>%
    dplyr::select(gene, sublineage_clean, meta_score)
  
  # For each gene, scale meta scores between 0 and 1 (across all sublineages)
  meta_score_lookup <- meta_score_lookup %>%
    group_by(gene) %>%
    mutate(meta_score_scaled = ifelse(max(meta_score, na.rm = TRUE) == min(meta_score, na.rm = TRUE),
                                      0,
                                      (meta_score - min(meta_score, na.rm = TRUE)) / (max(meta_score, na.rm = TRUE) - min(meta_score, na.rm = TRUE)))
    ) %>%
    ungroup()
  
  # Vectorized: Create a lookup matrix for (gene, sublineage)
  gene_ids <- rownames(expression_matrix)
  sublineages <- unique(meta_score_lookup$sublineage_clean)
  meta_score_mat <- matrix(0, nrow = length(gene_ids), ncol = length(sublineages),
                           dimnames = list(gene_ids, sublineages))
  meta_score_idx <- match(meta_score_lookup$gene, gene_ids)
  sublineage_idx <- match(meta_score_lookup$sublineage_clean, sublineages)
  valid <- !is.na(meta_score_idx) & !is.na(sublineage_idx)
  meta_score_mat[cbind(meta_score_idx[valid], sublineage_idx[valid])] <- meta_score_lookup$meta_score_scaled[valid]
  
  # For each cell, get its sublineage
  cell_sublineage <- cell_metadata_with_mapping$sublineage_match
  names(cell_sublineage) <- cell_metadata_with_mapping$cell_ID
  cell_ids <- colnames(expression_matrix)
  cell_sublineage_vec <- cell_sublineage[cell_ids]
  sublineage_col_idx <- match(cell_sublineage_vec, sublineages)

  # Build meta score matrix for all genes x cells (vectorized)
  meta_score_scaled_mat <- matrix(0, nrow = length(gene_ids), ncol = length(cell_ids),
                                   dimnames = list(gene_ids, cell_ids))
  valid_cells <- which(!is.na(sublineage_col_idx))
  if (length(valid_cells) > 0) {
    meta_score_scaled_mat[, valid_cells] <- meta_score_mat[, sublineage_col_idx[valid_cells], drop = FALSE]
  }

  # Group cells by cluster_name and tissue for group-wise median/deviation
  cell_groups <- cell_metadata_with_mapping %>%
    dplyr::select(cell_ID, cluster_name, tissue) %>%
    mutate(group_id = paste(cluster_name, tissue, sep = "|"))
  group_ids <- cell_groups$group_id
  names(group_ids) <- cell_groups$cell_ID

  # Prepare output vector for meta_score_expression
  meta_score_expression_vec <- setNames(rep(0, length(cell_ids)), cell_ids)

  expr_mat <- as.matrix(expression_matrix)

  # For each group, calculate gene medians and deviations within group
  unique_groups <- unique(group_ids[cell_ids])
  for (grp in unique_groups) {
    grp_cells <- names(group_ids)[group_ids == grp]
    grp_cells <- intersect(grp_cells, cell_ids)
    if (length(grp_cells) == 0) next
    grp_expr <- expr_mat[, grp_cells, drop = FALSE]
    # Median per gene within group
    grp_gene_medians <- matrixStats::rowMedians(grp_expr)
    # Squared deviation per gene per cell
    grp_deviation_sq <- sweep(grp_expr, 1, grp_gene_medians, FUN = function(x, med) (x - med)^2)
    # Scale squared deviations per gene between 0 and 1
    min_dev <- matrixStats::rowMins(grp_deviation_sq)
    max_dev <- matrixStats::rowMaxs(grp_deviation_sq)
    range_dev <- max_dev - min_dev
    range_dev[range_dev == 0] <- 1
    grp_deviation_sq_scaled <- sweep(grp_deviation_sq, 1, min_dev, "-")
    grp_deviation_sq_scaled <- sweep(grp_deviation_sq_scaled, 1, range_dev, "/")
    grp_deviation_sq_scaled[is.na(grp_deviation_sq_scaled)] <- 0
    # Sum (scaled deviation + meta score) for each cell in group
    grp_meta_score_scaled <- meta_score_scaled_mat[, grp_cells, drop = FALSE]
    grp_sum <- colSums(grp_deviation_sq_scaled + grp_meta_score_scaled, na.rm = TRUE)
    meta_score_expression_vec[grp_cells] <- grp_sum
  }

  # Scale meta_score_expression_vec between 0 and 1 (after summing)
  min_score <- min(meta_score_expression_vec, na.rm = TRUE)
  max_score <- max(meta_score_expression_vec, na.rm = TRUE)
  range_score <- max_score - min_score
  if (range_score == 0) {
    meta_score_expression_vec[] <- 0
  } else {
    meta_score_expression_vec <- (meta_score_expression_vec - min_score) / range_score
  }

  # Add meta scores to cell metadata
  cell_metadata_with_meta_scores <- cell_metadata %>%
    mutate(
      meta_score_expression = ifelse(cell_ID %in% names(meta_score_expression_vec), 
                                   meta_score_expression_vec[cell_ID], 0)
    )

  message(sprintf("Applied meta scores to %d cells", sum(cell_metadata_with_meta_scores$meta_score_expression > 0)))

  return(cell_metadata_with_meta_scores)
}

# Run Seurat UMAP
run_seurat_umap <- function(lung_ldm, cell_metadata, max_cells = 50000) {
  message("Running Seurat UMAP analysis using ds (downsampled matrix)...")
  
  # Always use all cells (no sampling)
  # Get expression matrix from ds (downsampled matrix)
  # Assumes load_umap_datasets has prepared lung_ldm$dataset$ds as a Matrix
  ds_matrix_obj <- lung_ldm$dataset$ds

  if (is.null(ds_matrix_obj) || !inherits(ds_matrix_obj, "Matrix")) {
    warning("ds matrix is not available or not a Matrix object in lung_ldm$dataset. Cannot proceed with UMAP.")
    return(data.frame(UMAP1=numeric(), UMAP2=numeric(), cell_ID=character()) %>% 
             left_join(cell_metadata, by="cell_ID"))
  }
  
  if (is.null(colnames(ds_matrix_obj))) {
    warning("ds (downsampled matrix) colnames are not available. Cannot proceed with UMAP.")
    return(data.frame(UMAP1=numeric(), UMAP2=numeric(), cell_ID=character()) %>% 
             left_join(cell_metadata, by="cell_ID"))
  }
  cells_to_use <- intersect(cell_metadata$cell_ID, colnames(ds_matrix_obj))
  
  if (length(cells_to_use) == 0) {
    warning("No cells to use after intersecting metadata with ds (downsampled matrix) colnames.")
    return(data.frame(UMAP1=numeric(), UMAP2=numeric(), cell_ID=character()) %>% 
             left_join(cell_metadata, by="cell_ID"))
  }
  expression_matrix <- ds_matrix_obj[, cells_to_use, drop = FALSE]
  
  # Filter cell metadata to cells in expression data
  cell_metadata <- cell_metadata %>% filter(cell_ID %in% cells_to_use)
  
  if (nrow(cell_metadata) == 0 || ncol(expression_matrix) == 0) {
    warning("No cells or features remaining for Seurat object creation.")
    return(data.frame(UMAP1=numeric(), UMAP2=numeric(), cell_ID=character()) %>% 
             left_join(cell_metadata, by="cell_ID"))
  }

  # Create Seurat object and run workflow
  seurat_obj <- CreateSeuratObject(
    counts = expression_matrix,
    meta.data = cell_metadata %>% 
      tibble::column_to_rownames("cell_ID") %>%
      dplyr::select(cluster_ID, cluster_name, tissue, meta_score_expression),
    min.cells = 3,
    min.features = 200
  ) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>%
    ScaleData(features = VariableFeatures(.), verbose = FALSE) %>%
    RunPCA(features = VariableFeatures(.), npcs = 50, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
  
  # Extract UMAP coordinates
  umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$cell_ID <- rownames(umap_coords)
  
  # Combine with metadata
  plot_data <- umap_coords %>% left_join(cell_metadata, by = "cell_ID")
  
  message("Seurat UMAP analysis complete!")
  
  return(plot_data)
}

# Create normal vs tumor plots
create_normal_vs_tumor_plots <- function(plot_data, output_dir) {
  message("Creating normal vs tumor comparison plots...")

  # Check if plot_data is valid for plotting
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning("Plot data is empty. Skipping normal vs tumor plot generation.")
    return(list(celltype_plot = NULL, metascore_plot = NULL))
  }
  
  required_cols <- c("UMAP1", "UMAP2", "cluster_name", "tissue", "meta_score_expression")
  missing_cols <- setdiff(required_cols, colnames(plot_data))
  if (length(missing_cols) > 0) {
    warning(sprintf("Plot data is missing required columns: %s. Skipping plot generation.", 
                    paste(missing_cols, collapse = ", ")))
    return(list(celltype_plot = NULL, metascore_plot = NULL))
  }
  
  if (all(is.na(plot_data$tissue)) || length(unique(na.omit(plot_data$tissue))) == 0) {
    warning("Tissue column in plot_data has no valid values for faceting. Skipping plot generation.")
    return(list(celltype_plot = NULL, metascore_plot = NULL))
  }

  # Create color palette for cell types
  unique_clusters <- sort(unique(plot_data$cluster_name))
  n_clusters <- length(unique_clusters)
  
  if (n_clusters <= 12) {
    colors <- brewer.pal(min(max(n_clusters, 3), 12), "Set3")
  } else {
    colors <- rainbow(n_clusters)
  }
  names(colors) <- unique_clusters
  
  # Calculate cluster centroids for labels
  cluster_centroids <- plot_data %>%
    group_by(cluster_name, tissue) %>%
    summarise(
      UMAP1_center = median(UMAP1, na.rm = TRUE),
      UMAP2_center = median(UMAP2, na.rm = TRUE),
      n_cells = n(),
      .groups = 'drop'
    ) %>%
    filter(n_cells >= 10)  # Only label clusters with enough cells
  
  # Plot 1: Cell types by tissue
  p1 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = cluster_name)) +
    geom_point(size = 0.4, alpha = 0.7) +
    geom_text(data = cluster_centroids, 
              aes(x = UMAP1_center, y = UMAP2_center, label = cluster_name),
              color = "black", size = 2.5, fontface = "bold",
              inherit.aes = FALSE) +
    scale_color_manual(values = colors) +
    facet_wrap(~ tissue, ncol = 2) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "Single Cell UMAP: Normal vs Tumor Tissue Comparison",
      subtitle = "Cell types colored consistently across tissues",
      color = "Cell Type",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1))
  
  # Plot 2: Meta score expression by tissue
  # Ensure meta_score_expression has non-NA values before calculating quantile
  if (all(is.na(plot_data$meta_score_expression))) {
    percentile_90 <- 0 # Default if all are NA
    warning("meta_score_expression column contains all NA values.")
  } else {
    percentile_90 <- quantile(plot_data$meta_score_expression, probs = 0.9, na.rm = TRUE)
  }
  min_score <- 0
  
  p2 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = meta_score_expression)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_viridis_c(
      name = "Meta Score\nWeighted\nExpression", 
      option = "plasma",
      limits = c(min_score, percentile_90),
      oob = scales::squish,
      labels = scales::comma_format(accuracy = 0.01)
    ) +
    facet_wrap(~ tissue, ncol = 2) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "Single Cell UMAP: Meta Score Expression by Tissue",
      subtitle = sprintf("Sum of (expression × meta score). Scale: 0 - %.2f (90th percentile)", percentile_90),
      x = "UMAP 1",
      y = "UMAP 2"
    )
  
  # Save plots
  ggsave(file.path(output_dir, "normal_vs_tumor_celltype.png"), p1, width = 16, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "normal_vs_tumor_celltype.pdf"), p1, width = 16, height = 8)
  
  ggsave(file.path(output_dir, "normal_vs_tumor_metascores.png"), p2, width = 16, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "normal_vs_tumor_metascores.pdf"), p2, width = 16, height = 8)
  
  return(list(celltype_plot = p1, metascore_plot = p2))
}

# Main function
generate_umap_plots <- function(output_dir = "output/6. plots/UMAP", max_cells = 50000) {
  dir_create(output_dir, recurse = TRUE)
  
  # Load data
  datasets <- load_umap_datasets()
  meta_scores <- load_meta_scores()
  
  # Prepare cell metadata
  cell_metadata <- prepare_cell_metadata(datasets$lung_ldm, datasets$table_s1, datasets$annots_list)
  
  # Apply meta scores
  cell_metadata <- apply_meta_scores(datasets$lung_ldm, cell_metadata, meta_scores)
  
  # Run UMAP (all cells, no sampling)
  plot_data <- run_seurat_umap(datasets$lung_ldm, cell_metadata, max_cells)
  
  # Create plots
  plots <- create_normal_vs_tumor_plots(plot_data, output_dir)
  
  # Save data
  write_csv(plot_data, file.path(output_dir, "umap_plot_data.csv"))
  
  message("UMAP analysis complete!")
  message(sprintf("Plots saved to: %s", output_dir))
  
  return(list(plot_data = plot_data, plots = plots))
}

# Run the analysis
if (interactive() || !exists(".umap_sourced")) {
  .umap_sourced <- TRUE
  umap_results <- generate_umap_plots(max_cells = 50000)
}
