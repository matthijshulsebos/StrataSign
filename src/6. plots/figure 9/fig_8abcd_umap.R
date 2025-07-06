library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(fs)
library(matrixStats) # Added for rowMedians

# Load utils
source("src/0. utils/format_utils.R")

# Load cp10k-normalized data and metadata
load_cp10k_datasets <- function() {
  message("Loading cp10k-normalized data and metadata for UMAP analysis...")

  # Load cp10k-normalized matrix from RDS
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

  # Ensure cluster_name exists, creating it if necessary
  if (!"cluster_name" %in% colnames(cell_metadata) || all(is.na(cell_metadata$cluster_name))) {
    message("Column 'cluster_name' not found or is empty in metadata, creating it from lineage/sub_lineage...")
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
  
  if (is.null(lung_ldm) || !inherits(lung_ldm, "Matrix")) {
    warning("ds (downsampled matrix) is not available or not a Matrix object. Cannot apply meta scores.")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  if (is.null(rownames(lung_ldm))) {
    warning("ds (downsampled matrix) rownames are not available. Cannot apply meta scores.")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  genes_in_expression_data <- intersect(genes_with_meta_scores, rownames(lung_ldm))
  if (length(genes_in_expression_data) == 0) {
    warning("No genes with meta scores found in ds (downsampled matrix) expression data")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  if (is.null(colnames(lung_ldm))) {
    warning("ds (downsampled matrix) colnames are not available. Cannot apply meta scores.")
    return(cell_metadata %>% mutate(meta_score_expression = 0))
  }
  cells_to_use <- intersect(cell_metadata$cell_ID, colnames(lung_ldm))
  expression_matrix <- lung_ldm[genes_in_expression_data, cells_to_use, drop = FALSE]
  
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
                                     meta_score_expression_vec[cell_ID], 
                                     0)
    )

  message(sprintf("Applied meta scores to %d cells.", nrow(cell_metadata_with_meta_scores)))

  return(cell_metadata_with_meta_scores)
}

# Run Seurat UMAP on cp10k-normalized data
run_cp10k_umap <- function(cp10k_matrix, cell_metadata, max_cells = 50000) {
  message("Running Seurat UMAP analysis using cp10k-normalized matrix...")
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

# Create normal vs tumor cell type plot and UMAPs for top 3 genes
create_cp10k_umap_plots <- function(plot_data, cp10k_matrix, cell_metadata, intersector_path, output_dir) {
  message("Creating UMAP plots: cell type and top 3 intersector genes...")
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning("Plot data is empty. Skipping plot generation.")
    return(list(celltype_plot = NULL, gene_umap_plots = list()))
  }
  # Cell type UMAP (Normal vs Tumor)
  unique_clusters <- sort(unique(plot_data$cluster_name))
  n_clusters <- length(unique_clusters)
  colors <- if (n_clusters <= 12) brewer.pal(min(max(n_clusters, 3), 12), "Set3") else rainbow(n_clusters)
  names(colors) <- unique_clusters
  cluster_centroids <- plot_data %>%
    group_by(cluster_name, tissue) %>%
    summarise(
      UMAP1_center = median(UMAP1, na.rm = TRUE),
      UMAP2_center = median(UMAP2, na.rm = TRUE),
      n_cells = n(),
      .groups = 'drop'
    ) %>%
    filter(n_cells >= 10)
  p_celltype <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = cluster_name)) +
    geom_point(size = 0.4, alpha = 0.7) +
    geom_text(data = cluster_centroids, aes(x = UMAP1_center, y = UMAP2_center, label = cluster_name),
              color = "black", size = 2.5, fontface = "bold", inherit.aes = FALSE) +
    scale_color_manual(values = colors) +
    facet_wrap(~ tissue, ncol = 2) +
    theme_minimal() +
    theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8),
          axis.title = element_text(size = 12), plot.title = element_text(size = 14), strip.text = element_text(size = 12, face = "bold")) +
    labs(title = "Single Cell UMAP: Normal vs Tumor Tissue Comparison", subtitle = "Cell types colored consistently across tissues", color = "Cell Type", x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1))
  ggsave(file.path(output_dir, "normal_vs_tumor_celltype.png"), p_celltype, width = 16, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "normal_vs_tumor_celltype.pdf"), p_celltype, width = 16, height = 8)
  # Top 3 genes from intersector (by meta score)
  meta_score_files <- list.files(intersector_path, pattern = "meta_scores.csv", recursive = TRUE, full.names = TRUE)
  ctnorm_global_metabolic_files <- meta_score_files[grepl("ctnorm_global/all_clusters/metabolic", meta_score_files)]
  all_meta_scores <- map_dfr(ctnorm_global_metabolic_files, ~ {
    tryCatch({
      read_csv(.x, show_col_types = FALSE) %>% mutate(source_file = basename(dirname(.x)))
    }, error = function(e) {
      message("Error reading ", .x, ": ", e$message)
      NULL
    })
  })
  if (nrow(all_meta_scores) == 0) {
    warning("No meta scores loaded from intersector for top gene selection.")
    return(list(celltype_plot = p_celltype, gene_umap_plots = list()))
  }
  # Pick top 3 genes by meta_score (highest absolute value, or as appropriate)
  top3 <- all_meta_scores %>% arrange(desc(abs(meta_score))) %>% filter(!is.na(gene)) %>% distinct(gene, .keep_all = TRUE) %>% head(3)
  gene_umap_plots <- list()
  for (i in seq_len(nrow(top3))) {
    gene <- top3$gene[i]
    if (!(gene %in% rownames(cp10k_matrix))) {
      warning(sprintf("Gene %s not found in cp10k matrix, skipping.", gene))
      next
    }
    gene_expr <- as.numeric(cp10k_matrix[gene, plot_data$cell_ID])
    plot_data$gene_expr <- gene_expr
    for (tissue_type in c("Normal", "Tumor")) {
      pd <- plot_data %>% filter(tissue == tissue_type)
      if (nrow(pd) == 0) next
      p_gene <- ggplot(pd, aes(x = UMAP1, y = UMAP2, color = gene_expr)) +
        geom_point(size = 0.4, alpha = 0.7) +
        scale_color_viridis_c(name = sprintf("%s expression", gene), option = "plasma") +
        theme_minimal() +
        labs(title = sprintf("UMAP: %s expression in %s", gene, tissue_type), x = "UMAP 1", y = "UMAP 2") +
        theme(legend.position = "right", plot.title = element_text(size = 13), axis.title = element_text(size = 11))
      fname <- file.path(output_dir, sprintf("umap_%s_%s.png", gene, tolower(tissue_type)))
      ggsave(fname, p_gene, width = 8, height = 7, dpi = 300)
      fname_pdf <- file.path(output_dir, sprintf("umap_%s_%s.pdf", gene, tolower(tissue_type)))
      ggsave(fname_pdf, p_gene, width = 8, height = 7)
      gene_umap_plots[[paste0(gene, "_", tissue_type)]] <- p_gene
    }
  }
  return(list(celltype_plot = p_celltype, gene_umap_plots = gene_umap_plots))
}

# Main function (cp10k version)
generate_cp10k_umap_plots <- function(output_dir = "output/6. plots/figure 8", max_cells = 50000, intersector_path = "output/3. intersector") {
  dir_create(output_dir, recurse = TRUE)
  # Load cp10k-normalized data and metadata
  datasets <- load_cp10k_datasets()
  cell_metadata <- datasets$cell_metadata
  # Run UMAP (all cells, no sampling)
  plot_data <- run_cp10k_umap(datasets$cp10k, cell_metadata, max_cells)
  # Create plots (cell type + top 3 genes)
  plots <- create_cp10k_umap_plots(plot_data, datasets$cp10k, cell_metadata, intersector_path, output_dir)
  # Save data
  write_csv(plot_data, file.path(output_dir, "umap_plot_data.csv"))
  message("UMAP analysis complete!")
  message(sprintf("Plots saved to: %s", output_dir))
  return(list(plot_data = plot_data, plots = plots))
}

# Run the analysis (cp10k version)
if (interactive() || !exists(".umap_sourced")) {
  .umap_sourced <- TRUE
  umap_results <- generate_cp10k_umap_plots(output_dir = "output/6. plots/figure 8", max_cells = 50000, intersector_path = "output/3. intersector")
}
