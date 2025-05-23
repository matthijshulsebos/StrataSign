library(data.table)
library(dplyr)
library(tidyr)
library(tibble)   # ← needed for column_to_rownames()

# Load fold‐change table once
load_fold_changes <- function(path="output/4. differential expression/feature_fold_changes.csv") {
  if (!file.exists(path)) stop("FC file not found: ", path)
  fc <- fread(path, colClasses=c(Feature="character",Value="numeric"))
  list(fold_changes=fc, has_clusters=any(grepl("@", fc$Feature)))
}

# Build FI matrix for top N genes
prepare_heatmap_data <- function(data_path, n_top_genes=50) {
  df <- fread(data_path, colClasses=c(Feature="character", Value="numeric")) %>% as_tibble()
  long <- df %>%
    separate(Feature, into=c("gene","cluster"), sep="@") %>%
    rename(value=Value) %>%
    filter(!is.na(cluster))
  top_genes <- long %>%
    group_by(gene) %>%
    summarise(total = sum(abs(value), na.rm=TRUE)) %>%
    arrange(desc(total)) %>%
    slice_head(n=n_top_genes) %>%
    pull(gene)
  mat <- long %>%
    filter(gene %in% top_genes) %>%
    pivot_wider(names_from=cluster, values_from=value, values_fill=0) %>%
    column_to_rownames("gene")
  as.matrix(mat)
}

# Align FC table to FI matrix
build_fc_matrix <- function(fi_mat, fc_data_table) { # fc_data_table is your main fold_changes data.table
  # Create a matrix initialized with the default value (0.1 from original script)
  output_matrix <- matrix(0.1, 
                          nrow = nrow(fi_mat), 
                          ncol = ncol(fi_mat),
                          dimnames = list(rownames(fi_mat), colnames(fi_mat)))

  # Work on a copy for splitting 'Feature' column if not already split
  fc_copy <- data.table::copy(fc_data_table)
  if (!all(c("gene", "cluster") %in% names(fc_copy))) {
      fc_copy[, c("gene", "cluster") := tstrsplit(Feature, "@", fixed=TRUE)]
  }

  # Filter the fold change data for relevant genes and clusters
  # This is efficient in data.table
  relevant_fc <- fc_copy[gene %in% rownames(fi_mat) & cluster %in% colnames(fi_mat)]

  if (nrow(relevant_fc) > 0) {
    # Update the output_matrix with values from relevant_fc
    # This loop is on a potentially much smaller set (relevant_fc)
    for (i in 1:nrow(relevant_fc)) {
      gene_name <- relevant_fc$gene[i]
      cluster_name <- relevant_fc$cluster[i]
      value <- relevant_fc$Value[i]
      
      # Check if the gene_name and cluster_name are valid for the matrix
      # (they should be, due to the filtering, but good for safety)
      if (gene_name %in% rownames(output_matrix) && cluster_name %in% colnames(output_matrix)) {
        output_matrix[gene_name, cluster_name] <- value
      }
    }
  }
  return(output_matrix)
}

# Mask FI where |log2FC| >= threshold (i.e., keep FI where |log2FC| < threshold)
build_filtered_fi_matrix <- function(fi_mat, fc_mat, threshold=1) { # RENAMED FUNCTION
  mask <- abs(fc_mat) < threshold
  fi_mat * mask # Element-wise multiplication; keeps fi_mat values where mask is TRUE, otherwise becomes 0
}
