library(fs)
library(purrr)
# library(tibble) # tibble is loaded in heatmap_data_utils.R

source("src/6. plots/heatmap/heatmap_data_utils.R")
source("src/6. plots/heatmap/complex_heatmap_wrapper.R")

create_all_heatplots <- function(
  models_base_dir  = "output/2. models",
  output_base_dir  = "output/6. plots/heatmap",
  n_top_genes      = 50 # Consistent with original
) {
  dir_create(output_base_dir, recurse = TRUE)

  fc_list <- tryCatch(load_fold_changes(), error = function(e){
    message("No FC data: ", e$message); return(NULL)
  })
  fold_changes <- if (!is.null(fc_list)) fc_list$fold_changes else NULL

  dataset_types <- dir_ls(models_base_dir, type="directory") %>% path_file()
  if (length(dataset_types) == 0) {
    message(paste0("No dataset type directories found in ", models_base_dir))
    return()
  }

  for (current_dataset_type in dataset_types) {
    message(paste0("\nProcessing dataset type: ", current_dataset_type))
    dataset_type_input_path <- file.path(models_base_dir, current_dataset_type)
    
    all_fi_csv_files <- dir_ls(dataset_type_input_path,
                               recurse = TRUE,
                               regexp = "feature_importance.*\\.csv$",
                               type = "file",
                               fail = FALSE)
    if (length(all_fi_csv_files) == 0) {
      message(paste0("  No feature importance CSV files found for dataset type: ", current_dataset_type))
      next
    }

    walk(all_fi_csv_files, function(csv_file_path){
      message(paste0("  Processing file: ", csv_file_path))

      # Parse model name, cluster info, and gene set from the file path (ADAPTED FROM ORIGINAL)
      relative_path_from_dataset_type <- path_rel(csv_file_path, dataset_type_input_path)
      path_parts <- path_split(relative_path_from_dataset_type)[[1]]
      
      model_algorithm <- "unknown_model"
      if (length(path_parts) > 1) {
        model_algorithm <- path_parts[1]
      } else {
         message(paste0("Could not determine model algorithm from path structure for: ", csv_file_path))
      }

      file_basename_no_ext <- path_file(csv_file_path) %>% path_ext_remove()
      parse_string <- sub("^feature_importance_", "", file_basename_no_ext)
      
      cluster_info_parsed <- "all_clusters"; gene_set_parsed <- "metabolic"; # Defaults

      # Simplified parsing logic based on original - adjust if more complex patterns are needed
      if (grepl("lcam_hi", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "lcam_hi"
      } else if (grepl("lcam_lo", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "lcam_lo"
      } else if (grepl("lcam_both", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "lcam_both"
      } else if (grepl("macrophages", parse_string, ignore.case = TRUE)) {
        cluster_info_parsed <- "macrophages"
      } # all_clusters is default if none of these match more specific parts of the name

      if (grepl("nonmetabolic", parse_string, ignore.case = TRUE)) {
        gene_set_parsed <- "nonmetabolic"
      } else if (grepl("random", parse_string, ignore.case = TRUE)) {
        gene_set_parsed <- "random"
      } # metabolic is default

      message(paste0("    Parsed Info: model_algo=", model_algorithm, 
                     ", cluster_info=", cluster_info_parsed, 
                     ", gene_set=", gene_set_parsed))

      specific_output_dir <- file.path(output_base_dir, current_dataset_type, model_algorithm, gene_set_parsed, cluster_info_parsed)
      dir_create(specific_output_dir, recurse = TRUE)
      
      plot_file_identifier <- file_basename_no_ext # Used for naming plot output files

      out_fi    <- file.path(specific_output_dir, paste0("fi_",    plot_file_identifier,".png"))
      out_fc    <- file.path(specific_output_dir, paste0("fc_",    plot_file_identifier,".png"))
      out_filtered_fi <- file.path(specific_output_dir, paste0("filtered_fi_", plot_file_identifier,".png")) # CHANGED

      tryCatch({
        mat <- prepare_heatmap_data(csv_file_path, n_top_genes)
        if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
          message("    Skipping heatmap generation due to empty or invalid feature matrix.")
          return() # Use return() inside walk's function
        }

        row_hc <- hclust(dist(mat, method="manhattan"), method="complete")
        col_hc <- hclust(dist(t(mat), method="manhattan"), method="complete")

        # Dynamic max clusters (FROM ORIGINAL)
        max_r_clust <- min(5, if(nrow(mat) <= 1) 1 else nrow(mat) -1, if(nrow(mat) < 2) 1 else 5)
        max_c_clust <- min(5, if(ncol(mat) <= 1) 1 else ncol(mat) -1, if(ncol(mat) < 2) 1 else 5)
        
        # Ensure min_clusters is at least 1 and not greater than max_clusters
        min_r_clust <- if (nrow(mat) > 1) 1 else 1 # At least 1 cluster
        min_c_clust <- if (ncol(mat) > 1) 1 else 1

        row_k <- if (nrow(mat) > 1) cut_dendrogram(row_hc, min_clusters = min_r_clust, max_clusters = max_r_clust, return_factor = FALSE) else 1
        col_k <- if (ncol(mat) > 1) cut_dendrogram(col_hc, min_clusters = min_c_clust, max_clusters = max_c_clust, return_factor = FALSE) else 1
        
        fi_max <- max(abs(mat), na.rm = TRUE)
        if (!is.finite(fi_max)) fi_max <- 1 # Fallback from original

        plot_matrix_heatmap(mat, "Feature\nImportance", row_hc, col_hc, row_k, col_k, out_fi, fixed_max = fi_max)
        
        if (!is.null(fold_changes)) {
          fc_mat    <- build_fc_matrix(mat, fold_changes) # from heatmap_data_utils.R
          plot_matrix_heatmap(fc_mat, "Log2 Fold\nChange", row_hc, col_hc, row_k, col_k, out_fc)
          
          filtered_fi_mat <- build_ratio_matrix(mat, fc_mat, threshold=1) # from heatmap_data_utils.R # VARIABLE NAME CHANGED
          plot_matrix_heatmap(filtered_fi_mat, "Filtered FI\n(|log2FC|<1)", row_hc, col_hc, row_k, col_k, out_filtered_fi, fixed_max = fi_max) # TITLE AND OUTPUT PATH CHANGED
        } else {
            message("    Fold change data not available. Skipping FC and Filtered FI heatmaps.") # MESSAGE UPDATED
        }
      }, error = function(e) {
        message(paste0("    Error processing heatmaps for ", csv_file_path, ": ", e$message))
        # try(print(sys.calls())) # Optional for deeper debugging
      })
    })
  }
  message("\nHeatplot generation complete.")
}

create_all_heatplots()
