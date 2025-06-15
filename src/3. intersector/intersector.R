library(dplyr)
library(tidyr)
library(readr)
library(fs)
library(RColorBrewer)
library(glue)
library(purrr)

# Source utility functions
source("src/0. utils/feature_name_utils.R")

# Create empty final table
create_empty_final_table <- function(model_names) {
  score_cols <- paste0(model_names, "_score")
  cols <- c("feature_id","gene","sublineage","meta_score","n_models_occur","fold_change", score_cols)
  
  df <- tibble()
  for (col in cols) {
    if (col %in% c("fold_change")) {
      df[[col]] <- numeric(0)
    } else {
      df[[col]] <- numeric(0)
    }
  }

  return(df[cols])
}

# Color map for sublineages to file
create_color_map <- function(items, output_path) {
  dir_create(dirname(output_path), recurse = TRUE)
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category=='qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (length(col_vector)==0) col_vector <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF")
  col_vector <- rep(col_vector, length.out = length(items))
  color_map  <- setNames(col_vector[seq_along(items)], items)

  saveRDS(color_map, file=output_path)
  message(glue("Color map for {length(items)} items saved to: {output_path}"))

  return(color_map)
}

# Load one file
load_and_extract_features <- function(file_path, model_name) {
  df <- tryCatch(read_csv(file_path, show_col_types=FALSE), error = function(e){
    message("Error reading ",file_path,": ",e$message); NULL
  })
  if (is.null(df) || !"Feature"%in%names(df)||!"Value"%in%names(df)||nrow(df)==0) return(NULL)
  
  df %>%
    transmute(
      feature_id = as.character(Feature),
      gene        = get_gene_from_feature(feature_id),
      sublineage  = get_simplified_sublineage(feature_id),
      abs_value   = abs(Value),
      model_type  = model_name
    ) %>%
    filter(!is.na(gene))
}

# Find and filter files
get_filtered_feature_files <- function(models_dir, target_gene_type, target_cell_type) {
  
  # Look for cell type directories (e.g., all_clusters, macrophages)
  cell_type_dirs <- dir_ls(models_dir, type = "directory")
  
  if (length(cell_type_dirs) == 0) {
    return(list())
  }

  all_found_files_list <- list() # Initialize an empty list to collect file info

  for (cell_type_dir_path in cell_type_dirs) {
    cell_type_name <- path_file(cell_type_dir_path) # Get 'all_clusters' from '.../cp10k/all_clusters'
    
    # Skip if this is not the target cell type
    if (cell_type_name != target_cell_type) {
      next
    }
    
    # Look for gene type directories (e.g., metabolic, nonmetabolic, random)
    gene_type_dirs <- dir_ls(cell_type_dir_path, type = "directory")
    
    for (gene_type_dir_path in gene_type_dirs) {
      gene_type_name <- path_file(gene_type_dir_path) # Get 'metabolic' from '.../all_clusters/metabolic'
      
      # Skip if this is not the target gene type
      if (gene_type_name != target_gene_type) {
        next
      }
      
      # Look for model directories
      model_dirs <- dir_ls(gene_type_dir_path, type = "directory")
      
      for (model_dir_path in model_dirs) {
        model_name_from_path <- path_file(model_dir_path)
        
        if (dir_exists(model_dir_path)) {
          # Define the glob pattern and its regex equivalent
          file_name_glob_pattern <- "feature_importance*.csv"
          regex_pattern_for_list_files <- utils::glob2rx(file_name_glob_pattern)
          
          # Use base R's list.files()
          files_in_model_dir <- list.files(path = model_dir_path, 
                                            pattern = regex_pattern_for_list_files, 
                                            full.names = TRUE, 
                                            recursive = FALSE,
                                            ignore.case = TRUE)
          
          if (length(files_in_model_dir) > 0) {
            for (f_path in files_in_model_dir) {
              all_found_files_list[[length(all_found_files_list) + 1]] <- list(
                file_path = as.character(f_path), 
                model_name = model_name_from_path
              )
            }
          }
        }
      }
    }
  }
  
  return(all_found_files_list) 
}


# Extract all into one data frame
extract_all_features <- function(file_model_list) {
  file_model_list %>% 
    map(~ load_and_extract_features(.x$file_path, .x$model_name)) %>% 
    compact() %>% 
    bind_rows()
}


# Compute model density (if >100 non zero features)
compute_model_density <- function(all_features, model_names, eps=1e-10) {
  map_lgl(model_names, function(m) {
    df <- filter(all_features, model_type==m)
    sum(df$abs_value>eps, na.rm=TRUE) > 100
  }) %>% set_names(model_names)
}


# Scale single model score table
scale_scores_for_model <- function(df, is_dense, threshold=0.7, eps=1e-10) {
  base_name <- paste0(unique(df$model_type), "_score")
  df0 <- select(df, feature_id, gene, sublineage, abs_value) %>%
         mutate(!!base_name := 0)
  cand <- filter(df0, abs_value>eps)
  if (nrow(cand)>0) {
    to_scale <- if (is_dense) {
      ord <- arrange(cand, desc(abs_value))
      cum <- cumsum(ord$abs_value)/sum(ord$abs_value)
      sel <- which(cum<=threshold)
      if (length(sel)==0) sel <- 1
      slice(ord, sel)
    } else cand
    mn <- min(to_scale$abs_value); mx <- max(to_scale$abs_value)
    scaled <- if (mx==mn) rep(0.75, nrow(to_scale)) else 0.5+0.5*(to_scale$abs_value-mn)/(mx-mn)
    upd <- tibble(feature_id=to_scale$feature_id, !!base_name:=scaled)
    df0 <- left_join(df0, upd, by="feature_id", suffix=c(".old",".new")) %>%
           mutate(!!base_name:=coalesce(.data[[paste0(base_name,".new")]], .data[[paste0(base_name,".old")]])) %>%
           select(-ends_with(".old"), -ends_with(".new"))
  }
  df0 %>% select(feature_id, gene, sublineage, !!base_name)
}


# Build model score table
build_model_score_tables <- function(all_feat, density_map) {
  models <- names(density_map)
  map(models, function(m) {
    dfm <- filter(all_feat, model_type==m)
    scale_scores_for_model(dfm, density_map[[m]])
  })
}


# Merge and compute meta score with fold change integration
merge_and_score <- function(score_tables, fold_changes_data = NULL) {
  merged <- reduce(score_tables, full_join, by=c("feature_id","gene","sublineage"))
  score_cols <- grep("_score$", names(merged), value=TRUE)
  
  # Add fold change information
  if (!is.null(fold_changes_data)) {
    merged <- merged %>%
      left_join(fold_changes_data %>% select(Feature, Value), 
                by = c("feature_id" = "Feature")) %>%
      rename(fold_change = Value)
  } else {
    merged <- merged %>%
      mutate(fold_change = NA_real_)
  }
  
  merged %>%
    mutate(across(all_of(score_cols), ~coalesce(.x,0))) %>%
    mutate(
      meta_score    = rowSums(select(., all_of(score_cols))),
      n_models_occur= rowSums(select(., all_of(score_cols))>0)
    )
}


# Filter meta scores and write to file
filter_and_write <- function(df, min_models, all_models, out_dir) {
  df_f <- filter(df, n_models_occur>=min_models) %>% arrange(desc(meta_score))
  if (nrow(df_f)==0) {
    message("No features met min_models_occurrence of ",min_models)
    empty <- create_empty_final_table(all_models)
    write_csv(empty, file.path(out_dir,"meta_scores.csv"))
    return(empty)
  }
  score_cols <- grep("_score$", names(df_f), value=TRUE) 
  cols <- c("feature_id","gene","sublineage","meta_score","n_models_occur","fold_change", sort(score_cols))
  write_csv(select(df_f, all_of(cols)), file.path(out_dir,"meta_scores.csv"))
  df_f[cols]
}


# Main intersector function
find_feature_intersection <- function(models_dir, output_dir, min_models_occurrence, cell_type_filter_val, gene_type_filter_val, fold_changes_data = NULL) {
  dir_create(output_dir, recurse=TRUE)

  # Pass target_gene_type and target_cell_type to the updated function
  fm_list <- get_filtered_feature_files(models_dir, gene_type_filter_val, cell_type_filter_val)
  
  if (length(fm_list)==0) {
    message("No files for gene=",gene_type_filter_val," cell=",cell_type_filter_val)
    return(NULL)
  }
  all_feat <- extract_all_features(fm_list)
  
  models   <- sort(unique(all_feat$model_type))

  if (nrow(all_feat)==0) {
    empty <- create_empty_final_table(models)
    write_csv(empty, file.path(output_dir,"meta_scores.csv"))
    return(empty)
  }

  density_map <- compute_model_density(all_feat, models)
  score_tables<- build_model_score_tables(all_feat, density_map)
  merged <- merge_and_score(score_tables, fold_changes_data)
  result <- filter_and_write(merged, min_models_occurrence, models, output_dir)

  return(result)
}


# Process a single dataset type, cell type, and gene type combination
process_dataset_cell_gene_combination <- function(type, cell_type, gene_type, base_models_dir, base_output_dir, min_models_param, fold_changes_data = NULL) {
  models_dir <- file.path(base_models_dir, type)
  output_dir <- file.path(base_output_dir, type, cell_type, gene_type)
  
  if (!dir.exists(models_dir)) {
    warning(glue("Models directory not found: {models_dir}."))
    return(NULL)
  }
  
  message(glue("Processing dataset type: {type}, cell type: {cell_type}, gene type: {gene_type}"))
  result <- find_feature_intersection(
    models_dir = models_dir,
    output_dir = output_dir,
    min_models_occurrence = min_models_param,
    cell_type_filter_val = cell_type,
    gene_type_filter_val = gene_type,
    fold_changes_data = fold_changes_data
  )
  
  if (is.null(result) || nrow(result) == 0) {
    message(glue("No results for dataset type: {type}, cell type: {cell_type}, gene type: {gene_type}"))
    return(NULL)
  }
  
  return(file.path(output_dir, "meta_scores.csv"))
}


# Generate sublineage color map
generate_sublineage_color_map <- function(meta_score_files, output_path) {
  sublineages <- meta_score_files %>%
    keep(file.exists) %>%
    map_dfr(~ read_csv(.x, show_col_types = FALSE) %>% select(sublineage)) %>%
    pull(sublineage) %>%
    unique() %>%
    sort() %>%
    discard(is.na)
  
  if (length(sublineages) > 0) {
    create_color_map(sublineages, output_path)
    message(glue("Sublineage color map created for {length(sublineages)} unique sublineages."))
  } else {
    create_color_map(character(0), output_path)
    message("No unique sublineages found. Saved an empty color map.")
  }
}


# Main function to run intersector for all combinations
run_intersector <- function(
  base_models_dir = "output/2. models",
  base_output_dir = "output/3. intersector",
  min_models_param = 2,
  dataset_types_to_process = c("ctnorm_global", "ctnorm_relative", "read_depth"),
  cell_types_to_process = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
  gene_types_to_process = c("metabolic", "nonmetabolic", "random"),
  sublineage_color_map_path = "output/3. intersector/sublineage_colors.rds",
  fold_changes_file = "output/4. fold changes/feature_fold_changes.csv"
) {
  
  # Load fold changes data
  fold_changes_data <- read_csv(fold_changes_file, show_col_types = FALSE)
  
  meta_score_files <- expand.grid(
    type = dataset_types_to_process,
    cell_type = cell_types_to_process,
    gene_type = gene_types_to_process,
    stringsAsFactors = FALSE
  ) %>%
    pmap(~ process_dataset_cell_gene_combination(..1, ..2, ..3, base_models_dir, base_output_dir, min_models_param, fold_changes_data)) %>%
    compact()
  
  generate_sublineage_color_map(meta_score_files, sublineage_color_map_path)
  message("Finished intersector analysis.")
}


# Run the intersector
run_intersector(
  base_models_dir = "output/2. models",
  base_output_dir = "output/3. intersector",
  min_models_param = 2,
  cell_types_to_process = c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both"),
  gene_types_to_process = c("metabolic", "nonmetabolic", "random"),
  sublineage_color_map_path = "output/3. intersector/sublineage_colors.rds"
)
