library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(viridis)
library(stringr)
library(patchwork)

# Helper function to safely save plots
safe_save_plot <- function(plot, filename, width = 8, height = 6) {
  if (is.null(plot)) {
    message("Error in safe_save_plot: Plot object is NULL. Cannot save.")
    return()
  }
  if (!inherits(plot, "ggplot")) {
    message("Error in safe_save_plot: Object to save is not a ggplot object. Cannot save.")
    return()
  }

  full_path <- normalizePath(filename, mustWork = FALSE)
  dir_creation_success <- dir.create(dirname(full_path), recursive = TRUE, showWarnings = FALSE)
  
  if (dir_creation_success) {
    message("Directory ensured/created: ", dirname(full_path))
  } else {
    # This case might not be hit if showWarnings is FALSE and dir already exists,
    # but good to have a conceptual check.
    # A more robust check is to see if the directory now exists.
    if (!dir.exists(dirname(full_path))) {
        message("Error in safe_save_plot: Failed to create directory: ", dirname(full_path))
        return()
    } else {
        message("Directory already existed or was created: ", dirname(full_path))
    }
  }
  
  png_file_path <- paste0(full_path, ".png")
  pdf_file_path <- paste0(full_path, ".pdf")
  
  message("Attempting to save PNG to: ", png_file_path)
  message("Attempting to save PDF to: ", pdf_file_path)
  
  # Save PNG
  tryCatch({
    png(png_file_path, width = width * 300, height = height * 300, res = 300)
    print(plot)
    dev.off()
    message("Successfully saved PNG: ", png_file_path)
  }, error = function(e_png) {
    message("Error saving PNG to ", png_file_path, ": ", e_png$message)
    if (exists("dev.off") && !is.null(dev.list())) dev.off() # Ensure device is closed on error
  })
  
  # Save PDF
  tryCatch({
    pdf(pdf_file_path, width = width, height = height)
    print(plot)
    dev.off()
    message("Successfully saved PDF: ", pdf_file_path)
  }, error = function(e_pdf) {
    message("Error saving PDF to ", pdf_file_path, ": ", e_pdf$message)
    if (exists("dev.off") && !is.null(dev.list())) dev.off() # Ensure device is closed on error
  })
}

# Load and clean model performance data
load_model_performance <- function(file_path = "output/2. models/model_performance.csv") {
  message("Loading model performance data from: ", file_path)
  
  if (!file.exists(file_path)) {
    message("Error in load_model_performance: File not found at: ", file_path)
    return(NULL) # Return NULL if file doesn't exist
  }
  
  perf_data <- read_csv(file_path, show_col_types = FALSE)
  
  if (nrow(perf_data) == 0) {
    message("Warning in load_model_performance: Loaded data has 0 rows from: ", file_path)
    # Still return the empty dataframe for downstream checks if needed
  }
  
  # Clean and filter the data
  clean_data <- perf_data %>%
    filter(Model != "unknown", Dataset != "unknown", Cell_Type != "unknown", Gene_Set != "unknown") %>%
    mutate(
      # Clean model names
      Model = case_when(
        str_to_lower(Model) == "lightgbm" ~ "LightGBM",
        str_to_lower(Model) == "xgboost" ~ "XGBoost",
        str_to_lower(Model) == "randomforest" ~ "Random Forest",
        str_to_lower(Model) == "elasticnet" ~ "Elastic Net",
        str_to_lower(Model) == "spls" ~ "sPLS",
        str_to_lower(Model) == "svm_rbf" ~ "SVM (RBF)",
        str_to_lower(Model) == "svm_linear" ~ "SVM (Linear)",
        TRUE ~ str_to_title(Model)
      ),
      # Clean dataset names for normalization comparison
      Dataset = case_when(
        Dataset == "cp10k" ~ "CP10K (Raw)",
        Dataset == "cp10k_ctnorm" ~ "CP10K (CT-norm)",
        Dataset == "cp10k_ctnorm_relative" ~ "CP10K (CT-norm Rel)",
        Dataset == "proportion_ctnorm" ~ "Proportion (CT-norm)",
        Dataset == "proportion_ctnorm_relative" ~ "Proportion (CT-norm Rel)",
        TRUE ~ Dataset
      ),
      # Clean cell type names
      Cell_Type = case_when(
        Cell_Type == "all_clusters" ~ "All Cell Types",
        Cell_Type == "lcam_both" ~ "LCAM (Both)",
        Cell_Type == "lcam_lo" ~ "LCAM (Low)",
        Cell_Type == "lcam_hi" ~ "LCAM (High)",
        Cell_Type == "macrophages" ~ "Macrophages",
        TRUE ~ str_to_title(Cell_Type)
      ),
      # Clean gene set names
      Gene_Set = case_when(
        Gene_Set == "metabolic" ~ "Metabolic",
        Gene_Set == "nonmetabolic" ~ "Non-metabolic",
        Gene_Set == "random" ~ "Random",
        TRUE ~ str_to_title(Gene_Set)
      )
    )
  
  message(sprintf("Loaded %d model performance records", nrow(clean_data)))
  return(clean_data)
}

# Create normalization method comparison plots
create_normalization_comparison_plots <- function(perf_data, output_dir) {
  
  # Filter for metabolic genes and all clusters
  filtered_data <- perf_data %>%
    filter(Gene_Set == "Metabolic", Cell_Type == "All Cell Types")
  
  if (nrow(filtered_data) == 0) {
    message("No data found for Metabolic genes with All Cell Types for normalization comparison")
    return(list())
  }
  
  # Single plot: AUC by normalization method with barplot and points
  p1 <- filtered_data %>%
    ggplot(aes(x = Dataset, y = AUC, fill = Dataset)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 2) +
    scale_fill_viridis_d() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = "Model AUC by Normalization Method",
      subtitle = "Metabolic genes, All cell types",
      x = "Normalization Method",
      y = "AUC"
    )
  
  safe_save_plot(p1, file.path(output_dir, "normalization_auc_comparison"), width = 10, height = 6)
  
  message("Created normalization comparison plot for metabolic genes with all cell types")
  
  return(list(auc_comparison = p1))
}

# Create model performance summary plots
create_model_performance_summary_plot <- function(perf_data, output_dir, target_dataset = "ctnorm_global") { # Changed default
  # Filter for the specified dataset, all cell types, and metabolic genes
  message(sprintf("Filtering model performance summary for Dataset: %s, Cell Type: All Cell Types, Gene Set: Metabolic", target_dataset))
  
  # Diagnostic: Print unique values before filtering
  if (!is.null(perf_data) && nrow(perf_data) > 0) {
    message("Unique 'Dataset' values in perf_data before filtering:")
    print(unique(perf_data$Dataset))
    message("Unique 'Cell_Type' values in perf_data before filtering:")
    print(unique(perf_data$Cell_Type))
    message("Unique 'Gene_Set' values in perf_data before filtering:")
    print(unique(perf_data$Gene_Set))
  } else {
    message("perf_data is NULL or empty before filtering in create_model_performance_summary_plot.")
    return(list(plot = NULL, data = NULL))
  }

  filtered_data <- perf_data %>%
    filter(
      Dataset == target_dataset, 
      Cell_Type == "All Cell Types", 
      Gene_Set == "Metabolic"
    )
  
  if (nrow(filtered_data) == 0) {
    message(sprintf("No data found for the specified filters for model performance summary (Dataset: %s, Cell_Type: All Cell Types, Gene_Set: Metabolic)", target_dataset))
    return(list(plot = NULL, data = NULL)) # Return NULL for data as well
  }
  
  # Order models by median AUC (or another primary metric like F1_Score)
  # Ensure AUC column exists and has non-NA values for ordering
  if ("AUC" %in% names(filtered_data) && sum(!is.na(filtered_data$AUC)) > 0) {
    model_order <- filtered_data %>%
      group_by(Model) %>%
      summarise(median_auc = median(AUC, na.rm = TRUE)) %>%
      arrange(desc(median_auc)) %>%
      pull(Model)
    
    filtered_data <- filtered_data %>%
      mutate(Model = factor(Model, levels = model_order))
  } else {
    message("AUC not available for ordering models or all AUC values are NA. Models will be ordered alphabetically.")
    # Fallback to alphabetical order if AUC is not suitable for ordering
    model_order <- sort(unique(filtered_data$Model))
    filtered_data <- filtered_data %>%
      mutate(Model = factor(Model, levels = model_order))
  }

  # Select and pivot longer for key performance metrics
  metrics_to_plot <- c("AUC", "Accuracy", "F1_Score", "Precision", "Recall", 
                       "Balanced_Accuracy", "Specificity", "MCC")
  
  # Check if all metric columns exist
  existing_metrics <- metrics_to_plot[metrics_to_plot %in% names(filtered_data)]
  
  if(length(existing_metrics) == 0) {
    message(sprintf("None of the specified metrics (%s) found in the data for model performance summary.", paste(metrics_to_plot, collapse=", ")))
    return(list(plot = NULL, data = filtered_data)) # Return filtered_data for inspection
  }
  
  if(length(existing_metrics) < length(metrics_to_plot)) {
    message(sprintf("Warning: Not all specified metrics found. Plotting available metrics: %s", paste(existing_metrics, collapse=", ")))
  }

  long_data <- filtered_data %>%
    dplyr::select(Model, all_of(existing_metrics)) %>% # Explicitly use dplyr::select
    pivot_longer(cols = all_of(existing_metrics), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = existing_metrics)) # Maintain order

  # Determine if coordinates will be flipped
  will_flip_coords <- length(unique(long_data$Model)) > 5

  # Create faceted bar plot
  p_summary <- long_data %>%
    ggplot(aes(x = Model, y = Value, fill = Model)) +
    geom_col(alpha = 0.8, position = position_dodge(width = 0.9)) 

  if (will_flip_coords) {
    p_summary <- p_summary + 
      geom_text(aes(label = round(Value, 3)), 
                position = position_dodge(width = 0.9), 
                hjust = -0.1, # Adjust for horizontal bars (text to the right)
                size = 2.5) +
      coord_flip()
  } else {
    p_summary <- p_summary +
      geom_text(aes(label = round(Value, 3)), 
                position = position_dodge(width = 0.9), 
                vjust = -0.25, # Adjust for vertical bars (text above)
                size = 2.5)
  }

  p_summary <- p_summary +
    facet_wrap(~ Metric, scales = "free_y", ncol = 3) + 
    scale_fill_viridis_d() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + # Expands y-axis (or new x-axis if flipped) for text
    theme_bw() + 
    theme(
      axis.text.x = if (will_flip_coords) element_text(angle = 0, hjust = 0.5) else element_text(angle = 45, hjust = 1),
      axis.text.y = if (will_flip_coords) element_text(angle = 0, hjust = 1) else waiver(), 
      legend.position = "none", 
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(), 
      panel.grid.major.x = if (will_flip_coords) element_line() else element_blank(), 
      panel.grid.major.y = if (will_flip_coords) element_blank() else element_line()  
    ) +
    labs( # Removed title and subtitle
      x = "Model",
      y = "Metric Value"
    )
  
  # Adjust width based on number of metrics; e.g. 3-4 metrics per row in facet_wrap
  num_metrics <- length(existing_metrics)
  plot_width <- ifelse(num_metrics > 5, 15, 12) # Increase width if many metrics
  plot_height <- ifelse(num_metrics > 3, 10, 8) # Increase height if many metrics in multiple rows

  safe_save_plot(p_summary, file.path(output_dir, sprintf("model_performance_summary_%s_ordered", gsub("[^a-zA-Z0-9_]", "_", target_dataset))), width = plot_width, height = plot_height)
  
  message(sprintf("Created model performance summary plot for Dataset: %s (Models ordered).", target_dataset))
  
  return(list(plot = p_summary, data = filtered_data)) # Return both plot and data
}


# Main function
generate_performance_plots <- function(performance_file = "output/2. models/model_performance.csv",
                                       focused_dataset = "ctnorm_global") { # Changed default
  # Define project root explicitly to ensure correct path construction
  # Adjust this path if your project root is different.
  project_root <- "c:/Users/mchul/Documents/StrataSign" 
  
  # Setup output directory using the explicit project root
  output_dir <- file.path(project_root, "output/6. plots/figure 2")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("Ensuring output directory exists: %s", normalizePath(output_dir, mustWork = FALSE)))

  # Test directory writability
  dummy_file_path <- file.path(output_dir, "write_test.tmp")
  can_write <- tryCatch({
    writeLines("test", dummy_file_path)
    if (file.exists(dummy_file_path)) {
      file.remove(dummy_file_path)
      TRUE
    } else {
      FALSE
    }
  }, error = function(e) {
    FALSE
  })
  
  if (can_write) {
    message("Successfully tested write access to output directory: ", output_dir)
  } else {
    message("Error: Failed to write to output directory: ", output_dir, ". Check permissions.")
    # Optionally, you might want to stop execution if writing is not possible
    # return(NULL) 
  }

  # Adjust performance_file path to be relative to project_root if it's not absolute
  # Assuming performance_file is relative to project_root if not an absolute path
  if (!file.exists(performance_file) && !startsWith(performance_file, "/") && !grepl("^[A-Za-z]:", performance_file)) {
    performance_file_path <- file.path(project_root, performance_file)
  } else {
    performance_file_path <- performance_file
  }
  
  if (!file.exists(performance_file_path)) {
    message(sprintf("Critical Error: Performance data file not found at resolved path: %s. Exiting.", normalizePath(performance_file_path, mustWork = FALSE)))
    return(NULL)
  }
  message(sprintf("Using performance file: %s", normalizePath(performance_file_path, mustWork = FALSE)))
  
  # Load and process data
  perf_data <- load_model_performance(performance_file_path)
  
  if (is.null(perf_data) || nrow(perf_data) == 0) {
    message("No performance data loaded. Exiting.")
    return(NULL)
  }
  
  # Create normalization comparison plots (Commented out as per user focus)
  # message("Creating normalization method comparison plots...")
  # norm_plots <- create_normalization_comparison_plots(perf_data, output_dir)
  
  # Create model performance summary plot for the focused dataset
  message(sprintf("Creating model performance summary plot for dataset: %s, All Cell Types, Metabolic gene set...", focused_dataset))
  summary_output <- create_model_performance_summary_plot(perf_data, output_dir, target_dataset = focused_dataset)
  
  model_summary_plot <- NULL
  if (!is.null(summary_output) && !is.null(summary_output$plot)) {
    model_summary_plot <- summary_output$plot
    
    # Save the filtered data used for this plot
    if (!is.null(summary_output$data) && nrow(summary_output$data) > 0) {
      focused_data_filename <- file.path(output_dir, sprintf("focused_plot_data_%s.csv", gsub("[^a-zA-Z0-9_]", "_", focused_dataset)))
      tryCatch({
        write_csv(summary_output$data, focused_data_filename)
        message(sprintf("Saved focused data for summary plot to: %s", normalizePath(focused_data_filename, mustWork=FALSE)))
      }, error = function(e) {
        message(sprintf("Error saving focused data to CSV: %s", e$message))
      })
    }
  } else {
    message("Model performance summary plot or its data was not generated.")
  }
  
  message("Performance plot generation complete!")
  message(sprintf("Output artifacts (plots, specific data subsets) saved to: %s", output_dir))
  
  # Return list, adjusted to reflect focus
  return(list(data = perf_data, model_summary_plot = model_summary_plot))
}

# Run the analysis
if (interactive() || !exists(".performance_analysis_sourced")) {
  .performance_analysis_sourced <- TRUE # Changed variable name
  # You can change "ctnorm_global" to any other dataset name present in your cleaned data
  # (e.g., "ctnorm_relative", "read_depth") based on the diagnostic output.
  performance_results <- generate_performance_plots(focused_dataset = "ctnorm_global") # Changed to an available dataset
}
