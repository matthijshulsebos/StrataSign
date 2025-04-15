# Script for investigating expression levels across tissue types

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

#==============================================================================
# CORE FUNCTIONS
#==============================================================================

# Creates a histogram to visualize expression distribution
create_histogram <- function(data, x, fill = "#3498db", title, subtitle, xlabel) {
  actual_column <- deparse(substitute(x))
  message("Histogram data summary for ", subtitle, ":")
  message("  Total sample count: ", nrow(data))
  message("  Number of zero values: ", sum(data[[actual_column]] == 0))
  message("  Percentage of zeros: ", round(sum(data[[actual_column]] == 0)/nrow(data)*100, 1), "%")
  
  p <- ggplot(data, aes(x = {{x}})) +
    geom_histogram(bins = 30, fill = fill) +
    geom_rug(alpha = 0.4) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlabel,
      y = "Count"
    ) +
    theme_minimal()
  
  # Apply log scale if the data is highly skewed
  if (max(data[[actual_column]], na.rm = TRUE) > 20 && 
      median(data[[actual_column]], na.rm = TRUE) < max(data[[actual_column]], na.rm = TRUE)/10) {
    message("  Using log scale since data is highly skewed")
    p <- p + scale_x_continuous(trans = "log1p")
  }
  
  return(p)
}

# Creates a barplot to compare expression between tissue types
create_barplot <- function(data, x, y, fill, title, subtitle, ylabel) {
  p <- ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    stat_summary(fun = mean, geom = "bar") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.2) +
    scale_fill_manual(values = c("NORMAL" = "#3498db", "TUMOR" = "#e74c3c")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Tissue Type",
      y = ylabel
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Saves a plot to a file in the output directory
save_plot <- function(p, filename, width = 7, height = 5, dpi = 300) {
  dir.create("output/figures/differential_expression", recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0("output/figures/differential_expression/", filename),
         p, width = width, height = height, dpi = dpi)
}

#==============================================================================
# MAIN FUNCTION
#==============================================================================

main <- function() {
  # Set up the output directory
  dir.create("output/figures/differential_expression", recursive = TRUE, showWarnings = FALSE)
  
  # Define the features we want to investigate
  feature_info <- data.frame(
    Feature = c("SARDH@T_activated_44", "VKORC1L1@MoMac-II_42"),
    gene = c("SARDH", "VKORC1L1"),
    cluster_name = c("T_activated", "MoMac-II"),
    cluster_ID = c(44, 42)
  )
  
  # Load the Leader dataset
  message("Loading Raw Leader data...")
  load("base/data/lung_ldm.rd")
  table_s1 <- read_csv("base/input_tables/table_s1_sample_table.csv") %>%
    mutate(tissue = toupper(tissue), sample_ID = as.character(sample_ID)) %>%
    filter(Use.in.Clustering.Model. == "Yes")
  
  # Load preprocessed metabolic training data
  message("Loading metabolic training data...")
  metabolic_path <- "output/training_data/raw/all_clusters/metabolic/X_train_all_clusters_metabolic.csv"
  y_train_path <- "output/training_data/raw/all_clusters/metabolic/y_train_all_clusters_metabolic.csv"
  
  if (file.exists(metabolic_path) && file.exists(y_train_path)) {
    metabolic_data <- read_csv(metabolic_path)
    tissue_data <- read_csv(y_train_path)
    metabolic_data$tissue <- toupper(tissue_data$x)
    message("Added tissue information from y_train file")
  } else {
    message("Metabolic data files not found")
    metabolic_data <- NULL
  }
  
  # Process each feature one by one
  for (i in 1:nrow(feature_info)) {
    gene_name <- feature_info$gene[i]
    cluster_id <- feature_info$cluster_ID[i]
    feature <- feature_info$Feature[i]
    cluster_name <- feature_info$cluster_name[i]
    feature_clean <- gsub("@", "_", feature)
    
    message(paste("\nProcessing feature:", feature))
    
    # Get raw counts for this gene and cluster
    counts <- lung_ldm$dataset$counts[, gene_name, as.character(cluster_id)]
    
    # Skip if no data is found
    if (is.null(counts)) {
      message("No data found for this feature")
      next
    }
    
    # Create a filtered dataset with only samples used in clustering
    raw_data <- data.frame(
      sample_ID = names(counts),
      count = as.numeric(counts)
    ) %>%
      inner_join(table_s1 %>% select(sample_ID, tissue), by = "sample_ID") %>%
      rename(condition = tissue)
    
    # Show how many samples were kept after filtering
    message("Samples after filtering: ", nrow(raw_data), " (out of ", length(counts), " total)")
    message("Raw data has ", round(mean(raw_data$count == 0)*100, 1), "% zero values")
    
    # Calculate statistics for raw data by tissue type
    tissue_stats <- raw_data %>%
      group_by(condition) %>%
      summarize(
        samples = n(),
        zeros = sum(count == 0),
        percent_zeros = round(100 * zeros/samples, 1),
        mean_all = round(mean(count), 3),
        mean_nonzero = round(mean(count[count > 0]), 3)
      )
    
    message("\nRaw data statistics by tissue type:")
    print(tissue_stats)
    
    # Show all raw count values by sample
    message("\nRAW COUNT VALUES FOR ", gene_name, " IN CLUSTER ", cluster_id, ":")
    message("==========================================================")
    
    # Sort values by tissue type and count
    sample_values <- raw_data %>%
      arrange(condition, desc(count))
    
    # Show normal samples first
    normal_samples <- sample_values %>% filter(condition == "NORMAL")
    message("\nNORMAL SAMPLES (", nrow(normal_samples), "):")
    if(nrow(normal_samples) > 0) {
      for(j in 1:nrow(normal_samples)) {
        message(sprintf("Sample ID: %-15s Value: %s", 
                       normal_samples$sample_ID[j], 
                       normal_samples$count[j]))
      }
    } else {
      message("No NORMAL samples found")
    }
    
    # Then show tumor samples
    tumor_samples <- sample_values %>% filter(condition == "TUMOR")
    message("\nTUMOR SAMPLES (", nrow(tumor_samples), "):")
    if(nrow(tumor_samples) > 0) {
      for(j in 1:nrow(tumor_samples)) {
        message(sprintf("Sample ID: %-15s Value: %s", 
                       tumor_samples$sample_ID[j], 
                       tumor_samples$count[j]))
      }
    } else {
      message("No TUMOR samples found")
    }
    
    message("==========================================================")
    
    # Create plots for raw expression data
    raw_barplot <- create_barplot(
      data = raw_data, 
      x = condition, 
      y = count, 
      fill = condition,
      title = paste("Raw Expression of", gene_name, "in", cluster_name, "cells"),
      subtitle = "From Leader et al. dataset",
      ylabel = "Mean Raw Expression"
    )
    
    raw_hist <- create_histogram(
      data = raw_data,
      x = count,
      title = paste("Distribution of Raw Expression Values for", gene_name),
      subtitle = paste("In", cluster_name, "cells"),
      xlabel = "Raw Expression"
    )
    
    # Save the raw data plots
    save_plot(raw_barplot, paste0(feature_clean, "_raw_expression_barplot.png"))
    save_plot(raw_hist, paste0(feature_clean, "_raw_expression_hist.png"))
    
    # Process corresponding metabolic data if available
    if (!is.null(metabolic_data) && feature %in% colnames(metabolic_data)) {
      message("Creating plots for metabolic training data...")
      
      # Extract relevant columns from metabolic dataset
      meta_data <- metabolic_data %>%
        select(tissue, !!feature) %>%
        rename(expression = !!feature)
      
      # Calculate statistics for metabolic data by tissue type
      meta_tissue_stats <- meta_data %>%
        group_by(tissue) %>%
        summarize(
          samples = n(),
          zeros = sum(expression == 0),
          percent_zeros = round(100 * zeros/samples, 1),
          mean_all = round(mean(expression), 3),
          mean_nonzero = round(mean(expression[expression > 0]), 3)
        )

      message("\nMetabolic data statistics by tissue type:")
      print(meta_tissue_stats)
      
      # Create plots for metabolic data
      meta_barplot <- create_barplot(
        data = meta_data, 
        x = tissue, 
        y = expression, 
        fill = tissue,
        title = paste("Expression of", feature, "in Metabolic Features Dataset"),
        subtitle = "From X_train_all_clusters_metabolic.csv",
        ylabel = "Mean Expression Value"
      )
      
      meta_hist <- create_histogram(
        data = meta_data,
        x = expression,
        title = paste("Distribution of", feature, "in Metabolic Features Dataset"),
        subtitle = "From X_train_all_clusters_metabolic.csv", 
        xlabel = "Expression Value"
      )
      
      # Save the metabolic data plots
      save_plot(meta_barplot, paste0("metabolic_", feature_clean, "_barplot.png"))
      save_plot(meta_hist, paste0("metabolic_", feature_clean, "_hist.png"))
    }
  }
  
  message("Processing complete. Results saved to output/figures/differential_expression/")
}

# Run main function
main()