library(ggplot2)
library(dplyr)
library(readr)
source("src/0. utils/feature_name_utils.R")

# Define combinations
norms <- c("ctnorm_global", "ctnorm_relative", "read_depth")
cell_types <- c("all_clusters", "macrophages", "lcam_hi", "lcam_lo", "lcam_both")
gene_sets <- c("metabolic", "nonmetabolic", "random")

# Load sublineage color mapping
sublineage_color_map_path <- "output/3. intersector/sublineage_colors.rds"
sublineage_colors <- readRDS(sublineage_color_map_path)
sublineage_colors_plot <- sublineage_colors

for (norm in norms) {
  for (cell_type in cell_types) {
    for (gene_set in gene_sets) {
      meta_scores_file <- file.path("output/3. intersector", norm, cell_type, gene_set, "meta_scores.csv")
      if (!file.exists(meta_scores_file)) next
      all_results <- read_csv(meta_scores_file, show_col_types = FALSE)

      # Select top 60 features by meta_score
      top_features_data <- all_results %>%
        arrange(desc(meta_score)) %>%
        head(60)

      # Use utils to extract gene and cluster_id
      parsed_features <- parse_feature_identifier(top_features_data$feature_id)
      top_features_data$gene <- parsed_features$gene

      # Use get_simplified_sublineage from utils
      top_features_data$sublineage <- get_simplified_sublineage(top_features_data$feature_id)

      # Also extract cluster_id for display
      top_features_data$cluster_id <- parsed_features$cluster_id

      # Display name with fold change direction
      top_features_data <- top_features_data %>%
        mutate(
          display_name = paste0(
            gene, "@", sublineage, " (", cluster_id, ")",
            ifelse(fold_change > 0, " \u2191", " \u2193")
          ),
          n_models_occur_f = as.factor(n_models_occur)
        )

      # Plot in chunks of 15 or less for the last chunk
      chunk_size <- 15
      n_chunks <- ceiling(nrow(top_features_data) / chunk_size)
      for (chunk in seq_len(n_chunks)) {
        start_idx <- (chunk - 1) * chunk_size + 1
        end_idx <- min(chunk * chunk_size, nrow(top_features_data))
        chunk_data <- top_features_data[start_idx:end_idx, ]
        p <- ggplot(chunk_data, aes(
            x = reorder(display_name, meta_score), 
            y = meta_score, 
            fill = sublineage
          )) +
          geom_col(color = "black", linewidth = 0.7) +
          geom_text(aes(label = n_models_occur, y = meta_score + max(meta_score) * 0.02), 
                    color = "black", fontface = "plain", size = 3.5, hjust = 0) +
          scale_fill_manual(values = sublineage_colors_plot, name = "Sublineage", na.value = "grey60") +
          guides(
            fill = guide_legend(title = "Sublineage")
          ) +
          labs(x = "Feature", y = "Meta score") +
          coord_flip() +
          theme_bw(base_size = 16, base_family = "sans") +
          theme(
            plot.title = element_text(face = "plain", hjust = 0.5, size = 18, color = "black"),
            plot.subtitle = element_text(face = "plain", hjust = 0.5, size = 15, color = "black"),
            axis.title = element_text(face = "plain", size = 16, color = "black"),
            axis.text = element_text(size = 13, color = "black"),
            axis.text.y = element_text(size = 11, hjust = 1, color = "black"),
            legend.position = "right",
            legend.background = element_rect(fill = NA, color = NA),
            legend.title = element_text(size = 13, face = "plain", color = "black"),
            legend.text = element_text(size = 12, color = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", size = 0.8),
            axis.line.x = element_line(color = "black", size = 0.8),
            axis.line.y = element_line(color = "black", size = 0.8),
            plot.margin = margin(10, 15, 10, 10, "pt"),
            strip.text = element_text(color = "black", face = "bold"),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)
          )

        # Safe names
        safe_norm <- gsub("[^a-zA-Z0-9_]+", "_", norm)
        safe_celltype <- gsub("[^a-zA-Z0-9_]+", "_", cell_type)
        safe_geneset <- gsub("[^a-zA-Z0-9_]+", "_", gene_set)

        # Save to figure
        top_scores_root <- file.path("output/6. plots/figure 7/top_scores")
        combo_dir <- file.path(top_scores_root, safe_norm, safe_celltype, safe_geneset)
        dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
        fname <- sprintf("barplot_chunk%d.png", chunk)
        ggsave(file.path(combo_dir, fname), plot = p, width = 8, height = 6, dpi = 300, units = "in", bg = "transparent", limitsize = FALSE)
      }
    }
  }
}
