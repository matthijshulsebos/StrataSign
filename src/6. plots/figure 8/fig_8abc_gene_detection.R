# Minimal detection plotting script
library(ggplot2)
library(dplyr)
library(readr)
library(stringr) # For string manipulation


# Load detection rates
detection_df <- read_csv("output/6. plots/data/detection/detection_rates.csv", show_col_types = FALSE)

# Ensure correct types and no hidden whitespace
detection_df <- detection_df %>%
  mutate(
    cell_type = as.character(cell_type),
    cluster_ID = as.character(cluster_ID)
  )


# Get all unique (cell_type, cluster_ID) pairs
all_clusters <- detection_df %>% distinct(cell_type, cluster_ID) %>% arrange(cell_type, cluster_ID)

# Get max samples for consistent axes
max_normal <- max(detection_df$n_normal, na.rm = TRUE)
max_tumor <- max(detection_df$n_tumor, na.rm = TRUE)

# Plot and save detection plot
plot_and_save_detection <- function(df, outdir, prefix, ct, clust, suffix, max_normal, max_tumor) {
  if (nrow(df) == 0) return()
  
  # Robustly calculate bin counts and max_count for the color scale
  plot_data <- df %>%
    group_by(n_detected_normal, n_detected_tumor) %>%
    summarise(count = n(), .groups = 'drop')
  
  max_count <- max(plot_data$count, na.rm = TRUE)
  
  # Determine breaks for the color scale legend
  if (max_count > 1) {
    breaks <- scales::log_breaks(n = 5)(c(1, max_count))
    breaks <- round(breaks[breaks >= 1 & breaks <= max_count])
    if (length(breaks) < 2) {
      breaks <- waiver()
    }
  } else {
    breaks <- waiver()
  }

  p <- ggplot(plot_data, aes(x = n_detected_normal, y = n_detected_tumor, fill = count)) +
    geom_tile(color = NA) +
    scale_fill_viridis_c(
      option = "C", 
      trans = "log1p", 
      name = "# genes",
      breaks = breaks,
      labels = scales::label_number(scale_cut = scales::cut_si("")),
      guide = guide_colorbar(background = element_rect(fill = "transparent", color = NA))
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    coord_fixed(xlim = c(0, max_normal), ylim = c(0, max_tumor)) +
    theme_bw(base_size = 14, base_family = "sans") +
    theme(
        plot.subtitle = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        strip.text = element_text(color = "black", face = "bold"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "# normal samples",
      y = "# tumor samples",
      subtitle = paste0("Sublineage: ", ct, " (", clust, ")"),
      fill = "# genes"
    )
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outname <- file.path(outdir, paste0(prefix, "_gene_detection_bin2d_", gsub("[ /]", "_", tolower(ct)), "_", clust, "_", suffix, ".png"))
  ggsave(outname, p, width = 6, height = 6, dpi = 300, bg = "transparent")
}

# Main plotting loop
for (i in seq_len(nrow(all_clusters))) {
  ct <- all_clusters$cell_type[i]
  clust <- all_clusters$cluster_ID[i]
  
  df_sub <- detection_df[which(detection_df$cell_type == ct & detection_df$cluster_ID == clust), ]

  # Remove rows with NA counts before plotting to prevent warnings
  df_sub <- df_sub[complete.cases(df_sub[, c("n_detected_normal", "n_detected_tumor")]), ]

  prefix <- paste0("fig_s6_", gsub("[ /]", "_", tolower(ct)), "_", clust)
  
  # Define output directories for unfiltered and important genes
  detection_root <- "output/6. plots/figure 8/detection"
  outdir_unfiltered <- file.path(detection_root, "unfiltered")
  outdir_important <- file.path(detection_root, "important")

  # All genes
  plot_and_save_detection(df_sub, outdir_unfiltered, prefix, ct, clust, "all", max_normal, max_tumor)

  # Important genes
  norm_method <- "ctnorm_global"
  gene_type_set <- "metabolic"
  cell_type_dirname <- "all_clusters"

  meta_scores_path <- file.path("output/3. intersector", norm_method, cell_type_dirname, gene_type_set, "meta_scores.csv")

  if (file.exists(meta_scores_path)) {
    meta_scores <- read_csv(meta_scores_path, show_col_types = FALSE)

    # Ensure required columns exist before proceeding
    if (all(c("feature_id", "sublineage") %in% names(meta_scores))) {
      # Extract cluster_ID from the feature_id string but we should use utils
      meta_scores <- meta_scores %>%
        mutate(cluster_ID_extracted = as.character(str_extract(feature_id, "\\d+$")))

      # Filter using the correct sublineage and the extracted cluster_ID
      important_features <- meta_scores[which(meta_scores$sublineage == ct & meta_scores$cluster_ID_extracted == clust), "feature_id", drop = TRUE]

      if (length(important_features) > 0) {
        df_filt <- df_sub[which(df_sub$feature_id %in% important_features), ]
        if (nrow(df_filt) > 0) {
          plot_and_save_detection(df_filt, outdir_important, prefix, ct, clust, "important", max_normal, max_tumor)
        }
      }
    }
  }
}

message("Completed writing figure 8 detection plots to file.")
