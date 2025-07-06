# Minimal detection plotting script for Sublineage Aggregation
library(ggplot2)
library(dplyr)
library(readr)
library(stringr) # For string manipulation


# Load detection rates
detection_df <- read_csv("output/6. plots/data/detection/detection_rates.csv", show_col_types = FALSE)

# Ensure correct types and no hidden whitespace
detection_df <- detection_df %>%
  mutate(
    cell_type = as.character(cell_type)
  )


# Get all unique cell_types
all_cell_types <- detection_df %>% distinct(cell_type) %>% arrange(cell_type)

# Hardcode axes to user-specified values
max_normal <- 29
max_tumor <- 35

# Helper: Plot and save detection plot for a sublineage
plot_and_save_detection <- function(df, outdir, prefix, ct, suffix, max_normal, max_tumor) {
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
      breaks <- waiver() # Let ggplot decide if we can't make good breaks
    }
  } else {
    breaks <- waiver() # Let ggplot handle cases with no overlapping points
  }

  p <- ggplot(plot_data, aes(x = n_detected_normal, y = n_detected_tumor, fill = count)) +
    geom_tile(color = NA) + # Use geom_tile with pre-counted data
    scale_fill_viridis_c(
      option = "C", 
      trans = "log1p", 
      name = "# genes",
      breaks = breaks,
      labels = scales::label_number(scale_cut = scales::cut_si(""))
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    coord_fixed(xlim = c(0, max_normal), ylim = c(0, max_tumor)) +
    theme_bw(base_size = 14, base_family = "sans") +
    labs(
      x = "# normal samples",
      y = "# tumor samples",
      subtitle = paste0("Sublineage: ", ct),
      fill = "# genes"
    ) +
    theme(
        plot.subtitle = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        strip.text = element_text(color = "black", face = "bold")
    )
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outname <- file.path(outdir, paste0(gsub("[ /]", "_", tolower(ct)), "_detection.png"))
  ggsave(outname, p, width = 6, height = 6, dpi = 300)
}

# Main plotting loop - one plot per cell_type
for (i in seq_len(nrow(all_cell_types))) {
  ct <- all_cell_types$cell_type[i]

  # Filter for the current cell_type. The data to plot is already in the correct columns.
  df_sub <- detection_df %>%
    filter(cell_type == ct) %>%
    filter(complete.cases(n_detected_normal, n_detected_tumor))

  prefix <- paste0("fig_s2_", gsub("[ /]", "_", tolower(ct)))
  outdir_unfiltered <- "output/6. plots/figure s2/detection/"

  # Plot detection for all features (genes) within the sublineage
  # Use the global max_normal and max_tumor for consistent axes across all plots
  if (nrow(df_sub) > 0) {
    plot_and_save_detection(df_sub, outdir_unfiltered, prefix, ct, "all_genes", max_normal, max_tumor)
  }
}

message("Completed writing figure S2 detection plots to file.")
