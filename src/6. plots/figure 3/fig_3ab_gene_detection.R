library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# Read the detection rates data
detection_df <- read_csv("output/6. plots/data/detection/detection_rates.csv", show_col_types = FALSE)

# We know the number of samples for normal and tumor
max_normal <- 29
max_tumor <- 35


# Plotting function for 2d density of detection rates
plot_and_save_detection <- function(df, outdir, prefix, ct, suffix, max_normal, max_tumor) {
  # Return immediately if no data
  if (nrow(df) == 0) return()
  
  # Plotting data to count genes that have these detection rates
  plot_data <- df %>%
    group_by(n_detected_normal, n_detected_tumor) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Check max counts for legend
  max_count <- max(plot_data$count, na.rm = TRUE)
  
  # Create breaks for the color scale
  if (max_count > 1) {
    breaks <- scales::log_breaks(n = 5)(c(1, max_count))
    breaks <- round(breaks[breaks >= 1 & breaks <= max_count])
    if (length(breaks) < 2) {
      breaks <- waiver()
    }
  } else {
    breaks <- waiver()
  }

  # Create the plot
  p <- ggplot(plot_data, aes(x = n_detected_normal, y = n_detected_tumor, fill = count)) +
    geom_tile(color = NA) +
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
      x = paste0("# normal samples (out of ", as.character(max_normal), ")"),
      y = paste0("# tumor samples (out of ", as.character(max_tumor), ")"),
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
  ggsave(outname, p, width = 6, height = 6, dpi = 300, bg = "transparent")
}


# Ensure cell_type is character for consistent filtering
detection_df <- detection_df %>% mutate(cell_type = as.character(cell_type))

# Get all unique cell types
all_cell_types <- detection_df %>% distinct(cell_type) %>% arrange(cell_type)

# Loop through each cell type and generate plots
for (i in seq_len(nrow(all_cell_types))) {
  # Current cell type
  ct <- all_cell_types$cell_type[i]

  # Filter the detection data for the current cell type
  df_sub <- detection_df %>%
    filter(cell_type == ct) %>%
    filter(complete.cases(n_detected_normal, n_detected_tumor))

  # Set output directory and prefix
  prefix <- paste0("fig_3_", gsub("[ /]", "_", tolower(ct)))
  outdir_unfiltered <- "output/6. plots/figure 3/detection/"

  # If there is data then plot
  if (nrow(df_sub) > 0) {
    plot_and_save_detection(df_sub, outdir_unfiltered, prefix, ct, "all_genes", max_normal, max_tumor)
  }
}

message("Completed writing figure 3 detection plots to file.")
