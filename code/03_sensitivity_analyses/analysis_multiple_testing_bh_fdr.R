#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  input_csv <- file.path(
    repo_root, "outputs", "sensitivity_analyses", "multiple_testing_bh_fdr", "inputs",
    "primary_climate_amr_tests.csv"
  )
  output_root <- file.path(
    repo_root, "outputs", "sensitivity_analyses", "multiple_testing_bh_fdr"
  )
} else if (length(args) == 2L) {
  input_csv <- args[[1]]
  output_root <- args[[2]]
} else {
  stop("Usage: Rscript analysis_multiple_testing_bh_fdr.R [input_csv output_root]")
}
input_csv <- normalizePath(input_csv, mustWork = TRUE)
dirs <- list(
  figure = file.path(output_root, "figures"),
  source = file.path(output_root, "source_data")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

raw <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "test_id", "phenotype", "climate", "lag_years", "EDF", "F_statistic",
  "p_raw", "p_BH_FDR"
)
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

alpha <- 0.05
m <- nrow(raw)
if (m != 24L) stop("Expected 24 prespecified primary tests; found ", m)
if (anyDuplicated(raw[c("phenotype", "climate")])) {
  stop("Phenotype-climate pairs are not unique")
}

raw <- raw %>%
  mutate(
    p_BH_recalculated = p.adjust(p_raw, method = "BH"),
    bh_absolute_difference = abs(p_BH_FDR - p_BH_recalculated)
  )
if (max(raw$bh_absolute_difference, na.rm = TRUE) > 1e-12) {
  stop("Stored BH-adjusted P values do not reproduce p.adjust(method = 'BH')")
}

phenotype_levels <- c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")
climate_levels <- c("Temperature", "Precipitation", "WetDays", "Humidity")
if (!setequal(raw$phenotype, phenotype_levels)) stop("Unexpected phenotype labels")
if (!setequal(raw$climate, climate_levels)) stop("Unexpected climate labels")

status_levels <- c("FDR supported", "Nominal only", "No nominal support")
classify_status <- function(p_raw, p_bh) {
  factor(
    case_when(
      p_bh <= alpha ~ "FDR supported",
      p_raw < alpha ~ "Nominal only",
      TRUE ~ "No nominal support"
    ),
    levels = status_levels
  )
}

palette_fill <- c(
  "FDR supported" = "#3B6F8F",
  "Nominal only" = "#C47A2C",
  "No nominal support" = "#D8D5CE"
)
palette_edge <- c(
  "FDR supported" = "#2C5872",
  "Nominal only" = "#99591D",
  "No nominal support" = "#A9A69E"
)
status_shapes <- c(
  "FDR supported" = 21,
  "Nominal only" = 23,
  "No nominal support" = 21
)
status_labels <- c(
  "FDR supported" = "BH–FDR supported",
  "Nominal only" = "Nominal P < 0.05 only",
  "No nominal support" = "No nominal support"
)

# Numerical underflow values are retained as zero in source data and shown at
# the declared plotting ceiling. No value is fitted, smoothed, or imputed.
score_cap <- 8
score_floor <- 10^(-score_cap)

ranked <- raw %>%
  arrange(p_BH_FDR, p_raw, test_id) %>%
  mutate(
    display_rank = row_number(),
    status = classify_status(p_raw, p_BH_FDR),
    minus_log10_q = pmin(-log10(pmax(p_BH_FDR, score_floor)), score_cap),
    q_display_capped = p_BH_FDR < score_floor,
    association = paste0(phenotype, " × ", climate)
  )

if (sum(ranked$p_raw < alpha) != 17L) stop("Expected 17 nominally supported tests")
if (sum(ranked$p_BH_FDR <= alpha) != 16L) stop("Expected 16 BH-FDR-supported tests")
if (sum(ranked$status == "Nominal only") != 1L) stop("Expected one nominal-only test")
if (sum(ranked$status == "No nominal support") != 7L) stop("Expected seven unsupported tests")
if (any(ranked$p_BH_FDR + 1e-15 < ranked$p_raw)) {
  stop("An adjusted P value is smaller than its raw P value")
}

boundary <- ranked %>% filter(status == "Nominal only")
if (nrow(boundary) != 1L || boundary$test_id != "CR-Ab__WetDays") {
  stop("The expected sole nominal-only association was not recovered")
}

theme_pub <- theme_classic(base_size = 8.0, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.38, colour = "#303030"),
    axis.ticks = element_line(linewidth = 0.34, colour = "#303030"),
    axis.ticks.length = unit(1.2, "mm"),
    axis.text = element_text(size = 7.2, colour = "#3D3D3D"),
    axis.title = element_text(size = 7.8, colour = "#222222"),
    plot.title = element_text(
      size = 8.4, face = "bold", colour = "#1F1F1F",
      margin = margin(b = 4.5)
    ),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7.1, colour = "#333333"),
    legend.key.width = unit(4.0, "mm"),
    legend.key.height = unit(4.0, "mm"),
    legend.spacing.x = unit(2.5, "mm"),
    legend.margin = margin(t = 0, unit = "pt"),
    plot.margin = margin(4, 6, 1, 3, unit = "pt")
  )

p_rank <- ggplot(ranked, aes(x = display_rank, y = minus_log10_q)) +
  geom_hline(
    yintercept = -log10(alpha), colour = "#767676",
    linewidth = 0.42, linetype = "22"
  ) +
  geom_segment(
    aes(xend = display_rank, y = 0, yend = minus_log10_q, colour = status),
    linewidth = 0.34, alpha = 0.36, show.legend = FALSE
  ) +
  geom_point(
    aes(fill = status, colour = status, shape = status),
    size = 3.0, stroke = 0.72
  ) +
  annotate(
    "text", x = 24.15, y = -log10(alpha) + 0.12,
    label = "q = 0.05", hjust = 1, vjust = 0,
    size = 2.45, family = "Arial", colour = "#616161"
  ) +
  annotate(
    "text", x = boundary$display_rank + 0.85, y = 1.84,
    label = "CR-Ab–wet days\nq = 0.0536", hjust = 0,
    size = 2.45, lineheight = 0.94, family = "Arial",
    fontface = "bold", colour = palette_edge[["Nominal only"]]
  ) +
  scale_fill_manual(values = palette_fill, labels = status_labels, drop = FALSE) +
  scale_colour_manual(values = palette_edge, labels = status_labels, drop = FALSE) +
  scale_shape_manual(values = status_shapes, labels = status_labels, drop = FALSE) +
  scale_x_continuous(
    limits = c(0.7, 24.35), breaks = c(1, 4, 8, 12, 16, 20, 24),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 8.18), breaks = c(0, 2, 4, 6, 8),
    expand = expansion(mult = c(0, 0))
  ) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(
        shape = unname(status_shapes), fill = unname(palette_fill),
        colour = unname(palette_edge), size = 3.1, stroke = 0.72
      )
    ),
    colour = "none", shape = "none"
  ) +
  labs(
    title = "16 of 24 associations retained FDR support",
    x = "Association rank by BH-adjusted P",
    y = expression(-log[10](q))
  ) +
  coord_cartesian(clip = "off") +
  theme_pub

matrix_phenotype_levels <- rev(phenotype_levels)
matrix_data <- raw %>%
  mutate(
    status = classify_status(p_raw, p_BH_FDR),
    phenotype = factor(phenotype, levels = matrix_phenotype_levels),
    climate = factor(climate, levels = climate_levels)
  )

support_counts <- matrix_data %>%
  group_by(climate) %>%
  summarise(n_supported = sum(p_BH_FDR <= alpha), .groups = "drop") %>%
  mutate(climate = factor(climate, levels = climate_levels)) %>%
  arrange(climate)

expected_counts <- c(6L, 5L, 4L, 1L)
if (!identical(support_counts$n_supported, expected_counts)) {
  stop("Unexpected support counts by climate: ", paste(support_counts$n_supported, collapse = ", "))
}

climate_axis_labels <- c(
  Temperature = "Temperature\n6/6",
  Precipitation = "Precipitation\n5/6",
  WetDays = "Wet days\n4/6",
  Humidity = "Relative humidity\n1/6"
)

p_matrix <- ggplot(matrix_data, aes(x = climate, y = phenotype)) +
  geom_hline(
    yintercept = seq(1, 6), colour = "#ECEBE8",
    linewidth = 0.34, show.legend = FALSE
  ) +
  geom_vline(
    xintercept = seq(1, 4), colour = "#F1F0ED",
    linewidth = 0.30, show.legend = FALSE
  ) +
  geom_point(
    aes(fill = status, colour = status, shape = status),
    size = 4.7, stroke = 0.78, show.legend = FALSE
  ) +
  scale_fill_manual(values = palette_fill, drop = FALSE) +
  scale_colour_manual(values = palette_edge, drop = FALSE) +
  scale_shape_manual(values = status_shapes, drop = FALSE) +
  scale_x_discrete(labels = climate_axis_labels, position = "top", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "Support was strongest for temperature",
    x = NULL, y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 8.0, base_family = "Arial") +
  theme(
    axis.text.x.top = element_text(
      size = 6.55, face = "bold", colour = "#333333",
      lineheight = 1.05, margin = margin(b = 7.0)
    ),
    axis.text.y = element_text(
      size = 7.2, colour = "#4A4A4A", margin = margin(r = 5.0)
    ),
    plot.title = element_text(
      size = 8.4, face = "bold", colour = "#1F1F1F",
      margin = margin(b = 9.0)
    ),
    plot.margin = margin(4, 3, 1, 4, unit = "pt")
  )

figure <- (p_rank | p_matrix) +
  plot_layout(widths = c(1.38, 1.10), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    plot.tag = element_text(size = 8.8, face = "bold", colour = "#111111"),
    plot.tag.position = c(0, 1)
  )

stem <- "bh_fdr_primary_associations"
width_mm <- 183
height_mm <- 89
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svg_path <- file.path(dirs$figure, paste0(stem, ".svg"))
pdf_path <- file.path(dirs$figure, paste0(stem, ".pdf"))
tiff_path <- file.path(dirs$figure, paste0(stem, ".tiff"))
png_path <- file.path(dirs$figure, paste0(stem, ".png"))

svglite::svglite(svg_path, width = width_in, height = height_in)
print(figure)
dev.off()

grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Arial")
print(figure)
dev.off()

ragg::agg_tiff(
  tiff_path, width = width_in, height = height_in, units = "in", res = 600,
  compression = "lzw", background = "white"
)
print(figure)
dev.off()

ragg::agg_png(
  png_path, width = width_in, height = height_in, units = "in", res = 300,
  background = "white"
)
print(figure)
dev.off()

source_data <- ranked %>%
  select(
    display_rank, test_id, phenotype, climate, lag_years, EDF, F_statistic,
    p_raw, p_BH_FDR, status, minus_log10_q, q_display_capped
  )
write.csv(
  source_data,
  file.path(dirs$source, "bh_fdr_primary_associations.csv"),
  row.names = FALSE, na = ""
)
message("BH-FDR outputs written to: ", normalizePath(output_root))
