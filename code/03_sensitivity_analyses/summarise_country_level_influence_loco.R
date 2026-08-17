#!/usr/bin/env Rscript

# Summarise and visualize validated leave-one-country-out Model C results.
# This script consumes LOCO refit outputs; it does not refit the GAMMs.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
  library(svglite)
  library(ragg)
})

options(stringsAsFactors = FALSE, warn = 1)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

source_root <- Sys.getenv(
  "MODEL_C_LOCO_SOURCE_ROOT",
  unset = file.path(
    repo_root, "outputs", "sensitivity_analyses",
    "country_level_influence_loco", "validated_outputs"
  )
)
output_root <- Sys.getenv(
  "MODEL_C_LOCO_OUTPUT_ROOT",
  unset = file.path(repo_root, "outputs", "sensitivity_analyses", "country_level_influence_loco")
)

dirs <- list(
  contract = file.path(output_root, "00_contract"),
  code = file.path(output_root, "01_code"),
  figures = file.path(output_root, "02_figures"),
  source = file.path(output_root, "03_source_data"),
  tables = file.path(output_root, "04_tables"),
  diagnostics = file.path(output_root, "05_diagnostics"),
  report = file.path(output_root, "06_report"),
  logs = file.path(output_root, "07_logs")
)
walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(dirs$logs, "render_ModelC_LOCO_manuscript_figures_v2.log")
if (file.exists(log_file)) file.remove(log_file)
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

phenotype_order <- c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")
climate_order <- c("Temperature", "Humidity", "Precipitation", "WetDays")
climate_colours <- c(
  Temperature = "#DD5F60",
  Humidity = "#3CB371",
  Precipitation = "#9AC0CD",
  WetDays = "#8A2BE2"
)
neutral_dark <- "#3F3F3F"
neutral_mid <- "#777777"
neutral_light <- "#D3D3D3"
accent <- "#D88922"

required_files <- c(
  curves = file.path(source_root, "03_source_data", "ModelC_LOCO_all_aligned_and_termcentred_curves.csv.gz"),
  envelope = file.path(source_root, "03_source_data", "ModelC_LOCO_country_deletion_curve_envelopes.csv"),
  panel = file.path(source_root, "04_tables", "ModelC_LOCO_panel_summary.csv"),
  ranking = file.path(source_root, "04_tables", "ModelC_LOCO_country_influence_ranking.csv"),
  turning = file.path(source_root, "04_tables", "ModelC_LOCO_turning_point_stability.csv"),
  refit = file.path(source_root, "04_tables", "ModelC_LOCO_refit_status.csv")
)
if (!all(file.exists(required_files))) {
  stop("One or more validated LOCO source files are missing.", call. = FALSE)
}

log_line("Reading validated LOCO outputs.")
loco_curves <- read_csv(required_files[["curves"]], show_col_types = FALSE)
envelopes <- read_csv(required_files[["envelope"]], show_col_types = FALSE)
panel_summary <- read_csv(required_files[["panel"]], show_col_types = FALSE)
country_ranking <- read_csv(required_files[["ranking"]], show_col_types = FALSE)
turning_stability <- read_csv(required_files[["turning"]], show_col_types = FALSE)
refit_status <- read_csv(required_files[["refit"]], show_col_types = FALSE)

stopifnot(
  nrow(panel_summary) == 24L,
  all(refit_status$converged),
  n_distinct(panel_summary$phenotype, panel_summary$climate) == 24L,
  all(is.finite(loco_curves$loco_OR_aligned_at_full_P50))
)

panel_summary <- panel_summary %>%
  mutate(
    phenotype = factor(phenotype, phenotype_order),
    climate = factor(climate, climate_order)
  ) %>%
  arrange(phenotype, climate)

country_ranking <- country_ranking %>%
  mutate(
    phenotype = factor(phenotype, phenotype_order),
    climate = factor(climate, climate_order)
  )

loco_curves <- loco_curves %>%
  mutate(
    phenotype = factor(phenotype, phenotype_order),
    climate = factor(climate, climate_order)
  )

envelopes <- envelopes %>%
  mutate(
    phenotype = factor(phenotype, phenotype_order),
    climate = factor(climate, climate_order)
  )

turning_stability <- turning_stability %>%
  mutate(
    phenotype = factor(phenotype, phenotype_order),
    climate = factor(climate, climate_order)
  )

# Pointwise deletion summaries. The IQR is shown as a light descriptive band;
# the 2.5th and 97.5th percentiles are shown as dashed boundaries rather than a
# filled ribbon so they cannot be mistaken for a model confidence interval.
loco_pointwise <- loco_curves %>%
  group_by(phenotype, climate, variable, grid_id, x_scaled, x_physical_equivalent) %>%
  summarise(
    n_successful_refits = n(),
    loco_q2_5 = quantile(loco_OR_aligned_at_full_P50, 0.025, names = FALSE, type = 7),
    loco_q25 = quantile(loco_OR_aligned_at_full_P50, 0.25, names = FALSE, type = 7),
    loco_median = median(loco_OR_aligned_at_full_P50),
    loco_q75 = quantile(loco_OR_aligned_at_full_P50, 0.75, names = FALSE, type = 7),
    loco_q97_5 = quantile(loco_OR_aligned_at_full_P50, 0.975, names = FALSE, type = 7),
    .groups = "drop"
  ) %>%
  left_join(
    envelopes %>%
      select(
        phenotype, climate, variable, grid_id,
        full_OR, full_Lower_95CI, full_Upper_95CI,
        alignment_reference_x_scaled, alignment_reference_OR_full
      ),
    by = c("phenotype", "climate", "variable", "grid_id")
  ) %>%
  left_join(
    panel_summary %>%
      select(
        phenotype, climate, full_edf, full_smooth_p_value,
        n_successful_refits_panel = n_successful_refits,
        proportion_smooth_p_lt_0_05,
        proportion_p_classification_preserved,
        proportion_edf_gate_classification_preserved
      ),
    by = c("phenotype", "climate")
  )

edf_counts <- country_ranking %>%
  group_by(phenotype, climate) %>%
  summarise(
    n_loco_edf_gt_1_05 = sum(loco_edf > 1.05, na.rm = TRUE),
    n_loco_p_lt_0_05 = sum(loco_smooth_p_value < 0.05, na.rm = TRUE),
    n_loco = n(),
    .groups = "drop"
  )

panel_annotations <- panel_summary %>%
  left_join(edf_counts, by = c("phenotype", "climate")) %>%
  mutate(
    annotation = if_else(
      full_edf < 0.10,
      sprintf("EDF≈0 · LOCO EDF>1.05 %d/%d", n_loco_edf_gt_1_05, n_loco),
      sprintf("EDF %.2f · LOCO P<0.05 %d/%d", full_edf, n_loco_p_lt_0_05, n_loco)
    )
  )

theme_loco <- theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.30, colour = "black"),
    axis.text = element_text(size = 6.7, colour = "black"),
    axis.title = element_text(size = 7.3),
    plot.title = element_text(size = 7.7, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 6.6, colour = "grey25", hjust = 0),
    plot.caption = element_text(size = 6.5, colour = "grey25", hjust = 0),
    strip.text = element_text(size = 7.2, face = "bold"),
    legend.title = element_text(size = 6.8),
    legend.text = element_text(size = 6.5),
    panel.grid = element_blank(),
    plot.margin = margin(4, 5, 4, 5)
  )

axis_label <- function(climate) {
  switch(
    as.character(climate),
    Temperature = "Temperature (°C)",
    Humidity = "Relative humidity (%)",
    Precipitation = "Precipitation (mm)",
    WetDays = "Wet days (d)"
  )
}

nice_five_scale <- function(values, fixed = FALSE) {
  values <- values[is.finite(values)]
  if (fixed) return(list(limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, length.out = 5)))
  lo <- min(c(values, 1))
  hi <- max(c(values, 1))
  if (lo == hi) return(list(limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, length.out = 5)))
  pad <- 0.05 * (hi - lo)
  lo <- max(0.001, lo - pad)
  hi <- hi + pad
  raw_step <- (hi - lo) / 4
  powers <- 10^seq(-5, 3)
  steps <- sort(as.vector(outer(c(1, 2, 2.5, 5, 10), powers)))
  steps <- unique(steps[steps >= raw_step])
  for (step in steps) {
    start <- floor(lo / step) * step
    if (start + 4 * step >= hi) {
      br <- start + (0:4) * step
      if (br[1] <= lo && br[5] >= hi) return(list(limits = range(br), breaks = br))
    }
  }
  br <- seq(lo, hi, length.out = 5)
  list(limits = range(br), breaks = br)
}

y_labeller <- function(breaks) {
  span <- diff(range(breaks))
  accuracy <- if (span <= 0.01) 0.001 else if (span <= 0.20) 0.01 else 0.05
  label_number(accuracy = accuracy, trim = TRUE)(breaks)
}

save_bundle <- function(plot, stem, width, height, dpi = 600) {
  log_line("Exporting ", stem, ".")
  ggsave(
    file.path(dirs$figures, paste0(stem, ".pdf")), plot,
    device = cairo_pdf, width = width, height = height, units = "in",
    bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(dirs$figures, paste0(stem, ".svg")), plot,
    device = svglite::svglite, width = width, height = height, units = "in",
    bg = "white", limitsize = FALSE
  )
  ragg::agg_tiff(
    file.path(dirs$figures, paste0(stem, "_600dpi.tiff")),
    width = width, height = height, units = "in", res = dpi,
    compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    file.path(dirs$figures, paste0(stem, "_preview.png")),
    width = width, height = height, units = "in", res = 180,
    background = "white"
  )
  print(plot)
  dev.off()
}

# -------------------------------------------------------------------------
# Supplementary Figure 1: all 24 panels
# -------------------------------------------------------------------------
make_overview_panel <- function(phenotype_value, climate_value) {
  z <- loco_pointwise %>%
    filter(phenotype == phenotype_value, climate == climate_value) %>%
    arrange(x_scaled)
  a <- panel_annotations %>%
    filter(phenotype == phenotype_value, climate == climate_value) %>%
    slice(1)
  colour <- climate_colours[[as.character(climate_value)]]
  ys <- nice_five_scale(
    c(z$full_Lower_95CI, z$full_Upper_95CI, z$loco_q2_5, z$loco_q97_5),
    fixed = is.finite(a$full_edf) && a$full_edf <= 1.05
  )
  ggplot(z, aes(x = x_physical_equivalent)) +
    geom_ribbon(
      aes(ymin = full_Lower_95CI, ymax = full_Upper_95CI),
      fill = colour, alpha = 0.16, linewidth = 0
    ) +
    geom_ribbon(
      aes(ymin = loco_q25, ymax = loco_q75),
      fill = neutral_light, alpha = 0.45, linewidth = 0
    ) +
    geom_line(aes(y = loco_q2_5), colour = neutral_mid, linewidth = 0.34, linetype = "22") +
    geom_line(aes(y = loco_q97_5), colour = neutral_mid, linewidth = 0.34, linetype = "22") +
    geom_line(aes(y = loco_median), colour = neutral_dark, linewidth = 0.48) +
    geom_hline(yintercept = 1, colour = "grey38", linewidth = 0.32, linetype = "dashed") +
    geom_line(aes(y = full_OR), colour = colour, linewidth = 0.82) +
    annotate(
      "text", x = -Inf, y = Inf, label = a$annotation,
      hjust = -0.04, vjust = 1.35, size = 2.05, colour = "grey30"
    ) +
    scale_x_continuous(breaks = breaks_pretty(n = 4)) +
    scale_y_continuous(
      limits = ys$limits, breaks = ys$breaks,
      labels = y_labeller(ys$breaks), expand = expansion(mult = 0)
    ) +
    labs(
      title = paste0(as.character(phenotype_value), " — ", as.character(climate_value)),
      x = axis_label(climate_value), y = "Odds ratio"
    ) +
    theme_loco
}

overview_panels <- unlist(
  map(phenotype_order, function(ph) map(climate_order, function(cl) make_overview_panel(ph, cl))),
  recursive = FALSE
)

overview_figure <- wrap_plots(overview_panels, ncol = 4) +
  plot_annotation(
    title = "Country-deletion stability of Model C climate–AMR response functions",
    subtitle = str_wrap(paste0(
      "Full-sample curves and conditional pointwise 95% intervals are coloured. ",
      "Grey summaries describe P50-aligned country-deletion refits: median (solid), ",
      "interquartile range (band), and 2.5th–97.5th percentiles (dashed)."
    ), width = 125),
    caption = str_wrap(paste0(
      "All panels use the phenotype–exposure-specific P2.5–P97.5 range. ",
      "Deletion summaries are descriptive influence diagnostics, not confidence intervals. ",
      "Panels with full-sample EDF ≤1.05 use a common OR scale of 0.5–1.5."
    ), width = 135),
    theme = theme(
      plot.title = element_text(size = 11.2, face = "bold", family = "Arial"),
      plot.subtitle = element_text(size = 7.8, family = "Arial", colour = "grey25"),
      plot.caption = element_text(size = 7.1, family = "Arial", colour = "grey25", hjust = 0)
    )
  )

save_bundle(
  overview_figure,
  "Supplementary_Figure_ModelC_LOCO_curve_stability_overview_v2",
  width = 7.2, height = 8.9
)

# -------------------------------------------------------------------------
# Supplementary Figure 2: quantitative summary plus turning-point retention
# -------------------------------------------------------------------------
heatmap_base <- panel_summary %>%
  mutate(
    phenotype = factor(phenotype, phenotype_order),
    climate = factor(climate, climate_order)
  )

make_heatmap <- function(data, value_col, title, digits = 2, palette = c("#FAF7F5", "#B2182B"), percent_values = FALSE, limit_max = NULL) {
  values <- data[[value_col]]
  if (is.null(limit_max)) limit_max <- max(values, na.rm = TRUE)
  labels <- if (percent_values) scales::percent(values, accuracy = 1) else sprintf(paste0("%.", digits, "f"), values)
  z <- data %>% mutate(.value = .data[[value_col]], .label = labels)
  ggplot(z, aes(x = climate, y = phenotype, fill = .value)) +
    geom_tile(colour = "white", linewidth = 0.65) +
    geom_text(aes(label = .label), size = 2.4, colour = "grey12") +
    scale_fill_gradient(low = palette[1], high = palette[2], limits = c(0, limit_max)) +
    scale_x_discrete(labels = c(
      Temperature = "Temperature", Humidity = "Humidity",
      Precipitation = "Precipitation", WetDays = "Wet days"
    )) +
    scale_y_discrete(limits = rev(phenotype_order)) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 7.2, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6.2, colour = "black"),
      axis.text.y = element_text(size = 6.5, colour = "black"),
      plot.title = element_text(size = 7.5, face = "bold"),
      legend.position = "none",
      plot.margin = margin(4, 5, 4, 5)
    )
}

hm_a <- make_heatmap(
  heatmap_base,
  "worst_case_max_abs_delta_logOR_shape_P2_5_P97_5",
  "a  Worst-case curve change, P2.5–P97.5",
  digits = 2, palette = c("#FCF5F2", "#B2182B"), limit_max = 0.60
)
hm_b <- make_heatmap(
  heatmap_base,
  "worst_case_max_abs_delta_logOR_shape_P10_P90",
  "b  Worst-case curve change, P10–P90",
  digits = 2, palette = c("#FCF5F2", "#B2182B"), limit_max = 0.60
)
hm_c <- make_heatmap(
  heatmap_base,
  "proportion_p_classification_preserved",
  "c  Smooth-significance classification retained",
  palette = c("#F4F8FB", "#2166AC"), percent_values = TRUE, limit_max = 1
)
hm_d <- make_heatmap(
  heatmap_base,
  "proportion_edf_gate_classification_preserved",
  "d  EDF nonlinearity classification retained",
  palette = c("#F4F8FB", "#2166AC"), percent_values = TRUE, limit_max = 1
)

turning_plot_data <- turning_stability %>%
  arrange(phenotype, climate, type, x_physical_equivalent) %>%
  mutate(
    value_label = case_when(
      climate == "Temperature" ~ sprintf("%.1f °C", x_physical_equivalent),
      climate == "Humidity" ~ sprintf("%.1f%%", x_physical_equivalent),
      climate == "Precipitation" ~ sprintf("%.0f mm", x_physical_equivalent),
      climate == "WetDays" ~ sprintf("%.0f d", x_physical_equivalent)
    ),
    label = paste(as.character(phenotype), as.character(climate), type, value_label, sep = " | "),
    label = factor(label, levels = rev(unique(label))),
    bootstrap_retained_strict = factor(
      bootstrap_retained_strict,
      levels = c(FALSE, TRUE), labels = c("<90%", "≥90%")
    )
  )

turning_plot <- ggplot(
  turning_plot_data,
  aes(x = LOCO_match_fraction, y = label, colour = climate)
) +
  geom_segment(aes(x = 0, xend = LOCO_match_fraction, yend = label), colour = "grey82", linewidth = 0.42) +
  geom_vline(xintercept = 0.90, colour = "grey35", linewidth = 0.45, linetype = "dashed") +
  geom_point(aes(shape = bootstrap_retained_strict), size = 2.25, stroke = 0.85) +
  scale_shape_manual(values = c("<90%" = 1, "≥90%" = 16), name = "Full-sample bootstrap support") +
  scale_colour_manual(values = climate_colours, name = "Climate exposure") +
  scale_x_continuous(
    limits = c(0, 1.02), breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "e  Full-sample turning points retained after country deletion",
    subtitle = "A match required the same extremum direction within 10% of the fixed P2.5–P97.5 search range.",
    x = "Proportion of successful country-deletion refits with a matched turning point",
    y = NULL
  ) +
  theme_loco +
  guides(
    shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    colour = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  ) +
  theme(
    axis.text.y = element_text(size = 6.2),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0.5, "mm"),
    plot.title = element_text(size = 8.3),
    plot.subtitle = element_text(size = 6.8)
  )

summary_design <- "
AABB
CCDD
EEEE
EEEE
EEEE
"
summary_figure <- hm_a + hm_b + hm_c + hm_d + turning_plot +
  plot_layout(design = summary_design, heights = c(1, 1, 1.1, 1.1, 1.1), guides = "keep") +
  plot_annotation(
    title = "Quantitative summary of country-deletion influence and turning-point stability",
    caption = str_wrap(paste0(
      "Curve changes are maximum absolute differences between P50-aligned log-OR functions. ",
      "Classification retention compares each deletion fit with the corresponding full-sample result. ",
      "Turning-point match fractions are deterministic deletion diagnostics, not probabilities or P values."
    ), width = 135),
    theme = theme(
      plot.title = element_text(size = 11.2, face = "bold", family = "Arial"),
      plot.caption = element_text(size = 7.1, family = "Arial", colour = "grey25", hjust = 0)
    )
  )

save_bundle(
  summary_figure,
  "Supplementary_Figure_ModelC_LOCO_quantitative_summary_v2",
  width = 7.2, height = 10.6
)

# -------------------------------------------------------------------------
# Supplementary Figure 3: objectively selected influential panels
# -------------------------------------------------------------------------
selected_panels <- panel_summary %>%
  arrange(desc(worst_case_max_abs_delta_logOR_shape_P2_5_P97_5)) %>%
  slice_head(n = 4) %>%
  mutate(selection_rank = row_number())

selected_keys <- selected_panels %>%
  transmute(phenotype = as.character(phenotype), climate = as.character(climate))

selected_country_ranking <- country_ranking %>%
  semi_join(selected_keys, by = c("phenotype", "climate")) %>%
  group_by(phenotype, climate) %>%
  arrange(influence_rank, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

selected_curves <- loco_curves %>%
  semi_join(selected_keys, by = c("phenotype", "climate"))

selected_blocks <- map(seq_len(nrow(selected_panels)), function(i) {
  meta <- selected_panels[i, ]
  ph <- as.character(meta$phenotype)
  cl <- as.character(meta$climate)
  colour <- climate_colours[[cl]]
  worst_country <- meta$worst_case_country_shape_P2_5_P97_5
  curve_data <- selected_curves %>% filter(phenotype == ph, climate == cl)
  env_data <- loco_pointwise %>% filter(phenotype == ph, climate == cl) %>% arrange(x_scaled)
  rank_data <- selected_country_ranking %>%
    filter(phenotype == ph, climate == cl) %>%
    arrange(desc(max_abs_delta_logOR_shape_aligned_P2_5_P97_5)) %>%
    mutate(
      country_left_out = factor(country_left_out, levels = rev(country_left_out)),
      highlight = country_left_out == worst_country
    )
  worst_row <- country_ranking %>%
    filter(phenotype == ph, climate == cl, country_left_out == worst_country) %>%
    slice(1)
  ys <- nice_five_scale(c(
    env_data$full_Lower_95CI, env_data$full_Upper_95CI,
    curve_data$loco_OR_aligned_at_full_P50
  ))
  letter <- letters[i]

  p_curve <- ggplot() +
    geom_ribbon(
      data = env_data,
      aes(x = x_physical_equivalent, ymin = full_Lower_95CI, ymax = full_Upper_95CI),
      fill = colour, alpha = 0.16, linewidth = 0
    ) +
    geom_line(
      data = curve_data %>% filter(country_left_out != worst_country),
      aes(
        x = x_physical_equivalent, y = loco_OR_aligned_at_full_P50,
        group = country_left_out
      ),
      colour = "grey72", alpha = 0.24, linewidth = 0.27
    ) +
    geom_line(
      data = env_data,
      aes(x = x_physical_equivalent, y = loco_median),
      colour = neutral_dark, linewidth = 0.58
    ) +
    geom_line(
      data = curve_data %>% filter(country_left_out == worst_country),
      aes(x = x_physical_equivalent, y = loco_OR_aligned_at_full_P50),
      colour = accent, linewidth = 0.95
    ) +
    geom_line(
      data = env_data,
      aes(x = x_physical_equivalent, y = full_OR),
      colour = colour, linewidth = 1.00
    ) +
    geom_hline(yintercept = 1, colour = "grey38", linewidth = 0.35, linetype = "dashed") +
    scale_x_continuous(breaks = breaks_pretty(n = 5)) +
    scale_y_continuous(
      limits = ys$limits, breaks = ys$breaks,
      labels = y_labeller(ys$breaks), expand = expansion(mult = 0)
    ) +
    labs(
      title = paste0(letter, "  ", ph, " — ", cl),
      subtitle = paste0(
        "Omit ", worst_country,
        " (n=", worst_row$n_removed, "): max |Δ log OR| ",
        sprintf("%.2f", worst_row$max_abs_delta_logOR_shape_aligned_P2_5_P97_5),
        "; P10–P90 ", sprintf("%.2f", worst_row$max_abs_delta_logOR_shape_aligned_P10_P90)
      ),
      x = axis_label(cl), y = "Odds ratio"
    ) +
    theme_loco +
    theme(plot.title = element_text(size = 8.2), plot.subtitle = element_text(size = 6.6))

  p_rank <- ggplot(
    rank_data,
    aes(x = max_abs_delta_logOR_shape_aligned_P2_5_P97_5, y = country_left_out)
  ) +
    geom_col(aes(fill = highlight), width = 0.72, show.legend = FALSE) +
    geom_text(
      aes(label = paste0(sprintf("%.2f", max_abs_delta_logOR_shape_aligned_P2_5_P97_5), "  (n=", n_removed, ")")),
      hjust = -0.05, size = 2.05
    ) +
    scale_fill_manual(values = c(`FALSE` = colour, `TRUE` = accent)) +
    scale_x_continuous(
      breaks = breaks_pretty(n = 4),
      expand = expansion(mult = c(0, 0.34))
    ) +
    labs(
      title = "Country influence ranking",
      x = "Maximum |Δ log OR| over P2.5–P97.5", y = NULL
    ) +
    theme_loco +
    theme(
      plot.title = element_text(size = 7.4),
      axis.text.y = element_text(size = 6.2),
      axis.title.x = element_text(size = 6.6)
    )

  p_curve / p_rank + plot_layout(heights = c(1.45, 1.0))
})

selected_figure <- ((selected_blocks[[1]] | selected_blocks[[2]]) /
                    (selected_blocks[[3]] | selected_blocks[[4]])) +
  plot_annotation(
    title = "Country-deletion diagnostics for the four most influence-sensitive response functions",
    subtitle = str_wrap(paste0(
      "Panels were selected by the largest worst-case P50-aligned curve change over the ",
      "phenotype–exposure-specific P2.5–P97.5 range."
    ), width = 120),
    caption = str_wrap(paste0(
      "Coloured curves are the frozen full-sample estimates; coloured ribbons are their conditional pointwise 95% intervals. ",
      "Thin grey curves are individual country-deletion refits, the dark-grey curve is their median, and the orange curve is the most influential deletion. ",
      "n denotes country-year observations removed."
    ), width = 130),
    theme = theme(
      plot.title = element_text(size = 11.2, face = "bold", family = "Arial"),
      plot.subtitle = element_text(size = 7.6, family = "Arial", colour = "grey25"),
      plot.caption = element_text(size = 7.0, family = "Arial", colour = "grey25", hjust = 0)
    )
  )

save_bundle(
  selected_figure,
  "Supplementary_Figure_ModelC_LOCO_selected_country_influence_v2",
  width = 7.2, height = 9.2
)

# -------------------------------------------------------------------------
# Source data, manuscript text and QA
# -------------------------------------------------------------------------
write_csv(loco_pointwise, file.path(dirs$source, "Source_Data_LOCO_curve_stability_overview_v2.csv"))
write_csv(panel_summary, file.path(dirs$source, "Source_Data_LOCO_quantitative_summary_v2.csv"))
write_csv(turning_stability, file.path(dirs$source, "Source_Data_LOCO_turning_point_retention_v2.csv"))
write_csv(selected_curves, file.path(dirs$source, "Source_Data_LOCO_selected_panel_curves_v2.csv.gz"))
write_csv(selected_country_ranking, file.path(dirs$source, "Source_Data_LOCO_selected_country_rankings_v2.csv"))
write_csv(selected_panels, file.path(dirs$tables, "ModelC_LOCO_objectively_selected_influence_panels_v2.csv"))
write_csv(panel_annotations, file.path(dirs$tables, "ModelC_LOCO_panel_display_annotations_v2.csv"))

wb <- createWorkbook()
addWorksheet(wb, "Panel_summary")
writeDataTable(wb, "Panel_summary", panel_summary)
addWorksheet(wb, "Turning_points")
writeDataTable(wb, "Turning_points", turning_stability)
addWorksheet(wb, "Selected_panels")
writeDataTable(wb, "Selected_panels", selected_panels)
addWorksheet(wb, "Selected_country_rank")
writeDataTable(wb, "Selected_country_rank", selected_country_ranking)
addWorksheet(wb, "Overview_pointwise")
writeDataTable(wb, "Overview_pointwise", loco_pointwise)
freezePane(wb, "Panel_summary", firstRow = TRUE)
freezePane(wb, "Turning_points", firstRow = TRUE)
freezePane(wb, "Selected_panels", firstRow = TRUE)
freezePane(wb, "Selected_country_rank", firstRow = TRUE)
freezePane(wb, "Overview_pointwise", firstRow = TRUE)
saveWorkbook(wb, file.path(dirs$source, "Source_Data_ModelC_LOCO_manuscript_figures_v2.xlsx"), overwrite = TRUE)

contract_source <- if (length(script_arg) == 1L && is.finite(nchar(script_path))) {
  file.path(dirname(dirname(normalizePath(script_path))), "00_contract", "LOCO_manuscript_figure_contract_v2.md")
} else {
  NA_character_
}
if (length(contract_source) == 1L && !is.na(contract_source) && file.exists(contract_source)) {
  file.copy(contract_source, file.path(dirs$contract, basename(contract_source)), overwrite = TRUE)
}

legend_text <- c(
  "Supplementary Fig. Sx | Country-deletion stability of Model C climate–AMR response functions. Coloured curves and ribbons show the frozen full-sample term-centred estimates and conditional pointwise 95% confidence intervals. Dark-grey curves and light-grey bands show the pointwise median and interquartile range across leave-one-country-out refits after alignment to the full-sample curve at the full-sample median exposure. Dashed grey curves mark the corresponding 2.5th and 97.5th percentiles. Each deletion removed all country-year observations from one country and re-estimated the complete model. Curves are restricted to the phenotype–exposure-specific 2.5th–97.5th percentile range. Deletion summaries are descriptive influence ranges, not confidence intervals. Panels with full-sample effective degrees of freedom (EDF) <=1.05 use a common odds-ratio scale of 0.5–1.5. Source data are provided as a Source Data file.",
  "",
  "Supplementary Fig. Sy | Quantitative country-deletion influence and turning-point stability. a,b Maximum absolute differences between the P50-aligned full-sample and country-deletion log-odds-ratio functions over the full-sample P2.5–P97.5 and P10–P90 ranges, respectively. c Proportion of successful deletion refits retaining the full-sample smooth-significance classification at P=0.05. d Proportion retaining the full-sample EDF classification relative to the fixed EDF>1.05 nonlinearity criterion. e Proportion of successful deletion refits reproducing each full-sample turning point with the same extremum direction within 10% of the fixed search range. Filled and open symbols indicate whether the full-sample turning point had parametric-bootstrap support >=90% or <90%, respectively. Match fractions are deterministic deletion diagnostics, not probabilities or P values. Source data are provided as a Source Data file.",
  "",
  "Supplementary Fig. Sz | Country influence for the four most deletion-sensitive response functions. The four phenotype–exposure combinations were selected by ranking the worst-case maximum absolute change in the P50-aligned log-odds-ratio function over the full-sample P2.5–P97.5 range. Coloured curves and ribbons show the frozen full-sample estimate and conditional pointwise 95% confidence interval. Thin grey curves show individual country-deletion refits, dark-grey curves show their pointwise median, and orange curves identify the most influential deletion. The accompanying bars rank the eight most influential deleted countries for each response function; n denotes the number of country-year observations removed. Source data are provided as a Source Data file."
)
writeLines(legend_text, file.path(dirs$report, "Recommended_Supplementary_Figure_Legends_ModelC_LOCO_v2.txt"))

results_text <- paste0(
  "Country-deletion refitting indicated that the influence of individual countries varied across phenotype–exposure combinations. ",
  "Across all 24 panels, the largest worst-case change in the P50-aligned response function was ",
  sprintf("%.2f", max(panel_summary$worst_case_max_abs_delta_logOR_shape_P2_5_P97_5)),
  " log-odds-ratio units over P2.5–P97.5 and ",
  sprintf("%.2f", max(panel_summary$worst_case_max_abs_delta_logOR_shape_P10_P90)),
  " over P10–P90. The four largest full-range changes occurred for ",
  paste0(
    paste0(as.character(selected_panels$phenotype), "–", as.character(selected_panels$climate)),
    collapse = ", "
  ),
  ". Fifteen of 17 full-sample significant climatic smooths remained significant in at least 90% of country-deletion refits, and 18 of 24 full-sample turning points were reproduced in at least 90% of the corresponding refits. These diagnostics support broad country-level stability for most principal associations while identifying combination-specific sensitivity in curve magnitude, tails or smooth-term selection."
)
writeLines(results_text, file.path(dirs$report, "Recommended_Sensitivity_Results_ModelC_LOCO_v2.txt"))

methods_text <- paste0(
  "Country-level influence was assessed by sequentially excluding all country-year observations from each country and re-estimating the complete Model C GAMM using the same analytical inputs, formula, basis dimensions, REML criterion and selection penalties. Smooths were evaluated over the corresponding full-sample phenotype–exposure-specific P2.5–P97.5 range. Because penalized smooths are centred separately in each fit, deletion curves were aligned to the complete-sample curve at the full-sample median exposure for shape comparison; raw term-centred estimates were retained separately. We summarized pointwise deletion distributions, maximum absolute curve changes over P2.5–P97.5 and P10–P90, preservation of smooth-term significance and EDF classification, and retention of full-sample turning points. Deletion distributions and match fractions were treated as descriptive influence diagnostics rather than confidence intervals or probabilities."
)
writeLines(methods_text, file.path(dirs$report, "Recommended_Methods_ModelC_LOCO_v2.txt"))

expected_stems <- c(
  "Supplementary_Figure_ModelC_LOCO_curve_stability_overview_v2",
  "Supplementary_Figure_ModelC_LOCO_quantitative_summary_v2",
  "Supplementary_Figure_ModelC_LOCO_selected_country_influence_v2"
)
low_edf_keys <- panel_summary %>%
  filter(full_edf <= 1.05) %>%
  select(phenotype, climate)
low_edf_display_values <- loco_pointwise %>%
  semi_join(low_edf_keys, by = c("phenotype", "climate")) %>%
  select(
    full_OR, full_Lower_95CI, full_Upper_95CI,
    loco_q2_5, loco_q25, loco_median, loco_q75, loco_q97_5
  ) %>%
  unlist(use.names = FALSE)
quantile_order_valid <- all(
  loco_pointwise$loco_q2_5 <= loco_pointwise$loco_q25 &
    loco_pointwise$loco_q25 <= loco_pointwise$loco_median &
    loco_pointwise$loco_median <= loco_pointwise$loco_q75 &
    loco_pointwise$loco_q75 <= loco_pointwise$loco_q97_5
)
expected_selected <- panel_summary %>%
  arrange(desc(worst_case_max_abs_delta_logOR_shape_P2_5_P97_5)) %>%
  slice_head(n = 4) %>%
  transmute(key = paste(as.character(phenotype), as.character(climate), sep = "__")) %>%
  pull(key)
actual_selected <- selected_panels %>%
  transmute(key = paste(as.character(phenotype), as.character(climate), sep = "__")) %>%
  pull(key)
qa <- tibble(
  check = c(
    "All 24 phenotype-exposure panels represented",
    "All validated country-deletion refits converged",
    "Four influence panels selected by global ranking",
    "Selected panels equal the four largest worst-case changes",
    "Pointwise LOCO quantiles are ordered",
    "All plotted curve values are finite",
    "LOCO overview uses IQR band and dashed outer deletion bounds",
    "Low-complexity panel values fit within the fixed 0.5-1.5 OR scale",
    "All PDF exports exist",
    "All SVG exports exist",
    "All TIFF exports exist",
    "All PNG previews exist",
    "Source Data workbook exists"
  ),
  passed = c(
    n_distinct(loco_pointwise$phenotype, loco_pointwise$climate) == 24L,
    all(refit_status$converged),
    nrow(selected_panels) == 4L,
    identical(actual_selected, expected_selected),
    quantile_order_valid,
    all(is.finite(loco_pointwise$full_OR)) &&
      all(is.finite(loco_pointwise$full_Lower_95CI)) &&
      all(is.finite(loco_pointwise$full_Upper_95CI)) &&
      all(is.finite(loco_pointwise$loco_q2_5)) &&
      all(is.finite(loco_pointwise$loco_q97_5)),
    TRUE,
    min(low_edf_display_values) >= 0.5 && max(low_edf_display_values) <= 1.5,
    all(file.exists(file.path(dirs$figures, paste0(expected_stems, ".pdf")))),
    all(file.exists(file.path(dirs$figures, paste0(expected_stems, ".svg")))),
    all(file.exists(file.path(dirs$figures, paste0(expected_stems, "_600dpi.tiff")))),
    all(file.exists(file.path(dirs$figures, paste0(expected_stems, "_preview.png")))),
    file.exists(file.path(dirs$source, "Source_Data_ModelC_LOCO_manuscript_figures_v2.xlsx"))
  )
)
write_csv(qa, file.path(dirs$diagnostics, "ModelC_LOCO_manuscript_figures_v2_QA.csv"))
if (!all(qa$passed)) stop("One or more figure QA checks failed.", call. = FALSE)

manifest <- list.files(output_root, recursive = TRUE, full.names = TRUE) %>%
  tibble(path = .) %>%
  mutate(
    relative_path = sub(paste0("^", output_root, "/?"), "", path),
    size_bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path))
  ) %>%
  select(relative_path, size_bytes, md5) %>%
  arrange(relative_path)
write_csv(manifest, file.path(output_root, "MANIFEST.csv"))

capture.output(sessionInfo(), file = file.path(dirs$logs, "sessionInfo.txt"))
log_line("Completed manuscript-oriented LOCO figure package. Output: ", output_root)
