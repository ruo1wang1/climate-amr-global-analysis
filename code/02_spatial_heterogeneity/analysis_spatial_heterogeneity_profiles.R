###############################################################################
# Spatial-heterogeneity profile analysis with optimal-k clustering
###############################################################################

suppressPackageStartupMessages({
  required_packages <- c("tidyverse", "cowplot", "viridis", "openxlsx", "scales")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)

input_csv <- file.path(
  revision_root,
  "outputs/historical_associations",
  "figure3_spatial_heterogeneity_optimal_k",
  "01_country_effect_tables",
  "ModelC_all_bacteria_country_level_climate_effects_with_or_ci.csv"
)

output_root <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "spatial_heterogeneity_profiles"
)

output_dirs <- list(
  figures = file.path(output_root, "01_figures"),
  tables = file.path(output_root, "02_tables"),
  workbook = file.path(output_root, "03_workbook")
)

invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

if (!file.exists(input_csv)) {
  stop("Country-level climate effect table not found: ", input_csv, call. = FALSE)
}

bacteria_order <- c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")
factor_order <- c("Temperature", "Relative Humidity", "Precipitation", "Wet Days")

factor_colors <- c(
  "Temperature" = "#DD5F60",
  "Relative Humidity" = "#3CB371",
  "Precipitation" = "#9AC0CD",
  "Wet Days" = "#8A2BE2"
)

country_effects <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
  mutate(
    Bacteria = factor(Bacteria, levels = bacteria_order),
    Climate_Factor = factor(Climate_Factor, levels = factor_order)
  )

heterogeneity_summary <- country_effects %>%
  group_by(Bacteria, Climate_Factor) %>%
  summarise(
    n_countries = dplyr::n(),
    Mean_OR = mean(OR, na.rm = TRUE),
    SD_OR = sd(OR, na.rm = TRUE),
    CV_OR = ifelse(abs(Mean_OR) > .Machine$double.eps, SD_OR / Mean_OR * 100, NA_real_),
    Median_OR = median(OR, na.rm = TRUE),
    Min_OR = min(OR, na.rm = TRUE),
    Max_OR = max(OR, na.rm = TRUE),
    Mean_P_Value = mean(P_Value, na.rm = TRUE),
    Significant_Countries = sum(P_Value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_OR = round(Mean_OR, 2),
    SD_OR = round(SD_OR, 2),
    CV_OR = round(CV_OR, 2),
    Median_OR = round(Median_OR, 2),
    Min_OR = round(Min_OR, 2),
    Max_OR = round(Max_OR, 2),
    Mean_P_Value = round(Mean_P_Value, 4)
  )

bacteria_summary <- heterogeneity_summary %>%
  group_by(Bacteria) %>%
  summarise(
    Mean_CV = round(mean(CV_OR, na.rm = TRUE), 2),
    Mean_OR = round(mean(Mean_OR, na.rm = TRUE), 2),
    Max_Heterogeneity_Factor = as.character(Climate_Factor[which.max(CV_OR)]),
    Max_CV = round(max(CV_OR, na.rm = TRUE), 2),
    .groups = "drop"
  )

factor_summary <- heterogeneity_summary %>%
  group_by(Climate_Factor) %>%
  summarise(
    Mean_CV = round(mean(CV_OR, na.rm = TRUE), 2),
    Mean_OR = round(mean(Mean_OR, na.rm = TRUE), 2),
    Most_Heterogeneous_Bacteria = as.character(Bacteria[which.max(CV_OR)]),
    Max_CV = round(max(CV_OR, na.rm = TRUE), 2),
    .groups = "drop"
  )

cv_bar_plot <- ggplot(heterogeneity_summary, aes(x = Bacteria, y = CV_OR, fill = Climate_Factor)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(
    aes(label = sprintf("%.2f", CV_OR)),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 2.8
  ) +
  scale_fill_manual(values = factor_colors) +
  labs(
    title = "Spatial Heterogeneity Across Countries",
    subtitle = "Coefficient of variation of country-level OR values",
    x = NULL,
    y = "CV of OR (%)",
    fill = "Climate Factor"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

mean_or_bar_plot <- ggplot(heterogeneity_summary, aes(x = Bacteria, y = Mean_OR, fill = Climate_Factor)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", Mean_OR)),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 2.8
  ) +
  scale_fill_manual(values = factor_colors) +
  labs(
    title = "Average Country-Level Effect Size",
    subtitle = "Mean OR across countries",
    x = NULL,
    y = "Mean OR",
    fill = "Climate Factor"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

cv_heatmap_plot <- ggplot(heterogeneity_summary, aes(x = Climate_Factor, y = Bacteria, fill = CV_OR)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", CV_OR)), size = 3.2, fontface = "bold") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(
    title = "Heatmap of Spatial Heterogeneity",
    subtitle = "Coefficient of variation of country-level OR values",
    x = NULL,
    y = NULL,
    fill = "CV (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

factor_summary_plot <- ggplot(factor_summary, aes(x = Climate_Factor, y = Mean_CV, fill = Climate_Factor)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", Mean_CV)), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = factor_colors, guide = "none") +
  labs(
    title = "Overall Heterogeneity Ranking",
    subtitle = "Average CV across six AMR phenotypes",
    x = NULL,
    y = "Average CV (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank()
  )

combined_plot <- cowplot::plot_grid(
  cv_bar_plot,
  mean_or_bar_plot,
  cv_heatmap_plot,
  factor_summary_plot,
  ncol = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  align = "hv"
)

ggsave(
  file.path(output_dirs$figures, "spatial_heterogeneity_profiles_summary.pdf"),
  combined_plot,
  width = 16,
  height = 12,
  dpi = 320
)

ggsave(
  file.path(output_dirs$figures, "spatial_heterogeneity_profiles_summary.png"),
  combined_plot,
  width = 16,
  height = 12,
  dpi = 320
)

ggsave(
  file.path(output_dirs$figures, "ModelC_SpatialHeterogeneity_CV_barplot.pdf"),
  cv_bar_plot,
  width = 11.5,
  height = 5.8,
  dpi = 320
)

ggsave(
  file.path(output_dirs$figures, "ModelC_SpatialHeterogeneity_mean_OR_barplot.pdf"),
  mean_or_bar_plot,
  width = 11.5,
  height = 5.8,
  dpi = 320
)

ggsave(
  file.path(output_dirs$figures, "ModelC_SpatialHeterogeneity_CV_heatmap.pdf"),
  cv_heatmap_plot,
  width = 7.5,
  height = 5.5,
  dpi = 320
)

write.csv(
  heterogeneity_summary,
  file.path(output_dirs$tables, "ModelC_spatial_heterogeneity_summary_by_bacteria_and_factor.csv"),
  row.names = FALSE
)

write.csv(
  bacteria_summary,
  file.path(output_dirs$tables, "ModelC_spatial_heterogeneity_summary_by_bacteria.csv"),
  row.names = FALSE
)

write.csv(
  factor_summary,
  file.path(output_dirs$tables, "ModelC_spatial_heterogeneity_summary_by_climate_factor.csv"),
  row.names = FALSE
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Country_Effects")
openxlsx::writeData(wb, "Country_Effects", country_effects)
openxlsx::addWorksheet(wb, "Bacteria_Factor_Summary")
openxlsx::writeData(wb, "Bacteria_Factor_Summary", heterogeneity_summary)
openxlsx::addWorksheet(wb, "Bacteria_Summary")
openxlsx::writeData(wb, "Bacteria_Summary", bacteria_summary)
openxlsx::addWorksheet(wb, "Factor_Summary")
openxlsx::writeData(wb, "Factor_Summary", factor_summary)

openxlsx::saveWorkbook(
  wb,
  file.path(output_dirs$workbook, "spatial_heterogeneity_profiles_source_data.xlsx"),
  overwrite = TRUE
)

message("Spatial heterogeneity profile analysis completed: ", output_root)
