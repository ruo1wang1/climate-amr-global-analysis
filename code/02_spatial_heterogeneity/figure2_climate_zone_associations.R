######## Model C Figure 2B Climate-Zone Plot (Legacy Style, Optimal-K, updated data) ########

suppressPackageStartupMessages({
  required_packages <- c("ggplot2", "dplyr", "cowplot", "png", "scales", "grid")
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
figure2_root <- file.path(
  revision_root,
  "outputs/historical_associations",
  "figure2_spatial_heterogeneity_optimal_k"
)

input_csv <- file.path(
  figure2_root,
  "03_Figure2B_climate_zone_associations",
  "ModelC_Figure2B_climate_zone_effects.csv"
)

figure2a_png <- file.path(
  figure2_root,
  "02_Figure2A_spatial_heatmaps",
  "ModelC_Figure2A_spatial_heterogeneity_heatmaps_optimal_k.png"
)

figure2c_png <- file.path(
  figure2_root,
  "04_Figure2C_gap_statistic_selection",
  "ModelC_Figure2C_gap_statistic_selection.png"
)

output_dir <- file.path(figure2_root, "03_Figure2B_climate_zone_associations_oldstyle")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_csv)) {
  stop("Figure 2B climate-zone effect file not found: ", input_csv, call. = FALSE)
}

climate_effects <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
  mutate(
    Bacteria = factor(
      Bacteria,
      levels = c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")
    ),
    Climate_Zone = factor(as.character(Climate_Zone), levels = c("Temperate", "Tropical")),
    OR = as.numeric(OR),
    OR_lower = as.numeric(CI_Lower),
    OR_upper = as.numeric(CI_Upper),
    significance = case_when(
      is.na(P_Value) ~ "",
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = sprintf("%.2f (%.2f-%.2f)%s", OR, OR_lower, OR_upper, significance)
  )

zone_colors <- c("Temperate" = "#F8766D", "Tropical" = "#00BFC4")

bar_plot <- ggplot(climate_effects, aes(x = Climate_Zone, y = OR, fill = Climate_Zone)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.48) +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.12, linewidth = 0.45) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray55", linewidth = 0.5) +
  geom_text(
    aes(
      label = label,
      y = ifelse(OR >= 1, OR_upper * 1.05, pmax(OR_lower * 0.95, 0.12))
    ),
    size = 2.6,
    hjust = 0.5,
    vjust = ifelse(climate_effects$OR >= 1, 0, 1)
  ) +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = scales::number_format(accuracy = 0.01),
    limits = c(0.20, 8)
  ) +
  scale_fill_manual(values = zone_colors) +
  facet_wrap(~ Bacteria, ncol = 6, strip.position = "top") +
  labs(
    x = "Climate Zone",
    y = "Odds Ratio (OR)",
    fill = "Climate Zone",
    title = "Figure 2B. Climate zone associations with AMR strains"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size = unit(0.8, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(0.15, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(4, 6, 4, 6, unit = "pt")
  )

forest_plot <- ggplot(climate_effects, aes(x = OR, y = Climate_Zone, color = Climate_Zone)) +
  geom_point(size = 2.6) +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.15, linewidth = 0.45) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray55", linewidth = 0.5) +
  geom_text(
    aes(
      label = sprintf("OR=%.2f%s", OR, significance),
      x = OR_upper * 1.08
    ),
    size = 2.6,
    hjust = 0
  ) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = scales::number_format(accuracy = 0.01),
    limits = c(0.20, 8),
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  scale_color_manual(values = zone_colors) +
  facet_wrap(~ Bacteria, nrow = 1, strip.position = "top") +
  labs(
    x = "Odds Ratio (OR)",
    y = NULL,
    color = "Climate Zone",
    title = "Figure 2B (Alternative forest-style view)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(4, 10, 4, 4, unit = "pt")
  )

ggsave(file.path(output_dir, "ModelC_Figure2B_climate_zone_barplot_oldstyle_optimal_k.pdf"), bar_plot, width = 12, height = 4.8, dpi = 320)
ggsave(file.path(output_dir, "ModelC_Figure2B_climate_zone_barplot_oldstyle_optimal_k.png"), bar_plot, width = 12, height = 4.8, dpi = 320)

ggsave(file.path(output_dir, "ModelC_Figure2B_climate_zone_forest_oldstyle_optimal_k.pdf"), forest_plot, width = 12, height = 4.8, dpi = 320)
ggsave(file.path(output_dir, "ModelC_Figure2B_climate_zone_forest_oldstyle_optimal_k.png"), forest_plot, width = 12, height = 4.8, dpi = 320)

if (file.exists(figure2a_png)) {
  fig2a_raster <- png::readPNG(figure2a_png)
  fig2a_grob <- grid::rasterGrob(fig2a_raster, interpolate = TRUE)

  panel_a <- cowplot::ggdraw() +
    cowplot::draw_grob(fig2a_grob) +
    cowplot::draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 26, fontface = "bold")

  panel_b <- cowplot::ggdraw() +
    cowplot::draw_plot(bar_plot) +
    cowplot::draw_label("B", x = 0.01, y = 0.98, hjust = 0, vjust = 1, size = 24, fontface = "bold")

  combined_ab <- cowplot::plot_grid(panel_a, panel_b, ncol = 1, rel_heights = c(1.9, 0.82))

  ggsave(file.path(output_dir, "ModelC_Figure2_spatial_heterogeneity_panels_oldstyleB_optimal_k.pdf"), combined_ab, width = 18, height = 24, dpi = 320)
  ggsave(file.path(output_dir, "ModelC_Figure2_spatial_heterogeneity_panels_oldstyleB_optimal_k.png"), combined_ab, width = 18, height = 24, dpi = 320)

  if (file.exists(figure2c_png)) {
    fig2c_raster <- png::readPNG(figure2c_png)
    fig2c_grob <- grid::rasterGrob(fig2c_raster, interpolate = TRUE)

    panel_c <- cowplot::ggdraw() +
      cowplot::draw_grob(fig2c_grob) +
      cowplot::draw_label("C", x = 0.01, y = 0.98, hjust = 0, vjust = 1, size = 24, fontface = "bold")

    combined_abc <- cowplot::plot_grid(panel_a, cowplot::plot_grid(panel_b, panel_c, ncol = 2, rel_widths = c(1.2, 1)), ncol = 1, rel_heights = c(1.9, 0.85))

    ggsave(file.path(output_dir, "ModelC_Figure2_spatial_heterogeneity_panels_oldstyleBC_optimal_k.pdf"), combined_abc, width = 18, height = 26, dpi = 320)
    ggsave(file.path(output_dir, "ModelC_Figure2_spatial_heterogeneity_panels_oldstyleBC_optimal_k.png"), combined_abc, width = 18, height = 26, dpi = 320)
  }
}

write.csv(climate_effects, file.path(output_dir, "ModelC_Figure2B_climate_zone_oldstyle_optimal_k_plot_ready_data.csv"), row.names = FALSE)

writeLines(
  c(
    "# Model C Figure 2B Old-Style Outputs (Optimal-K)",
    "",
    "- `ModelC_Figure2B_climate_zone_barplot_oldstyle_optimal_k.*`: legacy-style Figure 2B rendering aligned to the optimal-k clustering suite.",
    "- `ModelC_Figure2B_climate_zone_forest_oldstyle_optimal_k.*`: alternative forest-style rendering.",
    "- `ModelC_Figure2_spatial_heterogeneity_panels_oldstyleB_optimal_k.*`: combined Figure 2A+B using the optimal-k Figure 2A and old-style Figure 2B.",
    "- `ModelC_Figure2_spatial_heterogeneity_panels_oldstyleBC_optimal_k.*`: combined Figure 2A+B+C using the optimal-k clustering suite."
  ),
  con = file.path(output_dir, "README.md")
)

message("Old-style optimal-k Figure 2B outputs completed: ", output_dir)
