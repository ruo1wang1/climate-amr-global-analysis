######## Figure 1B Statistical Summary Plot for primary model ########

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(scales)
  library(cowplot)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
fig1b_input_csv <- file.path(
  revision_root,
  "data/source_data",
  "figure1_historical_associations",
  "01_csv",
  "Figure1_ModelC_Fig1B_summary.csv"
)

fig1b_output_dir <- file.path(
  revision_root,
  "outputs/historical_associations",
  "analysis_historical_associations_main_model",
  "05_figure1B_statistical_summary"
)

climate_colors <- c(
  "Temperature" = "#DD5F60",
  "Relative Humidity" = "#3CB371",
  "Precipitation" = "#9AC0CD",
  "Wet Days" = "#8A2BE2"
)

climate_shapes <- c(
  "Temperature" = 16,
  "Relative Humidity" = 17,
  "Precipitation" = 18,
  "Wet Days" = 15
)

prepare_fig1b_data <- function(input_csv) {
  if (!file.exists(input_csv)) {
    stop("Figure 1B summary CSV not found: ", input_csv, call. = FALSE)
  }

  read.csv(input_csv, stringsAsFactors = FALSE) %>%
    mutate(
      Factor = recode(
        Climate_Variable,
        "Temperature" = "Temperature",
        "Humidity" = "Relative Humidity",
        "Precipitation" = "Precipitation",
        "WetDays" = "Wet Days"
      ),
      significance_label = case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01 ~ "**",
        P_Value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Bacteria = factor(Bacteria, levels = c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")),
      Factor = factor(Factor, levels = rev(c("Temperature", "Relative Humidity", "Precipitation", "Wet Days"))),
      EDF_cat = cut(
        EDF,
        breaks = c(-0.001, 1, 2, 3, 4, 6),
        labels = c("1", "2", "3", "4", "5"),
        include.lowest = TRUE
      ),
      EDF_cat = factor(EDF_cat, levels = c("1", "2", "3", "4", "5")),
      label_text = sprintf("F=%.2f%s\nEDF=%.2f", F_statistic, significance_label, EDF)
    ) %>%
    arrange(Bacteria, Factor)
}

make_fig1b_plot <- function(plot_data) {
  max_f <- max(plot_data$F_statistic, na.rm = TRUE)
  x_limit <- max(25, ceiling(max_f / 5) * 5 + 10)
  x_breaks <- seq(0, x_limit, by = 5)
  label_nudge <- max(1.6, x_limit * 0.055)
  tick_color <- "gray60"

  base_plot <- ggplot(plot_data, aes(x = F_statistic, y = Factor)) +
    theme_bw(base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray96", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 11, color = "gray25"),
      axis.text.y = element_text(size = 12, face = "bold", color = "gray25"),
      strip.text = element_text(size = 13, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      legend.position = "right",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      panel.spacing = unit(0.25, "lines"),
      legend.box = "vertical",
      legend.key.size = unit(1.15, "lines"),
      legend.key = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", linewidth = 0.8),
      axis.ticks.length = unit(0.20, "cm"),
      axis.ticks = element_line(color = tick_color),
      axis.ticks.x = element_line(color = tick_color),
      axis.ticks.y = element_line(color = tick_color),
      plot.margin = margin(t = 8, r = 14, b = 8, l = 8)
    ) +
    geom_point(aes(size = EDF_cat, shape = Factor, color = Factor, alpha = EDF_cat)) +
    geom_text(
      aes(label = label_text),
      nudge_x = label_nudge,
      hjust = 0,
      size = 3.9,
      family = "serif",
      lineheight = 1.05
    ) +
    facet_wrap(~Bacteria, nrow = 1, strip.position = "top") +
    scale_color_manual(values = climate_colors, name = "Climate Factor") +
    scale_size_manual(values = c("1" = 2.6, "2" = 3.6, "3" = 4.6, "4" = 5.6, "5" = 6.6), name = "EDF Value") +
    scale_alpha_manual(values = c("1" = 0.35, "2" = 0.5, "3" = 0.7, "4" = 0.85, "5" = 1.0), name = "EDF Value") +
    scale_shape_manual(values = climate_shapes, name = "Climate Factor") +
    scale_x_continuous(
      limits = c(0, x_limit),
      breaks = x_breaks,
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    labs(x = "F value") +
    guides(
      shape = guide_legend(
        order = 1,
        override.aes = list(size = 5, alpha = 1, color = climate_colors)
      ),
      size = guide_legend(
        order = 2,
        override.aes = list(shape = 16, alpha = c(0.35, 0.5, 0.7, 0.85, 1.0), color = "black")
      ),
      alpha = "none",
      color = "none"
    )

  ggdraw() +
    draw_label("B", x = 0.018, y = 0.92, hjust = 0, vjust = 1, fontfamily = "serif", size = 34) +
    draw_plot(base_plot, x = 0.08, y = 0, width = 0.92, height = 1)
}

main <- function() {
  if (!dir.exists(fig1b_output_dir)) {
    dir.create(fig1b_output_dir, recursive = TRUE)
  }

  plot_data <- prepare_fig1b_data(fig1b_input_csv)
  plot_obj <- make_fig1b_plot(plot_data)

  plot_data_csv <- file.path(fig1b_output_dir, "ModelC_GAMM_Figure1B_plot_ready_data.csv")
  pdf_path <- file.path(fig1b_output_dir, "ModelC_GAMM_Figure1B_statistical_summary.pdf")
  png_path <- file.path(fig1b_output_dir, "ModelC_GAMM_Figure1B_statistical_summary.png")

  write.csv(plot_data, plot_data_csv, row.names = FALSE)

  ggsave(pdf_path, plot = plot_obj, width = 20, height = 5.8, device = cairo_pdf, dpi = 300)
  ggsave(png_path, plot = plot_obj, width = 20, height = 5.8, dpi = 300)

  cat("Saved Figure 1B to:\n")
  cat(pdf_path, "\n")
  cat(png_path, "\n")
  cat(plot_data_csv, "\n")
}

main()
