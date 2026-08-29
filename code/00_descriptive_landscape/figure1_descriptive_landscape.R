#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggh4x)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(rnaturalearth)
  library(scales)
  library(sf)
  library(svglite)
  library(tidyr)
})

options(stringsAsFactors = FALSE, scipen = 999)
sf::sf_use_s2(FALSE)
set.seed(20260828)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else getwd()
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

input_file <- Sys.getenv("FIGURE1_INPUT_FILE", unset = "")
if (!nzchar(input_file) || !file.exists(input_file)) {
  stop("Set FIGURE1_INPUT_FILE to the local country-year Figure 1 input CSV.")
}
output_dir <- Sys.getenv(
  "FIGURE1_OUTPUT_ROOT",
  unset = file.path(repo_root, "outputs", "figure1_descriptive_landscape")
)
figure_dir <- file.path(output_dir, "figures")
source_dir <- file.path(output_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

phenotypes <- c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")
zones <- c("Temperate", "Tropical", "Polar", "All countries")

dat <- read_csv(input_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    phenotype = factor(phenotype, levels = phenotypes),
    latitude_defined_zone = factor(
      latitude_defined_zone,
      levels = c("Temperate", "Tropical", "Polar")
    )
  ) %>%
  arrange(phenotype, iso3, year)

required_columns <- c(
  "phenotype", "year", "country", "iso3", "longitude", "latitude",
  "latitude_defined_zone", "reported_resistance_proportion",
  "resistance_pct", "tested_isolates", "estimated_resistant_isolate_records"
)
if (!all(required_columns %in% names(dat))) stop("Figure 1 input columns are incomplete.")
if (nrow(dat) != 4502L || anyDuplicated(dat[c("phenotype", "iso3", "year")])) {
  stop("Figure 1 input must contain 4,502 unique country-year-phenotype records.")
}
if (n_distinct(dat$phenotype) != 6L || n_distinct(dat$iso3) != 101L) {
  stop("Figure 1 input must contain six phenotypes and 101 countries or regions.")
}
if (!all(range(dat$year) == c(1999, 2022))) stop("Unexpected year range.")
if (any(!is.finite(dat$resistance_pct)) || any(dat$resistance_pct < 0 | dat$resistance_pct > 100)) {
  stop("Resistance percentages must be finite and within 0-100.")
}
if (any(!is.finite(dat$tested_isolates)) || any(dat$tested_isolates <= 0)) {
  stop("Tested-isolate counts must be positive.")
}

font_family <- if (
  requireNamespace("systemfonts", quietly = TRUE) &&
  nzchar(tryCatch(systemfonts::match_fonts("Arial")$path[[1]], error = function(e) ""))
) "Arial" else "sans"

theme_pub <- theme_classic(base_family = font_family, base_size = 6.2) +
  theme(
    axis.line = element_line(linewidth = 0.30, colour = "#626A70"),
    axis.ticks = element_line(linewidth = 0.28, colour = "#626A70"),
    axis.title = element_text(size = 6.2, face = "bold", colour = "#303438"),
    axis.text = element_text(size = 5.5, colour = "#444B50"),
    strip.text = element_text(size = 6.5, face = "bold", colour = "#282B2E"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E4E8EA", linewidth = 0.22),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.6),
    plot.tag = element_text(size = 8.3, face = "bold", colour = "#111111"),
    plot.tag.position = c(0.002, 0.997),
    plot.margin = margin(2, 3, 2, 3, unit = "pt")
  )

country_period <- dat %>%
  group_by(phenotype, iso3, country, longitude, latitude, latitude_defined_zone) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    reporting_years = n_distinct(year),
    observations = n(),
    tested_isolates = sum(tested_isolates),
    estimated_resistant_isolate_records = sum(estimated_resistant_isolate_records),
    resistance_pct = 100 * estimated_resistant_isolate_records / tested_isolates,
    .groups = "drop"
  ) %>%
  arrange(phenotype, country)

rate_breaks <- c(-Inf, 1, 5, 10, 25, 50, Inf)
rate_labels <- c("0-<1", "1-<5", "5-<10", "10-<25", "25-<50", ">=50")
map_colours <- c(
  "0-<1" = "#EFF5F7", "1-<5" = "#D5E5EC", "5-<10" = "#ACC6D6",
  "10-<25" = "#8997C4", "25-<50" = "#8B59AD", ">=50" = "#76206F",
  "No data" = "#ECEBE8"
)
country_period <- country_period %>%
  mutate(resistance_bin = cut(
    resistance_pct, breaks = rate_breaks, labels = rate_labels,
    right = FALSE, ordered_result = TRUE
  ))

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") %>%
  mutate(map_iso3 = adm0_a3) %>%
  st_make_valid()

map_sf <- map(phenotypes, function(ph) {
  values <- country_period %>%
    filter(as.character(phenotype) == ph) %>%
    select(iso3, resistance_bin)
  world %>%
    left_join(values, by = c("map_iso3" = "iso3")) %>%
    mutate(
      phenotype = factor(ph, levels = phenotypes),
      map_category = factor(
        if_else(is.na(resistance_bin), "No data", as.character(resistance_bin)),
        levels = c(rate_labels, "No data")
      )
    )
}) %>% bind_rows() %>% st_as_sf()

equal_earth <- "+proj=eqearth +datum=WGS84 +units=m +no_defs"
lon_line <- function(lon) st_linestring(cbind(rep(lon, 361), seq(-89.999, 89.999, length.out = 361)))
lat_line <- function(lat) st_linestring(cbind(seq(-179.999, 179.999, length.out = 721), rep(lat, 721)))
graticule <- st_sfc(
  c(map(seq(-150, 150, 30), lon_line), map(seq(-60, 60, 30), lat_line)), crs = 4326
) %>% st_transform(equal_earth)
latitude_values <- c(-66.5, -23.5, 23.5, 66.5)
latitude_lines <- st_sfc(map(latitude_values, lat_line), crs = 4326) %>%
  st_transform(equal_earth)
latitude_labels <- expand_grid(
  phenotype = c("3GCR-Ec", "CR-Ec"), latitude = latitude_values
) %>%
  mutate(
    phenotype = factor(phenotype, levels = phenotypes),
    longitude = -168,
    label = case_when(
      latitude == 66.5 ~ "66.5°N", latitude == 23.5 ~ "23.5°N",
      latitude == -23.5 ~ "23.5°S", TRUE ~ "66.5°S"
    )
  ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(equal_earth)
latitude_label_xy <- bind_cols(st_drop_geometry(latitude_labels), as_tibble(st_coordinates(latitude_labels)))

edge_lon <- 179.999
edge_lat <- 89.999
outline <- rbind(
  cbind(seq(-edge_lon, edge_lon, length.out = 721), -edge_lat),
  cbind(edge_lon, seq(-edge_lat, edge_lat, length.out = 361)),
  cbind(seq(edge_lon, -edge_lon, length.out = 721), edge_lat),
  cbind(-edge_lon, seq(edge_lat, -edge_lat, length.out = 361))
)
globe_outline <- st_sfc(st_linestring(outline), crs = 4326) %>% st_transform(equal_earth)

p_a <- ggplot() +
  geom_sf(data = graticule, colour = "#D9DDDF", linewidth = 0.045) +
  geom_sf(data = map_sf, aes(fill = map_category), colour = "#A6AAAD", linewidth = 0.095) +
  geom_sf(data = latitude_lines, colour = "#7C858A", linewidth = 0.072, linetype = "22") +
  geom_text(
    data = latitude_label_xy, aes(X, Y, label = label), hjust = 0,
    size = 1.32, family = font_family, colour = "#596066"
  ) +
  geom_sf(data = globe_outline, colour = "#555D62", linewidth = 0.092) +
  facet_wrap(~phenotype, ncol = 3) +
  scale_fill_manual(
    values = map_colours, limits = c(rate_labels, "No data"), drop = FALSE,
    name = "Aggregated prevalence of AMR (%)"
  ) +
  coord_sf(crs = equal_earth, datum = NA, expand = FALSE) +
  labs(tag = "a") +
  theme_void(base_family = font_family, base_size = 6.2) +
  theme(
    strip.text = element_text(size = 6.5, face = "bold", colour = "#282B2E"),
    panel.spacing = grid::unit(1.7, "mm"),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.6),
    legend.key.height = grid::unit(3.1, "mm"),
    legend.key.width = grid::unit(4.4, "mm"),
    legend.spacing.x = grid::unit(0.5, "mm"),
    legend.margin = margin(0), legend.box.margin = margin(-2, 0, 0, 0),
    plot.tag = element_text(size = 8.3, face = "bold"),
    plot.tag.position = c(0.002, 0.997),
    plot.margin = margin(1, 3, 1, 3, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"))

bootstrap_reps <- 2000L
bootstrap_median_ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(c(lower = NA_real_, upper = NA_real_))
  draws <- replicate(bootstrap_reps, median(sample(x, length(x), replace = TRUE)))
  setNames(quantile(draws, c(0.025, 0.975), names = FALSE), c("lower", "upper"))
}

zone_input <- bind_rows(
  country_period %>% mutate(zone = as.character(latitude_defined_zone)),
  country_period %>% mutate(zone = "All countries")
) %>% mutate(zone = factor(zone, levels = zones))

zone_summary <- zone_input %>%
  group_by(phenotype, zone) %>%
  group_modify(function(.x, .y) {
    interval <- bootstrap_median_ci(.x$resistance_pct)
    tibble(
      countries = n_distinct(.x$iso3),
      tested_isolates = sum(.x$tested_isolates),
      median = median(.x$resistance_pct),
      ci_lower = unname(interval[["lower"]]),
      ci_upper = unname(interval[["upper"]])
    )
  }) %>%
  ungroup() %>%
  mutate(
    zone = factor(zone, levels = zones),
    zone_x = recode(as.character(zone), Temperate = 1L, Tropical = 2L, Polar = 3L, `All countries` = 2L),
    show_ci = countries >= 3L & is.finite(ci_lower) & is.finite(ci_upper)
  )

zone_colours <- c(Temperate = "#C95D63", Tropical = "#2F78A8", Polar = "#7D8792", `All countries` = "#62568B")
zone_fills <- c(Temperate = "#F8E1E2", Tropical = "#DFEBF3", Polar = "#E8EBED", `All countries` = "#E8E4F0")
zone_shapes <- c(Temperate = 21, Tropical = 24, Polar = 23, `All countries` = 22)

size_limits <- range(zone_summary$tested_isolates)
tested_to_display <- function(x) {
  z <- (log10(x) - log10(size_limits[[1]])) / diff(log10(size_limits))
  pmin(1, pmax(0, z))^2.5
}
zone_summary <- zone_summary %>% mutate(size_display = tested_to_display(tested_isolates))
size_breaks <- c(1e3, 1e5, 1e6, 5e6)
size_labels <- c("1k", "100k", "1m", "5m")
size_positions <- tested_to_display(size_breaks)

b_breaks <- list(
  `3GCR-Ec` = c(0, 20, 40, 60, 80), `3GCR-Kp` = c(0, 20, 40, 60, 80),
  `CR-Ab` = c(0, 25, 50, 75, 100), `CR-Ec` = c(0, 1.5, 3, 4.5, 6),
  `CR-Kp` = c(0, 6, 12, 18, 24), `CR-Pa` = c(0, 10, 20, 30, 40)
)

b_scales <- function() list(
  scale_colour_manual(values = zone_colours, limits = zones, guide = "none"),
  scale_fill_manual(values = zone_fills, limits = zones, guide = "none"),
  scale_shape_manual(values = zone_shapes, limits = zones, guide = "none"),
  scale_radius(limits = c(0, 1), range = c(0.72, 2.45), guide = "none")
)

make_b_plot <- function(ph, section) {
  d <- zone_summary %>% filter(as.character(phenotype) == ph)
  d <- if (section == "zones") filter(d, as.character(zone) != "All countries") else filter(d, as.character(zone) == "All countries")
  brks <- b_breaks[[ph]]
  ggplot(d, aes(zone_x, colour = zone, fill = zone, shape = zone)) +
    geom_linerange(
      data = filter(d, show_ci), aes(ymin = ci_lower, ymax = ci_upper),
      linewidth = 0.56, lineend = "round", show.legend = FALSE
    ) +
    geom_point(aes(y = median, size = size_display), stroke = 0.46) +
    scale_x_continuous(limits = c(0.60, 3.40), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max(brks)), breaks = brks, labels = label_number(trim = TRUE), expand = c(0, 0)) +
    b_scales() +
    labs(title = if (section == "zones") ph else NULL, x = if (section == "all") ph else NULL, y = NULL) +
    theme_pub +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.line.x = element_line(colour = "#687177", linewidth = 0.28),
      axis.title.x = element_text(size = 6.3, face = "bold", margin = margin(2, 0, 0, 0)),
      axis.title.y = element_blank(), axis.text.y = element_text(size = 5.2),
      plot.title = element_text(size = 6.4, face = "bold", hjust = 0.5, margin = margin(0, 0, 1.5, 0)),
      legend.position = "none", plot.margin = margin(1, 1.6, 1, 1.6, unit = "pt")
    )
}

row_label <- function(label) ggplot() +
  annotate("text", 0, 0.5, label = label, hjust = 0, family = font_family, fontface = "bold", size = 2.15, colour = "#303438") +
  xlim(0, 1) + ylim(0, 1) + theme_void()

top_b <- wrap_plots(map(phenotypes, ~make_b_plot(.x, "zones")), nrow = 1)
bottom_b <- wrap_plots(map(phenotypes, ~make_b_plot(.x, "all")), nrow = 1)

zone_legend <- tibble(
  zone = factor(zones, levels = zones), y = c(9.0, 8.05, 7.10, 6.15),
  shape = unname(zone_shapes[zones]), colour = unname(zone_colours[zones]), fill = unname(zone_fills[zones])
)
size_legend <- tibble(
  label = size_labels, y = c(3.55, 2.75, 1.95, 1.15),
  point_size = 0.72 + (2.45 - 0.72) * size_positions
)

p_b_legend <- ggplot() +
  annotate("text", 0, 10.0, label = "Climate zones", hjust = 0, family = font_family, fontface = "bold", size = 2.05) +
  geom_point(data = zone_legend, aes(0.38, y), shape = zone_legend$shape, colour = zone_legend$colour, fill = zone_legend$fill, size = 2.05, stroke = 0.48) +
  geom_text(data = zone_legend, aes(0.88, y, label = zone), hjust = 0, family = font_family, size = 1.92, colour = "#303438") +
  annotate("text", 0, 4.55, label = "Tested isolates", hjust = 0, family = font_family, fontface = "bold", size = 2.05) +
  geom_point(data = size_legend, aes(0.38, y), shape = 21, colour = "#6F777C", fill = "white", size = size_legend$point_size, stroke = 0.42) +
  geom_text(data = size_legend, aes(0.88, y, label = label), hjust = 0, family = font_family, size = 1.90, colour = "#303438") +
  coord_cartesian(xlim = c(0, 3.2), ylim = c(0.5, 10.4), clip = "off") + theme_void()

p_b_body <- (row_label("Climate zones") / top_b / row_label("All countries") / bottom_b) +
  plot_layout(heights = c(0.13, 1.12, 0.13, 0.78))
p_b_y <- wrap_elements(full = grid::grobTree(
  grid::textGrob("b", x = grid::unit(0.12, "npc"), y = grid::unit(0.99, "npc"), just = c("left", "top"), gp = grid::gpar(fontfamily = font_family, fontsize = 8.3, fontface = "bold")),
  grid::textGrob("Prevalence of AMR (%)", rot = 90, gp = grid::gpar(fontfamily = font_family, fontsize = 6.2, fontface = "bold", col = "#303438"))
))
p_b <- (p_b_y | p_b_body | p_b_legend) + plot_layout(widths = c(0.045, 0.81, 0.145))

trajectory_points <- dat %>% select(phenotype, iso3, country, year, resistance_pct)
trajectory_lines <- dat %>%
  arrange(phenotype, iso3, year) %>%
  group_by(phenotype, iso3) %>%
  mutate(run_id = cumsum(is.na(lag(year)) | year - lag(year) > 1L)) %>%
  group_by(phenotype, iso3, country, run_id) %>%
  filter(n() >= 2L) %>%
  ungroup()

annual_summary <- dat %>%
  group_by(phenotype, year) %>%
  summarise(
    countries = n_distinct(iso3),
    estimate = mean(resistance_pct),
    annual_tested_isolates = sum(tested_isolates),
    .groups = "drop"
  )

format_axis <- function(x) vapply(x, function(v) sub("\\.?0+$", "", sprintf("%.3f", v)), character(1))
bounded_axis <- function(x) {
  candidates <- c(1, 2, 4, 6, 10, 20, 40, 60, 80, 100)
  upper <- min(candidates[candidates >= min(100, x * 1.02)])
  list(upper = upper, breaks = seq(0, upper, length.out = 5))
}
c_maxima <- bind_rows(
  trajectory_points %>% transmute(phenotype, value = resistance_pct),
  annual_summary %>% transmute(phenotype, value = estimate)
) %>%
  group_by(phenotype) %>% summarise(raw_upper = max(value), .groups = "drop") %>%
  mutate(phenotype = factor(phenotype, levels = phenotypes)) %>% arrange(phenotype)
c_axes <- c_maxima %>%
  mutate(axis = map2(as.character(phenotype), raw_upper, function(ph, mx) {
    if (ph == "CR-Ec") {
      top <- max(40, ceiling(mx / 10) * 10)
      list(upper = max(top, mx * 1.03), breaks = c(0, 0.1, 1, 10, top))
    } else bounded_axis(mx)
  }))
c_scales <- map(seq_len(nrow(c_axes)), function(i) {
  ph <- as.character(c_axes$phenotype[[i]])
  brks <- c_axes$axis[[i]]$breaks
  args <- list(limits = c(0, c_axes$axis[[i]]$upper), breaks = brks, labels = format_axis(brks), expand = expansion(mult = c(0, 0.03)))
  if (ph == "CR-Ec") args$trans <- pseudo_log_trans(base = 10, sigma = 0.1)
  do.call(scale_y_continuous, args)
})

trajectory_colour <- "#AAB5BB"
summary_colour <- "#3979A7"
p_c <- ggplot() +
  geom_line(
    data = trajectory_lines,
    aes(year, resistance_pct, group = interaction(iso3, run_id, drop = TRUE)),
    colour = trajectory_colour, linewidth = 0.17, alpha = 0.21, lineend = "round"
  ) +
  geom_point(data = trajectory_points, aes(year, resistance_pct), colour = trajectory_colour, size = 0.20, alpha = 0.14, stroke = 0) +
  geom_line(data = annual_summary, aes(year, estimate), colour = summary_colour, linewidth = 0.74, lineend = "round") +
  geom_point(data = annual_summary, aes(year, estimate), shape = 21, fill = "white", colour = summary_colour, size = 0.88, stroke = 0.34) +
  facet_wrap(~phenotype, ncol = 3, scales = "free_y") +
  scale_x_continuous(limits = c(1999, 2022), breaks = c(2000, 2010, 2020), minor_breaks = NULL, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, length.out = 5), labels = label_number(trim = TRUE), expand = expansion(mult = c(0, 0.03))) +
  ggh4x::facetted_pos_scales(y = c_scales) +
  labs(x = "Year", y = "AMR of different phenotypes (%)", tag = "c") +
  theme_pub +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#626A70", linewidth = 0.30),
    panel.spacing = grid::unit(2.5, "mm"),
    axis.text.x = element_text(size = 5.5), axis.text.y = element_text(size = 5.4),
    axis.title = element_text(size = 6.2, face = "bold"), strip.text = element_text(size = 6.5),
    plot.margin = margin(1, 4, 1, 4, unit = "pt")
  )

figure1 <- p_a / p_b / p_c + plot_layout(heights = c(1.00, 0.90, 1.04))
stem <- file.path(figure_dir, "Figure1_global_geographical_temporal_AMR")
ggsave(paste0(stem, ".pdf"), figure1, device = grDevices::cairo_pdf, width = 183, height = 205, units = "mm", bg = "white")
svglite(paste0(stem, ".svg"), width = 183 / 25.4, height = 205 / 25.4, bg = "white")
print(figure1)
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = 183, height = 205, units = "mm", res = 600, compression = "lzw", background = "white")
print(figure1)
dev.off()
ragg::agg_png(paste0(stem, ".png"), width = 183, height = 205, units = "mm", res = 300, background = "white")
print(figure1)
dev.off()

write_csv(country_period %>% mutate(resistance_bin = as.character(resistance_bin)), file.path(source_dir, "figure1_panel_a_country_period_resistance.csv"))
write_csv(zone_input, file.path(source_dir, "figure1_panel_b_country_input.csv"))
write_csv(zone_summary, file.path(source_dir, "figure1_panel_b_zone_summary.csv"))
write_csv(trajectory_points, file.path(source_dir, "figure1_panel_c_country_year_trajectories.csv"))
write_csv(annual_summary, file.path(source_dir, "figure1_panel_c_annual_country_means.csv"))

unmatched <- country_period %>% distinct(iso3) %>% filter(!iso3 %in% world$map_iso3)
qa <- tibble(
  check = c("input rows", "unique keys", "phenotypes", "countries", "years", "map identifiers", "panel a records", "panel c annual means"),
  passed = c(
    nrow(dat) == 4502L,
    !anyDuplicated(dat[c("phenotype", "iso3", "year")]),
    n_distinct(dat$phenotype) == 6L,
    n_distinct(dat$iso3) == 101L,
    all(range(dat$year) == c(1999, 2022)),
    nrow(unmatched) == 0L,
    sum(country_period$observations) == nrow(dat),
    max(abs(
      annual_summary$estimate -
        dat %>% group_by(phenotype, year) %>% summarise(x = mean(resistance_pct), .groups = "drop") %>% pull(x)
    )) < 1e-12
  )
)
write_csv(qa, file.path(source_dir, "figure1_qa_checks.csv"))
if (any(!qa$passed)) stop("Figure 1 quality-control checks failed.")

message("Figure 1 outputs written to: ", output_dir)
