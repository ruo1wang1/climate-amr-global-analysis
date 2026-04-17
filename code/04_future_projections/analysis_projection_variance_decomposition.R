######## Projection preparation: variance decomposition ########

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(ggplot2)
  library(scales)
  library(writexl)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
input_data_dir <- file.path(revision_root, "outputs/historical_associations", "model_ready_inputs")
lag_summary_path <- file.path(
  revision_root,
  "data",
  "source_data",
  "lag_selection",
  "historical_model_c",
  "01_csv",
  "historical_lag_summary_model_c.csv"
)
output_base <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "projection_preparation",
  "01_variance_decomposition"
)

dir.create(file.path(output_base, "01_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "02_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "03_workbook"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "04_metadata"), recursive = TRUE, showWarnings = FALSE)

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = "3GCREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = "3GCRKP_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRAB", title = "CR-Ab", file_name = "CRAB_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CREC", title = "CR-Ec", file_name = "CREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRKP", title = "CR-Kp", file_name = "CRKP_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRPA", title = "CR-Pa", file_name = "CRPA_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv")
)

load_lag_settings <- function(summary_csv) {
  lag_summary <- read.csv(summary_csv, stringsAsFactors = FALSE)
  lag_settings <- setNames(vector("list", nrow(lag_summary)), lag_summary$Display_Name)
  for (i in seq_len(nrow(lag_summary))) {
    lag_settings[[lag_summary$Display_Name[i]]] <- list(
      temp_lag = as.integer(lag_summary$TMP_lag[i]),
      precip_lag = as.integer(lag_summary$PREC_lag[i]),
      humid_lag = as.integer(lag_summary$HUM_lag[i]),
      wetdays_lag = as.integer(lag_summary$WET_lag[i])
    )
  }
  lag_settings
}

get_available_pls_components <- function(data) {
  pls_candidates <- paste0("PLS_Comp", 1:4)
  present <- pls_candidates[pls_candidates %in% names(data)]
  present[sapply(present, function(x) !all(is.na(data[[x]])))]
}

prepare_data <- function(file_path, lag_config) {
  data <- read.csv(file_path)

  data_processed <- data %>%
    mutate(
      year = as.numeric(as.character(year)),
      Region = factor(Region),
      climate_zone = case_when(
        abs(lat) > 66.5 ~ "Polar Zone",
        abs(lat) > 23.5 ~ "Temperate Zone",
        TRUE ~ "Tropical Zone"
      ),
      climate_zone = factor(climate_zone)
    ) %>%
    group_by(NAME) %>%
    mutate(location_id = cur_group_id()) %>%
    ungroup() %>%
    mutate(HUM = pmin(HUM, 100))

  data_processed %>%
    mutate(
      TMP_orig = TMP,
      PREC_orig = PREC,
      HUM_orig = HUM,
      WET_orig = WET
    ) %>%
    group_by(climate_zone) %>%
    mutate(
      across(c(TMP, PREC, HUM, WET), \(x) as.vector(scale(x)), .names = "{.col}_scaled")
    ) %>%
    group_by(location_id) %>%
    arrange(year) %>%
    mutate(
      TMP_scaled_lag = lag(TMP_scaled, lag_config$temp_lag),
      PREC_scaled_lag = lag(PREC_scaled, lag_config$precip_lag),
      HUM_scaled_lag = lag(HUM_scaled, lag_config$humid_lag),
      WET_scaled_lag = lag(WET_scaled, lag_config$wetdays_lag)
    ) %>%
    filter(
      !is.na(TMP_scaled_lag),
      !is.na(PREC_scaled_lag),
      !is.na(HUM_scaled_lag),
      !is.na(WET_scaled_lag)
    ) %>%
    ungroup()
}

build_model <- function(data, include_pls = TRUE) {
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
  pls_terms <- character(0)
  if (include_pls) {
    available_pls <- get_available_pls_components(data)
    pls_terms <- paste0("s(", available_pls, ", k = 10, bs = 'cr')")
  }

  rhs_terms <- c(
    "s(TMP_scaled_lag, k = 5, bs = 'cr')",
    "s(PREC_scaled_lag, k = 10, bs = 'cr')",
    "s(HUM_scaled_lag, k = 10, bs = 'cr')",
    "s(WET_scaled_lag, k = 10, bs = 'cr')",
    pls_terms,
    "s(lat, lon, bs = 'sos', k = 20)",
    "s(year, bs = 'cr', k = 8)",
    "s(Region, bs = 're')",
    "climate_zone"
  )

  model_formula <- as.formula(paste("logit_R ~", paste(rhs_terms, collapse = " + ")))

  tryCatch(
    bam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE,
      control = ctrl
    ),
    error = function(e) {
      gam(
        model_formula,
        data = data,
        family = gaussian(),
        method = "REML",
        select = TRUE
      )
    }
  )
}

safe_dev_explained <- function(model) {
  dev <- summary(model)$dev.expl * 100
  ifelse(is.finite(dev), dev, NA_real_)
}

lag_settings <- load_lag_settings(lag_summary_path)
results <- list()

for (spec in bacteria_specs) {
  lag_config <- lag_settings[[spec$title]]
  data_path <- file.path(input_data_dir, spec$file_name)
  prepared <- prepare_data(data_path, lag_config)

  full_model <- build_model(prepared, include_pls = TRUE)
  climate_only_model <- build_model(prepared, include_pls = FALSE)

  full_dev <- safe_dev_explained(full_model)
  climate_dev <- safe_dev_explained(climate_only_model)

  climate_component <- min(climate_dev, full_dev, na.rm = TRUE)
  pls_component <- max(full_dev - climate_dev, 0, na.rm = TRUE)

  climate_pct <- ifelse(full_dev > 0, climate_component / full_dev * 100, NA_real_)
  pls_pct <- ifelse(full_dev > 0, pls_component / full_dev * 100, NA_real_)

  results[[length(results) + 1]] <- tibble(
    Bacteria = spec$title,
    N = nrow(prepared),
    Full_Model_DevExplained_pct = round(full_dev, 3),
    Climate_Only_DevExplained_pct = round(climate_dev, 3),
    Climate_Component_pct_of_total = round(climate_pct, 3),
    PLS_Component_pct_of_total = round(pls_pct, 3),
    Climate_Component_raw = round(climate_component, 3),
    PLS_Component_raw = round(pls_component, 3),
    TMP_lag = lag_config$temp_lag,
    PREC_lag = lag_config$precip_lag,
    HUM_lag = lag_config$humid_lag,
    WET_lag = lag_config$wetdays_lag,
    Full_AIC = round(AIC(full_model), 3),
    Climate_Only_AIC = round(AIC(climate_only_model), 3),
    PLS_Count = length(get_available_pls_components(prepared))
  )
}

variance_df <- bind_rows(results) %>%
  mutate(
    Bacteria = factor(Bacteria, levels = c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa"))
  ) %>%
  arrange(Bacteria)

weighted_climate_pct <- with(
  variance_df,
  weighted.mean(Climate_Component_pct_of_total, Full_Model_DevExplained_pct, na.rm = TRUE)
)
weighted_pls_pct <- with(
  variance_df,
  weighted.mean(PLS_Component_pct_of_total, Full_Model_DevExplained_pct, na.rm = TRUE)
)

plot_df <- variance_df %>%
  select(Bacteria, Climate_Component_pct_of_total, PLS_Component_pct_of_total) %>%
  pivot_longer(
    cols = c(Climate_Component_pct_of_total, PLS_Component_pct_of_total),
    names_to = "Source",
    values_to = "Percentage"
  ) %>%
  mutate(
    Source = recode(
      Source,
      Climate_Component_pct_of_total = "Climate Variables",
      PLS_Component_pct_of_total = "PLS Components"
    )
  )

color_map <- c("Climate Variables" = "#2E62B0", "PLS Components" = "#C5336C")

s14_plot <- ggplot(plot_df, aes(x = Bacteria, y = Percentage, fill = Source)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 4.2,
    family = "serif"
  ) +
  scale_fill_manual(values = color_map) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(
    title = "Contribution Proportion of Climatic Variables and PLS Components",
    subtitle = "Derived from Multidimensional Covariates to the Explained Variance in Model C",
    x = NULL,
    y = "Percentage of Total Deviance Explained (%)",
    fill = "Source"
  ) +
  theme_bw(base_family = "serif", base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 12)
  )

tables_dir <- file.path(output_base, "01_tables")
figures_dir <- file.path(output_base, "02_figures")
workbook_dir <- file.path(output_base, "03_workbook")
metadata_dir <- file.path(output_base, "04_metadata")

write.csv(
  variance_df,
  file.path(tables_dir, "projection_variance_decomposition_variance_decomposition_by_bacteria.csv"),
  row.names = FALSE
)
write.csv(
  plot_df,
  file.path(tables_dir, "projection_variance_decomposition_plot_ready_data.csv"),
  row.names = FALSE
)
write.csv(
  tibble(
    Weighted_Climate_Contribution_pct = round(weighted_climate_pct, 3),
    Weighted_PLS_Contribution_pct = round(weighted_pls_pct, 3)
  ),
  file.path(tables_dir, "projection_variance_decomposition_weighted_average_summary.csv"),
  row.names = FALSE
)

ggsave(
  file.path(figures_dir, "projection_variance_decomposition_variance_decomposition.pdf"),
  s14_plot,
  width = 13,
  height = 8.5,
  dpi = 300
)
ggsave(
  file.path(figures_dir, "projection_variance_decomposition_variance_decomposition.png"),
  s14_plot,
  width = 13,
  height = 8.5,
  dpi = 300
)

write_xlsx(
  list(
    variance_by_bacteria = variance_df,
    variance_plot_ready = plot_df,
    weighted_summary = tibble(
      Weighted_Climate_Contribution_pct = round(weighted_climate_pct, 3),
      Weighted_PLS_Contribution_pct = round(weighted_pls_pct, 3)
    )
  ),
  file.path(workbook_dir, "projection_variance_decomposition_variance_decomposition_source_data.xlsx")
)

readme_lines <- c(
  "# analysis figures",
  "",
  "This folder contains the variance decomposition analysis under the current Model C data system.",
  "The full-model specification here is intentionally aligned with the current main Model C used for historical threshold and response analyses, rather than the broader lag-search fitting template used during lag selection.",
  "",
  "Definitions:",
  "- `Full_Model_DevExplained_pct`: explained deviance of the full Model C with dynamic PLS components.",
  "- `Climate_Only_DevExplained_pct`: explained deviance of the matched climate-only model using the same lag structure and spatiotemporal terms but excluding all PLS components.",
  "- `Climate_Component_pct_of_total`: climate-only explained deviance expressed as a percentage of the full-model explained deviance, truncated at 100% when the climate-only model slightly exceeds the full-model deviance due to smoothing/penalization differences.",
  "- `PLS_Component_pct_of_total`: residual contribution attributed to PLS components, defined as the non-negative difference between full-model and climate-only explained deviance.",
  "",
  sprintf("Weighted average climate contribution: %.2f%%", weighted_climate_pct),
  sprintf("Weighted average PLS contribution: %.2f%%", weighted_pls_pct)
)
writeLines(readme_lines, file.path(metadata_dir, "README.md"))

session_txt <- capture.output(sessionInfo())
writeLines(session_txt, file.path(metadata_dir, "sessionInfo.txt"))

cat("Saved variance decomposition outputs to:\n")
cat(output_base, "\n")
