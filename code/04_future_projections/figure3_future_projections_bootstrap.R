suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(openxlsx)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(digest)
})

set.seed(20260331)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
input_data_dir <- file.path(revision_root, "outputs/historical_associations", "model_ready_inputs")
legacy_cmip6_root <- Sys.getenv("CLIMATE_AMR_CMIP6_ROOT", unset = file.path(revision_root, "data", "cmip6_archive"))
wetdays_yearly_path <- Sys.getenv("CLIMATE_AMR_WET_DAYS_INPUT", unset = file.path(revision_root, "data", "projection_inputs", "yearly_wet_days_panel_data.csv"))
country_metadata_path <- file.path(legacy_cmip6_root, "cleaned_results", "climate_data_yearly_simple.csv")

projection_input_dir <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "projection_preparation",
  "03_projection_simplified_gamm_ready_inputs"
)

climate_input_dir <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "projection_preparation",
  "01_climate_scenario_inputs"
)
climate_data_path <- file.path(
  climate_input_dir,
  "cmip6_four_model_country_year_climate_input.csv"
)
climate_manifest_path <- file.path(
  climate_input_dir,
  "cmip6_four_model_country_year_climate_manifest.csv"
)

final_lag_settings_path <- file.path(
  revision_root,
  "data",
  "source_data",
  "lag_selection",
  "projection_simplified_model_c",
  "01_csv",
  "projection_simplified_final_lag_settings.csv"
)

variance_decomposition_path <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "projection_preparation",
  "01_variance_decomposition",
  "01_tables",
  "projection_variance_decomposition_variance_decomposition_by_bacteria.csv"
)

lag_uncertainty_results_path <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "projection_preparation",
  "02_simplified_model_lag_selection_S18_S19_S20",
  "01_full_results",
  "projection_lag_validation_all_simplified_lag_results.csv"
)

projection_version_tag <- "bootstrap_constraints"

results_root <- file.path(
  revision_root,
  "outputs/historical_associations",
  paste0("10_Figure3_FutureProjection_Simplified_", projection_version_tag)
)
source_data_root <- file.path(
  revision_root,
  "data/source_data",
  paste0("Figure3_ModelC_Projection_", projection_version_tag)
)

dir.create(file.path(results_root, "01_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "02_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "03_diagnostics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "04_workbook"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "05_metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_data_root, "01_csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_data_root, "02_workbook"), recursive = TRUE, showWarnings = FALSE)

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = "3GCREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = "3GCRKP_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRAB", title = "CR-Ab", file_name = "CRAB_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CREC", title = "CR-Ec", file_name = "CREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRKP", title = "CR-Kp", file_name = "CRKP_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRPA", title = "CR-Pa", file_name = "CRPA_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv")
)

bacteria_levels <- vapply(bacteria_specs, `[[`, character(1), "title")
climate_zone_levels <- c("Polar Zone", "Temperate Zone", "Tropical Zone")
heatmap_zone_levels <- c("Tropical Zone", "Temperate Zone", "Polar Zone")
scenario_levels <- c("ssp126", "ssp245", "ssp370", "ssp585")
scenario_labels <- c(
  "ssp126" = "SSP1-2.6",
  "ssp245" = "SSP2-4.5",
  "ssp370" = "SSP3-7.0",
  "ssp585" = "SSP5-8.5"
)
scenario_full_labels <- c(
  "ssp126" = "SSP1-2.6 (Sustainability)",
  "ssp245" = "SSP2-4.5 (Middle of the Road)",
  "ssp370" = "SSP3-7.0 (Regional Rivalry)",
  "ssp585" = "SSP5-8.5 (Fossil-fueled Development)"
)
scenario_colors <- c(
  "ssp126" = "#2C7BB6",
  "ssp245" = "#7AAED4",
  "ssp370" = "#F77F4F",
  "ssp585" = "#D73027"
)
period_colors <- c(
  "2030s" = "#D9C5B4",
  "2050s" = "#84A6A5",
  "2090s" = "#2F5E6F"
)
climate_models_expected <- c("GFDL-ESM4", "MIROC6", "MRI-ESM2-0", "NorESM2-LM")
assessment_periods <- list(
  "2030s" = 2030:2039,
  "2050s" = 2050:2059,
  "2090s" = 2090:2099
)

climate_specs <- list(
  list(raw = "TMP", lag = "TMP_scaled_lag", factor = "Temperature", lower = NA_real_, upper = NA_real_),
  list(raw = "PREC", lag = "PREC_scaled_lag", factor = "Precipitation", lower = 0, upper = NA_real_),
  list(raw = "HUM", lag = "HUM_scaled_lag", factor = "Humidity", lower = 0, upper = 100),
  list(raw = "WET", lag = "WET_scaled_lag", factor = "Wet Days", lower = 0, upper = NA_real_)
)

baseline_years <- 2010:2019
projection_years <- 2020:2100

model_settings <- list(
  use_development_modifier = FALSE,
  use_antibiotic_modifier = FALSE,
  use_ssp_diffusion_modifier = TRUE,
  memory_factors = c("ssp126" = 0.75, "ssp245" = 0.72, "ssp370" = 0.70, "ssp585" = 0.68),
  yearly_variation_factors = c("ssp126" = 1.0, "ssp245" = 1.1, "ssp370" = 1.2, "ssp585" = 1.3),
  scenario_uncertainty_factors = c("ssp126" = 1.0, "ssp245" = 1.0, "ssp370" = 1.0, "ssp585" = 1.0),
  antibiotic_factors = c("ssp126" = 0.85, "ssp245" = 1.0, "ssp370" = 1.15, "ssp585" = 1.30),
  scenario_development_params = list(
    "ssp126" = list(socioeconomic_factor = 0.60),
    "ssp245" = list(socioeconomic_factor = 0.80),
    "ssp370" = list(socioeconomic_factor = 1.20),
    "ssp585" = list(socioeconomic_factor = 1.50)
  ),
  development_modifier = list(
    min_strength = 0.20,
    max_strength = 0.50,
    offset = 0.25,
    pls_scale = 2.20,
    max_modifier = 1.35,
    min_modifier = 0.80
  ),
  crkp_controls = list(
    max_effect_limit = 0.45,
    variability_reduction = 0.25,
    max_growth_rate = 0.014,
    min_memory_factor = 0.76,
    enhanced_smoothing = 5,
    end_correction = TRUE
  ),
  time_series_correlation = list(
    default_persistence = 0.55,
    default_innovation_sd = 0.018
  ),
  empirical_recursive_uncertainty = list(
    enabled = TRUE,
    residual_floor_pct = 0.5,
    min_series_years = 8L,
    persistence_bounds = c(0.05, 0.85),
    innovation_sd_bounds = c(0.006, 0.050),
    initial_sd_bounds = c(0.004, 0.030),
    global_weight = 0.65,
    zone_weight = 0.35
  ),
  response_uncertainty = list(
    include_parameter_uncertainty = TRUE,
    temporal_correlation = 0.85,
    inflate_by_scenario = FALSE,
    inflate_by_time = FALSE
  ),
  response_or_bounds = c(0.50, 2.00),
  log_se_regularization = list(
    enabled = TRUE,
    edge_fraction = 0.08,
    anchor_window = 11L,
    edge_cap_multiplier = 1.20
  ),
  outer_loop_uncertainty = list(
    include_weight_perturbation = FALSE,
    include_memory_perturbation = FALSE,
    weight_bounds = c(0.80, 1.20),
    memory_bounds = c(0.90, 1.10)
  ),
  climate_model_uncertainty = list(
    bootstrap_resample_models = TRUE,
    export_central_model_spread = TRUE
  ),
  lag_uncertainty = list(
    delta_aic_threshold = 2.0,
    min_candidates = 1L,
    fallback_top_n = 2L,
    max_candidates = 4L
  ),
  lookup_extrapolation = list(
    method = "linear_tail",
    slope_fraction = 0.05,
    min_points = 5L
  ),
  temperature_tail_constraint = list(
    enabled = TRUE,
    start_quantile = 0.50,
    stabilize_tail_log_se = TRUE,
    tail_se_window = 11L,
    directions = c(
      "CR-Ab" = "increasing",
      "3GCR-Kp" = "increasing",
      "CR-Ec" = "decreasing"
    )
  ),
  ssp_diffusion_modifier = list(
    "ssp126" = list(net_annual_rate = -0.0006, max_cumulative = 0.94, ramp_up_years = 20),
    "ssp245" = list(net_annual_rate =  0.0003, max_cumulative = 1.08, ramp_up_years = 15),
    "ssp370" = list(net_annual_rate =  0.0007, max_cumulative = 1.16, ramp_up_years = 12),
    "ssp585" = list(net_annual_rate =  0.0012, max_cumulative = 1.25, ramp_up_years = 10)
  ),
  historical_support_quantiles = c(0.01, 0.99),
  historical_effect_quantiles = NULL,
  plot_settings = list(
    fig_width = 12,
    fig_height = 8,
    fig_dpi = 320,
    errbar_width = 0.3,
    errbar_linewidth = 0.7,
    errbar_alpha = 0.9,
    confidence_interval_alpha = 0.14,
    line_width = 1.2,
    line_alpha = 1.0,
    border_size = 0.5,
    axis_line_width = 0.5,
    grid_color = "gray92",
    grid_size = 0.2,
    title_size = 14,
    axis_title_size = 12,
    axis_text_size = 10,
    strip_text_size = 12,
    legend_title_size = 11,
    legend_text_size = 10,
    annotation_text_size = 3.4,
    end_label_nudge = 3,
    label_size = 3.0
  )
)

monte_carlo_settings <- list(
  n_simulations = 1000L,
  confidence_level = 0.95,
  seed = 12345L
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

clamp <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

weighted_mean_safe <- function(x, w) {
  valid <- !is.na(x) & !is.na(w) & is.finite(x) & is.finite(w) & w > 0
  if (!any(valid)) {
    return(mean(x, na.rm = TRUE))
  }
  sum(x[valid] * w[valid], na.rm = TRUE) / sum(w[valid], na.rm = TRUE)
}

fill_series <- function(years, values) {
  valid <- !is.na(values)
  if (sum(valid) == 0) {
    return(rep(NA_real_, length(values)))
  }
  if (sum(valid) == 1) {
    return(rep(values[valid][1], length(values)))
  }
  approx(x = years[valid], y = values[valid], xout = years, rule = 2)$y
}

safe_annual_aggregate <- function(x, mode = c("mean", "sum")) {
  mode <- match.arg(mode)
  if (all(is.na(x))) {
    return(NA_real_)
  }
  if (identical(mode, "sum")) {
    return(sum(x, na.rm = TRUE))
  }
  mean(x, na.rm = TRUE)
}

create_seed_from_string <- function(text_string) {
  hash_value <- digest(text_string, algo = "crc32")
  digits_only <- gsub("[^0-9]", "", hash_value)
  if (nchar(digits_only) == 0) {
    return(42L)
  }
  if (nchar(digits_only) > 9) {
    digits_only <- substr(digits_only, 1, 9)
  }
  seed_value <- suppressWarnings(as.integer(digits_only) %% 2147483647L)
  if (is.na(seed_value)) {
    return(42L)
  }
  seed_value
}

sanitize_file_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

significance_adjustment <- function(p_value) {
  case_when(
    is.na(p_value) ~ 0.1,
    p_value < 0.01 ~ 1.0,
    p_value < 0.05 ~ 0.7,
    p_value < 0.10 ~ 0.3,
    TRUE ~ 0.1
  )
}

bound_and_normalize_weights <- function(weights, min_w = 0.05, max_w = 0.70) {
  if (sum(weights, na.rm = TRUE) <= 0) {
    out <- rep(0.25, length(weights))
    names(out) <- names(weights)
    return(out)
  }

  w <- weights / sum(weights, na.rm = TRUE)
  names(w) <- names(weights)

  for (iter in seq_len(20)) {
    lower_mask <- w < min_w
    upper_mask <- w > max_w
    if (!any(lower_mask | upper_mask)) {
      break
    }

    fixed <- rep(FALSE, length(w))
    fixed[lower_mask | upper_mask] <- TRUE
    w[lower_mask] <- min_w
    w[upper_mask] <- max_w

    remaining_total <- 1 - sum(w[fixed])
    if (remaining_total <= 0) {
      w <- w / sum(w)
      break
    }

    free_idx <- which(!fixed)
    if (length(free_idx) == 0) {
      w <- w / sum(w)
      break
    }

    raw_free <- weights[free_idx]
    if (sum(raw_free, na.rm = TRUE) <= 0) {
      w[free_idx] <- remaining_total / length(free_idx)
    } else {
      w[free_idx] <- raw_free / sum(raw_free, na.rm = TRUE) * remaining_total
    }
  }

  w / sum(w, na.rm = TRUE)
}

load_pls_modifier_strengths <- function(path) {
  if (!file.exists(path)) {
    warning("Missing variance decomposition file at ", path, ". Using default development modifier strengths.")
    return(tibble(
      Bacteria = bacteria_levels,
      PLS_Component_pct_of_total = rep(6, length(bacteria_levels)),
      development_modifier_strength = rep(0.35, length(bacteria_levels))
    ))
  }

  read.csv(path, stringsAsFactors = FALSE) %>%
    transmute(
      Bacteria,
      PLS_Component_pct_of_total = as.numeric(PLS_Component_pct_of_total),
      pls_share = PLS_Component_pct_of_total / 100,
      development_modifier_strength = clamp(
        model_settings$development_modifier$offset +
          model_settings$development_modifier$pls_scale * pls_share,
        model_settings$development_modifier$min_strength,
        model_settings$development_modifier$max_strength
      )
    )
}

load_lag_uncertainty_catalog <- function(path) {
  if (!file.exists(path)) {
    stop("Missing lag uncertainty results file: ", path, call. = FALSE)
  }

  lag_df <- read.csv(path, stringsAsFactors = FALSE) %>%
    rename(
      TMP_lag = temp_lag,
      PREC_lag = prec_lag,
      HUM_lag = hum_lag,
      WET_lag = wet_lag
    ) %>%
    mutate(
      delta_aic = AIC - ave(AIC, Bacteria, FUN = min),
      akaike_weight_full = exp(-0.5 * delta_aic)
    ) %>%
    group_by(Bacteria) %>%
    mutate(akaike_weight_full = akaike_weight_full / sum(akaike_weight_full, na.rm = TRUE)) %>%
    ungroup()

  candidate_df <- lag_df %>%
    group_by(Bacteria) %>%
    arrange(AIC, .by_group = TRUE) %>%
    mutate(
      use_candidate = delta_aic <= model_settings$lag_uncertainty$delta_aic_threshold,
      fallback_required = sum(use_candidate, na.rm = TRUE) < model_settings$lag_uncertainty$min_candidates,
      use_candidate = if (first(fallback_required)) Rank <= model_settings$lag_uncertainty$fallback_top_n else use_candidate
    ) %>%
    filter(use_candidate) %>%
    slice_head(n = model_settings$lag_uncertainty$max_candidates) %>%
    mutate(
      lag_sampling_weight = akaike_weight_full / sum(akaike_weight_full, na.rm = TRUE),
      is_primary = Rank == min(Rank, na.rm = TRUE)
    ) %>%
    ungroup()

  candidate_df
}

baseline_column_for_factor <- function(climate_factor_name) {
  switch(
    climate_factor_name,
    "Temperature" = "baseline_temp",
    "Precipitation" = "baseline_precip",
    "Humidity" = "baseline_humidity",
    "Wet Days" = "baseline_wetdays",
    stop("Unsupported climate factor: ", climate_factor_name, call. = FALSE)
  )
}

compute_log_se_from_ci <- function(lower_ci, upper_ci) {
  lower_ci <- pmax(lower_ci, 1e-8)
  upper_ci <- pmax(upper_ci, lower_ci + 1e-8)
  pmax((log(upper_ci) - log(lower_ci)) / (2 * 1.96), 0)
}

apply_optional_bounds <- function(x, bounds) {
  if (is.null(bounds) || length(bounds) != 2 || any(!is.finite(bounds))) {
    return(x)
  }
  clamp(x, bounds[1], bounds[2])
}

current_uncertainty_mode <- function() {
  if (isTRUE(model_settings$climate_model_uncertainty$bootstrap_resample_models)) {
    return("4_gcm_spread_empirical_recursive_residual_calibration_response_parameter_uncertainty")
  }
  "empirical_recursive_residual_calibration_response_parameter_uncertainty_gcm_spread_exported_separately"
}

temperature_tail_constraint_direction <- function(bacteria_name) {
  if (!isTRUE(model_settings$temperature_tail_constraint$enabled)) {
    return("none")
  }
  directions <- model_settings$temperature_tail_constraint$directions
  if (!bacteria_name %in% names(directions)) {
    return("none")
  }
  unname(directions[[bacteria_name]]) %||% "none"
}

compute_relative_residual <- function(observed_pct, fitted_pct, floor_pct) {
  denominator <- pmax(fitted_pct, floor_pct, 0.1)
  clamp((observed_pct - fitted_pct) / denominator, -0.95, 0.95)
}

estimate_ar1_residual_process <- function(series_values) {
  settings <- model_settings$empirical_recursive_uncertainty
  values <- series_values[is.finite(series_values)]
  n_values <- length(values)

  default_persistence <- model_settings$time_series_correlation$default_persistence
  default_innovation_sd <- model_settings$time_series_correlation$default_innovation_sd

  if (n_values < (settings$min_series_years %||% 8L)) {
    return(list(
      persistence = default_persistence,
      innovation_sd = default_innovation_sd,
      series_sd = default_innovation_sd,
      n_years = n_values,
      calibrated = FALSE
    ))
  }

  lagged <- values[-n_values]
  leaded <- values[-1]
  valid <- is.finite(lagged) & is.finite(leaded)

  if (sum(valid) < max(4L, (settings$min_series_years %||% 8L) - 1L)) {
    return(list(
      persistence = default_persistence,
      innovation_sd = default_innovation_sd,
      series_sd = stats::sd(values, na.rm = TRUE) %||% default_innovation_sd,
      n_years = n_values,
      calibrated = FALSE
    ))
  }

  rho_raw <- suppressWarnings(stats::cor(lagged[valid], leaded[valid], use = "complete.obs"))
  rho <- if (is.finite(rho_raw)) {
    clamp(rho_raw, settings$persistence_bounds[1], settings$persistence_bounds[2])
  } else {
    default_persistence
  }

  innovation_values <- leaded[valid] - rho * lagged[valid]
  innovation_sd <- stats::sd(innovation_values, na.rm = TRUE)
  series_sd <- stats::sd(values, na.rm = TRUE)

  if (!is.finite(innovation_sd) || innovation_sd <= 0) {
    innovation_sd <- series_sd * sqrt(pmax(1 - rho^2, 0.1))
  }
  if (!is.finite(series_sd) || series_sd <= 0) {
    series_sd <- innovation_sd
  }

  innovation_sd <- clamp(innovation_sd, settings$innovation_sd_bounds[1], settings$innovation_sd_bounds[2])
  series_sd <- clamp(series_sd, settings$initial_sd_bounds[1], settings$initial_sd_bounds[2])

  list(
    persistence = rho,
    innovation_sd = innovation_sd,
    series_sd = series_sd,
    n_years = n_values,
    calibrated = TRUE
  )
}

combine_empirical_metric <- function(global_value, zone_value, default_value) {
  settings <- model_settings$empirical_recursive_uncertainty
  parts <- c(global_value, zone_value)
  weights <- c(settings$global_weight %||% 0.65, settings$zone_weight %||% 0.35)
  valid <- is.finite(parts)
  if (!any(valid)) {
    return(default_value)
  }
  sum(parts[valid] * weights[valid], na.rm = TRUE) / sum(weights[valid], na.rm = TRUE)
}

build_recursive_uncertainty_calibration <- function(historical_global_df, historical_zone_df, baseline_global_df, baseline_zone_df) {
  settings <- model_settings$empirical_recursive_uncertainty

  global_residual_df <- historical_global_df %>%
    left_join(
      baseline_global_df %>%
        select(Bacteria, baseline_resistance_pct),
      by = "Bacteria"
    ) %>%
    mutate(
      residual_floor_pct = pmax(settings$residual_floor_pct %||% 0.5, baseline_resistance_pct * 0.05),
      raw_relative_residual = compute_relative_residual(observed_R_pct, fitted_R_pct, residual_floor_pct)
    ) %>%
    group_by(Bacteria) %>%
    mutate(
      relative_residual = clamp(raw_relative_residual - mean(raw_relative_residual, na.rm = TRUE), -0.95, 0.95)
    ) %>%
    ungroup()

  zone_residual_df <- historical_zone_df %>%
    left_join(
      baseline_zone_df %>%
        select(Bacteria, climate_zone, baseline_resistance_pct),
      by = c("Bacteria", "climate_zone")
    ) %>%
    mutate(
      residual_floor_pct = pmax(settings$residual_floor_pct %||% 0.5, baseline_resistance_pct * 0.05),
      raw_relative_residual = compute_relative_residual(observed_R_pct, fitted_R_pct, residual_floor_pct)
    ) %>%
    group_by(Bacteria, climate_zone) %>%
    mutate(
      relative_residual = clamp(raw_relative_residual - mean(raw_relative_residual, na.rm = TRUE), -0.95, 0.95)
    ) %>%
    ungroup()

  zone_process_summary <- zone_residual_df %>%
    group_by(Bacteria, climate_zone) %>%
    arrange(year, .by_group = TRUE) %>%
    group_modify(~ {
      est <- estimate_ar1_residual_process(.x$relative_residual)
      tibble(
        zone_persistence = est$persistence,
        zone_innovation_sd = est$innovation_sd,
        zone_series_sd = est$series_sd,
        zone_n_years = est$n_years,
        zone_calibrated = est$calibrated
      )
    }) %>%
    ungroup()

  global_process_summary <- global_residual_df %>%
    group_by(Bacteria) %>%
    arrange(year, .by_group = TRUE) %>%
    group_modify(~ {
      est <- estimate_ar1_residual_process(.x$relative_residual)
      tibble(
        global_persistence = est$persistence,
        global_innovation_sd = est$innovation_sd,
        global_series_sd = est$series_sd,
        global_n_years = est$n_years,
        global_calibrated = est$calibrated
      )
    }) %>%
    ungroup()

  zone_combined_summary <- zone_process_summary %>%
    group_by(Bacteria) %>%
    summarise(
      zone_persistence = weighted_mean_safe(zone_persistence, zone_n_years),
      zone_innovation_sd = weighted_mean_safe(zone_innovation_sd, zone_n_years),
      zone_series_sd = weighted_mean_safe(zone_series_sd, zone_n_years),
      zone_series_count = sum(zone_calibrated, na.rm = TRUE),
      .groups = "drop"
    )

  calibration_df <- tibble(Bacteria = bacteria_levels) %>%
    left_join(global_process_summary, by = "Bacteria") %>%
    left_join(zone_combined_summary, by = "Bacteria") %>%
    mutate(
      calibrated_persistence = purrr::pmap_dbl(
        list(global_persistence, zone_persistence),
        ~ combine_empirical_metric(..1, ..2, model_settings$time_series_correlation$default_persistence)
      ),
      calibrated_innovation_sd = purrr::pmap_dbl(
        list(global_innovation_sd, zone_innovation_sd),
        ~ combine_empirical_metric(..1, ..2, model_settings$time_series_correlation$default_innovation_sd)
      ),
      calibrated_initial_sd = purrr::pmap_dbl(
        list(global_series_sd, zone_series_sd),
        ~ combine_empirical_metric(..1, ..2, model_settings$time_series_correlation$default_innovation_sd)
      ),
      calibrated_persistence = clamp(calibrated_persistence, settings$persistence_bounds[1], settings$persistence_bounds[2]),
      calibrated_innovation_sd = clamp(calibrated_innovation_sd, settings$innovation_sd_bounds[1], settings$innovation_sd_bounds[2]),
      calibrated_initial_sd = clamp(calibrated_initial_sd, settings$initial_sd_bounds[1], settings$initial_sd_bounds[2]),
      calibration_source = case_when(
        isTRUE(global_calibrated) & zone_series_count > 0 ~ "global_and_zone_historical_residuals",
        isTRUE(global_calibrated) ~ "global_historical_residuals",
        zone_series_count > 0 ~ "zone_historical_residuals",
        TRUE ~ "default_fallback"
      )
    )

  median_sd <- stats::median(calibration_df$calibrated_innovation_sd, na.rm = TRUE)
  if (!is.finite(median_sd) || median_sd <= 0) {
    median_sd <- model_settings$time_series_correlation$default_innovation_sd
  }

  calibration_df <- calibration_df %>%
    mutate(
      empirical_variation_factor = clamp(calibrated_innovation_sd / median_sd, 0.75, 1.30)
    )

  list(
    calibration_df = calibration_df,
    historical_global_residual_df = global_residual_df,
    historical_zone_residual_df = zone_residual_df,
    zone_summary_df = zone_process_summary
  )
}

stabilize_lookup_edge_log_se <- function(log_se) {
  settings <- model_settings$log_se_regularization
  if (!isTRUE(settings$enabled)) {
    return(list(log_se = pmax(log_se, 0), stabilized = rep(FALSE, length(log_se))))
  }

  values <- pmax(log_se, 0)
  n_values <- length(values)
  if (n_values < 8) {
    return(list(log_se = values, stabilized = rep(FALSE, n_values)))
  }

  edge_n <- max(3L, floor(n_values * (settings$edge_fraction %||% 0.08)))
  edge_n <- min(edge_n, floor((n_values - 2L) / 2L))
  if (edge_n < 1) {
    return(list(log_se = values, stabilized = rep(FALSE, n_values)))
  }

  anchor_window <- max(3L, settings$anchor_window %||% 11L)
  lower_anchor_idx <- seq.int(edge_n + 1L, min(n_values, edge_n + anchor_window))
  upper_anchor_idx <- seq.int(max(1L, n_values - edge_n - anchor_window + 1L), max(1L, n_values - edge_n))

  lower_anchor <- stats::median(values[lower_anchor_idx], na.rm = TRUE)
  upper_anchor <- stats::median(values[upper_anchor_idx], na.rm = TRUE)
  if (!is.finite(lower_anchor)) {
    lower_anchor <- values[min(n_values, edge_n + 1L)]
  }
  if (!is.finite(upper_anchor)) {
    upper_anchor <- values[max(1L, n_values - edge_n)]
  }

  out <- values
  lower_idx <- seq_len(edge_n)
  upper_idx <- seq.int(n_values - edge_n + 1L, n_values)
  lower_cap <- lower_anchor * (settings$edge_cap_multiplier %||% 1.20)
  upper_cap <- upper_anchor * (settings$edge_cap_multiplier %||% 1.20)
  out[lower_idx] <- pmin(out[lower_idx], lower_cap)
  out[upper_idx] <- pmin(out[upper_idx], upper_cap)

  list(
    log_se = pmax(out, 0),
    stabilized = out < (values - 1e-12)
  )
}

apply_temperature_upper_tail_constraint <- function(effect_or, relative_se, x_raw, reference_raw, bacteria_name, climate_factor_name) {
  if (!identical(climate_factor_name, "Temperature")) {
    return(list(
      effect_or = effect_or,
      log_se = relative_se,
      lower_ci = effect_or * exp(-1.96 * relative_se),
      upper_ci = effect_or * exp(1.96 * relative_se),
      direction = "none",
      applied = FALSE,
      tail_start_raw = NA_real_
    ))
  }

  direction <- temperature_tail_constraint_direction(bacteria_name)
  if (identical(direction, "none") || length(effect_or) == 0) {
    return(list(
      effect_or = effect_or,
      log_se = relative_se,
      lower_ci = effect_or * exp(-1.96 * relative_se),
      upper_ci = effect_or * exp(1.96 * relative_se),
      direction = direction,
      applied = FALSE,
      tail_start_raw = NA_real_
    ))
  }

  tail_start_raw <- max(
    reference_raw,
    as.numeric(stats::quantile(x_raw, probs = model_settings$temperature_tail_constraint$start_quantile, na.rm = TRUE, names = FALSE))
  )
  tail_idx <- which(x_raw >= tail_start_raw)
  if (length(tail_idx) == 0) {
    return(list(
      effect_or = effect_or,
      log_se = relative_se,
      lower_ci = effect_or * exp(-1.96 * relative_se),
      upper_ci = effect_or * exp(1.96 * relative_se),
      direction = direction,
      applied = FALSE,
      tail_start_raw = tail_start_raw
    ))
  }

  constrained <- effect_or
  constrained_log_se <- relative_se
  applied <- FALSE

  if (identical(direction, "increasing")) {
    running_value <- constrained[min(tail_idx)]
    for (j in tail_idx) {
      running_value <- max(running_value, constrained[j])
      if (!isTRUE(all.equal(constrained[j], running_value))) {
        applied <- TRUE
      }
      constrained[j] <- running_value
    }
  } else if (identical(direction, "decreasing")) {
    running_value <- constrained[min(tail_idx)]
    for (j in tail_idx) {
      running_value <- min(running_value, constrained[j])
      if (!isTRUE(all.equal(constrained[j], running_value))) {
        applied <- TRUE
      }
      constrained[j] <- running_value
    }
  }

  if (applied && isTRUE(model_settings$temperature_tail_constraint$stabilize_tail_log_se)) {
    anchor_idx <- min(tail_idx)
    half_window <- max(0L, floor((model_settings$temperature_tail_constraint$tail_se_window %||% 11L) / 2))
    window_start <- max(1L, anchor_idx - half_window)
    window_end <- min(length(relative_se), anchor_idx + half_window)
    anchor_se <- stats::median(relative_se[window_start:window_end], na.rm = TRUE)
    if (!is.finite(anchor_se)) {
      anchor_se <- relative_se[anchor_idx]
    }
    constrained_log_se[tail_idx] <- pmin(relative_se[tail_idx], anchor_se)
  }

  list(
    effect_or = constrained,
    log_se = constrained_log_se,
    lower_ci = constrained * exp(-1.96 * constrained_log_se),
    upper_ci = constrained * exp(1.96 * constrained_log_se),
    direction = direction,
    applied = applied,
    tail_start_raw = tail_start_raw
  )
}

compute_ssp_diffusion_multiplier <- function(scenario_name, years_passed) {
  if (!isTRUE(model_settings$use_ssp_diffusion_modifier)) {
    return(1.0)
  }

  params <- model_settings$ssp_diffusion_modifier[[scenario_name]]
  if (is.null(params)) {
    return(1.0)
  }

  years_passed <- pmax(years_passed, 0)
  ramp_up_years <- max(1, params$ramp_up_years %||% 10)
  effective_years <- years_passed * pmin(1, years_passed / ramp_up_years)
  modifier <- 1 + (params$net_annual_rate %||% 0) * effective_years

  lower_bound <- min(1, params$max_cumulative %||% 1)
  upper_bound <- max(1, params$max_cumulative %||% 1)
  clamp(modifier, lower_bound, upper_bound)
}

interpolate_with_linear_tail <- function(x, y, xout, slope_fraction, min_points) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  inside <- xout >= x_min & xout <= x_max

  approx_vals <- approx(x = x, y = y, xout = clamp(xout, x_min, x_max), rule = 2)$y
  out <- approx_vals

  n <- length(x)
  tail_n <- max(min_points, round(n * slope_fraction))
  tail_n <- min(tail_n, n)

  if (tail_n >= 2) {
    lower_df <- tibble(x = x[seq_len(tail_n)], y = y[seq_len(tail_n)])
    upper_df <- tibble(x = x[(n - tail_n + 1):n], y = y[(n - tail_n + 1):n])

    lower_fit <- stats::lm(y ~ x, data = lower_df)
    upper_fit <- stats::lm(y ~ x, data = upper_df)
    lower_slope <- stats::coef(lower_fit)[["x"]] %||% 0
    upper_slope <- stats::coef(upper_fit)[["x"]] %||% 0

    lower_boundary <- y[1]
    upper_boundary <- y[n]

    lower_mask <- !inside & xout < x_min
    upper_mask <- !inside & xout > x_max

    out[lower_mask] <- lower_boundary + lower_slope * (xout[lower_mask] - x_min)
    out[upper_mask] <- upper_boundary + upper_slope * (xout[upper_mask] - x_max)
  }

  out
}

generate_correlated_normal <- function(n, rho, seed_prefix) {
  set.seed(create_seed_from_string(seed_prefix))
  innovations <- rnorm(n)
  out <- numeric(n)
  out[1] <- innovations[1]
  if (n > 1) {
    for (i in 2:n) {
      out[i] <- rho * out[i - 1] + sqrt(1 - rho^2) * innovations[i]
    }
  }
  out
}

compute_development_modifier <- function(bacteria_name, scenario_name, years, modifier_strength) {
  if (!isTRUE(model_settings$use_development_modifier)) {
    return(rep(1, length(years)))
  }

  time_factor <- pmax(0, pmin(1, (years - min(projection_years)) / (max(projection_years) - min(projection_years))))
  scenario_factor <- model_settings$scenario_development_params[[scenario_name]]$socioeconomic_factor %||% 1.0
  antibiotic_factor <- if (isTRUE(model_settings$use_antibiotic_modifier)) {
    model_settings$antibiotic_factors[[scenario_name]] %||% 1.0
  } else {
    1.0
  }

  modifier <- 1 + (scenario_factor - 1) * time_factor * antibiotic_factor * modifier_strength
  clamp(
    modifier,
    model_settings$development_modifier$min_modifier,
    model_settings$development_modifier$max_modifier
  )
}

combine_weighted_effect_log <- function(temp_effect, humid_effect, precip_effect, wetdays_effect, weights_vec, development_modifier = NULL) {
  log_effect <- rep(0, length(temp_effect))
  if (!is.na(weights_vec["Temperature"]) && weights_vec["Temperature"] > 0) {
    log_effect <- log_effect + weights_vec["Temperature"] * log(pmax(temp_effect, 1e-8))
  }
  if (!is.na(weights_vec["Humidity"]) && weights_vec["Humidity"] > 0) {
    log_effect <- log_effect + weights_vec["Humidity"] * log(pmax(humid_effect, 1e-8))
  }
  if (!is.na(weights_vec["Precipitation"]) && weights_vec["Precipitation"] > 0) {
    log_effect <- log_effect + weights_vec["Precipitation"] * log(pmax(precip_effect, 1e-8))
  }
  if (!is.na(weights_vec["Wet Days"]) && weights_vec["Wet Days"] > 0) {
    log_effect <- log_effect + weights_vec["Wet Days"] * log(pmax(wetdays_effect, 1e-8))
  }
  if (!is.null(development_modifier)) {
    log_effect <- log_effect + log(pmax(development_modifier, 1e-8))
  }
  exp(log_effect)
}

sample_perturbed_weights <- function(base_weights, seed_prefix) {
  if (!isTRUE(model_settings$outer_loop_uncertainty$include_weight_perturbation)) {
    return(base_weights)
  }

  set.seed(create_seed_from_string(seed_prefix))
  bounds <- model_settings$outer_loop_uncertainty$weight_bounds
  names(base_weights) <- sub("^.*\\.", "", names(base_weights))
  multipliers <- stats::runif(length(base_weights), min = bounds[1], max = bounds[2])
  perturbed <- base_weights * multipliers
  names(perturbed) <- names(base_weights)
  bound_and_normalize_weights(perturbed)
}

sample_memory_multiplier <- function(seed_prefix) {
  if (!isTRUE(model_settings$outer_loop_uncertainty$include_memory_perturbation)) {
    return(1.0)
  }

  set.seed(create_seed_from_string(seed_prefix))
  bounds <- model_settings$outer_loop_uncertainty$memory_bounds
  stats::runif(1, min = bounds[1], max = bounds[2])
}

parse_cmip6_file_info <- function(file_path) {
  file_name <- basename(file_path)
  model_hits <- Filter(function(x) grepl(x, file_name, fixed = TRUE), climate_models_expected)
  scenario_hits <- Filter(function(x) grepl(x, file_name, fixed = TRUE), scenario_levels)
  if (length(model_hits) != 1 || length(scenario_hits) != 1) {
    return(NULL)
  }
  list(model = model_hits[[1]], ssp = scenario_hits[[1]])
}

read_country_metadata <- function(path) {
  read.csv(path, stringsAsFactors = FALSE) %>%
    select(NAME, Country_Code, Region, lon, lat) %>%
    distinct() %>%
    group_by(NAME) %>%
    summarise(
      Country_Code = Country_Code[which.max(!is.na(Country_Code))],
      Region = Region[which.max(!is.na(Region))],
      lon = mean(lon, na.rm = TRUE),
      lat = mean(lat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(NAME), NAME != "")
}

read_monthly_climate_variable_from_xlsx <- function(dir_path, value_name, annual_mode, target_countries) {
  files <- list.files(dir_path, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No CMIP6 xlsx files found in ", dir_path, call. = FALSE)
  }

  purrr::map_dfr(files, function(file_path) {
    file_info <- parse_cmip6_file_info(file_path)
    if (is.null(file_info)) {
      return(tibble())
    }

    raw_df <- readxl::read_xlsx(file_path)
    if (!"time" %in% names(raw_df)) {
      return(tibble())
    }

    raw_df %>%
      mutate(year = as.integer(substr(as.character(time), 1, 4))) %>%
      select(-time) %>%
      pivot_longer(cols = -year, names_to = "NAME", values_to = value_name) %>%
      filter(NAME %in% target_countries, year >= 2015, year <= 2100) %>%
      mutate(model = file_info$model, ssp = file_info$ssp) %>%
      group_by(NAME, model, ssp, year) %>%
      summarise(!!value_name := safe_annual_aggregate(.data[[value_name]], annual_mode), .groups = "drop")
  })
}

build_four_model_future_climate_input <- function(output_path, manifest_path) {
  metadata <- read_country_metadata(country_metadata_path)
  target_countries <- metadata$NAME

  message("Rebuilding four-model CMIP6 climate input from raw files...")

  tas_annual <- read_monthly_climate_variable_from_xlsx(
    dir_path = file.path(legacy_cmip6_root, "tas"),
    value_name = "tas_annual_mean",
    annual_mode = "mean",
    target_countries = target_countries
  )
  hurs_annual <- read_monthly_climate_variable_from_xlsx(
    dir_path = file.path(legacy_cmip6_root, "hurs"),
    value_name = "hurs_annual_mean",
    annual_mode = "mean",
    target_countries = target_countries
  )
  pr_annual <- read_monthly_climate_variable_from_xlsx(
    dir_path = file.path(legacy_cmip6_root, "pr"),
    value_name = "pr_annual_total",
    annual_mode = "sum",
    target_countries = target_countries
  )

  wetdays_annual <- read.csv(wetdays_yearly_path, stringsAsFactors = FALSE) %>%
    transmute(
      NAME = country,
      model = as.character(model),
      ssp = tolower(as.character(scenario)),
      year = as.integer(year),
      wet_days_yearly0.1 = as.numeric(wet_days_yearly)
    ) %>%
    filter(
      NAME %in% target_countries,
      model %in% climate_models_expected,
      ssp %in% scenario_levels,
      year >= 2015,
      year <= 2100
    )

  combined <- tas_annual %>%
    inner_join(hurs_annual, by = c("NAME", "model", "ssp", "year")) %>%
    inner_join(pr_annual, by = c("NAME", "model", "ssp", "year")) %>%
    inner_join(wetdays_annual, by = c("NAME", "model", "ssp", "year")) %>%
    left_join(metadata, by = "NAME") %>%
    mutate(pm25_annual_mean = NA_real_) %>%
    select(
      NAME, Country_Code, Region, model, ssp, year,
      hurs_annual_mean, pm25_annual_mean, pr_annual_total,
      tas_annual_mean, wet_days_yearly0.1, lon, lat
    ) %>%
    arrange(NAME, model, ssp, year)

  manifest <- combined %>%
    group_by(model, ssp) %>%
    summarise(
      n_countries = n_distinct(NAME),
      year_min = min(year, na.rm = TRUE),
      year_max = max(year, na.rm = TRUE),
      missing_humidity = sum(is.na(hurs_annual_mean)),
      missing_precipitation = sum(is.na(pr_annual_total)),
      missing_temperature = sum(is.na(tas_annual_mean)),
      missing_wetdays = sum(is.na(wet_days_yearly0.1)),
      .groups = "drop"
    )

  write.csv(combined, output_path, row.names = FALSE)
  write.csv(manifest, manifest_path, row.names = FALSE)

  combined
}

ensure_future_climate_input <- function(output_path, manifest_path) {
  rebuild_needed <- TRUE

  if (file.exists(output_path)) {
    existing <- read.csv(output_path, stringsAsFactors = FALSE)
    required_cols <- c(
      "NAME", "Country_Code", "Region", "model", "ssp", "year",
      "hurs_annual_mean", "pr_annual_total", "tas_annual_mean",
      "wet_days_yearly0.1", "lon", "lat"
    )
    existing_models <- sort(unique(as.character(existing$model)))
    if (
      all(required_cols %in% names(existing)) &&
      identical(existing_models, sort(climate_models_expected)) &&
      min(existing$year, na.rm = TRUE) <= 2015 &&
      max(existing$year, na.rm = TRUE) >= 2100
    ) {
      rebuild_needed <- FALSE
    }
  }

  if (rebuild_needed) {
    return(build_four_model_future_climate_input(output_path, manifest_path))
  }

  read.csv(output_path, stringsAsFactors = FALSE)
}

prepare_projection_model_data <- function(file_path, lag_row) {
  data <- read.csv(file_path, stringsAsFactors = FALSE) %>%
    mutate(
      year = as.numeric(as.character(year)),
      Region = factor(Region),
      climate_zone = case_when(
        abs(lat) > 66.5 ~ "Polar Zone",
        abs(lat) > 23.5 ~ "Temperate Zone",
        TRUE ~ "Tropical Zone"
      ),
      climate_zone = factor(climate_zone, levels = climate_zone_levels),
      observed_R_prop = case_when(
        is.na(R) ~ NA_real_,
        R > 1 ~ R / 100,
        TRUE ~ R
      ),
      observed_R_pct = observed_R_prop * 100,
      HUM = pmin(HUM, 100)
    ) %>%
    group_by(NAME) %>%
    mutate(location_id = cur_group_id()) %>%
    ungroup()

  scaling_df <- tibble(
    Climate_Factor = c("Temperature", "Precipitation", "Humidity", "Wet Days"),
    variable = c("TMP", "PREC", "HUM", "WET"),
    mean_value = c(mean(data$TMP, na.rm = TRUE), mean(data$PREC, na.rm = TRUE), mean(data$HUM, na.rm = TRUE), mean(data$WET, na.rm = TRUE)),
    sd_value = c(sd(data$TMP, na.rm = TRUE), sd(data$PREC, na.rm = TRUE), sd(data$HUM, na.rm = TRUE), sd(data$WET, na.rm = TRUE))
  ) %>%
    mutate(sd_value = ifelse(is.na(sd_value) | sd_value == 0, 1, sd_value))

  scale_lookup <- setNames(
    split(scaling_df[, c("mean_value", "sd_value")], scaling_df$variable),
    scaling_df$variable
  )

  support_df <- purrr::map_dfr(climate_specs, function(spec) {
    x <- data[[spec$raw]]
    q <- quantile(x, probs = model_settings$historical_support_quantiles, na.rm = TRUE)
    tibble(
      Climate_Factor = spec$factor,
      variable = spec$raw,
      support_lower = as.numeric(q[[1]]),
      support_upper = as.numeric(q[[2]]),
      observed_min = min(x, na.rm = TRUE),
      observed_max = max(x, na.rm = TRUE)
    )
  })

  lagged_data <- data %>%
    mutate(
      TMP_scaled = (TMP - scale_lookup$TMP$mean_value) / scale_lookup$TMP$sd_value,
      PREC_scaled = (PREC - scale_lookup$PREC$mean_value) / scale_lookup$PREC$sd_value,
      HUM_scaled = (HUM - scale_lookup$HUM$mean_value) / scale_lookup$HUM$sd_value,
      WET_scaled = (WET - scale_lookup$WET$mean_value) / scale_lookup$WET$sd_value
    ) %>%
    group_by(location_id) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      TMP_scaled_lag = lag(TMP_scaled, lag_row$TMP_lag[1]),
      PREC_scaled_lag = lag(PREC_scaled, lag_row$PREC_lag[1]),
      HUM_scaled_lag = lag(HUM_scaled, lag_row$HUM_lag[1]),
      WET_scaled_lag = lag(WET_scaled, lag_row$WET_lag[1])
    ) %>%
    filter(!is.na(TMP_scaled_lag), !is.na(PREC_scaled_lag), !is.na(HUM_scaled_lag), !is.na(WET_scaled_lag)) %>%
    ungroup()

  list(
    raw_data = data,
    lagged_data = lagged_data,
    scaling_df = scaling_df,
    support_df = support_df
  )
}

build_simplified_model <- function(data) {
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
  model_formula <- as.formula(
    paste(
      "logit_R ~",
      "s(TMP_scaled_lag, k = 5, bs = 'cr') +",
      "s(PREC_scaled_lag, k = 10, bs = 'cr') +",
      "s(HUM_scaled_lag, k = 10, bs = 'cr') +",
      "s(WET_scaled_lag, k = 10, bs = 'cr') +",
      "s(lat, lon, bs = 'sos', k = 20) +",
      "s(year, bs = 'cr', k = 8) +",
      "s(Region, bs = 're') +",
      "climate_zone"
    )
  )

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

extract_model_term_statistics <- function(model, bacteria, lag_row) {
  s_table <- as.data.frame(summary(model)$s.table)
  s_table$Term <- rownames(s_table)
  rownames(s_table) <- NULL

  s_table %>%
    transmute(
      Bacteria = bacteria,
      Term,
      Climate_Factor = recode(
        Term,
        "s(TMP_scaled_lag)" = "Temperature",
        "s(PREC_scaled_lag)" = "Precipitation",
        "s(HUM_scaled_lag)" = "Humidity",
        "s(WET_scaled_lag)" = "Wet Days",
        .default = NA_character_
      ),
      Lag_Years = case_when(
        Climate_Factor == "Temperature" ~ lag_row$TMP_lag[1],
        Climate_Factor == "Precipitation" ~ lag_row$PREC_lag[1],
        Climate_Factor == "Humidity" ~ lag_row$HUM_lag[1],
        Climate_Factor == "Wet Days" ~ lag_row$WET_lag[1],
        TRUE ~ NA_real_
      ),
      Term_Group = ifelse(is.na(Climate_Factor), "Adjustment", "Climate"),
      EDF = edf,
      Ref_DF = Ref.df,
      F_stat = F,
      P_value = `p-value`
    )
}

derive_weights_from_term_stats <- function(term_stats, bacteria) {
  climate_terms <- term_stats %>%
    filter(Bacteria == bacteria, Term_Group == "Climate", !is.na(Climate_Factor)) %>%
    mutate(
      significance_factor = significance_adjustment(P_value),
      edf_complexity_factor = pmax(EDF, 0),
      raw_score = pmax(F_stat, 0) * edf_complexity_factor * significance_factor
    )

  weight_vector <- setNames(climate_terms$raw_score, climate_terms$Climate_Factor)
  normalized_weights <- bound_and_normalize_weights(weight_vector)

  weights_df <- tibble(
    Bacteria = bacteria,
    Temperature = normalized_weights["Temperature"] %||% 0.25,
    Precipitation = normalized_weights["Precipitation"] %||% 0.25,
    Humidity = normalized_weights["Humidity"] %||% 0.25,
    `Wet Days` = normalized_weights["Wet Days"] %||% 0.25
  )

  details_df <- climate_terms %>%
    mutate(
      raw_weight = raw_score / sum(raw_score, na.rm = TRUE),
      bounded_weight = normalized_weights[Climate_Factor]
    ) %>%
    select(Bacteria, Climate_Factor, Lag_Years, EDF, Ref_DF, F_stat, P_value, significance_factor, edf_complexity_factor, raw_score, raw_weight, bounded_weight)

  list(weights = weights_df, details = details_df)
}

make_reference_prediction_row <- function(data) {
  ref_region <- names(sort(table(data$Region), decreasing = TRUE))[1]
  ref_zone <- names(sort(table(data$climate_zone), decreasing = TRUE))[1]

  tibble(
    TMP_scaled_lag = mean(data$TMP_scaled_lag, na.rm = TRUE),
    PREC_scaled_lag = mean(data$PREC_scaled_lag, na.rm = TRUE),
    HUM_scaled_lag = mean(data$HUM_scaled_lag, na.rm = TRUE),
    WET_scaled_lag = mean(data$WET_scaled_lag, na.rm = TRUE),
    lat = mean(data$lat, na.rm = TRUE),
    lon = mean(data$lon, na.rm = TRUE),
    year = mean(data$year, na.rm = TRUE),
    climate_zone = factor(ref_zone, levels = levels(data$climate_zone)),
    Region = factor(ref_region, levels = levels(data$Region))
  )
}

make_response_lookup <- function(model, lagged_data, scaling_df, support_df, baseline_zone_df, bacteria, lag_row, n_points = 401L) {
  reference_row <- make_reference_prediction_row(lagged_data)
  model_coef <- stats::coef(model)
  model_vcov <- stats::vcov(model, unconditional = TRUE)

  lookup_rows <- purrr::map_dfr(climate_specs, function(spec) {
    support_row <- support_df %>% filter(Climate_Factor == spec$factor)
    scale_row <- scaling_df %>% filter(Climate_Factor == spec$factor)
    lower_raw <- support_row$support_lower[1]
    upper_raw <- support_row$support_upper[1]
    lower_scaled <- (lower_raw - scale_row$mean_value[1]) / scale_row$sd_value[1]
    upper_scaled <- (upper_raw - scale_row$mean_value[1]) / scale_row$sd_value[1]

    x_scaled <- seq(lower_scaled, upper_scaled, length.out = n_points)
    pred_data <- reference_row[rep(1, length(x_scaled)), ]
    pred_data[[spec$lag]] <- x_scaled
    pred_lp <- suppressWarnings(predict(model, pred_data, type = "lpmatrix"))
    term_cols <- grep(paste0("^s\\(", spec$lag, "\\)"), colnames(pred_lp))
    if (length(term_cols) == 0) {
      stop("Unable to identify lpmatrix columns for ", spec$lag, call. = FALSE)
    }
    pred_term_lp <- pred_lp[, term_cols, drop = FALSE]
    fit_values <- as.numeric(pred_term_lp %*% model_coef[term_cols])

    x_raw <- x_scaled * scale_row$sd_value[1] + scale_row$mean_value[1]
    if (!is.na(spec$lower)) {
      x_raw <- pmax(x_raw, spec$lower)
    }
    if (!is.na(spec$upper)) {
      x_raw <- pmin(x_raw, spec$upper)
    }

    purrr::map_dfr(climate_zone_levels, function(zone_name) {
      zone_baseline <- baseline_zone_df %>% filter(climate_zone == zone_name)
      if (nrow(zone_baseline) == 0) {
        return(tibble())
      }

      baseline_col <- baseline_column_for_factor(spec$factor)
      reference_raw <- zone_baseline[[baseline_col]][1]
      reference_scaled <- (reference_raw - scale_row$mean_value[1]) / scale_row$sd_value[1]

      ref_data <- reference_row
      ref_data$climate_zone <- factor(zone_name, levels = levels(lagged_data$climate_zone))
      ref_data[[spec$lag]] <- reference_scaled
      ref_data <- ref_data[rep(1, length(x_scaled)), ]
      ref_lp <- suppressWarnings(predict(model, ref_data, type = "lpmatrix"))
      ref_term_lp <- ref_lp[, term_cols, drop = FALSE]
      ref_fit <- as.numeric(ref_term_lp %*% model_coef[term_cols])

      relative_fit <- fit_values - ref_fit
      diff_lp <- pred_term_lp - ref_term_lp
      diff_var <- rowSums((diff_lp %*% model_vcov[term_cols, term_cols, drop = FALSE]) * diff_lp)
      relative_se <- sqrt(pmax(diff_var, 0))
      edge_stabilized <- stabilize_lookup_edge_log_se(relative_se)
      relative_se <- edge_stabilized$log_se
      effect_or <- exp(relative_fit)
      constrained_temp <- apply_temperature_upper_tail_constraint(
        effect_or = effect_or,
        relative_se = relative_se,
        x_raw = x_raw,
        reference_raw = reference_raw,
        bacteria_name = bacteria,
        climate_factor_name = spec$factor
      )
      effect_or <- constrained_temp$effect_or
      relative_se <- constrained_temp$log_se
      lower_ci <- constrained_temp$lower_ci
      upper_ci <- constrained_temp$upper_ci

      observed_zone_scaled <- lagged_data %>%
        filter(climate_zone == zone_name) %>%
        pull(!!sym(spec$lag))
      observed_zone_raw <- observed_zone_scaled * scale_row$sd_value[1] + scale_row$mean_value[1]
      observed_effect_or <- approx(
        x = x_raw,
        y = effect_or,
        xout = clamp(observed_zone_raw, min(x_raw, na.rm = TRUE), max(x_raw, na.rm = TRUE)),
        rule = 2
      )$y
      effect_cap_lower <- NA_real_
      effect_cap_upper <- NA_real_
      if (!is.null(model_settings$historical_effect_quantiles)) {
        effect_cap_lower <- quantile(observed_effect_or, probs = model_settings$historical_effect_quantiles[1], na.rm = TRUE, names = FALSE)
        effect_cap_upper <- quantile(observed_effect_or, probs = model_settings$historical_effect_quantiles[2], na.rm = TRUE, names = FALSE)
        effect_cap_lower <- pmax(effect_cap_lower, 1e-4)
        effect_cap_upper <- pmax(effect_cap_upper, effect_cap_lower + 1e-4)

        effect_or <- clamp(effect_or, effect_cap_lower, effect_cap_upper)
        lower_ci <- clamp(lower_ci, effect_cap_lower, effect_cap_upper)
        upper_ci <- clamp(upper_ci, effect_cap_lower, effect_cap_upper)
      }

      effect_or <- apply_optional_bounds(effect_or, model_settings$response_or_bounds)
      lower_ci <- apply_optional_bounds(lower_ci, model_settings$response_or_bounds)
      upper_ci <- apply_optional_bounds(upper_ci, model_settings$response_or_bounds)
      lower_ci <- pmin(lower_ci, effect_or)
      upper_ci <- pmax(upper_ci, effect_or)

      tibble(
        Bacteria = bacteria,
        Climate_Factor = spec$factor,
        climate_zone = zone_name,
        Lag_Years = case_when(
          spec$factor == "Temperature" ~ lag_row$TMP_lag[1],
          spec$factor == "Precipitation" ~ lag_row$PREC_lag[1],
          spec$factor == "Humidity" ~ lag_row$HUM_lag[1],
          TRUE ~ lag_row$WET_lag[1]
        ),
        x_scaled = x_scaled,
        x_raw = x_raw,
        Effect_OR = effect_or,
        Lower_CI = lower_ci,
        Upper_CI = upper_ci,
        log_effect = log(pmax(effect_or, 1e-8)),
        log_se = pmax(relative_se, 0),
        support_lower = lower_raw,
        support_upper = upper_raw,
        effect_support_lower = effect_cap_lower,
        effect_support_upper = effect_cap_upper,
        projection_constraint_direction = constrained_temp$direction,
        projection_constraint_applied = constrained_temp$applied,
        projection_constraint_tail_start_raw = constrained_temp$tail_start_raw,
        projection_log_se_edge_stabilized = edge_stabilized$stabilized,
        observed_min = support_row$observed_min[1],
        observed_max = support_row$observed_max[1],
        mean_value = scale_row$mean_value[1],
        sd_value = scale_row$sd_value[1],
        reference_raw = reference_raw,
        reference_scaled = reference_scaled
      )
    })
  })

  lookup_rows
}

make_historical_outputs <- function(model, lagged_data, raw_data, bacteria) {
  fitted_logit <- as.numeric(predict(model, newdata = lagged_data))
  lagged_with_fit <- lagged_data %>%
    mutate(
      Bacteria = bacteria,
      fitted_logit_R = fitted_logit,
      fitted_R_prop = plogis(fitted_logit_R),
      fitted_R_pct = fitted_R_prop * 100
    )

  historical_zone_year <- lagged_with_fit %>%
    group_by(Bacteria, climate_zone, year) %>%
    summarise(
      observed_R_pct = mean(observed_R_pct, na.rm = TRUE),
      fitted_R_pct = mean(fitted_R_pct, na.rm = TRUE),
      residual_R_pct = observed_R_pct - fitted_R_pct,
      n_records = n(),
      .groups = "drop"
    )

  historical_global_year <- lagged_with_fit %>%
    group_by(Bacteria, year) %>%
    summarise(
      climate_zone = "Global Average",
      observed_R_pct = mean(observed_R_pct, na.rm = TRUE),
      fitted_R_pct = mean(fitted_R_pct, na.rm = TRUE),
      residual_R_pct = observed_R_pct - fitted_R_pct,
      n_records = n(),
      .groups = "drop"
    )

  baseline_zone <- raw_data %>%
    filter(year %in% baseline_years) %>%
    group_by(climate_zone) %>%
    summarise(
      baseline_resistance_prop = mean(observed_R_prop, na.rm = TRUE),
      baseline_resistance_pct = mean(observed_R_pct, na.rm = TRUE),
      baseline_temp = mean(TMP, na.rm = TRUE),
      baseline_precip = mean(PREC, na.rm = TRUE),
      baseline_humidity = mean(HUM, na.rm = TRUE),
      baseline_wetdays = mean(WET, na.rm = TRUE),
      n_records = n(),
      n_countries = n_distinct(NAME),
      .groups = "drop"
    ) %>%
    mutate(Bacteria = bacteria, .before = 1)

  baseline_global <- raw_data %>%
    filter(year %in% baseline_years) %>%
    summarise(
      Bacteria = bacteria,
      climate_zone = "Global Average",
      baseline_resistance_prop = mean(observed_R_prop, na.rm = TRUE),
      baseline_resistance_pct = mean(observed_R_pct, na.rm = TRUE),
      baseline_temp = mean(TMP, na.rm = TRUE),
      baseline_precip = mean(PREC, na.rm = TRUE),
      baseline_humidity = mean(HUM, na.rm = TRUE),
      baseline_wetdays = mean(WET, na.rm = TRUE),
      n_records = n(),
      n_countries = n_distinct(NAME)
    )

  list(
    historical_zone_year = historical_zone_year,
    historical_global_year = historical_global_year,
    baseline_zone = baseline_zone,
    baseline_global = baseline_global
  )
}

fit_projection_models <- function(lag_candidates_df) {
  curve_lookup_all <- list()
  scaling_all <- list()
  support_all <- list()
  term_stats_all <- list()
  primary_weight_all <- list()
  candidate_weight_all <- list()
  weight_details_all <- list()
  baseline_zone_all <- list()
  baseline_global_all <- list()
  history_zone_all <- list()
  history_global_all <- list()
  manifest_rows <- list()
  modifier_strength_df <- load_pls_modifier_strengths(variance_decomposition_path)

  for (spec in bacteria_specs) {
    candidate_rows <- lag_candidates_df %>%
      filter(Bacteria == spec$title) %>%
      arrange(Rank)

    if (nrow(candidate_rows) == 0) {
      warning("Missing lag candidate settings for ", spec$title)
      next
    }

    for (candidate_idx in seq_len(nrow(candidate_rows))) {
      lag_row <- candidate_rows[candidate_idx, , drop = FALSE]
      lag_combo <- lag_row$Lag_Combination[1]
      message("Fitting ", projection_version_tag, " simplified projection GAMM for ", spec$title, " / lag ", lag_combo)

      prepared <- prepare_projection_model_data(file.path(input_data_dir, spec$file_name), lag_row)
      model <- build_simplified_model(prepared$lagged_data)
      term_stats <- extract_model_term_statistics(model, spec$title, lag_row) %>%
        mutate(
          Lag_Combination = lag_combo,
          delta_aic = lag_row$delta_aic[1],
          lag_sampling_weight = lag_row$lag_sampling_weight[1],
          is_primary = lag_row$is_primary[1]
        )
      derived_weights <- derive_weights_from_term_stats(term_stats, spec$title)
      historical_outputs <- make_historical_outputs(model, prepared$lagged_data, prepared$raw_data, spec$title)
      response_lookup <- make_response_lookup(
        model = model,
        lagged_data = prepared$lagged_data,
        scaling_df = prepared$scaling_df,
        support_df = prepared$support_df,
        baseline_zone_df = historical_outputs$baseline_zone,
        bacteria = spec$title,
        lag_row = lag_row
      ) %>%
        mutate(
          Lag_Combination = lag_combo,
          delta_aic = lag_row$delta_aic[1],
          lag_sampling_weight = lag_row$lag_sampling_weight[1],
          is_primary = lag_row$is_primary[1]
        )

      model_rds_path <- file.path(
        results_root,
        "03_diagnostics",
        paste0("Projection_Simplified_ModelC_", sanitize_file_stem(spec$title), "_", sanitize_file_stem(lag_combo), ".rds")
      )
      saveRDS(
        list(
          bacteria = spec$title,
          lag_settings = lag_row,
          scaling = prepared$scaling_df,
          support = prepared$support_df,
          term_statistics = term_stats,
          weights = derived_weights$weights,
          model = model
        ),
        model_rds_path
      )

      curve_lookup_all[[paste(spec$title, lag_combo, sep = "__")]] <- response_lookup
      term_stats_all[[paste(spec$title, lag_combo, sep = "__")]] <- term_stats
      candidate_weight_all[[paste(spec$title, lag_combo, sep = "__")]] <- derived_weights$weights %>%
        mutate(
          Lag_Combination = lag_combo,
          delta_aic = lag_row$delta_aic[1],
          lag_sampling_weight = lag_row$lag_sampling_weight[1],
          is_primary = lag_row$is_primary[1]
        )
      weight_details_all[[paste(spec$title, lag_combo, sep = "__")]] <- derived_weights$details %>%
        mutate(
          Lag_Combination = lag_combo,
          delta_aic = lag_row$delta_aic[1],
          lag_sampling_weight = lag_row$lag_sampling_weight[1],
          is_primary = lag_row$is_primary[1]
        )

      if (isTRUE(lag_row$is_primary[1])) {
        scaling_all[[spec$title]] <- prepared$scaling_df %>% mutate(Bacteria = spec$title, .before = 1)
        support_all[[spec$title]] <- prepared$support_df %>% mutate(Bacteria = spec$title, .before = 1)
        primary_weight_all[[spec$title]] <- derived_weights$weights
        baseline_zone_all[[spec$title]] <- historical_outputs$baseline_zone
        baseline_global_all[[spec$title]] <- historical_outputs$baseline_global
        history_zone_all[[spec$title]] <- historical_outputs$historical_zone_year
        history_global_all[[spec$title]] <- historical_outputs$historical_global_year
      }

      manifest_rows[[paste(spec$title, lag_combo, sep = "__")]] <- tibble(
        Bacteria = spec$title,
        Lag_Combination = lag_combo,
        AIC = lag_row$AIC[1],
        BIC = lag_row$BIC[1],
        TMP_lag = lag_row$TMP_lag[1],
        PREC_lag = lag_row$PREC_lag[1],
        HUM_lag = lag_row$HUM_lag[1],
        WET_lag = lag_row$WET_lag[1],
        delta_aic = lag_row$delta_aic[1],
        lag_sampling_weight = lag_row$lag_sampling_weight[1],
        Rank = lag_row$Rank[1],
        is_primary = lag_row$is_primary[1],
        Data_File = spec$file_name,
        Model_RDS = model_rds_path,
        N_Rows_Fitted = nrow(prepared$lagged_data)
      )
    }
  }

  baseline_zone_bind <- bind_rows(baseline_zone_all)
  baseline_global_bind <- bind_rows(baseline_global_all)
  historical_zone_bind <- bind_rows(history_zone_all)
  historical_global_bind <- bind_rows(history_global_all)
  residual_calibration <- build_recursive_uncertainty_calibration(
    historical_global_df = historical_global_bind,
    historical_zone_df = historical_zone_bind,
    baseline_global_df = baseline_global_bind,
    baseline_zone_df = baseline_zone_bind
  )

  list(
    curve_lookup_df = bind_rows(curve_lookup_all),
    scaling_df = bind_rows(scaling_all),
    support_df = bind_rows(support_all),
    term_stats_df = bind_rows(term_stats_all),
    weights_df = bind_rows(primary_weight_all),
    candidate_weights_df = bind_rows(candidate_weight_all),
    weight_details_df = bind_rows(weight_details_all),
    baseline_zone_df = baseline_zone_bind,
    baseline_global_df = baseline_global_bind,
    historical_zone_df = historical_zone_bind,
    historical_global_df = historical_global_bind,
    recursive_uncertainty_df = residual_calibration$calibration_df,
    historical_global_residual_df = residual_calibration$historical_global_residual_df,
    historical_zone_residual_df = residual_calibration$historical_zone_residual_df,
    recursive_uncertainty_zone_summary_df = residual_calibration$zone_summary_df,
    manifest_df = bind_rows(manifest_rows),
    modifier_strength_df = modifier_strength_df
  )
}

prepare_climate_scenarios <- function(climate_data) {
  climate_long <- climate_data %>%
    transmute(
      country_name = NAME,
      climate_model = as.character(model),
      scenario = tolower(as.character(ssp)),
      year = as.integer(year),
      humidity = as.numeric(hurs_annual_mean),
      precipitation = as.numeric(pr_annual_total),
      temperature = as.numeric(tas_annual_mean),
      wetdays = as.numeric(wet_days_yearly0.1),
      lat = as.numeric(lat)
    ) %>%
    mutate(
      climate_zone = case_when(
        abs(lat) > 66.5 ~ "Polar Zone",
        abs(lat) > 23.5 ~ "Temperate Zone",
        TRUE ~ "Tropical Zone"
      )
    ) %>%
    filter(
      scenario %in% scenario_levels,
      climate_model %in% climate_models_expected,
      year >= min(projection_years),
      year <= max(projection_years)
    )

  zone_year_model_means <- climate_long %>%
    group_by(year, scenario, climate_model, climate_zone) %>%
    summarise(
      temperature = mean(temperature, na.rm = TRUE),
      humidity = mean(humidity, na.rm = TRUE),
      precipitation = mean(precipitation, na.rm = TRUE),
      wetdays = mean(wetdays, na.rm = TRUE),
      n_countries = n_distinct(country_name),
      .groups = "drop"
    )

  complete_grid <- expand.grid(
    year = projection_years,
    scenario = scenario_levels,
    climate_model = climate_models_expected,
    climate_zone = climate_zone_levels,
    stringsAsFactors = FALSE
  )

  filled_scenarios <- left_join(
    complete_grid,
    zone_year_model_means,
    by = c("year", "scenario", "climate_model", "climate_zone")
  ) %>%
    group_by(scenario, climate_model, climate_zone) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      temperature = fill_series(year, temperature),
      humidity = fill_series(year, humidity),
      precipitation = fill_series(year, precipitation),
      wetdays = fill_series(year, wetdays),
      n_countries = round(fill_series(year, n_countries))
    ) %>%
    ungroup()

  model_coverage <- filled_scenarios %>%
    group_by(scenario) %>%
    summarise(
      n_climate_models = n_distinct(climate_model[!is.na(temperature)]),
      .groups = "drop"
    )

  list(
    filled_scenarios = filled_scenarios,
    model_coverage = model_coverage,
    zone_year_model_means = filled_scenarios,
    raw_country_year = climate_long
  )
}

create_lagged_series <- function(values, lag_years, baseline_value) {
  out <- numeric(length(values))
  for (i in seq_along(values)) {
    out[i] <- if (i <= lag_years) baseline_value else values[i - lag_years]
  }
  out
}

interpolate_lookup_values <- function(lookup_df, bacteria_name, climate_factor_name, zone_name, raw_values, value_col, lag_combo = NULL) {
  rows <- lookup_df %>%
    filter(Bacteria == bacteria_name, Climate_Factor == climate_factor_name, climate_zone == zone_name)

  if (!is.null(lag_combo)) {
    rows <- rows %>% filter(Lag_Combination == lag_combo)
  }

  rows <- rows %>%
    arrange(x_raw)

  if (nrow(rows) == 0) {
    stop("No response lookup found for ", bacteria_name, " / ", climate_factor_name, call. = FALSE)
  }

  if (!value_col %in% names(rows)) {
    stop("Lookup column ", value_col, " not found for ", bacteria_name, " / ", climate_factor_name, call. = FALSE)
  }

  y_values <- rows[[value_col]]
  if (all(is.na(y_values))) {
    stop("Lookup column ", value_col, " is all NA for ", bacteria_name, " / ", climate_factor_name, call. = FALSE)
  }
  interp_df <- tibble(x_raw = rows$x_raw, y_value = y_values) %>%
    filter(is.finite(x_raw), is.finite(y_value)) %>%
    group_by(x_raw) %>%
    summarise(y_value = mean(y_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(x_raw)

  if (nrow(interp_df) < 2) {
    stop("Insufficient lookup support after deduplication for ", bacteria_name, " / ", climate_factor_name, call. = FALSE)
  }

  use_linear_tail <- identical(model_settings$lookup_extrapolation$method, "linear_tail") && !identical(value_col, "log_se")
  if (use_linear_tail) {
    return(interpolate_with_linear_tail(
      x = interp_df$x_raw,
      y = interp_df$y_value,
      xout = raw_values,
      slope_fraction = model_settings$lookup_extrapolation$slope_fraction %||% 0.05,
      min_points = model_settings$lookup_extrapolation$min_points %||% 5L
    ))
  }

  clamped_raw <- clamp(raw_values, min(interp_df$x_raw, na.rm = TRUE), max(interp_df$x_raw, na.rm = TRUE))
  approx(x = interp_df$x_raw, y = interp_df$y_value, xout = clamped_raw, rule = 2)$y
}

apply_crkp_controls <- function(effect_vec) {
  trend_component <- as.numeric(effect_vec)
  window_size <- model_settings$crkp_controls$enhanced_smoothing

  if (length(trend_component) > window_size) {
    weights <- rep(1 / window_size, window_size)
    ma_trend <- stats::filter(trend_component, weights, sides = 2)
    na_idx <- is.na(ma_trend)

    if (any(na_idx)) {
      if (sum(!na_idx) > 1) {
        ma_trend[na_idx] <- approx(
          x = which(!na_idx),
          y = ma_trend[!na_idx],
          xout = which(na_idx),
          rule = 2
        )$y
      } else {
        ma_trend[na_idx] <- trend_component[na_idx]
      }
    }

    fluctuation_component <- trend_component - ma_trend
    reduced_fluctuation <- fluctuation_component * (1 - model_settings$crkp_controls$variability_reduction)
    trend_component <- as.numeric(ma_trend + reduced_fluctuation)
  }

  if (model_settings$crkp_controls$end_correction && length(trend_component) > 3) {
    n_vals <- length(trend_component)
    trend_component[1] <- (3 * trend_component[1] + 2 * trend_component[2] + trend_component[3]) / 6
    trend_component[n_vals] <- (3 * trend_component[n_vals] + 2 * trend_component[n_vals - 1] + trend_component[n_vals - 2]) / 6
  }

  max_effect <- model_settings$crkp_controls$max_effect_limit * 1.5
  clamp(trend_component, 1 - max_effect, 1 + max_effect)
}

diffusion_base_rate_for_bacteria <- function(bacteria_name) {
  switch(
    bacteria_name,
    "CR-Ab" = 0.0018,
    "CR-Pa" = 0.0020,
    "3GCR-Kp" = 0.0022,
    "CR-Kp" = model_settings$crkp_controls$max_growth_rate,
    "3GCR-Ec" = 0.0018,
    "CR-Ec" = 0.0016,
    0.002
  )
}

simulate_recursive_path <- function(years, combined_effect, baseline_resistance, bacteria_name, scenario_name, seed_prefix, initial_variation = 0, memory_factor_multiplier = 1.0, recursive_persistence = NULL, recursive_innovation_sd = NULL) {
  is_crkp <- identical(bacteria_name, "CR-Kp")
  persistence <- recursive_persistence %||% model_settings$time_series_correlation$default_persistence
  innovation_sd <- recursive_innovation_sd %||% model_settings$time_series_correlation$default_innovation_sd

  combined_effect <- as.numeric(combined_effect)
  combined_effect[!is.finite(combined_effect)] <- 1
  combined_effect[1] <- 1.0 + initial_variation

  years_passed_vec <- years - min(years)
  scenario_variation <- model_settings$yearly_variation_factors[[scenario_name]]

  set.seed(create_seed_from_string(seed_prefix))
  innovations_vec <- rnorm(
    length(years),
    mean = 0,
    sd = innovation_sd * scenario_variation
  )

  current_persistence <- 0
  out <- numeric(length(years))

  for (i in seq_along(years)) {
    if (i == 1) {
      current_rate <- baseline_resistance * combined_effect[i]
    } else {
      prev_rate <- out[i - 1]
      if (!is.finite(prev_rate) || prev_rate <= 0) {
        prev_rate <- baseline_resistance
      }
      memory_factor <- clamp((model_settings$memory_factors[[scenario_name]] %||% 0.7) * memory_factor_multiplier, 0.50, 0.95)
      if (is_crkp) {
        memory_factor <- clamp(max(model_settings$crkp_controls$min_memory_factor, memory_factor + 0.10), 0.50, 0.97)
      }

      innovation_factor <- 1 - memory_factor
      years_passed <- years_passed_vec[i]
      diffusion_rate <- 1 + diffusion_base_rate_for_bacteria(bacteria_name) * min(years_passed, 40)
      diffusion_rate <- diffusion_rate * compute_ssp_diffusion_multiplier(scenario_name, years_passed)
      current_persistence <- persistence * current_persistence + innovations_vec[i]
      if (!is.finite(current_persistence)) {
        current_persistence <- 0
      }
      climate_effect <- baseline_resistance * combined_effect[i] * diffusion_rate * clamp(1 + current_persistence, 0.10, 2.50)
      if (!is.finite(climate_effect)) {
        climate_effect <- prev_rate
      }
      current_rate <- memory_factor * prev_rate + innovation_factor * climate_effect
      if (!is.finite(current_rate)) {
        current_rate <- prev_rate
      }
    }

    if (is_crkp && i > 1) {
      prev_rate <- if (is.finite(out[i - 1]) && out[i - 1] > 0) out[i - 1] else baseline_resistance
      growth_rate <- (current_rate - prev_rate) / prev_rate
      if (!is.finite(growth_rate)) {
        growth_rate <- 0
      }
      max_allowed_growth <- model_settings$crkp_controls$max_growth_rate * model_settings$yearly_variation_factors[[scenario_name]]
      if (abs(growth_rate) > max_allowed_growth) {
        current_rate <- prev_rate * (1 + sign(growth_rate) * max_allowed_growth)
      }
      current_rate <- min(max(current_rate, baseline_resistance * 0.45), baseline_resistance * 3.0)
    }

    out[i] <- clamp(current_rate, 0.01, 99.99)
  }

  out
}

prepare_projection_configs <- function(model_artifacts, lag_candidates_df, climate_inputs) {
  baseline_zone <- model_artifacts$baseline_zone_df
  lookup_df <- model_artifacts$curve_lookup_df
  candidate_weights_df <- model_artifacts$candidate_weights_df
  modifier_strength_df <- model_artifacts$modifier_strength_df
  recursive_uncertainty_df <- model_artifacts$recursive_uncertainty_df
  climate_scenarios <- climate_inputs$filled_scenarios

  config_list <- list()
  effect_rows <- list()
  recursive_rows <- list()
  idx <- 1L

  for (bacteria_name in bacteria_levels) {
    candidate_rows <- lag_candidates_df %>%
      filter(Bacteria == bacteria_name) %>%
      arrange(Rank)
    if (nrow(candidate_rows) == 0) {
      next
    }

    modifier_strength_match <- modifier_strength_df %>%
      filter(Bacteria == bacteria_name) %>%
      pull(development_modifier_strength)
    modifier_strength <- if (length(modifier_strength_match) == 0 || is.na(modifier_strength_match[1])) 0.35 else modifier_strength_match[1]
    recursive_uncertainty_row <- recursive_uncertainty_df %>%
      filter(Bacteria == bacteria_name)
    recursive_persistence <- recursive_uncertainty_row$calibrated_persistence[1] %||% model_settings$time_series_correlation$default_persistence
    recursive_innovation_sd <- recursive_uncertainty_row$calibrated_innovation_sd[1] %||% model_settings$time_series_correlation$default_innovation_sd
    initial_variation_sd <- recursive_uncertainty_row$calibrated_initial_sd[1] %||% model_settings$empirical_recursive_uncertainty$initial_sd_bounds[1]
    empirical_variation_factor <- recursive_uncertainty_row$empirical_variation_factor[1] %||% 1.0
    calibration_source <- recursive_uncertainty_row$calibration_source[1] %||% "default_fallback"

    for (candidate_idx in seq_len(nrow(candidate_rows))) {
      lag_row <- candidate_rows[candidate_idx, , drop = FALSE]
      lag_combo <- lag_row$Lag_Combination[1]
      weight_row <- candidate_weights_df %>%
        filter(Bacteria == bacteria_name, Lag_Combination == lag_combo)
      if (nrow(weight_row) == 0) {
        next
      }

      weights_vec <- c(
        "Temperature" = unname(weight_row$Temperature[1]),
        "Humidity" = unname(weight_row$Humidity[1]),
        "Precipitation" = unname(weight_row$Precipitation[1]),
        "Wet Days" = unname(weight_row$`Wet Days`[1])
      )

      for (scenario_name in scenario_levels) {
        for (climate_model_name in climate_models_expected) {
          for (zone_name in climate_zone_levels) {
            zone_baseline <- baseline_zone %>%
              filter(Bacteria == bacteria_name, climate_zone == zone_name)
            if (nrow(zone_baseline) == 0) {
              next
            }

            scenario_climate <- climate_scenarios %>%
              filter(
                climate_zone == zone_name,
                scenario == scenario_name,
                climate_model == climate_model_name,
                year >= min(projection_years),
                year <= max(projection_years)
              ) %>%
              arrange(year)
            if (nrow(scenario_climate) == 0) {
              next
            }

            years <- scenario_climate$year
            temp_lagged <- create_lagged_series(scenario_climate$temperature, lag_row$TMP_lag[1], zone_baseline$baseline_temp[1])
            humid_lagged <- create_lagged_series(scenario_climate$humidity, lag_row$HUM_lag[1], zone_baseline$baseline_humidity[1])
            precip_lagged <- create_lagged_series(scenario_climate$precipitation, lag_row$PREC_lag[1], zone_baseline$baseline_precip[1])
            wet_lagged <- create_lagged_series(scenario_climate$wetdays, lag_row$WET_lag[1], zone_baseline$baseline_wetdays[1])

            temp_effect <- pmax(interpolate_lookup_values(lookup_df, bacteria_name, "Temperature", zone_name, temp_lagged, "Effect_OR", lag_combo = lag_combo), 1e-06)
            humid_effect <- pmax(interpolate_lookup_values(lookup_df, bacteria_name, "Humidity", zone_name, humid_lagged, "Effect_OR", lag_combo = lag_combo), 1e-06)
            precip_effect <- pmax(interpolate_lookup_values(lookup_df, bacteria_name, "Precipitation", zone_name, precip_lagged, "Effect_OR", lag_combo = lag_combo), 1e-06)
            wet_effect <- pmax(interpolate_lookup_values(lookup_df, bacteria_name, "Wet Days", zone_name, wet_lagged, "Effect_OR", lag_combo = lag_combo), 1e-06)

            temp_log_se <- interpolate_lookup_values(lookup_df, bacteria_name, "Temperature", zone_name, temp_lagged, "log_se", lag_combo = lag_combo)
            humid_log_se <- interpolate_lookup_values(lookup_df, bacteria_name, "Humidity", zone_name, humid_lagged, "log_se", lag_combo = lag_combo)
            precip_log_se <- interpolate_lookup_values(lookup_df, bacteria_name, "Precipitation", zone_name, precip_lagged, "log_se", lag_combo = lag_combo)
            wet_log_se <- interpolate_lookup_values(lookup_df, bacteria_name, "Wet Days", zone_name, wet_lagged, "log_se", lag_combo = lag_combo)

            development_modifier <- compute_development_modifier(
              bacteria_name = bacteria_name,
              scenario_name = scenario_name,
              years = years,
              modifier_strength = modifier_strength
            )

            combined_effect <- combine_weighted_effect_log(
              temp_effect = temp_effect,
              humid_effect = humid_effect,
              precip_effect = precip_effect,
              wetdays_effect = wet_effect,
              weights_vec = weights_vec,
              development_modifier = development_modifier
            )
            if (identical(bacteria_name, "CR-Kp")) {
              combined_effect <- apply_crkp_controls(combined_effect)
            }

            set.seed(create_seed_from_string(paste(bacteria_name, lag_combo, zone_name, scenario_name, climate_model_name, "initial", sep = "_")))
            initial_variation <- rnorm(
              1,
              mean = 0,
              sd = initial_variation_sd * (model_settings$yearly_variation_factors[[scenario_name]] %||% 1.0)
            )

            config_list[[idx]] <- list(
              bacteria = bacteria_name,
              scenario = scenario_name,
              climate_model = climate_model_name,
              climate_zone = zone_name,
              Lag_Combination = lag_combo,
              lag_sampling_weight = lag_row$lag_sampling_weight[1],
              delta_aic = lag_row$delta_aic[1],
              is_primary = lag_row$is_primary[1],
              years = years,
              baseline_resistance = zone_baseline$baseline_resistance_pct[1],
              initial_variation = initial_variation,
              central_combined_effect = combined_effect,
              temp_effect = temp_effect,
              humid_effect = humid_effect,
              precip_effect = precip_effect,
              wet_effect = wet_effect,
              temp_log_se = temp_log_se,
              humid_log_se = humid_log_se,
              precip_log_se = precip_log_se,
              wet_log_se = wet_log_se,
              development_modifier = development_modifier,
              development_modifier_strength = modifier_strength,
              recursive_persistence = recursive_persistence,
              recursive_innovation_sd = recursive_innovation_sd,
              recursive_initial_sd = initial_variation_sd,
              empirical_variation_factor = empirical_variation_factor,
              recursive_uncertainty_source = calibration_source,
              weight_temperature = unname(weights_vec["Temperature"]),
              weight_humidity = unname(weights_vec["Humidity"]),
              weight_precipitation = unname(weights_vec["Precipitation"]),
              weight_wetdays = unname(weights_vec["Wet Days"]),
              zone_weight = zone_baseline$n_records[1],
              n_countries = round(mean(scenario_climate$n_countries, na.rm = TRUE))
            )

            effect_rows[[idx]] <- tibble(
              bacteria = bacteria_name,
              scenario = scenario_name,
              climate_model = climate_model_name,
              climate_zone = zone_name,
              Lag_Combination = lag_combo,
              lag_sampling_weight = lag_row$lag_sampling_weight[1],
              delta_aic = lag_row$delta_aic[1],
              is_primary = lag_row$is_primary[1],
              year = years,
              temperature = scenario_climate$temperature,
              humidity = scenario_climate$humidity,
              precipitation = scenario_climate$precipitation,
              wetdays = scenario_climate$wetdays,
              temperature_lagged = temp_lagged,
              humidity_lagged = humid_lagged,
              precipitation_lagged = precip_lagged,
              wetdays_lagged = wet_lagged,
              temp_effect_or = temp_effect,
              humid_effect_or = humid_effect,
              precip_effect_or = precip_effect,
              wetdays_effect_or = wet_effect,
              temp_log_se = temp_log_se,
              humid_log_se = humid_log_se,
              precip_log_se = precip_log_se,
              wetdays_log_se = wet_log_se,
              development_modifier = development_modifier,
              combined_effect_central = combined_effect,
              recursive_persistence = recursive_persistence,
              recursive_innovation_sd = recursive_innovation_sd,
              recursive_initial_sd = initial_variation_sd,
              empirical_variation_factor = empirical_variation_factor,
              recursive_uncertainty_source = calibration_source,
              weight_temperature = unname(weights_vec["Temperature"]),
              weight_humidity = unname(weights_vec["Humidity"]),
              weight_precipitation = unname(weights_vec["Precipitation"]),
              weight_wetdays = unname(weights_vec["Wet Days"]),
              development_modifier_strength = modifier_strength
            )

            recursive_rows[[idx]] <- tibble(
              bacteria = bacteria_name,
              scenario = scenario_name,
              climate_model = climate_model_name,
              climate_zone = zone_name,
              Lag_Combination = lag_combo,
              lag_sampling_weight = lag_row$lag_sampling_weight[1],
              delta_aic = lag_row$delta_aic[1],
              is_primary = lag_row$is_primary[1],
              zone_weight = zone_baseline$n_records[1],
              memory_factor = model_settings$memory_factors[[scenario_name]],
              innovation_factor = 1 - model_settings$memory_factors[[scenario_name]],
              persistence = recursive_persistence,
              innovation_sd = recursive_innovation_sd,
              initial_variation_sd = initial_variation_sd,
              diffusion_base_rate = diffusion_base_rate_for_bacteria(bacteria_name),
              ssp_diffusion_modifier_enabled = isTRUE(model_settings$use_ssp_diffusion_modifier),
              ssp_diffusion_modifier_net_rate = model_settings$ssp_diffusion_modifier[[scenario_name]]$net_annual_rate %||% 0,
              ssp_diffusion_modifier_max_cumulative = model_settings$ssp_diffusion_modifier[[scenario_name]]$max_cumulative %||% 1,
              ssp_diffusion_modifier_ramp_up_years = model_settings$ssp_diffusion_modifier[[scenario_name]]$ramp_up_years %||% NA_real_,
              bacteria_variation_factor = empirical_variation_factor,
              recursive_uncertainty_source = calibration_source,
              direct_socioeconomic_factor = model_settings$use_development_modifier,
              direct_antibiotic_factor = model_settings$use_antibiotic_modifier,
              development_modifier_strength = modifier_strength,
              uncertainty_mode = current_uncertainty_mode()
            )

            idx <- idx + 1L
          }
        }
      }
    }
  }

  list(
    configs = config_list,
    effect_components = bind_rows(effect_rows),
    recursive_parameters = bind_rows(recursive_rows),
    lag_candidates = lag_candidates_df
  )
}

build_central_prediction_diagnostics <- function(configs, baseline_global_df) {
  baseline_global_lookup <- baseline_global_df %>%
    transmute(bacteria = Bacteria, baseline_global_pct = baseline_resistance_pct)

  zone_predictions <- purrr::map_dfr(configs, function(cfg) {
    central_path <- simulate_recursive_path(
      years = cfg$years,
      combined_effect = cfg$central_combined_effect,
      baseline_resistance = cfg$baseline_resistance,
      bacteria_name = cfg$bacteria,
      scenario_name = cfg$scenario,
      seed_prefix = paste(cfg$bacteria, cfg$scenario, cfg$climate_model, cfg$climate_zone, "central", sep = "_"),
      initial_variation = cfg$initial_variation,
      recursive_persistence = cfg$recursive_persistence,
      recursive_innovation_sd = cfg$recursive_innovation_sd
    )

    tibble(
      year = cfg$years,
      bacteria = cfg$bacteria,
      scenario = cfg$scenario,
      climate_model = cfg$climate_model,
      climate_zone = cfg$climate_zone,
      zone_weight = cfg$zone_weight,
      resistance_rate = central_path
    )
  })

  global_model_predictions <- zone_predictions %>%
    group_by(bacteria, scenario, climate_model, year) %>%
    summarise(
      resistance_rate = weighted_mean_safe(resistance_rate, zone_weight),
      climate_zone = "Global Average",
      aggregation_weight_source = "baseline_zone_n_records",
      .groups = "drop"
    ) %>%
    group_by(bacteria, scenario, climate_model) %>%
    mutate(relative_change = 100 * (resistance_rate / first(resistance_rate) - 1)) %>%
    ungroup() %>%
    left_join(baseline_global_lookup, by = "bacteria") %>%
    mutate(relative_change_baseline = 100 * (resistance_rate / baseline_global_pct - 1)) %>%
    select(-baseline_global_pct)

  ensemble_predictions <- global_model_predictions %>%
    group_by(bacteria, scenario, year) %>%
    summarise(
      resistance_rate = mean(resistance_rate, na.rm = TRUE),
      climate_model = "Ensemble Mean",
      climate_zone = "Global Average",
      .groups = "drop"
    ) %>%
    group_by(bacteria, scenario, climate_model) %>%
    mutate(relative_change = 100 * (resistance_rate / first(resistance_rate) - 1)) %>%
    ungroup() %>%
    left_join(baseline_global_lookup, by = "bacteria") %>%
    mutate(relative_change_baseline = 100 * (resistance_rate / baseline_global_pct - 1)) %>%
    select(-baseline_global_pct)

  list(
    zone_predictions = zone_predictions,
    global_model_predictions = global_model_predictions,
    ensemble_predictions = ensemble_predictions
  )
}

sample_combined_effect_path <- function(cfg, draw_idx, sampled_weights_vec = NULL) {
  if (!isTRUE(model_settings$response_uncertainty$include_parameter_uncertainty)) {
    if (is.null(sampled_weights_vec)) {
      return(cfg$central_combined_effect)
    }
    combined_no_param <- combine_weighted_effect_log(
      temp_effect = cfg$temp_effect,
      humid_effect = cfg$humid_effect,
      precip_effect = cfg$precip_effect,
      wetdays_effect = cfg$wet_effect,
      weights_vec = sampled_weights_vec,
      development_modifier = cfg$development_modifier
    )
    if (identical(cfg$bacteria, "CR-Kp")) {
      combined_no_param <- apply_crkp_controls(combined_no_param)
    }
    return(combined_no_param)
  }

  years <- cfg$years
  rho <- model_settings$response_uncertainty$temporal_correlation
  scenario_scale <- if (isTRUE(model_settings$response_uncertainty$inflate_by_scenario)) {
    model_settings$scenario_uncertainty_factors[[cfg$scenario]] %||% 1.0
  } else {
    1.0
  }
  time_scale <- if (isTRUE(model_settings$response_uncertainty$inflate_by_time)) {
    1 + 0.35 * ((years - min(years)) / max(1, (max(years) - min(years))))
  } else {
    1.0
  }

  temp_noise <- generate_correlated_normal(length(years), rho, paste(cfg$bacteria, cfg$scenario, cfg$climate_model, cfg$climate_zone, draw_idx, "temp", sep = "_"))
  humid_noise <- generate_correlated_normal(length(years), rho, paste(cfg$bacteria, cfg$scenario, cfg$climate_model, cfg$climate_zone, draw_idx, "humid", sep = "_"))
  precip_noise <- generate_correlated_normal(length(years), rho, paste(cfg$bacteria, cfg$scenario, cfg$climate_model, cfg$climate_zone, draw_idx, "precip", sep = "_"))
  wet_noise <- generate_correlated_normal(length(years), rho, paste(cfg$bacteria, cfg$scenario, cfg$climate_model, cfg$climate_zone, draw_idx, "wet", sep = "_"))

  temp_draw <- exp(log(pmax(cfg$temp_effect, 1e-8)) + temp_noise * cfg$temp_log_se * scenario_scale * time_scale)
  humid_draw <- exp(log(pmax(cfg$humid_effect, 1e-8)) + humid_noise * cfg$humid_log_se * scenario_scale * time_scale)
  precip_draw <- exp(log(pmax(cfg$precip_effect, 1e-8)) + precip_noise * cfg$precip_log_se * scenario_scale * time_scale)
  wet_draw <- exp(log(pmax(cfg$wet_effect, 1e-8)) + wet_noise * cfg$wet_log_se * scenario_scale * time_scale)
  temp_draw[!is.finite(temp_draw)] <- cfg$temp_effect[!is.finite(temp_draw)]
  humid_draw[!is.finite(humid_draw)] <- cfg$humid_effect[!is.finite(humid_draw)]
  precip_draw[!is.finite(precip_draw)] <- cfg$precip_effect[!is.finite(precip_draw)]
  wet_draw[!is.finite(wet_draw)] <- cfg$wet_effect[!is.finite(wet_draw)]
  temp_draw <- apply_optional_bounds(temp_draw, model_settings$response_or_bounds)
  humid_draw <- apply_optional_bounds(humid_draw, model_settings$response_or_bounds)
  precip_draw <- apply_optional_bounds(precip_draw, model_settings$response_or_bounds)
  wet_draw <- apply_optional_bounds(wet_draw, model_settings$response_or_bounds)

  weights_vec <- sampled_weights_vec %||% c(
    "Temperature" = unname(cfg$weight_temperature),
    "Humidity" = unname(cfg$weight_humidity),
    "Precipitation" = unname(cfg$weight_precipitation),
    "Wet Days" = unname(cfg$weight_wetdays)
  )

  combined_draw <- combine_weighted_effect_log(
    temp_effect = temp_draw,
    humid_effect = humid_draw,
    precip_effect = precip_draw,
    wetdays_effect = wet_draw,
    weights_vec = weights_vec,
    development_modifier = cfg$development_modifier
  )

  if (identical(cfg$bacteria, "CR-Kp")) {
    combined_draw <- apply_crkp_controls(combined_draw)
  }

  combined_draw
}

run_projection_monte_carlo <- function(configs, baseline_global_df, climate_inputs, collect_zone_summary = TRUE) {
  grouped_configs <- split(configs, vapply(configs, function(x) paste(x$bacteria, x$scenario, sep = "__"), character(1)))
  baseline_global_lookup <- baseline_global_df %>%
    transmute(bacteria = Bacteria, baseline_global_pct = baseline_resistance_pct)

  annual_rows <- list()
  period_rows <- list()
  zone_annual_rows <- list()
  zone_period_rows <- list()
  draw_parameter_rows <- list()
  annual_idx <- 1L
  period_idx <- 1L
  zone_annual_idx <- 1L
  zone_period_idx <- 1L
  draw_parameter_idx <- 1L

  lower_prob <- (1 - monte_carlo_settings$confidence_level) / 2
  upper_prob <- 1 - lower_prob

  for (group_key in names(grouped_configs)) {
    config_group <- grouped_configs[[group_key]]
    bacteria_name <- config_group[[1]]$bacteria
    scenario_name <- config_group[[1]]$scenario
    message("Monte Carlo summary for ", bacteria_name, " / ", scenario_name)
    years <- config_group[[1]]$years
    n_years <- length(years)
    climate_models <- unique(vapply(config_group, `[[`, character(1), "climate_model"))
    zone_names <- climate_zone_levels[climate_zone_levels %in% unique(vapply(config_group, `[[`, character(1), "climate_zone"))]

    baseline_global <- baseline_global_lookup %>%
      filter(bacteria == bacteria_name) %>%
      pull(baseline_global_pct) %>%
      .[[1]]

    lag_candidate_table <- config_group %>%
      purrr::map_dfr(function(cfg) {
        tibble(
          Lag_Combination = cfg$Lag_Combination,
          lag_sampling_weight = cfg$lag_sampling_weight,
          delta_aic = cfg$delta_aic,
          is_primary = cfg$is_primary
        )
      }) %>%
      distinct() %>%
      arrange(delta_aic, Lag_Combination)

    draw_matrix <- matrix(NA_real_, nrow = n_years, ncol = monte_carlo_settings$n_simulations)
    zone_draw_arrays <- setNames(
      lapply(zone_names, function(x) matrix(NA_real_, nrow = n_years, ncol = monte_carlo_settings$n_simulations)),
      zone_names
    )
    zone_baseline_lookup <- setNames(
      lapply(zone_names, function(zone_name) {
        model_configs <- config_group[vapply(config_group, function(x) identical(x$climate_zone, zone_name), logical(1))]
        model_configs[[1]]$baseline_resistance
      }),
      zone_names
    )

    for (draw_idx in seq_len(monte_carlo_settings$n_simulations)) {
      set.seed(create_seed_from_string(paste(bacteria_name, scenario_name, draw_idx, "lag_combo", sep = "_")))
      sampled_lag_combo <- sample(
        x = lag_candidate_table$Lag_Combination,
        size = 1,
        prob = lag_candidate_table$lag_sampling_weight
      )
      selected_configs <- config_group[vapply(config_group, function(x) identical(x$Lag_Combination, sampled_lag_combo), logical(1))]
      if (length(selected_configs) == 0) {
        next
      }

      base_weights_vec <- c(
        "Temperature" = unname(selected_configs[[1]]$weight_temperature),
        "Humidity" = unname(selected_configs[[1]]$weight_humidity),
        "Precipitation" = unname(selected_configs[[1]]$weight_precipitation),
        "Wet Days" = unname(selected_configs[[1]]$weight_wetdays)
      )
      sampled_weights_vec <- sample_perturbed_weights(
        base_weights_vec,
        seed_prefix = paste(bacteria_name, scenario_name, sampled_lag_combo, draw_idx, "weights", sep = "_")
      )
      sampled_memory_multiplier <- sample_memory_multiplier(
        seed_prefix = paste(bacteria_name, scenario_name, sampled_lag_combo, draw_idx, "memory", sep = "_")
      )

      draw_parameter_rows[[draw_parameter_idx]] <- tibble(
        Bacteria = bacteria_name,
        scenario = scenario_name,
        draw = draw_idx,
        Lag_Combination = sampled_lag_combo,
        lag_sampling_weight = lag_candidate_table$lag_sampling_weight[match(sampled_lag_combo, lag_candidate_table$Lag_Combination)],
        memory_factor_multiplier = sampled_memory_multiplier,
        weight_temperature = unname(sampled_weights_vec["Temperature"]),
        weight_humidity = unname(sampled_weights_vec["Humidity"]),
        weight_precipitation = unname(sampled_weights_vec["Precipitation"]),
        weight_wetdays = unname(sampled_weights_vec["Wet Days"])
      )
      draw_parameter_idx <- draw_parameter_idx + 1L

      if (isTRUE(model_settings$climate_model_uncertainty$bootstrap_resample_models)) {
        set.seed(create_seed_from_string(paste(bacteria_name, scenario_name, draw_idx, "gcm_bootstrap", sep = "_")))
        model_resample_idx <- sample(seq_along(climate_models), size = length(climate_models), replace = TRUE)
      } else {
        model_resample_idx <- seq_along(climate_models)
      }

      model_matrix <- matrix(NA_real_, nrow = n_years, ncol = length(climate_models))
      zone_model_arrays <- setNames(
        lapply(zone_names, function(x) matrix(NA_real_, nrow = n_years, ncol = length(climate_models))),
        zone_names
      )

      for (model_idx in seq_along(climate_models)) {
        climate_model_name <- climate_models[[model_idx]]
        model_configs <- selected_configs[vapply(selected_configs, function(x) identical(x$climate_model, climate_model_name), logical(1))]
        if (length(model_configs) == 0) {
          next
        }

        zone_matrix <- matrix(NA_real_, nrow = n_years, ncol = length(model_configs))
        zone_weights <- vapply(model_configs, function(x) x$zone_weight %||% 1, numeric(1))

        for (zone_idx in seq_along(model_configs)) {
          cfg <- model_configs[[zone_idx]]
          combined_effect_draw <- sample_combined_effect_path(cfg, draw_idx, sampled_weights_vec = sampled_weights_vec)
          zone_path <- simulate_recursive_path(
            years = cfg$years,
            combined_effect = combined_effect_draw,
            baseline_resistance = cfg$baseline_resistance,
            bacteria_name = cfg$bacteria,
            scenario_name = cfg$scenario,
            seed_prefix = paste(cfg$bacteria, cfg$scenario, cfg$climate_model, cfg$climate_zone, draw_idx, sep = "_"),
            initial_variation = cfg$initial_variation,
            memory_factor_multiplier = sampled_memory_multiplier,
            recursive_persistence = cfg$recursive_persistence,
            recursive_innovation_sd = cfg$recursive_innovation_sd
          )
          zone_matrix[, zone_idx] <- zone_path
          zone_model_arrays[[cfg$climate_zone]][, model_idx] <- zone_path
        }

        model_mean <- vapply(
          seq_len(n_years),
          function(i) weighted_mean_safe(zone_matrix[i, ], zone_weights),
          numeric(1)
        )
        model_mean[is.nan(model_mean)] <- NA_real_
        model_matrix[, model_idx] <- model_mean
      }

      draw_mean <- rowMeans(model_matrix[, model_resample_idx, drop = FALSE], na.rm = TRUE)
      draw_mean[is.nan(draw_mean)] <- NA_real_
      draw_matrix[, draw_idx] <- draw_mean
      if (collect_zone_summary) {
        for (zone_name in zone_names) {
          zone_draw <- rowMeans(zone_model_arrays[[zone_name]][, model_resample_idx, drop = FALSE], na.rm = TRUE)
          zone_draw[is.nan(zone_draw)] <- NA_real_
          zone_draw_arrays[[zone_name]][, draw_idx] <- zone_draw
        }
      }
    }

    annual_rows[[annual_idx]] <- tibble(
      Bacteria = bacteria_name,
      scenario = scenario_name,
      year = years,
      resistance_rate_pct = rowMeans(draw_matrix, na.rm = TRUE),
      lower_95ci_pct = apply(draw_matrix, 1, quantile, probs = lower_prob, na.rm = TRUE),
      upper_95ci_pct = apply(draw_matrix, 1, quantile, probs = upper_prob, na.rm = TRUE),
      n_draws = monte_carlo_settings$n_simulations
    )
    annual_idx <- annual_idx + 1L

    for (period_name in names(assessment_periods)) {
      period_years <- assessment_periods[[period_name]]
      year_idx <- years %in% period_years
      period_draw_means <- colMeans(draw_matrix[year_idx, , drop = FALSE], na.rm = TRUE)
      relative_change_draws <- 100 * (period_draw_means / baseline_global - 1)

      period_rows[[period_idx]] <- tibble(
        Bacteria = bacteria_name,
        scenario = scenario_name,
        period = period_name,
        relative_change_pct = mean(relative_change_draws, na.rm = TRUE),
        lower_95ci_pct = quantile(relative_change_draws, probs = lower_prob, na.rm = TRUE),
        upper_95ci_pct = quantile(relative_change_draws, probs = upper_prob, na.rm = TRUE),
        n_years = sum(year_idx),
        n_draws = monte_carlo_settings$n_simulations
      )
      period_idx <- period_idx + 1L
    }

    if (collect_zone_summary) {
      for (zone_name in zone_names) {
        zone_draw_matrix <- zone_draw_arrays[[zone_name]]
        zone_annual_rows[[zone_annual_idx]] <- tibble(
          Bacteria = bacteria_name,
          scenario = scenario_name,
          climate_zone = zone_name,
          year = years,
          resistance_rate_pct = rowMeans(zone_draw_matrix, na.rm = TRUE),
          lower_95ci_pct = apply(zone_draw_matrix, 1, quantile, probs = lower_prob, na.rm = TRUE),
          upper_95ci_pct = apply(zone_draw_matrix, 1, quantile, probs = upper_prob, na.rm = TRUE),
          n_draws = monte_carlo_settings$n_simulations
        )
        zone_annual_idx <- zone_annual_idx + 1L

        for (period_name in names(assessment_periods)) {
          period_years <- assessment_periods[[period_name]]
          year_idx <- years %in% period_years
          period_draw_means <- colMeans(zone_draw_matrix[year_idx, , drop = FALSE], na.rm = TRUE)
          zone_baseline <- zone_baseline_lookup[[zone_name]]
          relative_change_draws <- 100 * (period_draw_means / zone_baseline - 1)

          zone_period_rows[[zone_period_idx]] <- tibble(
            Bacteria = bacteria_name,
            scenario = scenario_name,
            climate_zone = zone_name,
            period = period_name,
            relative_change_pct = mean(relative_change_draws, na.rm = TRUE),
            lower_95ci_pct = quantile(relative_change_draws, probs = lower_prob, na.rm = TRUE),
            upper_95ci_pct = quantile(relative_change_draws, probs = upper_prob, na.rm = TRUE),
            n_years = sum(year_idx),
            n_draws = monte_carlo_settings$n_simulations
          )
          zone_period_idx <- zone_period_idx + 1L
        }
      }
    }
  }

  annual_summary <- bind_rows(annual_rows) %>%
    left_join(climate_inputs$model_coverage, by = "scenario") %>%
    mutate(
      resistance_rate_prop = resistance_rate_pct / 100,
      lower_95ci_prop = lower_95ci_pct / 100,
      upper_95ci_prop = upper_95ci_pct / 100,
      Scenario_Label = unname(scenario_labels[scenario]),
      Scenario_Full_Label = unname(scenario_full_labels[scenario]),
      Bacteria = factor(Bacteria, levels = bacteria_levels)
    ) %>%
    arrange(Bacteria, scenario, year)

  period_summary <- bind_rows(period_rows) %>%
    mutate(
      Scenario_Label = unname(scenario_labels[scenario]),
      Scenario_Full_Label = unname(scenario_full_labels[scenario]),
      Bacteria = factor(Bacteria, levels = bacteria_levels),
      period = factor(period, levels = names(assessment_periods))
    ) %>%
    arrange(Bacteria, scenario, period)

  zone_annual_summary <- bind_rows(zone_annual_rows) %>%
    left_join(climate_inputs$model_coverage, by = "scenario") %>%
    mutate(
      Bacteria = factor(Bacteria, levels = bacteria_levels),
      climate_zone = factor(climate_zone, levels = climate_zone_levels),
      Scenario_Label = unname(scenario_labels[scenario])
    ) %>%
    arrange(climate_zone, Bacteria, scenario, year)

  zone_period_summary <- bind_rows(zone_period_rows) %>%
    mutate(
      Bacteria = factor(Bacteria, levels = bacteria_levels),
      climate_zone = factor(climate_zone, levels = heatmap_zone_levels),
      period = factor(period, levels = names(assessment_periods)),
      Scenario_Label = unname(scenario_labels[scenario])
    ) %>%
    arrange(Bacteria, period, scenario, climate_zone)

  list(
    annual_summary = annual_summary,
    period_summary = period_summary,
    zone_annual_summary = zone_annual_summary,
    zone_period_summary = zone_period_summary,
    draw_parameter_summary = bind_rows(draw_parameter_rows),
    uncertainty_mode = current_uncertainty_mode()
  )
}

build_endpoint_labels <- function(annual_summary) {
  annual_summary %>%
    group_by(Bacteria, scenario) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    transmute(
      Bacteria,
      scenario,
      Scenario_Label,
      Scenario_Full_Label,
      year,
      resistance_rate_pct,
      lower_95ci_pct,
      upper_95ci_pct,
      endpoint_label = sprintf("%.2f", round(resistance_rate_pct, 2)),
      endpoint_label_with_ci = sprintf("%.2f [%.2f, %.2f]", resistance_rate_pct, lower_95ci_pct, upper_95ci_pct)
    )
}

build_gcm_spread_summaries <- function(global_model_predictions) {
  model_level <- global_model_predictions %>%
    filter(climate_model != "Ensemble Mean")

  annual_gcm_summary <- model_level %>%
    group_by(bacteria, scenario, year) %>%
    summarise(
      gcm_mean_resistance_pct = mean(resistance_rate, na.rm = TRUE),
      gcm_min_resistance_pct = min(resistance_rate, na.rm = TRUE),
      gcm_max_resistance_pct = max(resistance_rate, na.rm = TRUE),
      gcm_sd_resistance_pct = sd(resistance_rate, na.rm = TRUE),
      n_climate_models = n_distinct(climate_model),
      .groups = "drop"
    ) %>%
    transmute(
      Bacteria = factor(bacteria, levels = bacteria_levels),
      scenario,
      year,
      gcm_mean_resistance_pct,
      gcm_min_resistance_pct,
      gcm_max_resistance_pct,
      gcm_sd_resistance_pct,
      gcm_range_width_pct = gcm_max_resistance_pct - gcm_min_resistance_pct,
      n_climate_models,
      Scenario_Label = unname(scenario_labels[scenario]),
      Scenario_Full_Label = unname(scenario_full_labels[scenario])
    ) %>%
    arrange(Bacteria, scenario, year)

  period_gcm_summary <- model_level %>%
    mutate(period = purrr::map_chr(year, function(y) {
      hit <- names(assessment_periods)[vapply(assessment_periods, function(yrs) y %in% yrs, logical(1))]
      if (length(hit) == 0) {
        return(NA_character_)
      }
      hit[[1]]
    })) %>%
    filter(!is.na(period)) %>%
    group_by(bacteria, scenario, climate_model, period) %>%
    summarise(
      period_mean_change_pct = mean(relative_change_baseline, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(bacteria, scenario, period) %>%
    summarise(
      gcm_mean_change_pct = mean(period_mean_change_pct, na.rm = TRUE),
      gcm_min_change_pct = min(period_mean_change_pct, na.rm = TRUE),
      gcm_max_change_pct = max(period_mean_change_pct, na.rm = TRUE),
      gcm_sd_change_pct = sd(period_mean_change_pct, na.rm = TRUE),
      n_climate_models = n_distinct(climate_model),
      .groups = "drop"
    ) %>%
    transmute(
      Bacteria = factor(bacteria, levels = bacteria_levels),
      scenario,
      period = factor(period, levels = names(assessment_periods)),
      gcm_mean_change_pct,
      gcm_min_change_pct,
      gcm_max_change_pct,
      gcm_sd_change_pct,
      gcm_range_width_pct = gcm_max_change_pct - gcm_min_change_pct,
      n_climate_models,
      Scenario_Label = unname(scenario_labels[scenario]),
      Scenario_Full_Label = unname(scenario_full_labels[scenario])
    ) %>%
    arrange(Bacteria, scenario, period)

  list(
    annual = annual_gcm_summary,
    period = period_gcm_summary
  )
}

plot_resistance_time_trends <- function(annual_summary) {
  settings <- model_settings$plot_settings
  plot_data <- annual_summary %>%
    mutate(
      scenario = factor(scenario, levels = scenario_levels),
      Bacteria = factor(Bacteria, levels = bacteria_levels)
    )

  final_values <- plot_data %>%
    group_by(Bacteria, scenario) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    mutate(label = sprintf("%.2f", round(resistance_rate_pct, 2)), x_position = year + settings$end_label_nudge)

  ggplot(plot_data, aes(x = year, y = resistance_rate_pct, color = scenario, group = scenario)) +
    geom_ribbon(aes(ymin = lower_95ci_pct, ymax = upper_95ci_pct, fill = scenario), alpha = settings$confidence_interval_alpha, color = NA) +
    geom_line(linewidth = settings$line_width, alpha = settings$line_alpha) +
    ggrepel::geom_text_repel(
      data = final_values,
      aes(x = x_position, y = resistance_rate_pct, label = label, color = scenario),
      size = settings$label_size,
      fontface = "bold",
      min.segment.length = 0,
      segment.color = "grey50",
      segment.size = 0.3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.15, "lines"),
      force = 1.8,
      nudge_x = 2.5,
      direction = "y",
      hjust = 0,
      max.overlaps = Inf,
      seed = monte_carlo_settings$seed,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_color_manual(values = scenario_colors, labels = scenario_labels, name = "Climate Scenario") +
    scale_fill_manual(values = scenario_colors, guide = "none") +
    scale_x_continuous(breaks = c(2020, 2040, 2060, 2080, 2100), limits = c(2020, 2110), expand = c(0.01, 0)) +
    coord_cartesian(clip = "off") +
    facet_wrap(~ Bacteria, scales = "free_y", ncol = 3) +
    labs(x = "Year", y = "Resistance Rate (%)") +
    theme_minimal() +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.title = element_text(size = settings$axis_title_size, face = "bold"),
      axis.text = element_text(size = settings$axis_text_size),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.line = element_line(color = "gray60", linewidth = settings$axis_line_width),
      axis.ticks = element_line(color = "gray60", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid.major = element_line(color = settings$grid_color, linewidth = settings$grid_size),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray60", fill = NA, linewidth = settings$border_size),
      panel.spacing = unit(1, "lines"),
      strip.background = element_rect(fill = "gray95", color = "gray60"),
      strip.text = element_text(size = settings$strip_text_size, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = settings$legend_title_size, face = "bold"),
      legend.text = element_text(size = settings$legend_text_size),
      plot.margin = margin(10, 40, 10, 10)
    )
}

plot_relative_changes_grid_with_errorbars <- function(period_summary) {
  settings <- model_settings$plot_settings
  plot_data <- period_summary %>%
    mutate(
      scenario = factor(scenario, levels = scenario_levels),
      Bacteria = factor(Bacteria, levels = bacteria_levels)
    )

  label_data <- plot_data %>%
    mutate(label = sprintf("%.2f", round(relative_change_pct, 2)), label_vjust = ifelse(relative_change_pct >= 0, -0.4, 1.5))

  ggplot(plot_data, aes(x = scenario, y = relative_change_pct, fill = period)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "gray40", linewidth = 0.2) +
    geom_errorbar(
      aes(ymin = lower_95ci_pct, ymax = upper_95ci_pct),
      position = position_dodge(width = 0.8),
      width = settings$errbar_width,
      linewidth = settings$errbar_linewidth,
      alpha = settings$errbar_alpha,
      color = "gray50"
    ) +
    geom_text(
      data = label_data,
      aes(label = label, vjust = label_vjust),
      position = position_dodge(width = 0.8),
      size = settings$annotation_text_size * 0.8,
      color = "black",
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_x_discrete(labels = function(x) scenario_labels[x]) +
    scale_fill_manual(values = period_colors, name = "Time Period") +
    facet_wrap(~ Bacteria, scales = "free_y", ncol = 3) +
    labs(x = "Climate Scenario", y = "Relative Change (%)") +
    theme_minimal() +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.title = element_text(size = settings$axis_title_size, face = "bold"),
      axis.text = element_text(size = settings$axis_text_size),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.line = element_line(color = "gray60", linewidth = settings$axis_line_width),
      axis.ticks = element_line(color = "gray60", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid.major.y = element_line(color = settings$grid_color, linewidth = settings$grid_size),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray60", fill = NA, linewidth = settings$border_size),
      panel.spacing = unit(1, "lines"),
      strip.background = element_rect(fill = "gray95", color = "gray60"),
      strip.text = element_text(size = settings$strip_text_size, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = settings$legend_title_size, face = "bold"),
      legend.text = element_text(size = settings$legend_text_size),
      plot.margin = margin(10, 10, 10, 10)
    )
}

plot_climate_zone_bar_panels <- function(zone_period_summary) {
  settings <- model_settings$plot_settings
  plot_data <- zone_period_summary %>%
    mutate(
      Bacteria = factor(Bacteria, levels = bacteria_levels),
      climate_zone = factor(climate_zone, levels = heatmap_zone_levels),
      scenario_label = factor(unname(scenario_labels[scenario]), levels = unname(scenario_labels[scenario_levels])),
      period = factor(period, levels = names(assessment_periods))
    )

  ggplot(plot_data, aes(x = climate_zone, y = relative_change_pct, fill = scenario)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_col(
      position = position_dodge(width = 0.82),
      width = 0.72,
      color = "gray35",
      linewidth = 0.2
    ) +
    geom_errorbar(
      aes(ymin = lower_95ci_pct, ymax = upper_95ci_pct),
      position = position_dodge(width = 0.82),
      width = 0.22,
      linewidth = settings$errbar_linewidth * 0.8,
      alpha = 0.9,
      color = "gray40"
    ) +
    coord_flip() +
    scale_fill_manual(values = scenario_colors, labels = scenario_labels, name = "Climate Scenario") +
    facet_grid(period ~ Bacteria, scales = "free_y") +
    labs(x = NULL, y = "Relative Change (%)") +
    theme_minimal() +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.title = element_text(size = settings$axis_title_size, face = "bold"),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray93", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "gray96", color = "gray70"),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = settings$legend_title_size, face = "bold"),
      legend.text = element_text(size = settings$legend_text_size),
      plot.margin = margin(10, 10, 10, 10)
    )
}

build_parameter_manifest <- function(model_artifacts, projection_artifacts) {
  lag_counts <- projection_artifacts$lag_candidates %>%
    group_by(Bacteria) %>%
    summarise(
      n_lag_candidates_used = n_distinct(Lag_Combination),
      min_delta_aic = min(delta_aic, na.rm = TRUE),
      max_delta_aic = max(delta_aic, na.rm = TRUE),
      .groups = "drop"
    )

  recursive_uncertainty_manifest <- model_artifacts$recursive_uncertainty_df %>%
    select(
      Bacteria,
      calibrated_persistence,
      calibrated_innovation_sd,
      calibrated_initial_sd,
      empirical_variation_factor,
      calibration_source
    )

  tibble(Bacteria = bacteria_levels) %>%
    left_join(lag_counts, by = "Bacteria") %>%
    left_join(recursive_uncertainty_manifest, by = "Bacteria") %>%
    mutate(
      diffusion_base_rate = vapply(Bacteria, diffusion_base_rate_for_bacteria, numeric(1)),
      uses_temperature_tail_constraint = Bacteria %in% names(model_settings$temperature_tail_constraint$directions),
      temperature_tail_constraint_direction = unname(model_settings$temperature_tail_constraint$directions[Bacteria]) %||% "none",
      uses_crkp_special_controls = Bacteria == "CR-Kp",
      response_or_lower = model_settings$response_or_bounds[1],
      response_or_upper = model_settings$response_or_bounds[2],
      lag_uncertainty_included = !is.na(n_lag_candidates_used) & n_lag_candidates_used > 1,
      methods_alignment_flag = case_when(
        Bacteria == "CR-Kp" ~ "Needs explicit methods description for phenotype-specific stabilization controls and bespoke diffusion settings",
        uses_temperature_tail_constraint ~ "Needs explicit methods description for phenotype-specific temperature-direction constraints and tail log-SE stabilization",
        lag_uncertainty_included ~ "Needs explicit methods description for Akaike-weighted lag-uncertainty propagation",
        TRUE ~ "Core projection settings align with shared simplified-model framework"
      )
    ) %>%
    select(
      Bacteria,
      diffusion_base_rate,
      calibrated_persistence,
      calibrated_innovation_sd,
      calibrated_initial_sd,
      empirical_variation_factor,
      calibration_source,
      uses_temperature_tail_constraint,
      temperature_tail_constraint_direction,
      uses_crkp_special_controls,
      n_lag_candidates_used,
      min_delta_aic,
      max_delta_aic,
      lag_uncertainty_included,
      response_or_lower,
      response_or_upper,
      methods_alignment_flag
    )
}

build_method_alignment_audit <- function(parameter_manifest) {
  tibble(
    item = c(
      "weights_from_simplified_model_c",
      "zone_specific_baseline_anchoring",
      "historical_support_clamping_and_linear_tail",
      "ssp_linked_diffusion_modifier",
      "lag_uncertainty_propagation",
      "historical_residual_calibrated_recursive_uncertainty",
      "crkp_special_stabilization_controls",
      "directional_temperature_tail_constraints",
      "response_or_bounds",
      "support_edge_log_se_regularization"
    ),
    aligned_with_original_main_methods = c(
      "Yes",
      "Yes",
      "Yes",
      "Partial",
      "Partial",
      "Partial",
      "Partial",
      "No",
      "Partial",
      "Partial"
    ),
    note = c(
      "Consistent with projection strategy based on simplified Model C smooth-term statistics.",
      "Consistent with the projection workflow using climate-zone-specific baselines.",
      "Consistent with the projection-safe response lookup strategy already used in the revised workflow.",
      "Scientifically defensible, but should be described as scenario parameterization rather than data-estimated socioeconomic effect.",
      "Now propagated by Akaike-weighted sampling of near-optimal lag combinations; this is stronger than the original best-lag-only description and must be added to Methods documentation.",
      "Recursive stochasticity is now calibrated from historical observed-versus-fitted residual series and should be described explicitly as an empirical uncertainty calibration step.",
      "CR-Kp-specific smoothing and growth-limit controls are reasonable as a phenotype-specific stabilization treatment, but must be transparently described in Methods documentation.",
      "Temperature-direction constraints are now applied only to phenotypes with independent Figure 2B support (CR-Ab increasing, 3GCR-Kp increasing, CR-Ec decreasing) and must be explicitly stated if retained in the main result.",
      "Global OR clipping improves robustness but should be described as a stricter numerical guardrail in Methods documentation.",
      "Support-edge log-SE regularization reduces boundary spikes without altering central response estimates and should be documented as a projection-stability safeguard."
    )
  ) %>%
    mutate(
      high_priority = item %in% c(
        "lag_uncertainty_propagation",
        "historical_residual_calibrated_recursive_uncertainty",
        "crkp_special_stabilization_controls",
        "directional_temperature_tail_constraints"
      )
    )
}

build_uncertainty_audit <- function(annual_summary, period_summary, lookup_df, lag_candidates_df) {
  annual_2100 <- annual_summary %>%
    filter(year == 2100) %>%
    mutate(annual_ci_width_2100 = upper_95ci_pct - lower_95ci_pct) %>%
    select(Bacteria, scenario, annual_ci_width_2100)

  period_2090s <- period_summary %>%
    filter(period == "2090s") %>%
    mutate(period_ci_width_2090s = upper_95ci_pct - lower_95ci_pct) %>%
    select(Bacteria, scenario, period_ci_width_2090s)

  temp_lookup <- lookup_df %>%
    filter(Climate_Factor == "Temperature") %>%
    group_by(Bacteria) %>%
    summarise(
      temperature_curve_mean_log_se = mean(log_se, na.rm = TRUE),
      temperature_curve_max_log_se = max(log_se, na.rm = TRUE),
      temperature_curve_mean_ci_width = mean(Upper_CI - Lower_CI, na.rm = TRUE),
      .groups = "drop"
    )

  log_se_regularization_audit <- lookup_df %>%
    group_by(Bacteria) %>%
    summarise(
      lookup_edge_stabilized_share = mean(projection_log_se_edge_stabilized, na.rm = TRUE),
      .groups = "drop"
    )

  lag_audit <- lag_candidates_df %>%
    group_by(Bacteria) %>%
    summarise(
      n_lag_candidates_used = n_distinct(Lag_Combination),
      n_lag_candidates_delta_aic_le_2 = sum(delta_aic <= 2, na.rm = TRUE),
      .groups = "drop"
    )

  annual_2100 %>%
    left_join(period_2090s, by = c("Bacteria", "scenario")) %>%
    left_join(temp_lookup, by = "Bacteria") %>%
    left_join(log_se_regularization_audit, by = "Bacteria") %>%
    left_join(lag_audit, by = "Bacteria") %>%
    arrange(Bacteria, scenario)
}

write_result_tables <- function(
  annual_summary,
  endpoint_labels,
  period_summary,
  zone_annual_summary,
  zone_period_summary,
  gcm_annual_summary,
  gcm_period_summary,
  central_predictions,
  model_artifacts,
  climate_inputs,
  projection_artifacts
) {
  parameter_manifest <- build_parameter_manifest(model_artifacts, projection_artifacts)
  method_alignment_audit <- build_method_alignment_audit(parameter_manifest)
  uncertainty_audit <- build_uncertainty_audit(
    annual_summary = annual_summary,
    period_summary = period_summary,
    lookup_df = model_artifacts$curve_lookup_df,
    lag_candidates_df = projection_artifacts$lag_candidates
  )

  annual_csv <- file.path(results_root, "01_tables", "Figure3A_annual_global_trajectories.csv")
  endpoint_csv <- file.path(results_root, "01_tables", "Figure3A_2100_endpoint_labels.csv")
  period_csv <- file.path(results_root, "01_tables", "Figure3B_period_relative_changes.csv")
  zone_annual_csv <- file.path(results_root, "01_tables", "climate_zone_annual_trajectories.csv")
  zone_period_csv <- file.path(results_root, "01_tables", "climate_zone_period_relative_changes.csv")
  gcm_annual_csv <- file.path(results_root, "01_tables", "gcm_spread_annual_trajectories.csv")
  gcm_period_csv <- file.path(results_root, "01_tables", "gcm_spread_period_changes.csv")

  write.csv(annual_summary, annual_csv, row.names = FALSE)
  write.csv(endpoint_labels, endpoint_csv, row.names = FALSE)
  write.csv(period_summary, period_csv, row.names = FALSE)
  write.csv(zone_annual_summary, zone_annual_csv, row.names = FALSE)
  write.csv(zone_period_summary, zone_period_csv, row.names = FALSE)
  write.csv(gcm_annual_summary, gcm_annual_csv, row.names = FALSE)
  write.csv(gcm_period_summary, gcm_period_csv, row.names = FALSE)

  write.csv(bind_rows(
    central_predictions$global_model_predictions %>% mutate(trajectory_type = "central_model"),
    central_predictions$ensemble_predictions %>% mutate(trajectory_type = "central_ensemble")
  ), file.path(results_root, "03_diagnostics", "Figure3_global_projection_predictions.csv"), row.names = FALSE)
  write.csv(central_predictions$zone_predictions, file.path(results_root, "03_diagnostics", "Figure3_zone_projection_predictions.csv"), row.names = FALSE)
  write.csv(projection_artifacts$effect_components, file.path(results_root, "03_diagnostics", "projection_effect_components_by_zone.csv"), row.names = FALSE)
  write.csv(projection_artifacts$recursive_parameters, file.path(results_root, "03_diagnostics", "projection_recursive_parameters.csv"), row.names = FALSE)
  write.csv(model_artifacts$weights_df, file.path(results_root, "03_diagnostics", "projection_climate_weights_used.csv"), row.names = FALSE)
  write.csv(model_artifacts$candidate_weights_df, file.path(results_root, "03_diagnostics", "projection_candidate_climate_weights_used.csv"), row.names = FALSE)
  write.csv(model_artifacts$modifier_strength_df, file.path(results_root, "03_diagnostics", "projection_development_modifier_strengths.csv"), row.names = FALSE)
  write.csv(model_artifacts$term_stats_df, file.path(results_root, "03_diagnostics", "projection_model_term_statistics_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$weight_details_df, file.path(results_root, "03_diagnostics", "projection_climate_weight_derivation_details_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$curve_lookup_df, file.path(results_root, "03_diagnostics", "projection_climate_response_lookup_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$scaling_df, file.path(results_root, "03_diagnostics", "projection_climate_scaling_parameters_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$support_df, file.path(results_root, "03_diagnostics", "projection_climate_support_ranges_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$baseline_zone_df, file.path(results_root, "03_diagnostics", "projection_baseline_by_zone_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$baseline_global_df, file.path(results_root, "03_diagnostics", "projection_baseline_global_refined.csv"), row.names = FALSE)
  write.csv(climate_inputs$zone_year_model_means, file.path(results_root, "03_diagnostics", "projection_climate_zone_year_model_means.csv"), row.names = FALSE)
  write.csv(model_artifacts$historical_global_df, file.path(results_root, "03_diagnostics", "projection_historical_fit_global_year_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$historical_zone_df, file.path(results_root, "03_diagnostics", "projection_historical_fit_zone_year_refined.csv"), row.names = FALSE)
  write.csv(model_artifacts$recursive_uncertainty_df, file.path(results_root, "03_diagnostics", "projection_recursive_uncertainty_calibration.csv"), row.names = FALSE)
  write.csv(model_artifacts$historical_global_residual_df, file.path(results_root, "03_diagnostics", "projection_historical_global_residual_series.csv"), row.names = FALSE)
  write.csv(model_artifacts$historical_zone_residual_df, file.path(results_root, "03_diagnostics", "projection_historical_zone_residual_series.csv"), row.names = FALSE)
  write.csv(model_artifacts$recursive_uncertainty_zone_summary_df, file.path(results_root, "03_diagnostics", "projection_recursive_uncertainty_zone_summary.csv"), row.names = FALSE)
  write.csv(model_artifacts$manifest_df, file.path(results_root, "05_metadata", "projection_model_manifest_refined.csv"), row.names = FALSE)
  write.csv(projection_artifacts$lag_candidates, file.path(results_root, "03_diagnostics", "projection_lag_uncertainty_catalog.csv"), row.names = FALSE)
  write.csv(parameter_manifest, file.path(results_root, "03_diagnostics", "projection_parameter_manifest_by_bacteria.csv"), row.names = FALSE)
  write.csv(method_alignment_audit, file.path(results_root, "05_metadata", "projection_method_alignment_audit.csv"), row.names = FALSE)
  write.csv(uncertainty_audit, file.path(results_root, "03_diagnostics", "projection_uncertainty_audit_by_bacteria.csv"), row.names = FALSE)
  if (!is.null(projection_artifacts$draw_parameter_summary)) {
    write.csv(projection_artifacts$draw_parameter_summary, file.path(results_root, "03_diagnostics", "projection_uncertainty_draw_parameters.csv"), row.names = FALSE)
  }

  metadata_df <- tibble(
    item = c(
      "projection_version",
      "projection_mode",
      "climate_zone_rule",
      "weight_derivation",
      "response_lookup_scaling",
      "response_lookup_anchor",
      "response_lookup_support",
      "response_lookup_extrapolation",
      "temperature_tail_constraint",
      "recursive_uncertainty_calibration",
      "support_edge_log_se_regularization",
      "response_or_bounds",
      "future_climate_models",
      "uncertainty_mode",
      "gcm_spread_export",
      "n_monte_carlo_simulations",
      "outer_loop_weight_uncertainty",
      "outer_loop_memory_uncertainty",
      "lag_uncertainty",
      "development_modifier",
      "modifier_calibration",
      "direct_antibiotic_factor",
      "ssp_diffusion_modifier",
      "recursive_uncertainty_export",
      "parameter_manifest_export",
      "method_alignment_audit_export",
      "uncertainty_audit_export"
    ),
    value = c(
      "Anchored projection",
      "Climate-response projection with phenotype-specific temperature-direction constraints, lag-selection uncertainty, empirical recursive-residual uncertainty calibration, and four-GCM bootstrap uncertainty",
      "Polar if |lat| > 66.5; Temperate if 23.5 < |lat| <= 66.5; Tropical otherwise",
      "Weights derived from refined simplified GAMM climate terms using F_stat x EDF x significance factor, then bounded to [0.05, 0.70] and normalized",
      "Global scaling across each phenotype-specific model-ready dataset",
      "Zone-specific baseline anchoring: Effect_OR = exp(term(x) - term(baseline climate in 2010-2019 for that zone))",
      paste0("Exposure support capped to empirical ", paste0(model_settings$historical_support_quantiles * 100, collapse = "th-"), "th percentiles"),
      "Future exposures outside historical support use linear tail extrapolation of the response lookup, while log-SE remains support-clamped and support-edge regularized",
      "Monotonic directional constraints applied to temperature response above the 50th percentile of support (and above zone baseline) for CR-Ab increasing, 3GCR-Kp increasing, and CR-Ec decreasing, with tail log-SE stabilized to the local boundary neighborhood",
      "Recursive stochasticity calibrated from historical observed-versus-fitted residual series using phenotype-specific AR(1)-style persistence and innovation SD estimates, combined across global and climate-zone annual fits",
      "All climate-response lookup curves use support-edge log-SE regularization to cap boundary spikes against the adjacent interior median without altering the central fitted effect",
      paste(model_settings$response_or_bounds, collapse = ", "),
      paste(climate_models_expected, collapse = ", "),
      if (isTRUE(model_settings$climate_model_uncertainty$bootstrap_resample_models)) {
        "95% CI from four-GCM spread, empirically calibrated recursive stochasticity, and response-function parameter uncertainty"
      } else {
        "95% CI from empirically calibrated recursive stochasticity and response-function parameter uncertainty; central four-GCM spread exported separately"
      },
      as.character(model_settings$climate_model_uncertainty$export_central_model_spread),
      as.character(monte_carlo_settings$n_simulations),
      as.character(model_settings$outer_loop_uncertainty$include_weight_perturbation),
      as.character(model_settings$outer_loop_uncertainty$include_memory_perturbation),
      "Lag uncertainty propagated by sampling simplified-model lag combinations with delta AIC <= 2.0 (up to four candidates per bacteria) using normalized Akaike weights",
      as.character(model_settings$use_development_modifier),
      "Direct development modifier retained in codebase but disabled for this projection",
      as.character(model_settings$use_antibiotic_modifier),
      "TRUE; applied multiplicatively to diffusion_rate with scenario-specific net annual rates and ramp-up bounds",
      "TRUE; phenotype-specific residual calibration exported to projection_recursive_uncertainty_calibration.csv",
      "TRUE; phenotype-specific parameters exported to projection_parameter_manifest_by_bacteria.csv",
      "TRUE; manuscript-method consistency checkpoints exported to projection_method_alignment_audit.csv",
      "TRUE; CI-width and lag-ambiguity diagnostics exported to projection_uncertainty_audit_by_bacteria.csv"
    )
  )
  write.csv(metadata_df, file.path(results_root, "05_metadata", paste0("Figure3_", tolower(projection_version_tag), "_metadata.csv")), row.names = FALSE)

  source_csv_map <- list(
    Figure3_ModelC_Fig3A_annual_trajectories = annual_summary,
    Figure3_ModelC_Fig3A_2100_labels = endpoint_labels,
    Figure3_ModelC_Fig3B_period_changes = period_summary,
    climate_zone_annual_trajectories = zone_annual_summary,
    climate_zone_period_bars = zone_period_summary,
    gcm_spread_annual_trajectories = gcm_annual_summary,
    gcm_spread_period_changes = gcm_period_summary,
    Figure3_ModelC_metadata = metadata_df,
    Figure3_ModelC_recursive_uncertainty_calibration = model_artifacts$recursive_uncertainty_df,
    Figure3_ModelC_parameter_manifest = parameter_manifest,
    Figure3_ModelC_method_alignment_audit = method_alignment_audit,
    Figure3_ModelC_uncertainty_audit = uncertainty_audit
  )

  purrr::iwalk(source_csv_map, function(df, nm) {
    write.csv(df, file.path(source_data_root, "01_csv", paste0(nm, ".csv")), row.names = FALSE)
  })

  wb <- createWorkbook()
  addWorksheet(wb, "Figure3A_Annual")
  writeData(wb, "Figure3A_Annual", annual_summary)
  addWorksheet(wb, "Figure3A_2100")
  writeData(wb, "Figure3A_2100", endpoint_labels)
  addWorksheet(wb, "Figure3B_Period")
  writeData(wb, "Figure3B_Period", period_summary)
  addWorksheet(wb, "Zone_Annual")
  writeData(wb, "Zone_Annual", zone_annual_summary)
  addWorksheet(wb, "Zone_Bars")
  writeData(wb, "Zone_Bars", zone_period_summary)
  addWorksheet(wb, "GCM_Annual_Spread")
  writeData(wb, "GCM_Annual_Spread", gcm_annual_summary)
  addWorksheet(wb, "GCM_Period_Spread")
  writeData(wb, "GCM_Period_Spread", gcm_period_summary)
  addWorksheet(wb, "Weights")
  writeData(wb, "Weights", model_artifacts$weights_df)
  addWorksheet(wb, "Candidate_Weights")
  writeData(wb, "Candidate_Weights", model_artifacts$candidate_weights_df)
  addWorksheet(wb, "Modifier_Strength")
  writeData(wb, "Modifier_Strength", model_artifacts$modifier_strength_df)
  addWorksheet(wb, "Recursive_UQ")
  writeData(wb, "Recursive_UQ", model_artifacts$recursive_uncertainty_df)
  addWorksheet(wb, "Weight_Details")
  writeData(wb, "Weight_Details", model_artifacts$weight_details_df)
  addWorksheet(wb, "Lag_Candidates")
  writeData(wb, "Lag_Candidates", projection_artifacts$lag_candidates)
  addWorksheet(wb, "Parameter_Manifest")
  writeData(wb, "Parameter_Manifest", parameter_manifest)
  addWorksheet(wb, "Method_Audit")
  writeData(wb, "Method_Audit", method_alignment_audit)
  addWorksheet(wb, "Uncertainty_Audit")
  writeData(wb, "Uncertainty_Audit", uncertainty_audit)
  addWorksheet(wb, "Metadata")
  writeData(wb, "Metadata", metadata_df)

  saveWorkbook(
    wb,
    file.path(results_root, "04_workbook", paste0("Figure3_ModelC_Projection_", projection_version_tag, "_source_tables.xlsx")),
    overwrite = TRUE
  )
  saveWorkbook(
    wb,
    file.path(source_data_root, "02_workbook", paste0("SourceData_Figure3_ModelC_Projection_", projection_version_tag, ".xlsx")),
    overwrite = TRUE
  )
}

lag_candidates_df <- load_lag_uncertainty_catalog(lag_uncertainty_results_path)
future_climate_data <- ensure_future_climate_input(climate_data_path, climate_manifest_path)
model_artifacts <- fit_projection_models(lag_candidates_df)
climate_inputs <- prepare_climate_scenarios(future_climate_data)
projection_artifacts <- prepare_projection_configs(model_artifacts, lag_candidates_df, climate_inputs)
primary_configs <- projection_artifacts$configs[vapply(projection_artifacts$configs, function(x) isTRUE(x$is_primary), logical(1))]
central_predictions <- build_central_prediction_diagnostics(primary_configs, model_artifacts$baseline_global_df)
gcm_spread_outputs <- build_gcm_spread_summaries(central_predictions$global_model_predictions)
mc_outputs <- run_projection_monte_carlo(
  configs = projection_artifacts$configs,
  baseline_global_df = model_artifacts$baseline_global_df,
  climate_inputs = climate_inputs,
  collect_zone_summary = TRUE
)

endpoint_labels <- build_endpoint_labels(mc_outputs$annual_summary)

fig3a <- plot_resistance_time_trends(mc_outputs$annual_summary)
fig3b <- plot_relative_changes_grid_with_errorbars(mc_outputs$period_summary)
zone_barplot <- plot_climate_zone_bar_panels(mc_outputs$zone_period_summary)
combined_figure <- fig3a / fig3b + plot_layout(heights = c(1.15, 1.0))

ggsave(
  file.path(results_root, "02_figures", "Figure3A_annual_global_trajectories.pdf"),
  fig3a,
  width = model_settings$plot_settings$fig_width,
  height = model_settings$plot_settings$fig_height,
  dpi = model_settings$plot_settings$fig_dpi,
  device = cairo_pdf
)
ggsave(
  file.path(results_root, "02_figures", "Figure3B_period_relative_changes.pdf"),
  fig3b,
  width = model_settings$plot_settings$fig_width,
  height = model_settings$plot_settings$fig_height,
  dpi = model_settings$plot_settings$fig_dpi,
  device = cairo_pdf
)
ggsave(
  file.path(results_root, "02_figures", "Figure3_combined.pdf"),
  combined_figure,
  width = 12,
  height = 14,
  dpi = model_settings$plot_settings$fig_dpi,
  device = cairo_pdf
)
ggsave(
  file.path(results_root, "02_figures", "Figure3_combined.png"),
  combined_figure,
  width = 12,
  height = 14,
  dpi = model_settings$plot_settings$fig_dpi
)
ggsave(
  file.path(results_root, "02_figures", "climate_zone_period_bars.pdf"),
  zone_barplot,
  width = 13,
  height = 8.5,
  dpi = model_settings$plot_settings$fig_dpi,
  device = cairo_pdf
)
ggsave(
  file.path(results_root, "02_figures", "climate_zone_period_bars.png"),
  zone_barplot,
  width = 13,
  height = 8.5,
  dpi = model_settings$plot_settings$fig_dpi
)

write_result_tables(
  annual_summary = mc_outputs$annual_summary,
  endpoint_labels = endpoint_labels,
  period_summary = mc_outputs$period_summary,
  zone_annual_summary = mc_outputs$zone_annual_summary,
  zone_period_summary = mc_outputs$zone_period_summary,
  gcm_annual_summary = gcm_spread_outputs$annual,
  gcm_period_summary = gcm_spread_outputs$period,
  central_predictions = central_predictions,
  model_artifacts = model_artifacts,
  climate_inputs = climate_inputs,
  projection_artifacts = c(projection_artifacts, list(draw_parameter_summary = mc_outputs$draw_parameter_summary))
)

message(projection_version_tag, " Figure 3 projection pipeline completed: ", results_root)
