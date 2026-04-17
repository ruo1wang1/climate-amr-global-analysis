######## Model C Projection Preparation: Final Simplified GAMM Inputs ########

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(writexl)
})

set.seed(20260330)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
input_data_dir <- file.path(revision_root, "outputs/historical_associations", "model_ready_inputs")
projection_prep_root <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "projection_preparation"
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
weights_path <- file.path(
  revision_root,
  "data",
  "source_data",
  "lag_selection",
  "projection_simplified_model_c",
  "01_csv",
  "projection_simplified_climate_weights.csv"
)
weight_details_path <- file.path(
  revision_root,
  "data",
  "source_data",
  "lag_selection",
  "projection_simplified_model_c",
  "01_csv",
  "projection_simplified_climate_weight_details.csv"
)
output_base <- file.path(
  projection_prep_root,
  "03_projection_simplified_gamm_ready_inputs"
)

dir.create(file.path(output_base, "01_model_objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "02_curve_lookup_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "03_baseline_and_history"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "04_parameter_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "05_response_function_helper"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "06_workbook"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "07_metadata"), recursive = TRUE, showWarnings = FALSE)

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = "3GCREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = "3GCRKP_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRAB", title = "CR-Ab", file_name = "CRAB_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CREC", title = "CR-Ec", file_name = "CREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRKP", title = "CR-Kp", file_name = "CRKP_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRPA", title = "CR-Pa", file_name = "CRPA_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv")
)

climate_specs <- list(
  list(term = "TMP_scaled_lag", factor = "Temperature", lower_bound = NA_real_, upper_bound = NA_real_),
  list(term = "PREC_scaled_lag", factor = "Precipitation", lower_bound = 0, upper_bound = NA_real_),
  list(term = "HUM_scaled_lag", factor = "Humidity", lower_bound = 0, upper_bound = 100),
  list(term = "WET_scaled_lag", factor = "Wet Days", lower_bound = 0, upper_bound = NA_real_)
)

sanitize_file_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

normalize_climate_factor <- function(x) {
  case_when(
    x %in% c("Temperature", "TMP", "temp") ~ "Temperature",
    x %in% c("Precipitation", "PREC", "precip", "precipitation") ~ "Precipitation",
    x %in% c("Humidity", "Relative Humidity", "HUM", "humid") ~ "Humidity",
    x %in% c("Wet Days", "WetDays", "WET", "wetdays") ~ "Wet Days",
    TRUE ~ x
  )
}

load_final_lag_settings <- function(path) {
  if (!file.exists(path)) {
    stop("Final simplified lag settings file not found: ", path, call. = FALSE)
  }

  lag_df <- read.csv(path, stringsAsFactors = FALSE)
  required_cols <- c("Bacteria", "TMP_lag", "PREC_lag", "HUM_lag", "WET_lag")
  missing_cols <- setdiff(required_cols, names(lag_df))
  if (length(missing_cols) > 0) {
    stop("Lag settings file is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  lag_df
}

load_projection_weights <- function(path) {
  if (!file.exists(path)) {
    stop("Projection weight table not found: ", path, call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

prepare_projection_data <- function(file_path, lag_row) {
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
      climate_zone = factor(climate_zone),
      observed_R_prop = case_when(
        is.na(R) ~ NA_real_,
        R > 1 ~ R / 100,
        TRUE ~ R
      ),
      observed_R_pct = observed_R_prop * 100
    ) %>%
    group_by(NAME) %>%
    mutate(location_id = cur_group_id()) %>%
    ungroup() %>%
    mutate(
      HUM = pmin(HUM, 100),
      TMP_orig = TMP,
      PREC_orig = PREC,
      HUM_orig = HUM,
      WET_orig = WET
    )

  scaling_lookup <- data_processed %>%
    group_by(climate_zone) %>%
    summarise(
      TMP_mean = mean(TMP, na.rm = TRUE),
      TMP_sd = sd(TMP, na.rm = TRUE),
      PREC_mean = mean(PREC, na.rm = TRUE),
      PREC_sd = sd(PREC, na.rm = TRUE),
      HUM_mean = mean(HUM, na.rm = TRUE),
      HUM_sd = sd(HUM, na.rm = TRUE),
      WET_mean = mean(WET, na.rm = TRUE),
      WET_sd = sd(WET, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -climate_zone,
      names_to = c("variable", ".value"),
      names_pattern = "([A-Z]+)_(mean|sd)"
    ) %>%
    mutate(
      Climate_Factor = recode(
        variable,
        TMP = "Temperature",
        PREC = "Precipitation",
        HUM = "Humidity",
        WET = "Wet Days"
      )
    ) %>%
    select(climate_zone, Climate_Factor, mean_value = mean, sd_value = sd)

  lagged_data <- data_processed %>%
    group_by(climate_zone) %>%
    mutate(
      across(c(TMP, PREC, HUM, WET), \(x) as.vector(scale(x)), .names = "{.col}_scaled")
    ) %>%
    ungroup() %>%
    group_by(location_id) %>%
    arrange(year) %>%
    mutate(
      TMP_scaled_lag = lag(TMP_scaled, lag_row$TMP_lag),
      PREC_scaled_lag = lag(PREC_scaled, lag_row$PREC_lag),
      HUM_scaled_lag = lag(HUM_scaled, lag_row$HUM_lag),
      WET_scaled_lag = lag(WET_scaled, lag_row$WET_lag)
    ) %>%
    filter(
      !is.na(TMP_scaled_lag),
      !is.na(PREC_scaled_lag),
      !is.na(HUM_scaled_lag),
      !is.na(WET_scaled_lag)
    ) %>%
    ungroup()

  list(
    raw_data = data_processed,
    lagged_data = lagged_data,
    scaling_lookup = scaling_lookup
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
        Climate_Factor == "Temperature" ~ lag_row$TMP_lag,
        Climate_Factor == "Precipitation" ~ lag_row$PREC_lag,
        Climate_Factor == "Humidity" ~ lag_row$HUM_lag,
        Climate_Factor == "Wet Days" ~ lag_row$WET_lag,
        TRUE ~ NA_real_
      ),
      Term_Group = ifelse(is.na(Climate_Factor), "Adjustment", "Climate"),
      EDF = edf,
      Ref_DF = Ref.df,
      F_stat = F,
      P_value = `p-value`
    )
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

make_climate_response_lookup <- function(model, lagged_data, scaling_lookup, bacteria, lag_row, n_points = 401) {
  reference_row <- make_reference_prediction_row(lagged_data)

  scaled_lookup_list <- list()
  zone_lookup_list <- list()

  for (spec in climate_specs) {
    x_scaled <- seq(
      min(lagged_data[[spec$term]], na.rm = TRUE),
      max(lagged_data[[spec$term]], na.rm = TRUE),
      length.out = n_points
    )

    pred_data <- reference_row[rep(1, length(x_scaled)), ]
    pred_data[[spec$term]] <- x_scaled

    pred <- suppressWarnings(predict(model, pred_data, type = "terms", se.fit = TRUE))
    term_col <- grep(spec$term, colnames(pred$fit))
    fit_values <- pred$fit[, term_col]
    se_values <- pred$se.fit[, term_col]

    scaled_lookup <- tibble(
      Bacteria = bacteria,
      Climate_Factor = spec$factor,
      Lag_Years = case_when(
        spec$factor == "Temperature" ~ lag_row$TMP_lag,
        spec$factor == "Precipitation" ~ lag_row$PREC_lag,
        spec$factor == "Humidity" ~ lag_row$HUM_lag,
        TRUE ~ lag_row$WET_lag
      ),
      x_scaled = x_scaled,
      Effect_OR = exp(fit_values),
      Lower_CI = exp(fit_values - 1.96 * se_values),
      Upper_CI = exp(fit_values + 1.96 * se_values)
    )

    zone_lookup <- scaled_lookup %>%
      crossing(
        scaling_lookup %>%
          filter(Climate_Factor == spec$factor) %>%
          select(climate_zone, mean_value, sd_value)
      ) %>%
      mutate(
        x_raw = x_scaled * sd_value + mean_value,
        x_raw = case_when(
          spec$factor == "Humidity" ~ pmin(pmax(x_raw, 0), 100),
          spec$factor %in% c("Precipitation", "Wet Days") ~ pmax(x_raw, 0),
          TRUE ~ x_raw
        )
      ) %>%
      select(
        Bacteria,
        climate_zone,
        Climate_Factor,
        Lag_Years,
        x_raw,
        x_scaled,
        Effect_OR,
        Lower_CI,
        Upper_CI,
        mean_value,
        sd_value
      )

    scaled_lookup_list[[spec$factor]] <- scaled_lookup
    zone_lookup_list[[spec$factor]] <- zone_lookup
  }

  list(
    scaled = bind_rows(scaled_lookup_list),
    by_zone = bind_rows(zone_lookup_list)
  )
}

make_historical_outputs <- function(model, lagged_data, raw_data, bacteria) {
  fitted_logit <- as.numeric(predict(model, newdata = lagged_data))
  lagged_with_fit <- lagged_data %>%
    mutate(
      Bacteria = bacteria,
      fitted_logit_R = fitted_logit,
      fitted_R_prop = plogis(fitted_logit_R),
      fitted_R_pct = fitted_R_prop * 100,
      observed_R_prop = coalesce(observed_R_prop, plogis(logit_R)),
      observed_R_pct = observed_R_prop * 100
    )

  historical_zone_year <- lagged_with_fit %>%
    group_by(Bacteria, climate_zone, year) %>%
    summarise(
      observed_R_pct = mean(observed_R_pct, na.rm = TRUE),
      fitted_R_pct = mean(fitted_R_pct, na.rm = TRUE),
      n_records = n(),
      .groups = "drop"
    )

  historical_global_year <- historical_zone_year %>%
    group_by(Bacteria, year) %>%
    summarise(
      climate_zone = "Global Average",
      observed_R_pct = mean(observed_R_pct, na.rm = TRUE),
      fitted_R_pct = mean(fitted_R_pct, na.rm = TRUE),
      n_records = sum(n_records, na.rm = TRUE),
      .groups = "drop"
    )

  historical_observed_global <- raw_data %>%
    mutate(Bacteria = bacteria) %>%
    group_by(Bacteria, year) %>%
    summarise(
      climate_zone = "Global Average",
      observed_R_prop = mean(observed_R_prop, na.rm = TRUE),
      observed_R_pct = mean(observed_R_pct, na.rm = TRUE),
      n_records = n(),
      .groups = "drop"
    )

  baseline_zone <- raw_data %>%
    filter(year >= 2010, year <= 2019) %>%
    group_by(climate_zone) %>%
    summarise(
      baseline_resistance_prop = mean(observed_R_prop, na.rm = TRUE),
      baseline_resistance_pct = mean(observed_R_pct, na.rm = TRUE),
      baseline_temp = mean(TMP, na.rm = TRUE),
      baseline_precip = mean(PREC, na.rm = TRUE),
      baseline_humidity = mean(HUM, na.rm = TRUE),
      baseline_wetdays = mean(WET, na.rm = TRUE),
      n_records = n(),
      .groups = "drop"
    ) %>%
    mutate(Bacteria = bacteria) %>%
    relocate(Bacteria)

  baseline_global <- baseline_zone %>%
    summarise(
      Bacteria = first(Bacteria),
      climate_zone = "Global Average",
      baseline_resistance_prop = mean(baseline_resistance_prop, na.rm = TRUE),
      baseline_resistance_pct = mean(baseline_resistance_pct, na.rm = TRUE),
      baseline_temp = mean(baseline_temp, na.rm = TRUE),
      baseline_precip = mean(baseline_precip, na.rm = TRUE),
      baseline_humidity = mean(baseline_humidity, na.rm = TRUE),
      baseline_wetdays = mean(baseline_wetdays, na.rm = TRUE),
      n_records = sum(n_records, na.rm = TRUE)
    )

  list(
    lagged_with_fit = lagged_with_fit,
    historical_zone_year = historical_zone_year,
    historical_global_year = historical_global_year,
    historical_observed_global = historical_observed_global,
    baseline_zone = baseline_zone,
    baseline_global = baseline_global
  )
}

write_response_helper <- function(output_path, helper_dir) {
  helper_lines <- c(
    "# Auto-generated helper for primary model projection simplified GAMM inputs",
    sprintf("projection_input_dir <- '%s'", output_path),
    "default_lag_path <- file.path(projection_input_dir, '04_parameter_tables', 'projection_simplified_final_lag_settings.csv')",
    "default_weight_path <- file.path(projection_input_dir, '04_parameter_tables', 'projection_simplified_climate_weights.csv')",
    "default_scale_path <- file.path(projection_input_dir, '03_baseline_and_history', 'projection_climate_scaling_parameters.csv')",
    "default_lookup_scaled_path <- file.path(projection_input_dir, '02_curve_lookup_tables', 'projection_climate_response_lookup_scaled.csv')",
    "default_lookup_zone_path <- file.path(projection_input_dir, '02_curve_lookup_tables', 'projection_climate_response_lookup_by_zone.csv')",
    "",
    "normalize_climate_factor <- function(x) {",
    "  if (x %in% c('Temperature', 'TMP', 'temp')) return('Temperature')",
    "  if (x %in% c('Precipitation', 'PREC', 'precip', 'precipitation')) return('Precipitation')",
    "  if (x %in% c('Humidity', 'Relative Humidity', 'HUM', 'humid')) return('Humidity')",
    "  if (x %in% c('Wet Days', 'WetDays', 'WET', 'wetdays')) return('Wet Days')",
    "  x",
    "}",
    "",
    "load_projection_parameter_tables <- function(",
    "  lag_path = default_lag_path,",
    "  weight_path = default_weight_path,",
    "  scale_path = default_scale_path,",
    "  lookup_scaled_path = default_lookup_scaled_path,",
    "  lookup_zone_path = default_lookup_zone_path",
    ") {",
    "  list(",
    "    lags = read.csv(lag_path, stringsAsFactors = FALSE),",
    "    weights = read.csv(weight_path, stringsAsFactors = FALSE, check.names = FALSE),",
    "    scaling = read.csv(scale_path, stringsAsFactors = FALSE),",
    "    lookup_scaled = read.csv(lookup_scaled_path, stringsAsFactors = FALSE),",
    "    lookup_zone = read.csv(lookup_zone_path, stringsAsFactors = FALSE)",
    "  )",
    "}",
    "",
    "get_projection_lag_settings <- function(bacteria, lag_df = NULL, lag_path = default_lag_path) {",
    "  if (is.null(lag_df)) lag_df <- read.csv(lag_path, stringsAsFactors = FALSE)",
    "  out <- lag_df[lag_df$Bacteria == bacteria, , drop = FALSE]",
    "  if (nrow(out) == 0) stop('No lag settings found for ', bacteria, call. = FALSE)",
    "  out[1, ]",
    "}",
    "",
    "standardize_climate_value <- function(bacteria, climate_zone, climate_variable, raw_value, scaling_df = NULL, scale_path = default_scale_path) {",
    "  climate_variable <- normalize_climate_factor(climate_variable)",
    "  if (is.null(scaling_df)) scaling_df <- read.csv(scale_path, stringsAsFactors = FALSE)",
    "  row <- scaling_df[scaling_df$Bacteria == bacteria & scaling_df$climate_zone == climate_zone & scaling_df$Climate_Factor == climate_variable, , drop = FALSE]",
    "  if (nrow(row) == 0) stop('No scaling row found for ', bacteria, ' / ', climate_zone, ' / ', climate_variable, call. = FALSE)",
    "  (raw_value - row$mean_value[1]) / row$sd_value[1]",
    "}",
    "",
    "predict_climate_effect_scaled <- function(bacteria, climate_variable, scaled_value, lookup_df = NULL, lookup_path = default_lookup_scaled_path) {",
    "  climate_variable <- normalize_climate_factor(climate_variable)",
    "  if (is.null(lookup_df)) lookup_df <- read.csv(lookup_path, stringsAsFactors = FALSE)",
    "  rows <- lookup_df[lookup_df$Bacteria == bacteria & lookup_df$Climate_Factor == climate_variable, , drop = FALSE]",
    "  if (nrow(rows) == 0) stop('No scaled response lookup found for ', bacteria, ' / ', climate_variable, call. = FALSE)",
    "  stats::approx(rows$x_scaled, rows$Effect_OR, xout = scaled_value, rule = 2)$y",
    "}",
    "",
    "predict_climate_effect <- function(",
    "  bacteria,",
    "  climate_zone,",
    "  climate_variable,",
    "  raw_value,",
    "  lookup_df = NULL,",
    "  scaling_df = NULL,",
    "  lookup_path = default_lookup_zone_path,",
    "  scale_path = default_scale_path",
    ") {",
    "  climate_variable <- normalize_climate_factor(climate_variable)",
    "  if (is.null(lookup_df)) lookup_df <- read.csv(lookup_path, stringsAsFactors = FALSE)",
    "  if (is.null(scaling_df)) scaling_df <- read.csv(scale_path, stringsAsFactors = FALSE)",
    "  scaled_value <- standardize_climate_value(bacteria, climate_zone, climate_variable, raw_value, scaling_df = scaling_df)",
    "  rows <- lookup_df[lookup_df$Bacteria == bacteria & lookup_df$climate_zone == climate_zone & lookup_df$Climate_Factor == climate_variable, , drop = FALSE]",
    "  if (nrow(rows) == 0) stop('No zone-specific response lookup found for ', bacteria, ' / ', climate_zone, ' / ', climate_variable, call. = FALSE)",
    "  stats::approx(rows$x_scaled, rows$Effect_OR, xout = scaled_value, rule = 2)$y",
    "}"
  )

  writeLines(helper_lines, file.path(helper_dir, "projection_simplified_response_functions.R"))
}

lag_settings_df <- load_final_lag_settings(final_lag_settings_path)
weights_df <- load_projection_weights(weights_path)
weight_details_df <- if (file.exists(weight_details_path)) {
  read.csv(weight_details_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  tibble()
}

curve_scaled_all <- list()
curve_zone_all <- list()
scaling_all <- list()
term_stats_all <- list()
baseline_zone_all <- list()
baseline_global_all <- list()
history_zone_all <- list()
history_global_all <- list()
history_observed_all <- list()
manifest_rows <- list()

for (spec in bacteria_specs) {
  message("Preparing final simplified projection GAMM inputs for ", spec$title)

  lag_row <- lag_settings_df %>% filter(Bacteria == spec$title)
  if (nrow(lag_row) == 0) {
    warning("No final lag settings found for ", spec$title, "; skipping.")
    next
  }

  prepared <- prepare_projection_data(file.path(input_data_dir, spec$file_name), lag_row)
  model <- build_simplified_model(prepared$lagged_data)
  term_stats <- extract_model_term_statistics(model, spec$title, lag_row)
  response_lookup <- make_climate_response_lookup(model, prepared$lagged_data, prepared$scaling_lookup, spec$title, lag_row)
  historical_outputs <- make_historical_outputs(model, prepared$lagged_data, prepared$raw_data, spec$title)

  stem <- sanitize_file_stem(spec$title)
  model_rds_path <- file.path(output_base, "01_model_objects", paste0("Projection_Simplified_ModelC_", stem, ".rds"))
  model_summary_path <- file.path(output_base, "01_model_objects", paste0("Projection_Simplified_ModelC_", stem, "_summary.txt"))

  saveRDS(
    list(
      bacteria = spec$title,
      lag_settings = lag_row,
      model = model,
      term_statistics = term_stats,
      scaling_lookup = prepared$scaling_lookup
    ),
    model_rds_path
  )

  sink(model_summary_path)
  cat("Projection simplified GAMM:", spec$title, "\n\n")
  cat("Final lag settings used:\n")
  print(lag_row)
  cat("\nModel summary:\n")
  print(summary(model))
  cat("\n")
  print(gam.check(model))
  sink()

  curve_scaled_all[[spec$title]] <- response_lookup$scaled
  curve_zone_all[[spec$title]] <- response_lookup$by_zone
  scaling_all[[spec$title]] <- prepared$scaling_lookup %>% mutate(Bacteria = spec$title, .before = 1)
  term_stats_all[[spec$title]] <- term_stats
  baseline_zone_all[[spec$title]] <- historical_outputs$baseline_zone
  baseline_global_all[[spec$title]] <- historical_outputs$baseline_global
  history_zone_all[[spec$title]] <- historical_outputs$historical_zone_year
  history_global_all[[spec$title]] <- historical_outputs$historical_global_year
  history_observed_all[[spec$title]] <- historical_outputs$historical_observed_global

  manifest_rows[[length(manifest_rows) + 1]] <- tibble(
    Bacteria = spec$title,
    Model_RDS = model_rds_path,
    Model_Summary = model_summary_path
  )
}

curve_scaled_df <- bind_rows(curve_scaled_all)
curve_zone_df <- bind_rows(curve_zone_all)
scaling_df <- bind_rows(scaling_all) %>%
  mutate(Climate_Factor = normalize_climate_factor(Climate_Factor))
term_stats_df <- bind_rows(term_stats_all)
baseline_zone_df <- bind_rows(baseline_zone_all)
baseline_global_df <- bind_rows(baseline_global_all)
history_zone_df <- bind_rows(history_zone_all)
history_global_df <- bind_rows(history_global_all)
history_observed_df <- bind_rows(history_observed_all)
manifest_df <- bind_rows(manifest_rows)

write.csv(curve_scaled_df, file.path(output_base, "02_curve_lookup_tables", "projection_climate_response_lookup_scaled.csv"), row.names = FALSE)
write.csv(curve_zone_df, file.path(output_base, "02_curve_lookup_tables", "projection_climate_response_lookup_by_zone.csv"), row.names = FALSE)
write.csv(scaling_df, file.path(output_base, "03_baseline_and_history", "projection_climate_scaling_parameters.csv"), row.names = FALSE)
write.csv(baseline_zone_df, file.path(output_base, "03_baseline_and_history", "projection_baseline_by_zone.csv"), row.names = FALSE)
write.csv(baseline_global_df, file.path(output_base, "03_baseline_and_history", "projection_baseline_global.csv"), row.names = FALSE)
write.csv(history_zone_df, file.path(output_base, "03_baseline_and_history", "projection_historical_fit_zone_year.csv"), row.names = FALSE)
write.csv(history_global_df, file.path(output_base, "03_baseline_and_history", "projection_historical_fit_global_year.csv"), row.names = FALSE)
write.csv(history_observed_df, file.path(output_base, "03_baseline_and_history", "projection_historical_observed_global_year.csv"), row.names = FALSE)
write.csv(lag_settings_df, file.path(output_base, "04_parameter_tables", "projection_simplified_final_lag_settings.csv"), row.names = FALSE)
write.csv(weights_df, file.path(output_base, "04_parameter_tables", "projection_simplified_climate_weights.csv"), row.names = FALSE)
write.csv(weight_details_df, file.path(output_base, "04_parameter_tables", "projection_simplified_climate_weight_details.csv"), row.names = FALSE)
write.csv(term_stats_df, file.path(output_base, "04_parameter_tables", "projection_simplified_model_term_statistics.csv"), row.names = FALSE)
write.csv(manifest_df, file.path(output_base, "07_metadata", "projection_input_manifest.csv"), row.names = FALSE)

write_xlsx(
  list(
    Final_Lags = lag_settings_df,
    Climate_Weights = weights_df,
    Weight_Details = weight_details_df,
    Term_Statistics = term_stats_df,
    Scaling = scaling_df,
    Baseline_By_Zone = baseline_zone_df,
    Baseline_Global = baseline_global_df,
    Historical_Fit_Global = history_global_df,
    Historical_Observed_Global = history_observed_df,
    Response_Lookup_Scaled = curve_scaled_df,
    Response_Lookup_By_Zone = curve_zone_df
  ),
  file.path(output_base, "06_workbook", "Projection_Simplified_GAMM_Preparation_source_data.xlsx")
)

write_response_helper(output_base, file.path(output_base, "05_response_function_helper"))

readme_lines <- c(
  "# Projection Simplified GAMM Ready Inputs",
  "",
  "This folder contains machine-readable inputs prepared from the final simplified projection GAMMs under the primary model data system.",
  "",
  "Key principles:",
  "- Lag settings come from the BEST SIMPLIFIED MODEL LAG selected during the projection-specific lag search.",
  "- Climate weights are taken from projection climate-weight summary and correspond to the same best simplified GAMMs.",
  "- Response lookup tables are exported as multiplicative odds-ratio effects (OR, centered around 1) rather than direct resistance-rate predictions.",
  "- Scaling parameters are stored by bacteria and climate zone so future CMIP6 climate values can be standardized consistently before response lookup.",
  "",
  "Important files:",
  "- `02_curve_lookup_tables/projection_climate_response_lookup_scaled.csv`: scaled-domain response lookup for each bacteria and climate factor.",
  "- `02_curve_lookup_tables/projection_climate_response_lookup_by_zone.csv`: climate-zone-specific raw-value lookup derived from the scaled-domain response curves.",
  "- `03_baseline_and_history/projection_baseline_by_zone.csv`: 2010-2019 baseline resistance and climate values by bacteria and climate zone.",
  "- `04_parameter_tables/projection_simplified_final_lag_settings.csv`: final lag settings used in downstream projection modeling.",
  "- `04_parameter_tables/projection_simplified_climate_weights.csv`: final climatic-factor weights aligned with projection climate-weight summary.",
  "- `05_response_function_helper/projection_simplified_response_functions.R`: helper functions for later projection scripts."
)
writeLines(readme_lines, file.path(output_base, "07_metadata", "README.md"))
writeLines(capture.output(sessionInfo()), file.path(output_base, "07_metadata", "sessionInfo.txt"))

cat("Saved final simplified projection GAMM preparation outputs to:\n")
cat(output_base, "\n")
