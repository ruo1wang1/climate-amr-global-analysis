#!/usr/bin/env Rscript

# Figure 2 turning-point analysis and observation-level influence diagnostics.

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
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

input_root <- Sys.getenv(
  "FIGURE2_INPUT_ROOT",
  unset = file.path(repo_root, "outputs", "historical_associations", "figure2_inputs")
)
analysis_input_dir <- input_root
model_input_root <- input_root
output_root <- Sys.getenv(
  "FIGURE2_OUTPUT_ROOT",
  unset = file.path(repo_root, "outputs", "historical_associations", "figure2_current")
)

source(file.path(script_dir, "figure2_turning_point_core.R"), local = FALSE)

required_public_inputs <- c(
  "bacteria_four_factors_lag_summary.csv",
  "Figure2_v3_main_figure_strict_thresholds_max2.csv",
  "Figure2_v3_bootstrap_support_for_observed_turning_points.csv",
  "Figure2_reference_termcentred_curve_source_data.csv"
)
missing_public_inputs <- required_public_inputs[
  !file.exists(file.path(analysis_input_dir, required_public_inputs))
]
if (length(missing_public_inputs) > 0L) {
  stop("Missing Figure 2 analysis inputs: ", paste(missing_public_inputs, collapse = ", "), call. = FALSE)
}

required_model_inputs <- c(
  file.path("frozen_models", v3_specs$model_file),
  file.path("model_ready_inputs", v3_specs$input_file)
)
missing_model_inputs <- required_model_inputs[
  !file.exists(file.path(model_input_root, required_model_inputs))
]
if (length(missing_model_inputs) > 0L) {
  stop(
    "Model-ready panels or fitted Model C objects are missing. Set FIGURE2_INPUT_ROOT ",
    "to a directory containing frozen_models/, model_ready_inputs/ and the required CSV inputs.",
    call. = FALSE
  )
}

dirs <- list(
  root = output_root,
  contract = file.path(output_root, "00_figure_contract"),
  inputs = file.path(output_root, "00_inputs"),
  input_models = file.path(output_root, "00_inputs", "frozen_models"),
  input_data = file.path(output_root, "00_inputs", "model_ready_inputs"),
  code = file.path(output_root, "01_code"),
  figures = file.path(output_root, "02_figures"),
  panels = file.path(output_root, "02_figures", "single_panels"),
  source = file.path(output_root, "03_source_data"),
  tables = file.path(output_root, "04_tables"),
  diagnostics = file.path(output_root, "05_diagnostics"),
  sensitivity_models = file.path(output_root, "05_diagnostics", "influence_excluded_models"),
  report = file.path(output_root, "06_report"),
  logs = file.path(output_root, "07_logs")
)
walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(dirs$logs, "Figure2_threshold_ORCI_supplementary_features_v9_two_decimal_turning_point_labels_run.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Figure 2 zero-cascade redraw started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("Model input root:", model_input_root, "\n")
cat("Output:", output_root, "\n\n")

# Keep an auditable, self-contained copy of the small frozen inputs used here.
file.copy(
  file.path(analysis_input_dir, "bacteria_four_factors_lag_summary.csv"),
  file.path(dirs$inputs, "bacteria_four_factors_lag_summary.csv"),
  overwrite = TRUE
)
file.copy(
  file.path(script_dir, "figure2_turning_point_core.R"),
  file.path(dirs$code, "figure2_turning_point_core.R"),
  overwrite = TRUE
)
file.copy(
  file.path(analysis_input_dir, "Figure2_v3_main_figure_strict_thresholds_max2.csv"),
  file.path(dirs$inputs, "Figure2_v3_main_figure_strict_thresholds_max2.csv"),
  overwrite = TRUE
)
file.copy(
  file.path(analysis_input_dir, "Figure2_v3_bootstrap_support_for_observed_turning_points.csv"),
  file.path(dirs$inputs, "Figure2_v3_bootstrap_support_for_observed_turning_points.csv"),
  overwrite = TRUE
)

for (i in seq_len(nrow(v3_specs))) {
  file.copy(
    file.path(model_input_root, "frozen_models", v3_specs$model_file[i]),
    file.path(dirs$input_models, v3_specs$model_file[i]),
    overwrite = TRUE
  )
  file.copy(
    file.path(model_input_root, "model_ready_inputs", v3_specs$input_file[i]),
    file.path(dirs$input_data, v3_specs$input_file[i]),
    overwrite = TRUE
  )
}

v3_load_lags_cache <<- v3_load_lags(
  file.path(dirs$inputs, "bacteria_four_factors_lag_summary.csv")
)

models <- setNames(lapply(seq_len(nrow(v3_specs)), function(i) {
  readRDS(file.path(dirs$input_models, v3_specs$model_file[i]))
}), v3_specs$phenotype)

prepared <- setNames(lapply(seq_len(nrow(v3_specs)), function(i) {
  v3_prepare_data(
    file.path(dirs$input_data, v3_specs$input_file[i]),
    v3_specs$phenotype[i],
    v3_load_lags_cache
  )
}), v3_specs$phenotype)

model_data_qa <- map_dfr(v3_specs$phenotype, function(phenotype) {
  v3_validate_model_data(models[[phenotype]], prepared[[phenotype]], phenotype)
})
write.csv(
  model_data_qa,
  file.path(dirs$diagnostics, "Figure2_model_frame_exact_match_QA.csv"),
  row.names = FALSE
)
if (any(!model_data_qa$passed)) {
  stop("Frozen model-frame validation failed; no figure was produced.", call. = FALSE)
}

palette <- c(
  Temperature = "#DD5F60",
  Humidity = "#3CB371",
  Precipitation = "#9AC0CD",
  WetDays = "#8A2BE2",
  neutral_dark = "#2E2E2E",
  neutral_mid = "#7D7D7D",
  neutral_light = "#D6D6D6",
  tail_line = "#A7A7A7",
  tail_fill = "#D9D9D9",
  threshold = "#D73027",
  reference = "#666666"
)

physical_axis_specs <- tribble(
  ~climate,        ~unit_display, ~axis_title,
  "Temperature",   "°C",         "Temperature",
  "Humidity",      "%",          "Relative humidity",
  "Precipitation", "mm",         "Precipitation",
  "WetDays",       "d",          "Wet days"
)

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(quantile(x, p, type = 7, names = FALSE))
}

format_physical_value <- function(x, climate) {
  x <- ifelse(abs(x) < 0.005, 0, x)
  if (climate == "Precipitation") {
    formatC(x, format = "f", digits = 2, big.mark = ",")
  } else {
    formatC(x, format = "f", digits = 2)
  }
}

format_threshold_location <- function(x, climate) {
  value <- format_physical_value(x, climate)
  case_when(
    climate == "Temperature" ~ paste0(value, " °C"),
    climate == "Humidity" ~ paste0(value, "%"),
    climate == "Precipitation" ~ paste0(value, " mm"),
    climate == "WetDays" ~ paste0(value, " d"),
    TRUE ~ value
  )
}

# Return exactly five evenly spaced, human-readable ticks that cover the data.
# A five-tick contract is more reproducible than pretty(..., n=4), whose output
# count varies with the range. Zero is not permitted on an odds-ratio axis unless
# it is needed to cover the plotted interval.
five_nice_breaks <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) stop("Cannot derive y-axis breaks from empty values.")
  lo <- min(values)
  hi <- max(values)
  if (hi <= lo) {
    delta <- max(abs(lo) * 0.05, 0.01)
    lo <- lo - delta
    hi <- hi + delta
  }
  raw_step <- (hi - lo) / 4
  exponent <- floor(log10(raw_step))
  steps <- sort(unique(as.vector(outer(
    c(1, 2, 2.5, 5, 10),
    10 ^ seq(exponent - 2, exponent + 3),
    `*`
  ))))
  candidates <- map_dfr(steps[steps >= raw_step * (1 - 1e-12)], function(step) {
    anchors <- seq(floor(lo / step) - 4, floor(lo / step) + 1)
    map_dfr(anchors, function(anchor) {
      start <- anchor * step
      end <- start + 4 * step
      if (start <= lo + 1e-12 && end >= hi - 1e-12 && start >= 0) {
        tibble(
          step = step,
          start = start,
          end = end,
          padding = (lo - start) + (end - hi),
          includes_one = any(abs(seq(start, end, length.out = 5) - 1) < step * 1e-8)
        )
      } else tibble()
    })
  })
  if (nrow(candidates) == 0L) {
    return(seq(lo, hi, length.out = 5L))
  }
  choice <- candidates %>%
    arrange(padding, desc(includes_one), step) %>%
    slice(1)
  seq(choice$start, choice$end, length.out = 5L)
}

# Construct a panel-specific physical x axis without changing the fitted curve
# domain. Candidate major ticks use climate-appropriate round increments, while
# the display limits are allowed to extend only far enough to contain the
# observed min-max range. Except when a natural zero-bounded axis is available,
# the first and last labelled ticks have exactly the same distance to the panel
# boundary. This avoids the asymmetric terminal segments produced by
# pretty(xr, n = 4) combined with independently padded limits.
balanced_physical_x_axis <- function(values, climate) {
  xr <- range(values[is.finite(values)])
  if (length(xr) != 2L || !all(is.finite(xr)) || xr[2] <= xr[1]) {
    stop("Cannot derive balanced x axis from the supplied values.", call. = FALSE)
  }

  span <- diff(xr)
  data_pad <- 0.015 * span

  if (climate == "Temperature") {
    steps <- 5
    offsets <- 0
    target_ticks <- 5L
    coverage_weight <- 0.05
  } else if (climate == "Humidity") {
    steps <- 10
    offsets <- 0
    target_ticks <- 5L
    coverage_weight <- 0.05
  } else if (climate == "WetDays") {
    steps <- 50
    offsets <- 0
    target_ticks <- 5L
    coverage_weight <- 0.40
  } else if (climate == "Precipitation") {
    # A 400-mm shifted grid avoids a large unsupported 0-300 mm display gap in
    # panels such as CR-Pa, while the conventional 500-mm grid remains
    # preferred when it fits the observed range well.
    steps <- c(400, 500)
    offsets <- c(0, 0.25, 0.50, 0.75)
    target_ticks <- 5L
    coverage_weight <- 0.40
  } else {
    stop("Unknown climate axis: ", climate, call. = FALSE)
  }

  candidates <- list()
  candidate_index <- 1L
  for (step in steps) {
    for (offset_fraction in offsets) {
      offset <- offset_fraction * step
      grid_index <- seq(
        floor((xr[1] - offset) / step) - 2,
        ceiling((xr[2] - offset) / step) + 2
      )
      grid <- grid_index * step + offset

      for (n_ticks in 3:7) {
        if (length(grid) < n_ticks) next
        for (start_index in seq_len(length(grid) - n_ticks + 1L)) {
          breaks <- grid[start_index:(start_index + n_ticks - 1L)]

          natural_zero_axis <-
            climate %in% c("Precipitation", "WetDays") &&
            abs(breaks[1]) < 1e-10 &&
            breaks[1] <= xr[1] && breaks[n_ticks] >= xr[2]

          if (natural_zero_axis) {
            limits <- c(0, breaks[n_ticks])
            edge_gap <- 0
          } else {
            required_gap <- max(
              0.04 * step,
              breaks[1] - (xr[1] - data_pad),
              (xr[2] + data_pad) - breaks[n_ticks]
            )
            gap_unit <- step / 20
            edge_gap <- ceiling(required_gap / gap_unit - 1e-12) * gap_unit
            limits <- c(breaks[1] - edge_gap, breaks[n_ticks] + edge_gap)
          }

          if (limits[1] > xr[1] + 1e-9 || limits[2] < xr[2] - 1e-9) next

          outside_tick_fraction <- (
            max(0, breaks[1] - xr[1]) +
              max(0, xr[2] - breaks[n_ticks])
          ) / span
          expansion_fraction <- diff(limits) / span - 1
          score <-
            expansion_fraction +
            0.025 * abs(n_ticks - target_ticks) +
            coverage_weight * outside_tick_fraction +
            ifelse(step == 400, 0.02, 0) +
            ifelse(offset_fraction == 0, 0, 0.015)

          candidates[[candidate_index]] <- tibble(
            score = score,
            step = step,
            offset_fraction = offset_fraction,
            n_ticks = n_ticks,
            edge_gap = edge_gap,
            display_min = limits[1],
            display_max = limits[2],
            expansion_fraction = expansion_fraction,
            outside_tick_fraction = outside_tick_fraction,
            breaks = list(breaks),
            limits = list(limits)
          )
          candidate_index <- candidate_index + 1L
        }
      }
    }
  }

  if (length(candidates) == 0L) {
    stop("No valid balanced x-axis candidate for ", climate, call. = FALSE)
  }

  choice <- bind_rows(candidates) %>%
    arrange(score, expansion_fraction, abs(n_ticks - target_ticks)) %>%
    slice(1)

  list(
    observed_range = xr,
    limits = choice$limits[[1]],
    breaks = choice$breaks[[1]],
    step = choice$step,
    edge_gap = choice$edge_gap,
    expansion_fraction = choice$expansion_fraction,
    outside_tick_fraction = choice$outside_tick_fraction
  )
}

build_physical_mapping <- function(dat, phenotype, factor_row) {
  variable <- factor_row$variable
  raw_lag <- paste0(factor_row$raw, "_orig_lag")
  if (!raw_lag %in% names(dat)) {
    stop("Missing lagged physical exposure: ", raw_lag, call. = FALSE)
  }
  axis_spec <- physical_axis_specs %>% filter(climate == factor_row$climate)

  mapping <- dat %>%
    transmute(
      phenotype = phenotype,
      climate = factor_row$climate,
      variable = variable,
      climate_zone = as.character(climate_zone),
      x_scaled = .data[[variable]],
      x_physical = .data[[raw_lag]]
    ) %>%
    filter(is.finite(x_scaled), is.finite(x_physical)) %>%
    group_by(phenotype, climate, variable, climate_zone) %>%
    group_modify(~ {
      fit <- lm(x_physical ~ x_scaled, data = .x)
      fitted_values <- fitted(fit)
      total_ss <- sum((.x$x_physical - mean(.x$x_physical))^2)
      tibble(
        map_intercept = unname(coef(fit)[1]),
        map_slope = unname(coef(fit)[2]),
        mapping_r_squared = if (total_ss > 0) {
          1 - sum((.x$x_physical - fitted_values)^2) / total_ss
        } else 1,
        max_abs_mapping_residual = max(abs(.x$x_physical - fitted_values)),
        n_mapping_observations = nrow(.x),
        physical_observed_min = min(.x$x_physical),
        physical_p2_5 = safe_quantile(.x$x_physical, 0.025),
        physical_p97_5 = safe_quantile(.x$x_physical, 0.975),
        physical_observed_max = max(.x$x_physical)
      )
    }) %>%
    ungroup() %>%
    mutate(
      unit_display = axis_spec$unit_display,
      axis_title = axis_spec$axis_title
    )

  reference <- mapping %>% filter(climate_zone == "Temperate Zone")
  if (nrow(reference) != 1L) {
    reference <- mapping %>% slice_max(n_mapping_observations, n = 1, with_ties = FALSE)
  }
  list(all_zones = mapping, reference = reference)
}

compute_penalized_influence <- function(model, dat, phenotype) {
  X <- suppressWarnings(predict(model, type = "lpmatrix"))
  if (nrow(X) != nrow(model$model)) {
    stop("Lpmatrix row count mismatch for ", phenotype, call. = FALSE)
  }
  sigma2 <- as.numeric(model$sig2)
  prior_weights <- model$prior.weights
  if (is.null(prior_weights) || length(prior_weights) != nrow(X)) {
    prior_weights <- rep(1, nrow(X))
  }
  leverage <- prior_weights * rowSums((X %*% model$Vp) * X) / sigma2
  leverage <- pmin(pmax(leverage, 0), 1 - 1e-8)
  response_residual <- as.numeric(residuals(model, type = "response"))
  effective_parameters <- sum(leverage)
  studentized_residual <- response_residual /
    sqrt(sigma2 * pmax(1 - leverage, 1e-8))
  cooks_distance <- (response_residual^2 / (effective_parameters * sigma2)) *
    leverage / pmax((1 - leverage)^2, 1e-12)
  n_obs <- nrow(X)
  leverage_cutoff <- 2 * effective_parameters / n_obs
  cooks_cutoff <- 4 / n_obs

  tibble(
    phenotype = phenotype,
    observation_id = seq_len(n_obs),
    country = as.character(dat$NAME),
    year = dat$year,
    response = as.numeric(model$model[[all.vars(formula(model))[1]]]),
    fitted = as.numeric(fitted(model)),
    response_residual = response_residual,
    leverage = leverage,
    studentized_residual = studentized_residual,
    cooks_distance = cooks_distance,
    effective_parameters = effective_parameters,
    leverage_cutoff = leverage_cutoff,
    cooks_cutoff = cooks_cutoff,
    high_leverage = leverage > leverage_cutoff,
    large_studentized_residual = abs(studentized_residual) > 3,
    high_cooks_distance = cooks_distance > cooks_cutoff,
    influence_exclusion_flag = large_studentized_residual | high_cooks_distance
  )
}

refit_influence_excluded_model <- function(model, influence, phenotype) {
  keep <- !influence$influence_exclusion_flag
  if (sum(keep) < 0.80 * length(keep)) {
    stop("Influence rule would remove over 20% of observations for ", phenotype, call. = FALSE)
  }
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
  fit <- bam(
    formula(model), data = model$model[keep, , drop = FALSE],
    family = gaussian(), method = "REML", select = TRUE,
    discrete = FALSE, use.chol = FALSE, control = ctrl
  )
  attr(fit, "influence_exclusion_rule") <-
    "abs(internally studentized residual) > 3 OR approximate penalized Cook distance > 4/n"
  attr(fit, "n_excluded") <- sum(!keep)
  fit
}

make_quantile_bins <- function(x, n_bins = 10L) {
  br <- unique(as.numeric(quantile(
    x[is.finite(x)], seq(0, 1, length.out = n_bins + 1L),
    type = 7, names = FALSE
  )))
  if (length(br) < 3L) {
    br <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 3L)
  }
  cut(x, breaks = br, include.lowest = TRUE, labels = FALSE)
}

make_histogram <- function(x, support_limits, phenotype, climate, variable,
                           lag_years, n_breaks = 16L) {
  xr <- range(x, na.rm = TRUE)
  breaks <- seq(xr[1], xr[2], length.out = n_breaks + 1L)
  h <- hist(x, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = TRUE)
  tibble(
    phenotype = phenotype,
    climate = climate,
    variable = variable,
    lag_years = lag_years,
    bin_id = seq_along(h$counts),
    xmin = head(h$breaks, -1L),
    xmax = tail(h$breaks, -1L),
    xmid = h$mids,
    count = h$counts,
    support_class = if_else(
      h$mids >= support_limits[1] & h$mids <= support_limits[2],
      "P2.5-P97.5", "Observed tail"
    )
  )
}

make_partial_residual_data <- function(model, dat, phenotype, factor_row,
                                       support_limits, influence, physical_reference,
                                       n_bins = 10L) {
  variable <- factor_row$variable
  term_matrix <- suppressWarnings(predict(model, type = "terms"))
  idx <- v3_target_term_index(term_matrix, variable)
  term_eta <- as.numeric(term_matrix[, idx])
  response_residual <- as.numeric(residuals(model, type = "response"))
  partial_eta <- term_eta + response_residual
  x <- as.numeric(model$model[[variable]])
  raw_lag <- paste0(factor_row$raw, "_orig_lag")
  x_physical_actual <- as.numeric(dat[[raw_lag]])
  x_physical_equivalent <- physical_reference$map_intercept +
    physical_reference$map_slope * x

  raw <- tibble(
    phenotype = phenotype,
    climate = factor_row$climate,
    variable = variable,
    lag_years = as.integer(v3_load_lags_cache[[phenotype]][[factor_row$raw]]),
    observation_id = seq_along(x),
    country = as.character(dat$NAME),
    year = dat$year,
    x_scaled = x,
    x_physical_actual = x_physical_actual,
    x_physical_equivalent = x_physical_equivalent,
    term_eta = term_eta,
    response_residual = response_residual,
    partial_eta = partial_eta,
    partial_OR = exp(partial_eta),
    support_class = if_else(
      x >= support_limits[1] & x <= support_limits[2],
      "P2.5-P97.5", "Observed tail"
    ),
    quantile_bin = make_quantile_bins(x, n_bins)
  ) %>%
    left_join(
      influence %>%
        select(
          observation_id, leverage, studentized_residual, cooks_distance,
          high_leverage, large_studentized_residual, high_cooks_distance,
          influence_exclusion_flag
        ),
      by = "observation_id"
    )

  binned <- raw %>%
    filter(is.finite(quantile_bin)) %>%
    group_by(phenotype, climate, variable, lag_years, quantile_bin) %>%
    summarise(
      x_mean_scaled = mean(x_scaled),
      x_median_scaled = median(x_scaled),
      x_min_scaled = min(x_scaled),
      x_max_scaled = max(x_scaled),
      x_mean_physical_equivalent = mean(x_physical_equivalent),
      x_median_physical_equivalent = median(x_physical_equivalent),
      x_min_physical_equivalent = min(x_physical_equivalent),
      x_max_physical_equivalent = max(x_physical_equivalent),
      partial_eta_mean = mean(partial_eta),
      partial_eta_median = median(partial_eta),
      partial_eta_q25 = safe_quantile(partial_eta, 0.25),
      partial_eta_q75 = safe_quantile(partial_eta, 0.75),
      partial_OR_geometric_mean = exp(partial_eta_mean),
      partial_OR_median = exp(partial_eta_median),
      partial_OR_q25 = exp(partial_eta_q25),
      partial_OR_q75 = exp(partial_eta_q75),
      n_observations = n(),
      n_countries = n_distinct(country),
      n_influence_flagged = sum(influence_exclusion_flag),
      .groups = "drop"
    ) %>%
    mutate(
      support_class = if_else(
        x_mean_scaled >= support_limits[1] & x_mean_scaled <= support_limits[2],
        "P2.5-P97.5", "Observed tail"
      )
    )

  qa <- tibble(
    phenotype = phenotype,
    climate = factor_row$climate,
    variable = variable,
    n_model = nrow(model$model),
    n_partial_residual_rows = nrow(raw),
    n_rows_in_bins = sum(binned$n_observations),
    n_bins_requested = n_bins,
    n_bins_obtained = nrow(binned),
    max_partial_residual_identity_error = max(
      abs(raw$partial_eta - raw$term_eta - raw$response_residual),
      na.rm = TRUE
    )
  )
  list(raw = raw, binned = binned, qa = qa)
}

make_influence_excluded_sensitivity_curve <- function(
    sensitivity_model, main_curve, original_x, phenotype, factor_row) {
  x_grid <- main_curve$x_scaled
  sensitivity <- v3_term_curve(
    sensitivity_model, factor_row$variable, x_grid, with_se = TRUE
  )
  reference_x <- safe_quantile(original_x, 0.50)
  main_reference_eta <- approx(
    main_curve$x_scaled, main_curve$eta,
    xout = reference_x, ties = "ordered"
  )$y
  sensitivity_reference_eta <- approx(
    sensitivity$x_scaled, sensitivity$eta,
    xout = reference_x, ties = "ordered"
  )$y
  sensitivity %>%
    transmute(
      phenotype = phenotype,
      climate = factor_row$climate,
      variable = factor_row$variable,
      x_scaled,
      eta_unaligned = eta,
      OR_unaligned = OR,
      eta_aligned_to_main_at_P50 = eta - sensitivity_reference_eta + main_reference_eta,
      OR_aligned_to_main_at_P50 = exp(eta_aligned_to_main_at_P50),
      alignment_reference_x_scaled = reference_x,
      n_sensitivity_model = nrow(sensitivity_model$model),
      n_excluded = attr(sensitivity_model, "n_excluded")
    )
}

curve_list <- list()
support_list <- list()
histogram_list <- list()
rug_list <- list()
partial_raw_list <- list()
partial_bin_list <- list()
partial_qa_list <- list()
physical_mapping_list <- list()
sensitivity_curve_list <- list()
panel_key <- 1L

influence_by_phenotype <- setNames(map(v3_specs$phenotype, function(phenotype) {
  compute_penalized_influence(
    models[[phenotype]], prepared[[phenotype]], phenotype
  )
}), v3_specs$phenotype)

sensitivity_models <- setNames(map(v3_specs$phenotype, function(phenotype) {
  model_path <- file.path(
    dirs$sensitivity_models,
    paste0("ModelC_", v3_safe_stem(phenotype), "_influence_excluded_sensitivity.rds")
  )
  if (file.exists(model_path)) {
    cat("Loading cached influence-excluded sensitivity fit:", phenotype, "\n")
    fit <- readRDS(model_path)
  } else {
    cat("Influence-excluded sensitivity refit:", phenotype, "\n")
    fit <- refit_influence_excluded_model(
      models[[phenotype]], influence_by_phenotype[[phenotype]], phenotype
    )
    saveRDS(fit, model_path)
  }
  fit
}), v3_specs$phenotype)

for (i in seq_len(nrow(v3_specs))) {
  phenotype <- v3_specs$phenotype[i]
  model <- models[[phenotype]]
  dat <- prepared[[phenotype]]
  influence <- influence_by_phenotype[[phenotype]]
  sensitivity_model <- sensitivity_models[[phenotype]]
  cat("Processing", phenotype, "\n")

  for (j in seq_len(nrow(v3_factors))) {
    factor_row <- v3_factors[j, ]
    analysis <- v3_observed_analysis(
      model, dat, phenotype, factor_row,
      n_curve = 801L, n_search = 1001L
    )
    physical <- build_physical_mapping(dat, phenotype, factor_row)
    physical_reference <- physical$reference
    curve <- analysis$curve %>%
      mutate(
        reference_zone = physical_reference$climate_zone,
        map_intercept = physical_reference$map_intercept,
        map_slope = physical_reference$map_slope,
        x_physical_equivalent = map_intercept + map_slope * x_scaled,
        support_class = if_else(inside_P2_5_P97_5, "P2.5-P97.5", "Observed tail"),
        tail_group = case_when(
          x_scaled < analysis$support$p2_5_scaled ~ "Lower tail",
          x_scaled > analysis$support$p97_5_scaled ~ "Upper tail",
          TRUE ~ "P2.5-P97.5"
        )
      )
    support_row <- analysis$support %>%
      mutate(
        reference_zone = physical_reference$climate_zone,
        map_intercept = physical_reference$map_intercept,
        map_slope = physical_reference$map_slope,
        unit_display = physical_reference$unit_display,
        axis_title = physical_reference$axis_title,
        observed_min_physical_equivalent = map_intercept + map_slope * observed_min_scaled,
        p2_5_physical_equivalent = map_intercept + map_slope * p2_5_scaled,
        p97_5_physical_equivalent = map_intercept + map_slope * p97_5_scaled,
        observed_max_physical_equivalent = map_intercept + map_slope * observed_max_scaled,
        actual_all_zone_physical_min = min(dat[[paste0(factor_row$raw, "_orig_lag")]], na.rm = TRUE),
        actual_all_zone_physical_p2_5 = safe_quantile(dat[[paste0(factor_row$raw, "_orig_lag")]], 0.025),
        actual_all_zone_physical_p97_5 = safe_quantile(dat[[paste0(factor_row$raw, "_orig_lag")]], 0.975),
        actual_all_zone_physical_max = max(dat[[paste0(factor_row$raw, "_orig_lag")]], na.rm = TRUE)
      )
    x <- as.numeric(model$model[[factor_row$variable]])

    histogram <- make_histogram(
      x,
      c(support_row$p2_5_scaled, support_row$p97_5_scaled),
      phenotype, factor_row$climate, factor_row$variable,
      support_row$lag_years
    ) %>%
      mutate(
        xmin_physical_equivalent = physical_reference$map_intercept +
          physical_reference$map_slope * xmin,
        xmax_physical_equivalent = physical_reference$map_intercept +
          physical_reference$map_slope * xmax,
        xmid_physical_equivalent = physical_reference$map_intercept +
          physical_reference$map_slope * xmid
      )
    rug <- tibble(
      phenotype = phenotype,
      climate = factor_row$climate,
      variable = factor_row$variable,
      lag_years = support_row$lag_years,
      observation_id = seq_along(x),
      country = as.character(dat$NAME),
      year = dat$year,
      x_scaled = x,
      x_physical_actual = dat[[paste0(factor_row$raw, "_orig_lag")]],
      x_physical_equivalent = physical_reference$map_intercept +
        physical_reference$map_slope * x,
      support_class = if_else(
        x >= support_row$p2_5_scaled & x <= support_row$p97_5_scaled,
        "P2.5-P97.5", "Observed tail"
      )
    ) %>%
      left_join(
        influence %>%
          select(
            observation_id, leverage, studentized_residual, cooks_distance,
            influence_exclusion_flag
          ),
        by = "observation_id"
      )
    partial <- make_partial_residual_data(
      model, dat, phenotype, factor_row,
      c(support_row$p2_5_scaled, support_row$p97_5_scaled),
      influence = influence,
      physical_reference = physical_reference,
      n_bins = 10L
    )
    sensitivity_curve <- make_influence_excluded_sensitivity_curve(
      sensitivity_model, curve, x, phenotype, factor_row
    ) %>%
      mutate(
        x_physical_equivalent = physical_reference$map_intercept +
          physical_reference$map_slope * x_scaled
      )

    curve_list[[panel_key]] <- curve
    support_list[[panel_key]] <- support_row
    histogram_list[[panel_key]] <- histogram
    rug_list[[panel_key]] <- rug
    partial_raw_list[[panel_key]] <- partial$raw
    partial_bin_list[[panel_key]] <- partial$binned
    partial_qa_list[[panel_key]] <- partial$qa
    physical_mapping_list[[panel_key]] <- physical$all_zones
    sensitivity_curve_list[[panel_key]] <- sensitivity_curve
    panel_key <- panel_key + 1L
  }
}

curves <- bind_rows(curve_list)
support <- bind_rows(support_list)
histograms <- bind_rows(histogram_list)
rug_data <- bind_rows(rug_list)
partial_raw <- bind_rows(partial_raw_list)
partial_bins <- bind_rows(partial_bin_list)
partial_qa <- bind_rows(partial_qa_list)
physical_mapping <- bind_rows(physical_mapping_list)
sensitivity_curves <- bind_rows(sensitivity_curve_list)
influence_observations <- bind_rows(influence_by_phenotype)

influence_summary <- influence_observations %>%
  group_by(phenotype) %>%
  summarise(
    n_model = n(),
    effective_parameters = first(effective_parameters),
    leverage_cutoff = first(leverage_cutoff),
    cooks_cutoff = first(cooks_cutoff),
    n_high_leverage = sum(high_leverage),
    n_large_studentized_residual = sum(large_studentized_residual),
    n_high_cooks_distance = sum(high_cooks_distance),
    n_influence_excluded = sum(influence_exclusion_flag),
    fraction_influence_excluded = mean(influence_exclusion_flag),
    max_leverage = max(leverage),
    max_abs_studentized_residual = max(abs(studentized_residual)),
    max_cooks_distance = max(cooks_distance),
    .groups = "drop"
  )

if (any(partial_qa$n_model != partial_qa$n_partial_residual_rows) ||
    any(partial_qa$n_model != partial_qa$n_rows_in_bins) ||
    any(partial_qa$max_partial_residual_identity_error > 1e-12)) {
  stop("Adjusted partial-residual QA failed.", call. = FALSE)
}

physical_mapping_qa <- physical_mapping %>%
  mutate(
    mapping_pass = mapping_r_squared > 0.999999999 &
      max_abs_mapping_residual < 1e-6 & map_slope > 0
  )
if (any(!physical_mapping_qa$mapping_pass)) {
  stop("Climate-zone-specific physical-axis mapping QA failed.", call. = FALSE)
}

sensitivity_curve_qa <- map_dfr(v3_specs$phenotype, function(phenotype_value) {
  map_dfr(v3_factors$climate, function(climate_value) {
    main <- curves %>%
      filter(phenotype == phenotype_value, climate == climate_value) %>%
      arrange(x_scaled)
    sensitivity <- sensitivity_curves %>%
      filter(phenotype == phenotype_value, climate == climate_value) %>%
      arrange(x_scaled)
    supp <- support %>%
      filter(phenotype == phenotype_value, climate == climate_value)
    inside <- main$x_scaled >= supp$p2_5_scaled & main$x_scaled <= supp$p97_5_scaled
    delta_eta <- sensitivity$eta_aligned_to_main_at_P50 - main$eta
    tibble(
      phenotype = phenotype_value,
      climate = climate_value,
      n_grid = nrow(main),
      n_excluded = unique(sensitivity$n_excluded),
      max_abs_log_OR_difference_full_minmax = max(abs(delta_eta)),
      max_abs_log_OR_difference_P2_5_P97_5 = max(abs(delta_eta[inside])),
      median_abs_log_OR_difference_P2_5_P97_5 = median(abs(delta_eta[inside])),
      max_OR_ratio_P2_5_P97_5 = max(exp(delta_eta[inside])),
      min_OR_ratio_P2_5_P97_5 = min(exp(delta_eta[inside]))
    )
  })
})

# Exact curve preservation against the archived source data.
v3_curve_source <- read.csv(
  file.path(
    analysis_input_dir, "Figure2_reference_termcentred_curve_source_data.csv"
  ),
  check.names = FALSE
)
curve_qa <- map_dfr(v3_specs$phenotype, function(phenotype_value) {
  map_dfr(v3_factors$climate, function(climate_value) {
    current <- curves %>%
      filter(phenotype == phenotype_value, climate == climate_value) %>%
      arrange(x_scaled)
    archived <- v3_curve_source %>%
      filter(phenotype == phenotype_value, climate == climate_value) %>%
      arrange(x_scaled)
    if (nrow(current) != nrow(archived)) {
      return(tibble(
        phenotype = phenotype_value, climate = climate_value,
        n_compared = min(nrow(current), nrow(archived)),
        max_abs_x_difference = Inf,
        max_abs_OR_difference = Inf,
        max_abs_lower_difference = Inf,
        max_abs_upper_difference = Inf
      ))
    }
    tibble(
      phenotype = phenotype_value,
      climate = climate_value,
      n_compared = nrow(current),
      max_abs_x_difference = max(abs(current$x_scaled - archived$x_scaled)),
      max_abs_OR_difference = max(abs(current$OR - archived$OR)),
      max_abs_lower_difference = max(abs(current$Lower_95CI - archived$Lower_95CI)),
      max_abs_upper_difference = max(abs(current$Upper_95CI - archived$Upper_95CI))
    )
  })
})
if (any(curve_qa$n_compared != 801L) ||
    any(curve_qa$max_abs_x_difference > 1e-12) ||
    any(curve_qa$max_abs_OR_difference > 1e-12) ||
    any(curve_qa$max_abs_lower_difference > 1e-12) ||
    any(curve_qa$max_abs_upper_difference > 1e-12)) {
  stop("Frozen-curve preservation QA failed.", call. = FALSE)
}

bootstrap_turning_points <- read.csv(
  file.path(dirs$inputs, "Figure2_v3_bootstrap_support_for_observed_turning_points.csv"),
  check.names = FALSE
) %>%
  left_join(
    support %>%
      select(
        phenotype, climate, reference_zone, map_intercept, map_slope,
        unit_display, axis_title
      ),
    by = c("phenotype", "climate")
  ) %>%
  mutate(
    x_physical_equivalent = map_intercept + map_slope * x_scaled,
    location_median_physical_equivalent = map_intercept +
      map_slope * location_median_scaled,
    location_lower_95_physical_equivalent = map_intercept +
      map_slope * location_lower_95_scaled,
    location_upper_95_physical_equivalent = map_intercept +
      map_slope * location_upper_95_scaled,
    evidence_class = if_else(
      support_attempted >= 0.90,
      "Bootstrap-supported turning point (BS >=90%)",
      "Estimated turning point (BS <90%)"
    )
  )

strict_thresholds <- bootstrap_turning_points %>%
  filter(support_attempted >= 0.90)

display_thresholds <- bootstrap_turning_points %>%
  # Do not impose an unplanned 50% display threshold. Every geometrically
  # detected observed turning point is shown with its exact bootstrap support.
  arrange(phenotype, climate, desc(support_attempted), desc(prominence_eta))

# Display-only audit: all main-figure turning-point locations use two decimal
# places. Raw candidate locations remain unrounded in Source Data.
turning_point_label_precision_qa <- display_thresholds %>%
  mutate(
    location_display = mapply(
      format_threshold_location,
      x_physical_equivalent,
      climate,
      USE.NAMES = FALSE
    ),
    location_has_two_decimals = grepl(
      "\\.[0-9]{2} (°C|mm|d)$|\\.[0-9]{2}%$",
      location_display
    ),
    or_ci_display = sprintf(
      "%.2f (%.2f–%.2f)",
      observed_OR,
      observed_Lower_95CI,
      observed_Upper_95CI
    ),
    or_ci_has_two_decimals = grepl(
      "^-?[0-9]+\\.[0-9]{2} \\(-?[0-9]+\\.[0-9]{2}–-?[0-9]+\\.[0-9]{2}\\)$",
      or_ci_display
    )
  ) %>%
  select(
    phenotype, climate, type, x_physical_equivalent,
    location_display, location_has_two_decimals,
    observed_OR, observed_Lower_95CI, observed_Upper_95CI,
    or_ci_display, or_ci_has_two_decimals
  )

if (any(!turning_point_label_precision_qa$location_has_two_decimals) ||
    any(!turning_point_label_precision_qa$or_ci_has_two_decimals)) {
  stop("Turning-point label precision QA failed.", call. = FALSE)
}

# Panel-level audit explaining why a turning point is or is not available.
# The "without EDF gate" columns are diagnostic only and never feed the figure.
turning_point_panel_audit <- map_dfr(seq_len(nrow(v3_specs)), function(i) {
  phenotype_value <- v3_specs$phenotype[i]
  model <- models[[phenotype_value]]
  map_dfr(seq_len(nrow(v3_factors)), function(j) {
    factor_row <- v3_factors[j, ]
    variable <- factor_row$variable
    x <- model$model[[variable]]
    q <- as.numeric(quantile(x, c(0.025, 0.975), type = 7, names = FALSE))
    search_curve <- v3_term_curve(
      model, variable, seq(q[1], q[2], length.out = 1001L), with_se = FALSE
    )
    stat <- v3_smooth_statistics(model, variable)
    detected <- v3_detect_turning_points(search_curve, stat$edf)
    geometry_without_edf_gate <- v3_detect_turning_points(
      search_curve, edf = max(2, stat$edf), edf_cutoff = 1.05
    )
    d_eta <- diff(search_curve$eta) / diff(search_curve$x_scaled)
    derivative_tolerance <- max(1e-10, max(abs(d_eta), na.rm = TRUE) * 1e-5)
    monotonic_direction <- case_when(
      all(d_eta >= -derivative_tolerance) ~ "monotone increasing",
      all(d_eta <= derivative_tolerance) ~ "monotone decreasing",
      TRUE ~ "non-monotone"
    )
    panel_candidates <- bootstrap_turning_points %>%
      filter(
        .data$phenotype == !!phenotype_value,
        .data$climate == !!factor_row$climate
      )
    support_text <- if (nrow(panel_candidates) == 0L) "" else paste0(
      panel_candidates$type, " ",
      sprintf("%.2f%%", 100 * panel_candidates$support_attempted),
      collapse = "; "
    )
    display_complexity_class <- case_when(
      stat$edf < 0.10 ~ "shrunk towards zero",
      stat$edf <= 1.05 & nrow(geometry_without_edf_gate) > 0L ~
        "low-complexity curved fit",
      stat$edf <= 1.05 ~ "approximately linear",
      TRUE ~ "flexible smooth"
    )
    explanation <- case_when(
      nrow(detected) > 0L ~
        "Objective internal turning point(s) detected and displayed with exact bootstrap support",
      stat$edf <= 1.05 & nrow(geometry_without_edf_gate) > 0L ~
        paste0(
          "A derivative sign change is visible, but the fitted term met the fixed ",
          "low-complexity gate (EDF <=1.05); it is labelled as a low-complexity ",
          "curved fit and no nonlinear turning point is declared"
        ),
      stat$edf <= 1.05 ~
        paste0(
          "The fitted term met the fixed low-complexity/near-zero gate ",
          "(EDF <=1.05); no nonlinear turning point was searched"
        ),
      nrow(geometry_without_edf_gate) == 0L ~
        paste0(
          "No first-derivative sign change occurred within P2.5-P97.5 (",
          monotonic_direction, ")"
        ),
      TRUE ~ "No objective turning point retained"
    )
    tibble(
      phenotype = phenotype_value,
      climate = factor_row$climate,
      variable = variable,
      edf = stat$edf,
      smooth_p_value = stat$p_value,
      display_complexity_class = display_complexity_class,
      log_odds_amplitude_P2_5_P97_5 = diff(range(search_curve$eta)),
      monotonic_direction = monotonic_direction,
      n_detected_under_reported_algorithm = nrow(detected),
      detected_types = paste(detected$type, collapse = "; "),
      bootstrap_support_for_detected_candidates = support_text,
      n_geometric_sign_changes_without_edf_gate = nrow(geometry_without_edf_gate),
      diagnostic_types_without_edf_gate = paste(
        geometry_without_edf_gate$type, collapse = "; "
      ),
      explanation = explanation
    )
  })
})

# Targeted audit requested for the two wet-days GMax estimates that round to
# 117 d. The search grids and standardized candidate locations are retained at
# full precision so that an apparent equality after rounding cannot conceal a
# shared-grid error.
wet_day_phenotypes <- c("3GCR-Ec", "CR-Ab")
wet_day_grids <- setNames(map(wet_day_phenotypes, function(phenotype_value) {
  supp <- support %>%
    filter(phenotype == phenotype_value, climate == "WetDays")
  seq(supp$p2_5_scaled, supp$p97_5_scaled, length.out = 1001L)
}), wet_day_phenotypes)
wet_day_grids_identical <- identical(
  wet_day_grids[["3GCR-Ec"]], wet_day_grids[["CR-Ab"]]
)
wet_days_gmax_grid_audit <- map_dfr(wet_day_phenotypes, function(phenotype_value) {
  supp <- support %>%
    filter(phenotype == phenotype_value, climate == "WetDays")
  candidate <- bootstrap_turning_points %>%
    filter(
      phenotype == phenotype_value,
      climate == "WetDays",
      type == "GMax"
    ) %>%
    slice(1)
  grid <- wet_day_grids[[phenotype_value]]
  tibble(
    phenotype = phenotype_value,
    n_model = supp$n_model,
    observed_min_scaled = supp$observed_min_scaled,
    observed_max_scaled = supp$observed_max_scaled,
    p2_5_scaled = supp$p2_5_scaled,
    p97_5_scaled = supp$p97_5_scaled,
    search_grid_first_scaled = grid[1],
    search_grid_last_scaled = grid[length(grid)],
    search_grid_step_scaled = grid[2] - grid[1],
    search_grid_length = length(grid),
    gmax_x_scaled = candidate$x_scaled,
    gmax_physical_equivalent_days = candidate$x_physical_equivalent,
    gmax_display_rounded_days = round(candidate$x_physical_equivalent),
    map_intercept = supp$map_intercept,
    map_slope = supp$map_slope,
    search_grids_identical_exactly = wet_day_grids_identical
  )
})

# Verify that the CR-Pa temperature EDF and the plotted curve are extracted
# from the same fitted smooth term. EDF is an effective-complexity measure and
# does not mathematically constrain a penalized smooth to be an exact line.
crpa_model <- models[["CR-Pa"]]
crpa_variable <- "TMP_scaled_lag"
crpa_stat <- v3_smooth_statistics(crpa_model, crpa_variable)
crpa_supp <- support %>%
  filter(phenotype == "CR-Pa", climate == "Temperature")
crpa_grid <- seq(crpa_supp$p2_5_scaled, crpa_supp$p97_5_scaled, length.out = 1001L)
crpa_curve <- v3_term_curve(crpa_model, crpa_variable, crpa_grid, with_se = TRUE)
crpa_nd <- v3_prediction_template(crpa_model, crpa_variable, crpa_grid)
crpa_prediction <- predict(
  crpa_model, crpa_nd, type = "terms", se.fit = TRUE,
  unconditional = FALSE
)
crpa_term_index <- v3_target_term_index(crpa_prediction$fit, crpa_variable)
crpa_geometry_without_gate <- v3_detect_turning_points(
  crpa_curve, edf = max(2, crpa_stat$edf), edf_cutoff = 1.05
)
crpa_linear_fit <- lm(eta ~ x_scaled, data = crpa_curve)
crpa_derivative <- diff(crpa_curve$eta) / diff(crpa_curve$x_scaled)
crpa_temperature_edf_curve_audit <- tibble(
  phenotype = "CR-Pa",
  climate = "Temperature",
  smooth_term_from_summary = paste0("s(", crpa_variable, ")"),
  smooth_term_from_prediction = colnames(crpa_prediction$fit)[crpa_term_index],
  same_term_label = smooth_term_from_summary == smooth_term_from_prediction,
  max_abs_curve_difference_from_selected_term = max(abs(
    crpa_curve$eta - crpa_prediction$fit[, crpa_term_index]
  )),
  edf = crpa_stat$edf,
  smooth_p_value = crpa_stat$p_value,
  log_odds_amplitude_P2_5_P97_5 = diff(range(crpa_curve$eta)),
  OR_min_P2_5_P97_5 = min(crpa_curve$OR),
  OR_max_P2_5_P97_5 = max(crpa_curve$OR),
  linear_fit_R_squared = summary(crpa_linear_fit)$r.squared,
  max_abs_log_odds_deviation_from_linear_fit = max(abs(residuals(crpa_linear_fit))),
  derivative_min = min(crpa_derivative),
  derivative_max = max(crpa_derivative),
  n_geometric_turning_points_without_edf_gate = nrow(crpa_geometry_without_gate),
  geometric_type_without_edf_gate = paste(crpa_geometry_without_gate$type, collapse = "; "),
  geometric_x_scaled_without_edf_gate = if (nrow(crpa_geometry_without_gate) > 0L) {
    crpa_geometry_without_gate$x_scaled[1]
  } else NA_real_,
  geometric_physical_equivalent_without_edf_gate = if (
    nrow(crpa_geometry_without_gate) > 0L
  ) {
    crpa_supp$map_intercept +
      crpa_supp$map_slope * crpa_geometry_without_gate$x_scaled[1]
  } else NA_real_,
  reported_display_class = "low-complexity curved fit",
  interpretation = paste0(
    "Curve and EDF are from the same penalized smooth. The derivative changes ",
    "sign, but the fixed EDF<=1.05 gate excludes nonlinear turning-point ",
    "declaration; the panel must not be labelled approximately linear."
  )
)

# Supplementary curve-feature catalogue. These descriptors deliberately reuse
# the geometry of the original workflow (extrema, rapid change and stability)
# while removing phenotype-specific overrides, target locations, LOESS and
# priority scores. Chg and Stb are descriptive curve features; the 1,000-refit
# bootstrap support applies only to internal extrema in the primary algorithm.
v6_pointwise_summary <- function(curve, x0) {
  eta0 <- approx(curve$x_scaled, curve$eta, xout = x0, ties = "ordered")$y
  se0 <- approx(curve$x_scaled, curve$se, xout = x0, ties = "ordered")$y
  z0 <- if (is.finite(se0) && se0 > 0) eta0 / se0 else NA_real_
  p0 <- if (is.finite(z0)) 2 * pnorm(-abs(z0)) else NA_real_
  tibble(
    eta = eta0,
    se = se0,
    OR = exp(eta0),
    Lower_95CI = exp(eta0 - 1.96 * se0),
    Upper_95CI = exp(eta0 + 1.96 * se0),
    pointwise_p_normal = p0,
    pointwise_evidence = case_when(
      is.na(p0) ~ "not estimable",
      p0 < 0.05 ~ "Ppointwise <0.05",
      p0 < 0.10 ~ "0.05<=Ppointwise<0.10",
      TRUE ~ "Ppointwise>=0.10"
    ),
    T_flag = is.finite(p0) & p0 >= 0.10
  )
}

v6_contiguous_runs <- function(flag) {
  if (length(flag) == 0L) return(tibble(start = integer(), end = integer()))
  rr <- rle(flag)
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1L
  tibble(value = rr$values, start = starts, end = ends) %>%
    filter(value) %>%
    select(start, end)
}

v6_detect_supplementary_features <- function(model, phenotype, factor_row, supp) {
  variable <- factor_row$variable
  x_grid <- seq(supp$p2_5_scaled, supp$p97_5_scaled, length.out = 1001L)
  curve <- v3_term_curve(model, variable, x_grid, with_se = TRUE)
  eta_amplitude <- diff(range(curve$eta))
  if (!is.finite(supp$edf) || supp$edf <= 1.05 || eta_amplitude < 1e-4) {
    return(tibble())
  }

  x <- curve$x_scaled
  eta <- curve$eta
  derivative <- c(
    (eta[2] - eta[1]) / (x[2] - x[1]),
    (eta[3:length(eta)] - eta[1:(length(eta) - 2L)]) /
      (x[3:length(x)] - x[1:(length(x) - 2L)]),
    (eta[length(eta)] - eta[length(eta) - 1L]) /
      (x[length(x)] - x[length(x) - 1L])
  )
  abs_derivative <- abs(derivative)
  q10 <- as.numeric(quantile(abs_derivative, 0.10, type = 7, names = FALSE))
  q90 <- as.numeric(quantile(abs_derivative, 0.90, type = 7, names = FALSE))
  width <- diff(range(x))

  # One Chg point is selected from each contiguous high-slope region.
  chg_runs <- v6_contiguous_runs(abs_derivative >= q90)
  chg <- map_dfr(seq_len(nrow(chg_runs)), function(k) {
    idx_set <- seq.int(chg_runs$start[k], chg_runs$end[k])
    idx <- idx_set[which.max(abs_derivative[idx_set])]
    # A boundary maximum is retained but explicitly flagged in Source Data.
    pointwise <- v6_pointwise_summary(curve, x[idx])
    bind_cols(
      tibble(
        phenotype = phenotype,
        climate = factor_row$climate,
        variable = variable,
        feature_family = "change characteristic",
        canonical_type = "Chg",
        x_start_scaled = x[idx],
        x_end_scaled = x[idx],
        x_scaled = x[idx],
        derivative = derivative[idx],
        absolute_derivative = abs_derivative[idx],
        derivative_quantile_rule = "one maximum per contiguous |f'(x)|>=Q90 region",
        near_support_boundary =
          x[idx] <= min(x) + 0.05 * width || x[idx] >= max(x) - 0.05 * width
      ),
      pointwise
    )
  })
  if (nrow(chg) > 0L) {
    chg <- chg %>% mutate(
      display_type = if_else(T_flag, "TChg", "Chg")
    )
  }

  # Stb is an interval, not a point. Only contiguous low-slope regions spanning
  # at least 5% of the supported x range are retained; the midpoint is used
  # solely to report a representative term-centred OR.
  stb_runs <- v6_contiguous_runs(abs_derivative <= q10)
  stb <- map_dfr(seq_len(nrow(stb_runs)), function(k) {
    i1 <- stb_runs$start[k]
    i2 <- stb_runs$end[k]
    interval_width <- x[i2] - x[i1]
    if (!is.finite(interval_width) || interval_width < 0.05 * width) {
      return(tibble())
    }
    x0 <- mean(c(x[i1], x[i2]))
    pointwise <- v6_pointwise_summary(curve, x0)
    bind_cols(
      tibble(
        phenotype = phenotype,
        climate = factor_row$climate,
        variable = variable,
        feature_family = "stability interval",
        canonical_type = "Stb",
        x_start_scaled = x[i1],
        x_end_scaled = x[i2],
        x_scaled = x0,
        derivative = approx(x, derivative, xout = x0)$y,
        absolute_derivative = approx(x, abs_derivative, xout = x0)$y,
        derivative_quantile_rule =
          "contiguous |f'(x)|<=Q10 interval with width>=5% of P2.5-P97.5",
        near_support_boundary =
          x0 <= min(x) + 0.05 * width || x0 >= max(x) - 0.05 * width
      ),
      pointwise
    )
  })
  if (nrow(stb) > 0L) {
    stb <- stb %>% mutate(
      display_type = if_else(T_flag, "TStb", "Stb")
    )
  }

  bind_rows(chg, stb) %>%
    mutate(
      edf = supp$edf,
      smooth_p_value = supp$p_value,
      p2_5_scaled = supp$p2_5_scaled,
      p97_5_scaled = supp$p97_5_scaled,
      x_start_physical_equivalent = supp$map_intercept +
        supp$map_slope * x_start_scaled,
      x_end_physical_equivalent = supp$map_intercept +
        supp$map_slope * x_end_scaled,
      x_physical_equivalent = supp$map_intercept + supp$map_slope * x_scaled,
      unit_display = supp$unit_display,
      bootstrap_support = NA_real_,
      bootstrap_assessment =
        "not assessed; descriptive supplementary curve feature"
    )
}

supplementary_curve_features <- map_dfr(seq_len(nrow(v3_specs)), function(i) {
  phenotype_value <- v3_specs$phenotype[i]
  map_dfr(seq_len(nrow(v3_factors)), function(j) {
    supp <- support %>%
      filter(
        .data$phenotype == !!phenotype_value,
        .data$climate == !!v3_factors$climate[j]
      )
    v6_detect_supplementary_features(
      models[[phenotype_value]], phenotype_value, v3_factors[j, ], supp
    )
  })
})

extrema_feature_catalog <- bootstrap_turning_points %>%
  mutate(
    feature_family = "internal extremum",
    canonical_type = type,
    pointwise_se =
      (log(observed_Upper_95CI) - log(observed_Lower_95CI)) / (2 * 1.96),
    pointwise_p_normal = if_else(
      is.finite(pointwise_se) & pointwise_se > 0,
      2 * pnorm(-abs(eta / pointwise_se)), NA_real_
    ),
    pointwise_evidence = case_when(
      is.na(pointwise_p_normal) ~ "not estimable",
      pointwise_p_normal < 0.05 ~ "Ppointwise <0.05",
      pointwise_p_normal < 0.10 ~ "0.05<=Ppointwise<0.10",
      TRUE ~ "Ppointwise>=0.10"
    ),
    T_flag = is.finite(pointwise_p_normal) & pointwise_p_normal >= 0.10,
    display_type = if_else(T_flag, paste0("T", type), type),
    x_start_scaled = x_scaled,
    x_end_scaled = x_scaled,
    x_start_physical_equivalent = x_physical_equivalent,
    x_end_physical_equivalent = x_physical_equivalent,
    derivative = NA_real_,
    absolute_derivative = NA_real_,
    derivative_quantile_rule = "first-derivative sign change",
    near_support_boundary = FALSE,
    bootstrap_support = support_attempted,
    bootstrap_assessment = if_else(
      support_attempted >= 0.90,
      "met 900/1,000 stability criterion",
      "did not meet 900/1,000 stability criterion"
    )
  )

supplementary_full_feature_catalog <- bind_rows(
  extrema_feature_catalog %>%
    transmute(
      phenotype, climate, variable, feature_family, canonical_type,
      display_type,
      x_start_scaled, x_end_scaled, x_scaled,
      x_start_physical_equivalent, x_end_physical_equivalent,
      x_physical_equivalent, unit_display,
      OR = observed_OR, Lower_95CI = observed_Lower_95CI,
      Upper_95CI = observed_Upper_95CI,
      pointwise_p_normal, pointwise_evidence, T_flag,
      derivative, absolute_derivative, derivative_quantile_rule,
      near_support_boundary, bootstrap_support, bootstrap_assessment,
      conditional_location_lower_scaled = location_lower_95_scaled,
      conditional_location_upper_scaled = location_upper_95_scaled,
      n_observations_within_neighborhood,
      n_countries_within_neighborhood
    ),
  supplementary_curve_features %>%
    transmute(
      phenotype, climate, variable, feature_family, canonical_type,
      display_type,
      x_start_scaled, x_end_scaled, x_scaled,
      x_start_physical_equivalent, x_end_physical_equivalent,
      x_physical_equivalent, unit_display,
      OR, Lower_95CI, Upper_95CI,
      pointwise_p_normal, pointwise_evidence, T_flag,
      derivative, absolute_derivative, derivative_quantile_rule,
      near_support_boundary, bootstrap_support, bootstrap_assessment,
      conditional_location_lower_scaled = NA_real_,
      conditional_location_upper_scaled = NA_real_,
      n_observations_within_neighborhood = NA_integer_,
      n_countries_within_neighborhood = NA_integer_
    )
) %>%
  arrange(phenotype, climate, feature_family, x_scaled)

theme_zero_cascade <- function(base_size = 6.0) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "#303030", linewidth = 0.35),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.30, colour = "#303030"),
      axis.ticks.length = grid::unit(1.0, "mm"),
      axis.text = element_text(size = base_size - 0.25, colour = "#202020"),
      axis.title = element_text(size = base_size + 0.15, colour = "#202020"),
      strip.background = element_rect(fill = "white", colour = "#303030", linewidth = 0.35),
      strip.text = element_text(size = base_size + 0.1, face = "bold"),
      plot.margin = margin(2.0, 2.0, 0.5, 2.0, unit = "pt")
    )
}

format_p <- function(p) {
  case_when(
    !is.finite(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

make_curve_panel <- function(phenotype, factor_row, show_y = FALSE,
                             add_thresholds = FALSE,
                             show_partial_bins = TRUE) {
  climate <- factor_row$climate
  curve <- curves %>% filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  supp <- support %>% filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  hist <- histograms %>% filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  rug <- rug_data %>% filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  bins <- partial_bins %>% filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  threshold <- display_thresholds %>%
    filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  panel_audit <- turning_point_panel_audit %>%
    filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  colour <- unname(palette[[climate]])

  central <- curve %>% filter(tail_group == "P2.5-P97.5")
  tails <- curve %>% filter(tail_group != "P2.5-P97.5")
  central_bins <- bins %>% filter(support_class == "P2.5-P97.5")
  tail_bins <- bins %>% filter(support_class == "Observed tail")

  xr <- c(
    supp$observed_min_physical_equivalent,
    supp$observed_max_physical_equivalent
  )
  x_axis <- balanced_physical_x_axis(xr, climate)
  xlim <- x_axis$limits
  xbreaks <- x_axis$breaks

  partial_y <- if (isTRUE(show_partial_bins)) bins$partial_OR_geometric_mean else numeric()
  low_complexity_panel <- is.finite(supp$edf) && supp$edf <= 1.05
  if (low_complexity_panel) {
    # A common scale prevents near-null or low-complexity terms from
    # appearing structured solely because of panel-specific axis expansion.
    ylim <- c(0.5, 1.5)
    ybreaks <- seq(0.5, 1.5, length.out = 5L)
  } else {
    yr <- range(
      c(curve$Lower_95CI, curve$Upper_95CI, partial_y, 1),
      finite = TRUE
    )
    ypad <- max(diff(yr) * 0.075, 0.012)
    ylim <- c(max(0, yr[1] - ypad), yr[2] + ypad)
    ybreaks <- five_nice_breaks(ylim)
    ylim <- range(ybreaks)
  }
  if (length(ybreaks) != 5L) stop("Every curve panel must have exactly five y ticks.")

  main <- ggplot() +
    geom_hline(
      yintercept = 1, colour = palette[["reference"]],
      linetype = "dashed", linewidth = 0.30
    ) +
    geom_ribbon(
      data = tails,
      aes(x = x_physical_equivalent, ymin = Lower_95CI, ymax = Upper_95CI,
          group = tail_group),
      fill = palette[["tail_fill"]], alpha = 0.42, colour = NA
    ) +
    geom_ribbon(
      data = central,
      aes(x = x_physical_equivalent, ymin = Lower_95CI, ymax = Upper_95CI),
      fill = colour, alpha = 0.17, colour = NA
    ) +
    geom_line(
      data = tails,
      aes(x = x_physical_equivalent, y = OR, group = tail_group),
      colour = palette[["tail_line"]], linetype = "22", linewidth = 0.48
    ) +
    geom_line(
      data = central,
      aes(x = x_physical_equivalent, y = OR),
      colour = colour, linewidth = 0.60, lineend = "round"
    ) +
    geom_vline(
      xintercept = c(
        supp$p2_5_physical_equivalent,
        supp$p97_5_physical_equivalent
      ),
      colour = palette[["neutral_mid"]], linetype = "dashed", linewidth = 0.30
    ) +
    geom_point(
      data = if (isTRUE(show_partial_bins)) central_bins else central_bins[0, ],
      aes(x = x_mean_physical_equivalent, y = partial_OR_geometric_mean),
      shape = 21, size = 1.25, stroke = 0.32,
      fill = "white", colour = palette[["neutral_dark"]]
    ) +
    geom_point(
      data = if (isTRUE(show_partial_bins)) tail_bins else tail_bins[0, ],
      aes(x = x_mean_physical_equivalent, y = partial_OR_geometric_mean),
      shape = 21, size = 1.20, stroke = 0.30,
      fill = "white", colour = palette[["neutral_mid"]], alpha = 0.72
    ) +
    annotate(
      "text",
      x = c(
        supp$p2_5_physical_equivalent,
        supp$p97_5_physical_equivalent
      ),
      y = ylim[2], label = c("P2.5", "P97.5"),
      colour = palette[["neutral_mid"]], size = 1.60,
      vjust = 1.05, hjust = c(-0.06, 1.06), family = "Arial"
    ) +
    scale_x_continuous(
      breaks = xbreaks, labels = NULL,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = ylim, breaks = ybreaks,
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = xlim, clip = "on") +
    labs(
      x = NULL,
      y = if (show_y) {
        paste0(
          phenotype, "\n(n = ", supp$n_model,
          "; ", supp$n_countries, " countries)\nOR (95% CI)"
        )
      } else NULL
    ) +
    theme_zero_cascade() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = if (show_y) element_text(lineheight = 0.90) else element_blank(),
      plot.margin = margin(2.0, 8.0, 0, 2.0, unit = "pt")
    )

  if (low_complexity_panel) {
    low_complexity_label <- if (supp$edf < 0.10) {
      "EDF ≈ 0\nshrunk towards zero"
    } else if (panel_audit$n_geometric_sign_changes_without_edf_gate > 0L) {
      paste0("low-complexity curved fit\nEDF = ", sprintf("%.2f", supp$edf))
    } else {
      paste0("approximately linear\nEDF = ", sprintf("%.2f", supp$edf))
    }
    main <- main + annotate(
      "text", x = mean(xr), y = 1.32, label = low_complexity_label,
      family = "Arial", size = 1.55, lineheight = 0.95,
      colour = palette[["neutral_mid"]]
    )
  }

  if (isTRUE(add_thresholds) && nrow(threshold) > 0L) {
    threshold <- threshold %>%
      mutate(
        label_short = paste0(
          type, " ",
          map_chr(
            seq_along(x_physical_equivalent),
            ~ format_threshold_location(
              x_physical_equivalent[.x], .env$climate
            )
          ), "\n",
          "OR ", sprintf("%.2f", observed_OR), " (",
          sprintf("%.2f", observed_Lower_95CI), "–",
          sprintf("%.2f", observed_Upper_95CI), ")"
        ),
        label_y = case_when(
          direction == "Max" & observed_OR <= ylim[2] - 0.20 * diff(ylim) ~
            observed_OR + 0.18 * diff(ylim),
          direction == "Max" ~
            observed_OR - 0.22 * diff(ylim),
          direction == "Min" & observed_OR >= ylim[1] + 0.18 * diff(ylim) ~
            observed_OR - 0.18 * diff(ylim),
          TRUE ~
            observed_OR + 0.18 * diff(ylim)
        )
      )
    threshold_robust <- threshold %>%
      filter(evidence_class == "Bootstrap-supported turning point (BS >=90%)")
    threshold_candidate <- threshold %>%
      filter(evidence_class != "Bootstrap-supported turning point (BS >=90%)")
    main <- main +
      geom_vline(
        data = threshold, aes(xintercept = x_physical_equivalent),
        colour = palette[["threshold"]], linetype = "longdash", linewidth = 0.34
      ) +
      geom_point(
        data = threshold_robust, aes(x = x_physical_equivalent, y = observed_OR),
        shape = 21, size = 1.35, stroke = 0.35,
        fill = palette[["threshold"]], colour = palette[["threshold"]]
      ) +
      geom_point(
        data = threshold_candidate,
        aes(x = x_physical_equivalent, y = observed_OR),
        shape = 21, size = 1.25, stroke = 0.35,
        fill = "white", colour = palette[["threshold"]]
      ) +
      ggrepel::geom_text_repel(
        data = threshold,
        aes(x = x_physical_equivalent, y = label_y, label = label_short),
        family = "Arial", size = 1.34, colour = palette[["neutral_dark"]],
        lineheight = 0.90, direction = "y", seed = 20260808,
        box.padding = 0.12, point.padding = 0.08,
        min.segment.length = 0, segment.size = 0.16,
        segment.colour = palette[["neutral_mid"]],
        force = 1.5, max.overlaps = Inf, max.time = 3,
        ylim = ylim
      )
  }

  histogram_panel <- ggplot(hist) +
    geom_rect(
      data = hist %>% filter(support_class == "Observed tail"),
      aes(xmin = xmin_physical_equivalent, xmax = xmax_physical_equivalent,
          ymin = 0, ymax = count),
      fill = palette[["tail_fill"]], colour = "white", linewidth = 0.08
    ) +
    geom_rect(
      data = hist %>% filter(support_class == "P2.5-P97.5"),
      aes(xmin = xmin_physical_equivalent, xmax = xmax_physical_equivalent,
          ymin = 0, ymax = count),
      fill = colour, alpha = 0.58, colour = "white", linewidth = 0.08
    ) +
    geom_rug(
      data = rug,
      aes(x = x_physical_equivalent),
      sides = "b", colour = palette[["neutral_dark"]],
      alpha = 0.24, linewidth = 0.12, length = grid::unit(0.7, "mm")
    ) +
    geom_rug(
      data = rug %>% filter(influence_exclusion_flag),
      aes(x = x_physical_equivalent),
      sides = "b", colour = palette[["threshold"]],
      alpha = 0.90, linewidth = 0.28, length = grid::unit(1.05, "mm")
    ) +
    geom_vline(
      xintercept = c(
        supp$p2_5_physical_equivalent,
        supp$p97_5_physical_equivalent
      ),
      colour = palette[["neutral_mid"]], linetype = "dashed", linewidth = 0.28
    ) +
    scale_x_continuous(
      breaks = xbreaks,
      labels = label_number(accuracy = 1, big.mark = ","),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(xlim = xlim, clip = "on") +
    labs(
      x = paste0(
        supp$axis_title, " (", supp$unit_display,
        "; lag ", supp$lag_years, " yr)"
      ),
      y = NULL
    ) +
    theme_zero_cascade() +
    theme(
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.30, colour = "#303030"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 5.45, margin = margin(t = 1.0, unit = "pt")),
      plot.margin = margin(0, 8.0, 2.0, 2.0, unit = "pt")
    )

  cowplot::plot_grid(
    main, histogram_panel,
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(0.69, 0.31)
  )
}

make_support_key <- function(show_partial_bins = TRUE) {
  key <- ggplot() +
    annotate("segment", x = 0.015, xend = 0.07, y = 0.73, yend = 0.73,
             colour = palette[["Temperature"]], linewidth = 0.62) +
    annotate("text", x = 0.08, y = 0.73, label = "P2.5–P97.5 model curve",
             hjust = 0, size = 1.72, family = "Arial") +
    annotate("segment", x = 0.28, xend = 0.335, y = 0.73, yend = 0.73,
             colour = palette[["tail_line"]], linetype = "22", linewidth = 0.52) +
    annotate("text", x = 0.345, y = 0.73, label = "Observed tail (limited support)",
             hjust = 0, size = 1.72, family = "Arial") +
    annotate("segment", x = 0.59, xend = 0.59, y = 0.57, yend = 0.88,
             colour = palette[["threshold"]], linewidth = 0.42) +
    annotate("text", x = 0.605, y = 0.73,
             label = "Red rug: influence flag",
             hjust = 0, size = 1.72, family = "Arial") +
    annotate("rect", xmin = 0.015, xmax = 0.035, ymin = 0.10, ymax = 0.44,
             fill = palette[["Precipitation"]], alpha = 0.58) +
    annotate("segment", x = 0.01, xend = 0.04, y = 0.08, yend = 0.08,
             colour = palette[["neutral_dark"]], linewidth = 0.18) +
    annotate("text", x = 0.045, y = 0.27, label = "Histogram + all-observation rug",
             hjust = 0, size = 1.72, family = "Arial") +
    annotate("segment", x = 0.36, xend = 0.36, y = 0.08, yend = 0.45,
             colour = palette[["threshold"]], linetype = "longdash", linewidth = 0.42) +
    annotate("text", x = 0.375, y = 0.27,
             label = "Estimated turning point",
             hjust = 0, size = 1.72, family = "Arial") +
    annotate("point", x = 0.59, y = 0.27, shape = 21, size = 1.35,
             stroke = 0.34, fill = palette[["threshold"]], colour = palette[["threshold"]]) +
    annotate("text", x = 0.605, y = 0.27,
             label = "Filled: met 900/1,000 stability criterion",
             hjust = 0, size = 1.72, family = "Arial") +
    annotate("point", x = 0.84, y = 0.27, shape = 21, size = 1.35,
             stroke = 0.34, fill = "white", colour = palette[["threshold"]]) +
    annotate("text", x = 0.855, y = 0.27,
             label = "Open: criterion not met",
             hjust = 0, size = 1.72, family = "Arial") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(0, 0, 0, 0))
  if (isTRUE(show_partial_bins)) {
    key <- key +
      annotate("point", x = 0.82, y = 0.73, shape = 21, size = 1.45,
               stroke = 0.32, fill = "white", colour = palette[["neutral_dark"]]) +
      annotate("text", x = 0.835, y = 0.73,
               label = "Adjusted partial-residual bin mean",
               hjust = 0, size = 1.72, family = "Arial")
  }
  key
}

# The first climate column previously carried the complete rotated phenotype
# label inside its own plot object. Because cowplot allocated equal outer widths
# to the four plot objects, that label reduced only the Temperature panel's data
# viewport. Move the shared row label into a dedicated 4%-width column so all
# four climate panels can use identical outer and inner dimensions.
row_label_fraction <- 0.04

make_external_row_label <- function(label_text, size = 6.15) {
  cowplot::ggdraw() +
    cowplot::draw_label(
      label_text,
      x = 0.5, y = 0.5,
      hjust = 0.5, vjust = 0.5,
      angle = 90,
      fontfamily = "Arial",
      size = size,
      lineheight = 0.90,
      colour = "#202020"
    )
}

make_equal_panel_title_row <- function(titles) {
  title_grid <- cowplot::plot_grid(
    plotlist = map(titles, function(x) {
      ggplot() + theme_void(base_family = "Arial") +
        ggtitle(x) +
        theme(plot.title = element_text(hjust = 0.5, size = 7.0, face = "bold"))
    }),
    ncol = 4
  )
  cowplot::plot_grid(
    ggplot() + theme_void(), title_grid,
    ncol = 2,
    rel_widths = c(row_label_fraction, 1 - row_label_fraction)
  )
}

make_equal_panel_row <- function(plot_list, row_label_text, label_size = 6.15) {
  panel_grid <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 4,
    align = "h",
    axis = "tb"
  )
  cowplot::plot_grid(
    make_external_row_label(row_label_text, size = label_size), panel_grid,
    ncol = 2,
    rel_widths = c(row_label_fraction, 1 - row_label_fraction)
  )
}

build_figure_a <- function(add_thresholds = FALSE,
                           show_partial_bins = TRUE) {
  panel_list <- list()
  for (i in seq_len(nrow(v3_specs))) {
    phenotype <- v3_specs$phenotype[i]
    for (j in seq_len(nrow(v3_factors))) {
      key <- paste(phenotype, v3_factors$climate[j], sep = "__")
      panel_list[[key]] <- make_curve_panel(
        phenotype, v3_factors[j, ], show_y = FALSE,
        add_thresholds = add_thresholds,
        show_partial_bins = show_partial_bins
      )
    }
  }

  title_row <- make_equal_panel_title_row(v3_factors$title)
  phenotype_rows <- map(v3_specs$phenotype, function(phenotype) {
    keys <- paste(phenotype, v3_factors$climate, sep = "__")
    panel_support <- support %>%
      filter(.data$phenotype == !!phenotype) %>%
      slice(1)
    row_label_text <- paste0(
      phenotype, "\n(n = ", panel_support$n_model,
      "; ", panel_support$n_countries, " countries)\nOR (95% CI)"
    )
    make_equal_panel_row(panel_list[keys], row_label_text)
  })
  body <- cowplot::plot_grid(
    plotlist = c(list(title_row), phenotype_rows),
    ncol = 1,
    rel_heights = c(0.10, rep(1, nrow(v3_specs)))
  )
  cowplot::plot_grid(
    make_support_key(show_partial_bins = show_partial_bins), body,
    ncol = 1, rel_heights = c(0.20, 3.80)
  )
}

stats_df <- support %>%
  mutate(
    phenotype = factor(phenotype, levels = v3_specs$phenotype),
    climate = factor(climate, levels = rev(v3_factors$climate)),
    p_symbol = format_p(p_value),
    label = paste0(
      "F=", formatC(statistic, format = "f", digits = 2), p_symbol,
      "\nEDF=", formatC(edf, format = "f", digits = 2)
    )
  )

f_max <- max(stats_df$statistic, na.rm = TRUE)
figure_b <- ggplot(
  stats_df,
  aes(x = statistic, y = climate, colour = climate, shape = climate, size = edf)
) +
  geom_segment(
    aes(x = 0, xend = statistic, yend = climate),
    colour = "#E2E2E2", linewidth = 0.32
  ) +
  geom_point(alpha = 0.95) +
  geom_text(
    aes(label = label), hjust = 0,
    nudge_x = 0.030 * f_max, size = 1.72,
    colour = palette[["neutral_dark"]], family = "Arial", lineheight = 0.90
  ) +
  facet_grid(. ~ phenotype) +
  scale_colour_manual(values = palette[v3_factors$climate], guide = "none") +
  scale_shape_manual(
    values = c(Temperature = 16, Humidity = 17, Precipitation = 18, WetDays = 15),
    guide = "none"
  ) +
  scale_size_continuous(range = c(1.1, 3.0), guide = "none") +
  scale_x_continuous(
    limits = c(0, f_max * 1.31),
    breaks = pretty(c(0, f_max), n = 4), expand = expansion(mult = 0)
  ) +
  labs(x = "F statistic", y = NULL) +
  theme_zero_cascade(base_size = 6.0) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.28),
    axis.text.y = element_text(size = 5.4),
    strip.text = element_text(size = 6.1, face = "bold"),
    legend.position = "none",
    plot.margin = margin(2, 4, 2, 4, unit = "pt")
  )

influence_plot_data <- influence_observations %>%
  mutate(
    phenotype = factor(phenotype, levels = v3_specs$phenotype),
    flag_label = if_else(influence_exclusion_flag, "Influence-exclusion flag", "Not flagged")
  )
influence_labels <- influence_plot_data %>%
  group_by(phenotype) %>%
  slice_max(cooks_distance, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = paste0(country, "–", year))

figure_influence <- ggplot(
  influence_plot_data,
  aes(x = leverage, y = studentized_residual)
) +
  geom_hline(yintercept = c(-3, 3), colour = palette[["neutral_mid"]],
             linetype = "dashed", linewidth = 0.32) +
  geom_vline(
    data = influence_summary,
    aes(xintercept = leverage_cutoff),
    colour = palette[["neutral_mid"]], linetype = "dotted", linewidth = 0.32
  ) +
  geom_point(
    aes(size = cooks_distance, colour = flag_label),
    alpha = 0.70, stroke = 0
  ) +
  ggrepel::geom_text_repel(
    data = influence_labels,
    aes(label = label),
    size = 2.05, family = "Arial", seed = 20260808,
    box.padding = 0.25, point.padding = 0.12,
    min.segment.length = 0, segment.size = 0.22,
    max.overlaps = Inf
  ) +
  facet_wrap(~ phenotype, ncol = 3, scales = "free_x") +
  scale_colour_manual(
    values = c(
      "Not flagged" = "#A8A8A8",
      "Influence-exclusion flag" = palette[["threshold"]]
    ),
    name = NULL
  ) +
  scale_size_continuous(range = c(0.7, 4.2), name = "Approx. Cook distance") +
  labs(
    x = "Penalized-smoother leverage",
    y = "Internally studentized response residual"
  ) +
  theme_zero_cascade(base_size = 7.0) +
  theme(
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.25),
    panel.border = element_rect(fill = NA, colour = "#303030", linewidth = 0.35),
    strip.text = element_text(size = 7.0, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

make_sensitivity_panel <- function(phenotype, factor_row, show_y = FALSE) {
  climate <- factor_row$climate
  colour <- unname(palette[[climate]])
  main <- curves %>%
    filter(.data$phenotype == !!phenotype, .data$climate == !!climate) %>%
    arrange(x_scaled)
  sensitivity <- sensitivity_curves %>%
    filter(.data$phenotype == !!phenotype, .data$climate == !!climate) %>%
    arrange(x_scaled)
  supp <- support %>%
    filter(.data$phenotype == !!phenotype, .data$climate == !!climate)
  qa <- sensitivity_curve_qa %>%
    filter(.data$phenotype == !!phenotype, .data$climate == !!climate)

  reference_x <- unique(sensitivity$alignment_reference_x_scaled)
  main_reference_eta <- approx(
    main$x_scaled, main$eta, xout = reference_x, ties = "ordered"
  )$y
  main <- main %>% mutate(OR_relative_P50 = exp(eta - main_reference_eta))
  sensitivity <- sensitivity %>%
    mutate(OR_relative_P50 = exp(eta_aligned_to_main_at_P50 - main_reference_eta))

  xr <- range(main$x_physical_equivalent)
  x_axis <- balanced_physical_x_axis(xr, climate)
  xlim <- x_axis$limits
  xbreaks <- x_axis$breaks
  yr <- range(c(main$OR_relative_P50, sensitivity$OR_relative_P50, 1), finite = TRUE)
  ypad <- max(diff(yr) * 0.07, 0.015)
  ylim <- c(max(0, yr[1] - ypad), yr[2] + ypad)
  ybreaks <- five_nice_breaks(ylim)
  ylim <- range(ybreaks)

  ggplot() +
    annotate(
      "rect", xmin = xr[1], xmax = supp$p2_5_physical_equivalent,
      ymin = -Inf, ymax = Inf, fill = "#EFEFEF", alpha = 0.55
    ) +
    annotate(
      "rect", xmin = supp$p97_5_physical_equivalent, xmax = xr[2],
      ymin = -Inf, ymax = Inf, fill = "#EFEFEF", alpha = 0.55
    ) +
    geom_hline(yintercept = 1, colour = palette[["reference"]],
               linetype = "dashed", linewidth = 0.28) +
    geom_line(
      data = main,
      aes(x = x_physical_equivalent, y = OR_relative_P50),
      colour = colour, linewidth = 0.58
    ) +
    geom_line(
      data = sensitivity,
      aes(x = x_physical_equivalent, y = OR_relative_P50),
      colour = palette[["neutral_dark"]], linetype = "dotdash", linewidth = 0.42
    ) +
    geom_vline(
      xintercept = c(supp$p2_5_physical_equivalent, supp$p97_5_physical_equivalent),
      colour = palette[["neutral_mid"]], linetype = "dashed", linewidth = 0.28
    ) +
    annotate(
      "text", x = mean(xr), y = ylim[2],
      label = paste0(
        "max |ΔlogOR| within P2.5–P97.5 = ",
        formatC(qa$max_abs_log_OR_difference_P2_5_P97_5,
                format = "f", digits = 2)
      ),
      size = 1.55, vjust = 1.05, colour = palette[["neutral_mid"]],
      family = "Arial"
    ) +
    scale_x_continuous(
      breaks = xbreaks,
      labels = label_number(accuracy = 1, big.mark = ","),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = ylim, breaks = ybreaks,
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = xlim, clip = "on") +
    labs(
      x = paste0(supp$axis_title, " (", supp$unit_display, ")"),
      y = if (show_y) paste0(phenotype, "\nOR relative to P50") else NULL
    ) +
    theme_zero_cascade(base_size = 5.7) +
    theme(
      axis.title.x = element_text(size = 5.1),
      axis.title.y = if (show_y) element_text(lineheight = 0.9) else element_blank(),
      plot.margin = margin(2, 8, 2, 2, unit = "pt")
    )
}

build_sensitivity_figure <- function() {
  panels <- list()
  for (i in seq_len(nrow(v3_specs))) {
    phenotype <- v3_specs$phenotype[i]
    for (j in seq_len(nrow(v3_factors))) {
      key <- paste(phenotype, v3_factors$climate[j], sep = "__")
      panels[[key]] <- make_sensitivity_panel(
        phenotype, v3_factors[j, ], show_y = FALSE
      )
    }
  }
  title_row <- make_equal_panel_title_row(v3_factors$title)
  rows <- map(v3_specs$phenotype, function(phenotype) {
    keys <- paste(phenotype, v3_factors$climate, sep = "__")
    make_equal_panel_row(
      panels[keys],
      paste0(phenotype, "\nOR relative to P50"),
      label_size = 5.85
    )
  })
  key <- ggplot() +
    annotate("segment", x = 0.18, xend = 0.28, y = 0.5, yend = 0.5,
             colour = palette[["Temperature"]], linewidth = 0.62) +
    annotate("text", x = 0.30, y = 0.5, label = "Frozen Model C",
             hjust = 0, size = 2.0, family = "Arial") +
    annotate("segment", x = 0.56, xend = 0.66, y = 0.5, yend = 0.5,
             colour = palette[["neutral_dark"]], linetype = "dotdash", linewidth = 0.45) +
    annotate("text", x = 0.68, y = 0.5,
             label = "Influence-excluded sensitivity (aligned at P50)",
             hjust = 0, size = 2.0, family = "Arial") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = "Arial")
  body <- cowplot::plot_grid(
    plotlist = c(list(title_row), rows), ncol = 1,
    rel_heights = c(0.10, rep(1, nrow(v3_specs)))
  )
  cowplot::plot_grid(key, body, ncol = 1, rel_heights = c(0.12, 3.88))
}

figure_sensitivity <- build_sensitivity_figure()

figure_a_clean <- build_figure_a(add_thresholds = FALSE)
figure_a_thresholds <- build_figure_a(add_thresholds = TRUE)
figure_a_clean_no_partial_bins <- build_figure_a(
  add_thresholds = FALSE, show_partial_bins = FALSE
)
figure_a_thresholds_no_partial_bins <- build_figure_a(
  add_thresholds = TRUE, show_partial_bins = FALSE
)

combined_clean <- cowplot::plot_grid(
  figure_a_clean, figure_b,
  labels = c("a", "b"),
  label_fontfamily = "Arial", label_fontface = "bold", label_size = 9,
  label_x = 0.003, ncol = 1, rel_heights = c(4.8, 1.0)
)
combined_thresholds <- cowplot::plot_grid(
  figure_a_thresholds, figure_b,
  labels = c("a", "b"),
  label_fontfamily = "Arial", label_fontface = "bold", label_size = 9,
  label_x = 0.003, ncol = 1, rel_heights = c(4.8, 1.0)
)
combined_clean_no_partial_bins <- cowplot::plot_grid(
  figure_a_clean_no_partial_bins, figure_b,
  labels = c("a", "b"),
  label_fontfamily = "Arial", label_fontface = "bold", label_size = 9,
  label_x = 0.003, ncol = 1, rel_heights = c(4.8, 1.0)
)
combined_thresholds_no_partial_bins <- cowplot::plot_grid(
  figure_a_thresholds_no_partial_bins, figure_b,
  labels = c("a", "b"),
  label_fontfamily = "Arial", label_fontface = "bold", label_size = 9,
  label_x = 0.003, ncol = 1, rel_heights = c(4.8, 1.0)
)

save_pub <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width, height = height,
    family = "Arial", onefile = TRUE
  )
  suppressWarnings(print(plot))
  dev.off()
  svglite::svglite(
    paste0(stem, ".svg"), width = width, height = height,
    system_fonts = list(sans = "Arial")
  )
  suppressWarnings(print(plot))
  dev.off()
  ragg::agg_tiff(
    paste0(stem, "_600dpi.tiff"), width = width, height = height,
    units = "in", res = dpi, background = "white", compression = "lzw"
  )
  suppressWarnings(print(plot))
  dev.off()
  ragg::agg_png(
    paste0(stem, ".png"), width = width, height = height,
    units = "in", res = 300, background = "white"
  )
  suppressWarnings(print(plot))
  dev.off()
}

stem_a_clean <- file.path(
  dirs$figures,
  "Figure2a_physical_axes_support_influence_clean"
)
stem_a_thresholds <- file.path(
  dirs$figures,
  "Figure2a_turning_points_termcentred_OR_95CI_with_support_symbols"
)
stem_combined_clean <- file.path(
  dirs$figures,
  "Figure2_physical_axes_support_influence_clean"
)
stem_combined_thresholds <- file.path(
  dirs$figures,
  "Figure2_turning_points_termcentred_OR_95CI_with_support_symbols"
)
stem_a_clean_no_partial_bins <- file.path(
  dirs$figures,
  "Figure2a_physical_axes_support_clean_without_partial_residual_bins"
)
stem_a_thresholds_no_partial_bins <- file.path(
  dirs$figures,
  "Figure2a_main_turning_points_ORCI_without_partial_residual_bins"
)
stem_combined_clean_no_partial_bins <- file.path(
  dirs$figures,
  "Figure2_physical_axes_support_clean_without_partial_residual_bins"
)
stem_combined_thresholds_no_partial_bins <- file.path(
  dirs$figures,
  "Figure2_main_turning_points_ORCI_without_partial_residual_bins"
)
stem_influence <- file.path(
  dirs$figures,
  "Supplementary_Figure_observation_leverage_and_influence"
)
stem_sensitivity <- file.path(
  dirs$figures,
  "Supplementary_Figure_influence_excluded_curve_sensitivity"
)

# v9 is a focused manuscript-figure update. Companion and diagnostic graphics
# remain available in v8; only the recommended Figure 2 is re-exported here.
save_pub(combined_thresholds_no_partial_bins, stem_combined_thresholds_no_partial_bins, 183, 245)

# Reproduce the recommended primary figure's y-axis calculation as a tabular
# QA record. This asserts the five-tick contract independently of the rendered
# artwork and makes every panel's limits auditable.
main_y_axis_five_tick_qa <- map_dfr(seq_len(nrow(v3_specs)), function(i) {
  phenotype_value <- v3_specs$phenotype[i]
  map_dfr(seq_len(nrow(v3_factors)), function(j) {
    climate_value <- v3_factors$climate[j]
    panel_curve <- curves %>%
      filter(phenotype == phenotype_value, climate == climate_value)
    panel_support <- support %>%
      filter(phenotype == phenotype_value, climate == climate_value)
    if (panel_support$edf <= 1.05) {
      breaks <- seq(0.5, 1.5, length.out = 5L)
    } else {
      yr <- range(
        c(panel_curve$Lower_95CI, panel_curve$Upper_95CI, 1),
        finite = TRUE
      )
      ypad <- max(diff(yr) * 0.075, 0.012)
      breaks <- five_nice_breaks(c(max(0, yr[1] - ypad), yr[2] + ypad))
    }
    tibble(
      phenotype = phenotype_value,
      climate = climate_value,
      edf = panel_support$edf,
      n_y_ticks = length(breaks),
      y_tick_1 = breaks[1],
      y_tick_2 = breaks[2],
      y_tick_3 = breaks[3],
      y_tick_4 = breaks[4],
      y_tick_5 = breaks[5],
      passed_exactly_five_ticks = length(breaks) == 5L
    )
  })
})
if (any(!main_y_axis_five_tick_qa$passed_exactly_five_ticks)) {
  stop("Five-tick y-axis QA failed.", call. = FALSE)
}

# Independently audit the panel-specific x-axis display contract. The observed
# curve domain must remain inside the display limits, and terminal axis segments
# must be symmetric unless both are exactly zero on a natural zero-bounded axis.
main_x_axis_balanced_qa <- map_dfr(seq_len(nrow(support)), function(i) {
  panel_support <- support[i, ]
  xr <- c(
    panel_support$observed_min_physical_equivalent,
    panel_support$observed_max_physical_equivalent
  )
  axis <- balanced_physical_x_axis(xr, panel_support$climate)
  left_terminal_gap <- axis$breaks[1] - axis$limits[1]
  right_terminal_gap <- axis$limits[2] - tail(axis$breaks, 1)
  covers_observed_range <-
    axis$limits[1] <= xr[1] + 1e-9 && axis$limits[2] >= xr[2] - 1e-9
  balanced_terminal_gaps <-
    abs(left_terminal_gap - right_terminal_gap) <= 1e-9

  tibble(
    phenotype = panel_support$phenotype,
    climate = panel_support$climate,
    observed_min = xr[1],
    observed_max = xr[2],
    display_min = axis$limits[1],
    display_max = axis$limits[2],
    tick_step = axis$step,
    n_x_ticks = length(axis$breaks),
    x_ticks = paste(
      format(axis$breaks, trim = TRUE, scientific = FALSE),
      collapse = " | "
    ),
    left_terminal_gap = left_terminal_gap,
    right_terminal_gap = right_terminal_gap,
    display_expansion_percent = 100 * axis$expansion_fraction,
    covers_observed_range = covers_observed_range,
    balanced_terminal_gaps = balanced_terminal_gaps,
    passed = covers_observed_range && balanced_terminal_gaps
  )
})
if (any(!main_x_axis_balanced_qa$passed)) {
  stop("Balanced x-axis QA failed.", call. = FALSE)
}

write.csv(curves, file.path(dirs$source, "Figure2_termcentred_curve_source_data.csv"), row.names = FALSE)
write.csv(histograms, file.path(dirs$source, "Figure2_histogram_source_data.csv"), row.names = FALSE)
write.csv(rug_data, file.path(dirs$source, "Figure2_rug_observations_source_data.csv"), row.names = FALSE)
write.csv(partial_raw, file.path(dirs$source, "Figure2_adjusted_partial_residual_observations.csv"), row.names = FALSE)
write.csv(partial_bins, file.path(dirs$source, "Figure2_adjusted_partial_residual_bins.csv"), row.names = FALSE)
write.csv(sensitivity_curves, file.path(dirs$source, "Figure2_influence_excluded_sensitivity_curves.csv"), row.names = FALSE)
write.csv(support, file.path(dirs$tables, "Figure2_panel_minmax_P2_5_P97_5_and_model_statistics.csv"), row.names = FALSE)
write.csv(strict_thresholds, file.path(dirs$tables, "Figure2_v3_strict_thresholds_with_physical_crosswalk_and_local_support.csv"), row.names = FALSE)
write.csv(display_thresholds, file.path(dirs$tables, "Figure2_displayed_all_objectively_detected_turning_points.csv"), row.names = FALSE)
write.csv(bootstrap_turning_points, file.path(dirs$tables, "Figure2_all_observed_turning_points_with_1000_refit_support.csv"), row.names = FALSE)
write.csv(turning_point_panel_audit, file.path(dirs$tables, "Figure2_turning_point_detection_status_for_all_24_panels.csv"), row.names = FALSE)
write.csv(wet_days_gmax_grid_audit, file.path(dirs$diagnostics, "Figure2_wet_days_117d_GMax_grid_and_scaled_location_audit.csv"), row.names = FALSE)
write.csv(crpa_temperature_edf_curve_audit, file.path(dirs$diagnostics, "Figure2_CR_Pa_temperature_EDF_curve_term_consistency_audit.csv"), row.names = FALSE)
write.csv(main_y_axis_five_tick_qa, file.path(dirs$diagnostics, "Figure2_main_y_axis_exactly_five_ticks_QA.csv"), row.names = FALSE)
write.csv(main_x_axis_balanced_qa, file.path(dirs$diagnostics, "Figure2_balanced_physical_x_axis_QA.csv"), row.names = FALSE)
write.csv(turning_point_label_precision_qa, file.path(dirs$diagnostics, "Figure2_turning_point_label_two_decimal_QA.csv"), row.names = FALSE)
write.csv(supplementary_curve_features, file.path(dirs$tables, "Supplementary_curve_features_Chg_Stb_and_T_flags.csv"), row.names = FALSE)
write.csv(supplementary_full_feature_catalog, file.path(dirs$tables, "Supplementary_full_threshold_and_curve_feature_catalog.csv"), row.names = FALSE)
write.csv(physical_mapping, file.path(dirs$tables, "Figure2_climate_zone_specific_physical_axis_crosswalk.csv"), row.names = FALSE)
write.csv(influence_observations, file.path(dirs$diagnostics, "Figure2_observation_level_leverage_and_influence.csv"), row.names = FALSE)
write.csv(influence_summary, file.path(dirs$diagnostics, "Figure2_influence_summary_by_phenotype.csv"), row.names = FALSE)
write.csv(sensitivity_curve_qa, file.path(dirs$diagnostics, "Figure2_influence_excluded_curve_comparison.csv"), row.names = FALSE)
write.csv(physical_mapping_qa, file.path(dirs$diagnostics, "Figure2_physical_axis_mapping_QA.csv"), row.names = FALSE)
write.csv(partial_qa, file.path(dirs$diagnostics, "Figure2_partial_residual_QA.csv"), row.names = FALSE)
write.csv(curve_qa, file.path(dirs$diagnostics, "Figure2_frozen_curve_preservation_QA.csv"), row.names = FALSE)

wb <- createWorkbook()
workbook_data <- list(
  Curves = curves,
  Panel_support = support,
  Partial_residual_bins = partial_bins,
  Partial_residual_raw = partial_raw,
  Histogram = histograms,
  Rug_observations = rug_data,
  Strict_v3_thresholds = strict_thresholds,
  Displayed_all_turning_points = display_thresholds,
  All_turning_points = bootstrap_turning_points,
  Turning_point_panel_audit = turning_point_panel_audit,
  Wet_days_117d_grid_audit = wet_days_gmax_grid_audit,
  CRPa_temperature_EDF_audit = crpa_temperature_edf_curve_audit,
  Main_y_axis_5tick_QA = main_y_axis_five_tick_qa,
  Balanced_physical_x_axis_QA = main_x_axis_balanced_qa,
  Turning_point_label_2dp_QA = turning_point_label_precision_qa,
  Supplementary_curve_features = supplementary_curve_features,
  Full_feature_catalog = supplementary_full_feature_catalog,
  Physical_axis_crosswalk = physical_mapping,
  Influence_observations = influence_observations,
  Influence_summary = influence_summary,
  Sensitivity_curves = sensitivity_curves,
  Sensitivity_comparison = sensitivity_curve_qa,
  Curve_preservation_QA = curve_qa,
  Physical_mapping_QA = physical_mapping_qa,
  Partial_residual_QA = partial_qa,
  Model_frame_QA = model_data_qa
)
for (nm in names(workbook_data)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, workbook_data[[nm]])
  freezePane(wb, nm, firstRow = TRUE)
  addFilter(wb, nm, rows = 1, cols = seq_len(ncol(workbook_data[[nm]])))
}
saveWorkbook(
  wb,
  file.path(dirs$source, "SourceData_Figure2_threshold_ORCI_and_supplementary_features_v9_two_decimal_turning_point_labels.xlsx"),
  overwrite = TRUE
)

contract_text <- c(
  "# Figure contract",
  "",
  "Core conclusion: the frozen fully adjusted Model C curves are presented on intuitive phenotype-specific physical-equivalent axes together with direct evidence on data support, adjusted observations and influence sensitivity.",
  "",
  "Figure archetype: quantitative grid (6 AMR phenotypes x 4 climatic exposures) plus the unchanged smooth-term summary panel.",
  "",
  "Backend: R only (mgcv, ggplot2 and cowplot).",
  "",
  "Primary display rules:",
  "- observed phenotype-specific min-max on the exact model-input scale, labelled using the phenotype-specific temperate-zone physical-unit crosswalk;",
  "- fitted curves remain restricted to the observed panel-specific min-max range, while display-only physical-axis limits use round ticks and equal terminal gaps (or natural zero-bounded endpoints);",
  "- the shared phenotype/sample-size/OR label occupies a dedicated row-label column, leaving all four climate-response viewports with equal physical dimensions;",
  "- P2.5-P97.5 solid coloured curve and band;",
  "- observed tails in pale grey with dashed curves;",
  "- 16-bin histogram plus observation rug;",
  "- adjusted partial-residual bin means are omitted from the uncluttered primary export and retained in a reviewer-evidence companion export and Source Data;",
  "- red rug marks for a fixed observation-level influence-screening rule used only in sensitivity analysis;",
  "- influence-excluded sensitivity curves aligned at the panel median exposure and reported in a separate Supplementary Figure;",
  "- every objectively detected internal extremum is shown by the same red long-dashed line and labelled with its term-centred OR and pointwise 95% interval;",
  "- turning-point locations are displayed to two decimal places in physical-equivalent units, while full-precision values are retained in Source Data;",
  "- filled versus open red markers indicate whether the original 900/1,000 bootstrap-stability criterion was met, without printing bootstrap-support values or imposing another cutoff;",
  "- Chg and Stb descriptors are generated by a phenotype-agnostic supplementary feature algorithm; a T prefix records pointwise P>=0.10 but is not a geometric category; Rap/Tre aliases are discontinued; none is overlaid on the main figure;",
  "- every exposure-response panel has exactly five labelled y-axis ticks; EDF<0.10 is labelled as shrunk towards zero, 0.10<=EDF<=1.05 is labelled approximately linear only when no derivative sign change is present, and a low-complexity term with a sign change is labelled as a low-complexity curved fit;",
  "- term-centred OR and the existing conditional pointwise 95% interval retained exactly;",
  "- sensitivity refits are diagnostic only and do not replace the primary models or feed Results 2 or 3.",
  "",
  "Reviewer risk controlled: tail estimates are visually distinguished; sample support is explicit; adjusted observations remain available in a dedicated evidence version; model-parameter uncertainty is not mislabelled as a prediction interval."
)
writeLines(contract_text, file.path(dirs$contract, "Figure2_physical_axes_thresholds_influence_contract.md"))

legend_text <- paste0(
  "Fig. 2 | Climate–AMR exposure–response associations and data support. ",
  "a, Term-centred odds-ratio (OR) curves from the frozen fully adjusted Model C generalized additive models for six AMR phenotypes and four lagged climatic factors. Fitted curves retain the observed minimum-to-maximum range of each exact model input; display-only limits extend slightly where needed to provide round physical-unit ticks and balanced terminal axis segments. Each panel uses a phenotype-specific temperate-zone physical-equivalent axis; crosswalks for all climate zones are supplied in Source Data. Solid coloured curves and pointwise model-based 95% intervals show P2.5–P97.5; pale grey dashed curves and bands show the limited-support tails. Lower subpanels show 16-bin histograms and country-year rugs; red rug marks identify observations meeting the fixed influence-screening rule. Red long-dashed lines denote internal turning points identified among terms with EDF>1.05 by first-derivative sign changes within P2.5–P97.5. Labels give the turning-point type, physical-equivalent location and term-centred pointwise OR (95% interval); this interval is not a confidence interval for turning-point location. Filled and open red markers indicate whether the fixed 900-of-1,000 bootstrap-stability criterion was or was not met, respectively. All objectively detected candidates are displayed; exact bootstrap support, conditional location intervals and local observation and country counts are reported in Source Data. Panels with EDF<=1.05 use a common OR scale of 0.5–1.5. Terms with EDF<0.10 are labelled as shrunk towards zero; other low-complexity terms are labelled approximately linear unless a derivative sign change remains visible. Shallow geometric candidates below the fixed complexity gate, including the CR-Pa temperature candidate, are not declared as turning points and are reported in Source Data. ",
  "b, F statistics and effective degrees of freedom (EDF) for the four climatic smooths; symbol size is proportional to EDF. n denotes country-year observations. P values are approximate overall smooth-term tests returned by mgcv; *P<0.05, **P<0.01 and ***P<0.001. Source data are provided as a Source Data file."
)
writeLines(legend_text, file.path(dirs$report, "Figure2_proposed_legend.txt"))

reviewer_evidence_legend <- paste0(
  legend_text,
  " In the reviewer-evidence companion, open circles are geometric means of adjusted partial residuals in ten equal-frequency exposure bins; for observation i and exposure term j, the partial residual was defined on the model's Gaussian identity scale as the fitted contribution of term j plus the full-model response residual. These descriptive summaries are not independent observations and their dispersion should not be compared with the shaded model-parameter confidence bands, which are not outcome prediction intervals."
)
writeLines(
  reviewer_evidence_legend,
  file.path(dirs$report, "Figure2_reviewer_evidence_companion_legend.txt")
)

supplementary_influence_legend <- paste0(
  "Supplementary Fig. X | Observation-level leverage and influence diagnostics for the six frozen Model C fits. ",
  "Each point represents one country-year model observation. The vertical dotted line is twice the mean penalized-smoother leverage (2p_eff/n), and horizontal dashed lines mark internally studentized response residuals of -3 and 3. ",
  "Point size is proportional to approximate penalized Cook distance. Red points meet the fixed influence-exclusion rule used for sensitivity analysis, defined as |studentized residual| > 3 or Cook distance > 4/n. The three observations with the largest Cook distances in each phenotype are labelled. High leverage or influence is not interpreted as a data error; the flags define a reproducible sensitivity set."
)
writeLines(
  supplementary_influence_legend,
  file.path(dirs$report, "Supplementary_Figure_influence_diagnostics_legend.txt")
)

supplementary_sensitivity_legend <- paste0(
  "Supplementary Fig. Y | Influence-excluded sensitivity of climate-AMR exposure-response shapes. ",
  "Coloured curves show the frozen Model C smooths and black dot-dashed curves show models refitted after excluding observations with |studentized residual| > 3 or approximate penalized Cook distance > 4/n. ",
  "Both curves are expressed relative to the panel-specific median exposure and aligned at OR = 1 at that point, so differences represent changes in curve shape rather than changes in the mgcv term-centering constraint. ",
  "Grey regions lie outside the phenotype-exposure-specific P2.5-P97.5 model-input range. Text reports the maximum absolute log-OR difference between curves within P2.5-P97.5. Physical axes are phenotype-specific temperate-zone equivalents; all climate-zone crosswalks are supplied in Source Data."
)
writeLines(
  supplementary_sensitivity_legend,
  file.path(dirs$report, "Supplementary_Figure_influence_excluded_sensitivity_legend.txt")
)

methods_text <- c(
  "Physical-axis, adjusted partial-residual and influence-sensitivity display for Figure 2",
  "",
  "The primary models were not refitted. Because climatic exposures entered Model C after standardization within climate zone, a unique pooled global conversion from the model-input scale to physical units does not exist. For display only, each phenotype-specific x axis was labelled using the exact intercept and slope linking the standardized and lagged physical exposure within the temperate-zone observations of that analytical sample; crosswalks for polar, temperate and tropical observations were retained in Source Data. For each frozen Model C fit and climatic smooth j, observation-level adjusted partial residuals were calculated on the Gaussian identity (logit-resistance) scale as r_ij = f_j(x_ij) + e_i, where f_j(x_ij) is the fitted term-centred smooth contribution and e_i is the response residual from the complete multivariable model. Observations were divided into ten equal-frequency bins using the empirical distribution of the exact standardized model input. Open circles show exp(mean(r_ij)) at the mean exposure in each bin. Observation-level leverage was calculated as diag[X Vp X' W]/sigma2 using the model lpmatrix and conditional penalized coefficient covariance. Internally studentized residuals and approximate penalized Cook distances were then calculated, and a fixed sensitivity rule excluded observations with |studentized residual| > 3 or Cook distance > 4/n. Six diagnostic sensitivity models were fitted with the original formula, Gaussian identity family, REML and select=TRUE; these fits were not used in the primary curve, threshold detection or downstream Results 2–3. Sensitivity curves were aligned to the frozen curve at the panel median exposure to compare shape. Primary curves retain the original term-centred OR and conditional pointwise interval calculations."
)
writeLines(methods_text, file.path(dirs$report, "Figure2_adjusted_partial_residual_methods.txt"))

result1_wording <- c(
  paste0(
    "Across the 24 phenotype–climate exposure–response functions, the objective derivative-based procedure identified 24 internal turning-point estimates in 16 panels within the phenotype-specific P2.5–P97.5 exposure ranges. Seven estimates met the fixed 900-of-1,000 parametric-bootstrap stability criterion. For transparent reporting, this criterion was applied as a labelling rather than exclusion rule: lower-support estimates were retained in Fig. 2 to represent fitted-curve geometry but were not interpreted as bootstrap-supported turning points. Labels report term-centred pointwise ORs and model-based 95% intervals, whereas exact bootstrap support and conditional location intervals are reported separately in Supplementary Table X. Panels with EDF<=1.05 were displayed on a common OR scale. Terms with EDF<0.10 were described as shrunk towards zero; terms with 0.10<=EDF<=1.05 were described as approximately linear only in the absence of an internal derivative sign change. The CR-Pa temperature term (EDF=0.996) contained a shallow derivative sign change and was therefore described as a low-complexity curved fit, but it was not assigned a nonlinear turning point under the fixed EDF gate."
  ),
  "",
  paste0(
    "Parametric-bootstrap support and influence-exclusion sensitivity were interpreted as complementary diagnostics. For example, the CR-Ab temperature GMax was reproduced in 90.1% of parametric-bootstrap refits, whereas excluding observations meeting the influence rule changed the fitted curve by a maximum absolute delta log(OR) of 0.548 within P2.5–P97.5. The former quantifies sampling stability conditional on the observed data configuration; the latter assesses sensitivity to individual influential observations. We therefore report both diagnostics and do not use bootstrap support alone as evidence that a turning point is insensitive to outliers."
  )
)
writeLines(
  result1_wording,
  file.path(dirs$report, "Recommended_Result1_turning_point_reporting.txt")
)

terminology_ledger <- c(
  "# Threshold and curve-feature terminology ledger",
  "",
  "- estimated turning point: an internal derivative sign change in the observed fitted smooth; not a causal or clinical threshold",
  "- bootstrap-stability criterion: reproduction in >=900 of 1,000 attempted full parametric-bootstrap refits; not a P value",
  "- bootstrap support (BS): conditional stability under model-based sampling noise around the fitted values; it does not assess sensitivity to the observed data configuration or to individual influential observations",
  "- term-centred pointwise OR: exp{f_j(x)} under the mgcv smooth-centering constraint; not an OR relative to P50 or another physical reference exposure",
  "- pointwise 95% interval: exp{f_j(x) +/- 1.96 SE[f_j(x)]}; does not quantify turning-point location uncertainty",
  "- Chg: supplementary rapid-change point defined by a high absolute first derivative",
  "- Stb: supplementary low-derivative interval, not a point",
  "- T prefix: descriptive pointwise P>=0.10 flag; never used for feature detection or retention",
  "- Rap/Tre: discontinued legacy aliases; they are not executed or reported as geometric classes in the audited analysis",
  "- low-complexity curved fit: a term with 0.10<=EDF<=1.05 that nevertheless contains a derivative sign change; the visual descriptor prevents EDF from being misread as proof of an exact straight line"
)
writeLines(
  terminology_ledger,
  file.path(dirs$report, "Threshold_and_curve_feature_terminology_ledger.md")
)

audit_text <- c(
  "# Figure 2 physical-axis and influence-sensitivity audit",
  "",
  "## Frozen inferential objects",
  "",
  "- Six existing Model C objects were read without refitting.",
  "- Model frames were reconstructed and checked variable-by-variable against the saved model frames.",
  "- All 24 curves and confidence limits were reproduced exactly against the validated v3 source data.",
  "- The term-centred OR definition and conditional pointwise 95% interval calculation were unchanged.",
  "",
  "## New display-only evidence",
  "",
  "- Each panel uses its own observed model-input min-max and P2.5-P97.5.",
  "- Estimates outside P2.5-P97.5 are retained but encoded as limited-support tails.",
  "- Histograms and rugs expose the exact country-year exposure support.",
  "- Adjusted partial-residual bin means expose how the observations underlying the fitted smooth align with the term curve after adjustment for the other model components.",
  "- Physical labels use an exact phenotype-specific temperate-zone crosswalk; no global mean/SD inverse was used.",
  "- Observation-level leverage, studentized residual and approximate penalized Cook distance are supplied for every model row.",
  "- Six influence-excluded refits are sensitivity models only; they do not replace the frozen primary models.",
  "",
  "## Threshold overlays",
  "",
  "The recommended export overlays every objectively detected internal extremum using the same red long-dashed line and reports its term-centred pointwise OR (95% interval). Filled red markers indicate extrema meeting the original 900/1,000 bootstrap-stability criterion and open red markers indicate that the criterion was not met; the numerical support is not printed in the main panel and every candidate remains visible. Pointwise OR significance is not a detection or retention criterion. Exact stability, conditional location intervals and local observation and country support remain in Source Data. No legacy hard-coded target was reintroduced.",
  "A phenotype-agnostic supplementary feature catalogue implements rapid-change points (Chg) as one maximum of |f'(x)| per contiguous region above its within-curve 90th percentile and stability regions (Stb) as contiguous intervals below the 10th percentile spanning at least 5% of the supported exposure range. A T prefix records pointwise P>=0.10 but does not define or retain a feature. Rap/Tre aliases are discontinued. Chg and Stb were not assigned post hoc bootstrap support because the archived refits implemented point matching for extrema only, whereas Chg is defined by relative derivative quantiles and Stb is an interval requiring a separately specified overlap rule.",
  "",
  "## Interpretation boundary",
  "",
  "The uncluttered primary export omits adjusted partial-residual circles. The reviewer-evidence companion retains them as descriptive adjusted response summaries; they are not validation confidence intervals and are not expected to lie inside the shaded term-wise confidence band. Parametric-bootstrap support is conditional on the original observed data configuration and is not a test of resistance to influential observations; the separately reported influence-exclusion refits address that distinct question. Formal claims about clustered uncertainty still require the separate country-cluster sensitivity analysis."
)
writeLines(audit_text, file.path(dirs$report, "Figure2_physical_axes_thresholds_influence_scientific_audit.md"))

threshold_axis_audit <- c(
  "# Figure 2 axis and threshold audit",
  "",
  "## Phenotype-exposure prediction domains",
  "",
  "All 24 panels use an independently constructed prediction grid. For each phenotype-exposure combination, the curve spans the minimum-to-maximum value of that exact lagged, within-climate-zone standardized variable in the corresponding frozen Model C model frame; P2.5 and P97.5 are also calculated separately for that combination. Similar-looking axes do not imply that a common grid was used. They mainly reflect overlapping country-year climate coverage across phenotypes and the use of phenotype-specific temperate-zone physical-equivalent labels.",
  "",
  "The 3GCR-Ec and CR-Ab wet-days GMax estimates round to 117 d, but their standardized search grids, grid steps and candidate z values differ. Their unrounded temperate-zone physical-equivalent locations differ by approximately 0.22 d. The equality at the displayed precision is therefore a rounding coincidence rather than shared-grid reuse; the full-precision audit is supplied in `05_diagnostics`.",
  "",
  "Because the model standardized climate exposures separately within polar, temperate and tropical zones, a standardized value has three different physical interpretations. A unique pooled raw-unit x axis therefore does not exist. The main figure uses the exact temperate-zone conversion as a labelled reference and supplies all three zone-specific crosswalks in Source Data.",
  "",
  "## Legacy threshold engine versus strict v3",
  "",
  "The legacy code is not an executable implementation of the manuscript S2 description. It predicts on fixed physical ranges shared across phenotypes, back-transforms climate-zone-standardized inputs with pooled mean and SD values, applies phenotype/exposure-specific overrides to linearity and LOESS span, uses many hard-coded target locations and priority bonuses, and does not execute the stated 1,000-replicate threshold bootstrap.",
  "",
  "Strict v3 is a deliberate replacement, not a cosmetic redraw. It searches only the combination-specific P2.5-P97.5 domain on the actual model-input scale, evaluates the fitted GAM smooth directly, defines internal extrema by first-derivative sign changes, removes phenotype-specific priorities and OR-significance filters, and uses 1,000 full Gaussian response simulations followed by identical REML bam refits. A turning point is retained only when matched in at least 900 of 1,000 attempted refits.",
  "",
  "The current manuscript S2 text must therefore be replaced by the packaged method. Linear-versus-nonlinear R-squared screening, adaptive LOESS, phenotype-specific target values and the maximum-of-five priority rule are not retained. Chg and Stb descriptors are executed only by the separately specified, phenotype-agnostic supplementary feature algorithm and are not used to select main-figure turning points. A T prefix records pointwise P>=0.10; Rap/Tre aliases are discontinued.",
  "",
  "## Main-figure display decision",
  "",
  "The recommended manuscript export shows all internal extrema detected under the reported objective algorithm. Each label gives the term-centred OR and pointwise 95% interval. All extrema use the same red long-dashed line; marker fill alone indicates whether the original >=90% stability criterion was met, and exact support is retained in Source Data. No candidate is removed for low support. Adjusted partial-residual bin means are removed from the primary figure to reduce clutter but retained in a reviewer-evidence companion, the Source Data workbook and the observation-level diagnostic outputs.",
  "",
  "The CR-Pa temperature curve and its EDF=0.996 are extracted from the same s(TMP_scaled_lag) term. EDF measures effective model complexity and does not impose an exact algebraic line. Because the fitted term contains a shallow derivative sign change, the panel is labelled `low-complexity curved fit` rather than `approximately linear`; the fixed EDF<=1.05 gate still prevents it from being declared a nonlinear turning point."
)
writeLines(
  threshold_axis_audit,
  file.path(dirs$report, "Figure2_axis_and_threshold_method_audit.md")
)

replacement_s2_path <- file.path(
  dirs$report, "Recommended_replacement_Supplementary_Method_S2.txt"
)
writeLines(
  c(
    "S2. Detection and stability assessment of turning points in climate–AMR associations",
    "",
    "We characterized internal turning points in the fitted climate–AMR exposure–response functions using a fixed, phenotype-agnostic procedure based on curve geometry. For each phenotype–exposure combination, the term-centred climatic smooth from Model C was evaluated directly on a regular grid of 1,001 points spanning the empirical 2.5th to 97.5th percentiles of the exact standardized model input. This produced 1,000 equal intervals, each corresponding to 0.1% of the searched exposure range. Fitted curves retained the observed minimum-to-maximum range for visualization; display-only physical-axis limits could extend slightly to provide round tick values and balanced terminal segments, but no model values were extrapolated into these margins. Values outside the empirical 2.5th to 97.5th percentile interval were not searched for turning points. Turning-point searches were restricted to smooth terms with effective degrees of freedom (EDF) >1.05, using a fixed operational margin above EDF=1 to distinguish nonlinear structure from approximately linear fits. For display, terms with EDF<0.10 were classified as shrunk towards zero; terms with 0.10<=EDF<=1.05 were classified as approximately linear when no internal derivative sign change was present and as low-complexity curved fits otherwise.",
    "",
    "First derivatives of the fitted smooths were calculated by finite differences on the linear-predictor scale. Derivative values with an absolute magnitude <=10^-5 of the maximum absolute derivative, subject to a minimum numerical tolerance of 10^-10, were treated as zero. A positive-to-negative derivative sign change defined a local maximum, whereas a negative-to-positive change defined a local minimum; zero-crossing locations were obtained by linear interpolation between adjacent grid points. A local extremum located within two grid steps of the highest or lowest fitted value in the searched interval was classified as a global maximum (GMax) or global minimum (GMin), respectively; other extrema were classified as Max or Min. Curves with a total log-odds amplitude <10^-4 and extrema with a local log-odds prominence <10^-4 were treated as numerically negligible. The same criteria were applied to all phenotypes, climatic exposures and bootstrap refits. Candidate detection was independent of pointwise P values and OR magnitude.",
    "",
    "The stability of each observed-data turning point was assessed using 1,000 parametric-bootstrap refits. For each phenotype, bootstrap responses were generated as y* = fitted(y) + epsilon, where epsilon was sampled independently from N(0, sigma^2) using the estimated residual variance. Each simulated dataset was fitted using the same Model C formula, Gaussian identity family, restricted maximum likelihood estimation and select=TRUE penalties, after which turning points were redetected using the procedure described above. Bootstrap extrema were matched one-to-one to observed-data extrema of the same direction within a location tolerance equal to 10% of the corresponding P2.5–P97.5 exposure range. Matching maximized the number of paired extrema and then minimized their total location difference. Bootstrap support was calculated as the number of matched refits divided by 1,000, and support >=90% met the fixed bootstrap-stability criterion. Conditional 95% location intervals were calculated from the 2.5th and 97.5th percentiles of matched locations. All detected observed-data turning points were reported, with marker fill indicating whether the 90% stability criterion was met; exact support values and location intervals are provided in Source Data.",
    "",
    "At each turning-point location x_t, the term-centred contribution f_j(x_t) and its conditional standard error were extracted from Model C. Pointwise ORs and 95% confidence intervals were calculated as exp{f_j(x_t)} and exp{f_j(x_t) +/- 1.96 SE[f_j(x_t)]}, respectively. These intervals quantify uncertainty in the fitted smooth contribution at the estimated location; turning-point location uncertainty was summarized separately by the conditional bootstrap location intervals. Pointwise ORs and confidence intervals were not used to detect or retain turning points.",
    "",
    "Rapid-change and stability features were characterized separately for supplementary reporting using the same fitted curves within P2.5–P97.5. A rapid-change point (Chg) was defined as the location of the maximum absolute first derivative within each contiguous region above the curve-specific 90th percentile of |f'(x)|. A stability region (Stb) was defined as a contiguous interval below the 10th percentile of |f'(x)| with a width of at least 5% of the searched exposure range. The midpoint of each Stb interval was used to report a representative pointwise OR. A T prefix was added when the representative pointwise P value was >=0.10; this descriptor did not affect feature detection or retention. Chg, Stb and T-prefixed features were reported in the Supplementary Information and were not used to select main-figure turning points. Bootstrap support was evaluated for extrema only. All detection and matching were conducted on the standardized model-input scale."
  ),
  replacement_s2_path
)

main_methods_threshold_text <- c(
  "Construction of GAMMs and identification of climate–AMR turning points",
  "",
  "On the basis of the selected lag structures, phenotype-specific generalized additive mixed models were fitted to the logit-transformed resistance rates. Climatic variables and partial least-squares components were represented using cubic regression splines (bs=\"cr\"). Geographic variation was modelled using a spherical spline for latitude and longitude (bs=\"sos\"), temporal variation using a cubic regression spline for year, and regional variation using a random-effect smooth (bs=\"re\"); climate zone was included as a categorical covariate.",
  "",
  "The basis dimension was set to k=5 for temperature, k=10 for precipitation, relative humidity, wet days and each partial least-squares component, k=20 for the spatial smooth and k=8 for the temporal smooth. Models were fitted using bam with a Gaussian identity family and restricted maximum likelihood estimation. Selection penalties (select=TRUE) were applied to permit weakly supported smooth terms to shrink towards zero.",
  "",
  "Climate–AMR exposure–response functions were subsequently evaluated for internal turning points using a fixed, phenotype-agnostic first-derivative procedure. Analyses were restricted to the phenotype–exposure-specific 2.5th to 97.5th percentile range of the standardized model input to limit inference in sparsely observed tails. Positive-to-negative and negative-to-positive derivative sign changes defined local maxima and minima, respectively. Candidate stability was assessed using 1,000 parametric-bootstrap refits, with extrema matched by direction using a location tolerance defined relative to the searched exposure range. Candidates reproduced in at least 90% of refits were designated bootstrap-supported turning points, and all detected candidates were reported. Operational definitions of nonlinear complexity, numerical differentiation, candidate matching and classification are provided in Supplementary Method S2."
)
writeLines(
  main_methods_threshold_text,
  file.path(dirs$report, "Recommended_replacement_main_Methods_GAMM_and_turning_points.txt")
)

s2_revision_rationale <- c(
  "# Supplementary Method S2 revision and code-alignment record",
  "",
  "## Overall conclusion",
  "",
  "The replacement S2 describes the executed audited algorithm. It supersedes the original three-stage text and must not be combined with the legacy R-squared screen, adaptive LOESS, phenotype-specific target priorities or maximum-of-five selection rule, none of which defines the current output.",
  "",
  "## Five wording decisions",
  "",
  "1. `Pre-specified` was replaced by `fixed, phenotype-agnostic` because the EDF<=1.05 complexity gate was introduced during revision rather than documented in the original analysis. The defensible property is uniform application across all phenotype-exposure combinations, not temporal pre-specification.",
  "2. The original manuscript described >=900/1,000 as a retention rule. The revised analysis applies the same numerical criterion as a labelling rule and displays all objectively detected candidates. This reporting change is stated explicitly and does not change the fitted models, candidate geometry or bootstrap calculations.",
  "3. The general algorithm permits three reasons for an empty panel, but all eight empty panels in the present result arose from the EDF gate. This observed outcome is stated explicitly.",
  "4. CR-Pa temperature has EDF=0.996 and a shallow diagnostic geometric minimum at approximately 6.69 °C when the EDF gate is ignored. The candidate is retained in Source Data but is not declared as a turning point under the fixed algorithm; the figure legend now states this boundary.",
  "5. Chg and Stb were not assigned post hoc bootstrap support. This is not because such resampling is impossible in principle. The archived refits implemented point matching for extrema only; Chg is defined by a relative derivative-quantile region and Stb is an interval, so a distinct fixed point-matching or interval-overlap rule would be required. No such rule was introduced retrospectively.",
  "",
  "## Terminology",
  "",
  "- turning point: an internal fitted-curve extremum meeting the fixed complexity and numerical criteria; not a causal, biological or clinical threshold",
  "- bootstrap support: conditional reproduction frequency under model-based response simulation and full refitting; not a P value",
  "- Chg: descriptive rapid-change feature based on a within-curve derivative quantile",
  "- Stb: descriptive low-derivative interval",
  "- T prefix: pointwise P>=0.10 descriptor only; not a geometric category or retention rule"
)
writeLines(
  s2_revision_rationale,
  file.path(dirs$report, "Supplementary_Method_S2_revision_and_code_alignment.md")
)
response_letter_disclosure <- c(
  "Recommended disclosure for the response letter or revision history",
  "",
  "In the original Supplementary Method S2, the >=900/1,000 bootstrap criterion was described as a retention rule. In the revised analysis, we retained the same numerical criterion but use it as a labelling rule: all turning-point candidates detected by the fixed phenotype-agnostic geometry algorithm are displayed, while marker fill indicates whether the criterion was met. We made this reporting change to avoid selectively suppressing model-derived candidates and now state it explicitly. The change does not affect the fitted Model C curves, observed-data candidate locations, the 1,000 bootstrap refits or the calculated bootstrap-support values."
)
writeLines(
  response_letter_disclosure,
  file.path(dirs$report, "S2_retention_to_labelling_change_for_response_letter.txt")
)
readme_text <- c(
  "# Figure 2 equal-panel, balanced physical-axis redraw",
  "",
  "## Recommended primary file",
  "",
  "`02_figures/Figure2_main_turning_points_ORCI_without_partial_residual_bins.pdf`",
  "",
  "This version is recommended for the manuscript main figure. It shows every objectively detected internal turning point using a uniform red long-dashed line and reports the physical-equivalent location to two decimal places together with the term-centred pointwise OR (95% interval). Filled versus open red markers indicate whether the original 900/1,000 stability criterion was met; numerical bootstrap support is retained in Source Data and no secondary display cutoff is imposed. Every curve panel has exactly five y-axis ticks. Panels with EDF<=1.05 use a common 0.5–1.5 OR scale, with separate annotations for terms shrunk towards zero, approximately linear terms and low-complexity curved fits.",
  "",
  "The fitted x range remains the observed panel-specific minimum-to-maximum range. Display-only limits use round physical-unit ticks and either equal terminal gaps or natural zero-bounded endpoints; no curve is extrapolated. The full 24-panel axis audit is saved in `05_diagnostics/Figure2_balanced_physical_x_axis_QA.csv`.",
  "",
  "The phenotype/sample-size/OR label is placed in a dedicated row-label column rather than inside the Temperature plot. Consequently, all four climate-response viewports have the same physical width and height within every phenotype row.",
  "",
  "## Companion clean file",
  "",
  "`02_figures/Figure2_turning_points_termcentred_OR_95CI_with_support_symbols.pdf`",
  "",
  "This reviewer-evidence companion retains turning-point labels and adjusted partial-residual bin means. The supplementary feature catalogue implements Chg and Stb using phenotype-agnostic rules; a T prefix records pointwise P>=0.10, while Rap/Tre aliases are discontinued.",
  "",
  "## Zero-cascade guarantee",
  "",
  "The six frozen primary Model C objects were not refitted, and all 24 primary curves and confidence limits reproduce v3 to numerical precision. Six separately identified influence-excluded models are diagnostic sensitivity fits only and do not feed threshold detection or Results 2–3.",
  "",
  "The influence-diagnostic and influence-excluded curve-comparison figures are supplied separately in `02_figures` so that substantial sensitivity in selected panels is visible without compressing or obscuring the primary Figure 2 curves.",
  "",
  "See `04_tables/Supplementary_full_threshold_and_curve_feature_catalog.csv` for the complete extrema/Chg/Stb/T catalogue, `05_diagnostics` for numerical QA (including the CR-Pa EDF/curve term-consistency check and the two 117-d wet-days grid audit), and `06_report` for the proposed legend, methods text and scientific audit."
)
writeLines(readme_text, file.path(output_root, "README.md"))

writeLines(capture.output(sessionInfo()), file.path(dirs$report, "R_sessionInfo.txt"))

script_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(file_arg) == 1L && file.exists(file_arg)) {
  file.copy(file_arg, file.path(dirs$code, basename(file_arg)), overwrite = TRUE)
}

input_files <- c(
  list.files(dirs$input_models, full.names = TRUE),
  list.files(dirs$input_data, full.names = TRUE),
  file.path(dirs$inputs, "bacteria_four_factors_lag_summary.csv"),
  file.path(dirs$inputs, "Figure2_v3_main_figure_strict_thresholds_max2.csv"),
  file.path(dirs$inputs, "Figure2_v3_bootstrap_support_for_observed_turning_points.csv")
)
input_fingerprint <- tibble(
  file = basename(input_files),
  relative_path = sub(paste0("^", output_root, "/"), "", input_files),
  md5 = unname(tools::md5sum(input_files)),
  size_bytes = file.info(input_files)$size
)
write.csv(input_fingerprint, file.path(dirs$diagnostics, "Figure2_frozen_input_fingerprints.csv"), row.names = FALSE)

manifest_paths <- setdiff(
  list.files(output_root, recursive = TRUE, all.files = FALSE),
  c("MANIFEST.csv", "07_logs/Figure2_threshold_ORCI_supplementary_features_v9_two_decimal_turning_point_labels_run.log")
)
manifest <- tibble(
  relative_path = manifest_paths,
  size_bytes = file.info(file.path(output_root, manifest_paths))$size,
  md5 = unname(tools::md5sum(file.path(output_root, manifest_paths)))
)
write.csv(manifest, file.path(output_root, "MANIFEST.csv"), row.names = FALSE)

cat("\nCompleted:", format(Sys.time()), "\n")
cat("Recommended threshold-preserving Figure 2:", paste0(stem_combined_thresholds_no_partial_bins, ".pdf"), "\n")
cat("Reviewer-evidence companion Figure 2:", paste0(stem_combined_thresholds, ".pdf"), "\n")
cat("Output bundle:", output_root, "\n")
