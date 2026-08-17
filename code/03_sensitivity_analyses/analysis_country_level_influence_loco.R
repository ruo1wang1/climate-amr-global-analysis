#!/usr/bin/env Rscript

# Leave-one-country-out refitting of the six Model C GAMMs.

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(readr)
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

analysis_root <- Sys.getenv(
  "MODEL_C_LOCO_INPUT_ROOT",
  unset = file.path(repo_root, "outputs", "historical_associations", "figure1_current")
)
output_root <- Sys.getenv(
  "MODEL_C_LOCO_REFIT_OUTPUT_ROOT",
  unset = file.path(
    repo_root, "outputs", "sensitivity_analyses",
    "country_level_influence_loco", "validated_outputs"
  )
)
n_cores <- suppressWarnings(as.integer(Sys.getenv("MODEL_C_LOCO_CORES", unset = "3")))
if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
if (.Platform$OS.type == "windows") n_cores <- 1L

dirs <- list(
  contract = file.path(output_root, "00_contract"),
  code = file.path(output_root, "01_code"),
  figures = file.path(output_root, "02_figures"),
  source = file.path(output_root, "03_source_data"),
  tables = file.path(output_root, "04_tables"),
  diagnostics = file.path(output_root, "05_diagnostics"),
  report = file.path(output_root, "06_report"),
  logs = file.path(output_root, "07_logs"),
  checkpoints = file.path(output_root, "08_checkpoints")
)
walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(dirs$logs, "ModelC_LOCO_sensitivity_v1_run.log")
if (file.exists(log_file)) file.remove(log_file)
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  # Rscript encodes spaces in --file paths as ~+~ on some macOS builds.
  script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
  if (file.exists(script_path)) {
    script_path <- normalizePath(script_path, mustWork = TRUE)
    file.copy(script_path, file.path(dirs$code, basename(script_path)), overwrite = TRUE)
  }
}

log_line("Starting Model C country-deletion analysis with ", n_cores, " worker(s).")
source(file.path(
  repo_root, "code", "01_historical_associations", "figure1_turning_point_core.R"
))
v3_load_lags_cache <- v3_load_lags(
  file.path(analysis_root, "00_inputs", "bacteria_four_factors_lag_summary.csv")
)

phenotype_order <- v3_specs$phenotype
climate_order <- v3_factors$climate
climate_colours <- setNames(v3_factors$colour, v3_factors$climate)

safe_quantile <- function(x, p) {
  as.numeric(quantile(x[is.finite(x)], p, type = 7, names = FALSE, na.rm = TRUE))
}

safe_curve_cor <- function(x, y, method = "spearman", tolerance = 1e-10) {
  sx <- sd(x)
  sy <- sd(y)
  if (is.finite(sx) && is.finite(sy) && sx <= tolerance && sy <= tolerance) return(1)
  if (!is.finite(sx) || !is.finite(sy) || sx <= tolerance || sy <= tolerance) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}

bind_nonempty <- function(x) {
  x <- keep(x, ~ !is.null(.x) && is.data.frame(.x) && nrow(.x) > 0L)
  if (length(x) == 0L) tibble() else bind_rows(x)
}

fit_model_c <- function(model, data) {
  warnings <- character()
  ctrl <- model$control
  ctrl$nthreads <- 1L
  started <- proc.time()[[3]]
  fit <- tryCatch(
    withCallingHandlers(
      bam(
        formula(model), data = droplevels(data), family = gaussian(),
        method = "REML", select = TRUE, discrete = FALSE, use.chol = FALSE,
        control = ctrl
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[[3]] - started
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL, elapsed_seconds = elapsed, converged = FALSE,
      warning_text = paste(unique(warnings), collapse = " | "),
      error_text = conditionMessage(fit)
    ))
  }
  converged <- isTRUE(fit$converged) && all(is.finite(coef(fit)))
  list(
    fit = fit, elapsed_seconds = elapsed, converged = converged,
    warning_text = paste(unique(warnings), collapse = " | "),
    error_text = if (converged) "" else "Non-converged fit or non-finite coefficient"
  )
}

models <- setNames(map(seq_len(nrow(v3_specs)), function(i) {
  readRDS(file.path(
    analysis_root, "00_inputs", "frozen_models", v3_specs$model_file[i]
  ))
}), phenotype_order)

prepared <- setNames(map(seq_len(nrow(v3_specs)), function(i) {
  v3_prepare_data(
    file.path(
      analysis_root, "00_inputs", "model_ready_inputs", v3_specs$input_file[i]
    ),
    v3_specs$phenotype[i],
    v3_load_lags_cache
  )
}), phenotype_order)

data_validation <- map_dfr(phenotype_order, function(phenotype) {
  v3_validate_model_data(models[[phenotype]], prepared[[phenotype]], phenotype)
})
write_csv(data_validation, file.path(dirs$diagnostics, "ModelC_LOCO_input_alignment_QA.csv"))
if (!all(data_validation$passed)) {
  stop("Prepared-data alignment with frozen Model C failed.", call. = FALSE)
}

archived_curves <- read_csv(
  file.path(analysis_root, "03_source_data", "Figure1_termcentred_curve_source_data.csv"),
  show_col_types = FALSE
)
full_candidates <- read_csv(
  file.path(
    analysis_root, "04_tables",
    "Figure1_all_observed_turning_points_with_1000_refit_support.csv"
  ),
  show_col_types = FALSE
)

# Refit the complete analytical sample before country deletion. Exact
# reproduction verifies that the refitting function preserves Model C.
full_refit_qa <- map_dfr(phenotype_order, function(phenotype) {
  original <- models[[phenotype]]
  result <- fit_model_c(original, original$model)
  if (is.null(result$fit)) {
    return(tibble(
      phenotype = phenotype, n_model = nrow(original$model),
      converged = FALSE, max_abs_coefficient_difference = Inf,
      max_abs_fitted_difference = Inf, max_abs_edf_difference = Inf,
      elapsed_seconds = result$elapsed_seconds,
      warning_text = result$warning_text, error_text = result$error_text
    ))
  }
  refit <- result$fit
  tibble(
    phenotype = phenotype,
    n_model = nrow(original$model),
    converged = result$converged,
    max_abs_coefficient_difference = max(abs(coef(refit) - coef(original))),
    max_abs_fitted_difference = max(abs(fitted(refit) - fitted(original))),
    max_abs_edf_difference = max(abs(refit$edf - original$edf)),
    elapsed_seconds = result$elapsed_seconds,
    warning_text = result$warning_text,
    error_text = result$error_text
  )
})
write_csv(
  full_refit_qa,
  file.path(dirs$diagnostics, "ModelC_LOCO_full_sample_refit_reproducibility.csv")
)
if (any(!full_refit_qa$converged) ||
    any(full_refit_qa$max_abs_coefficient_difference > 1e-10) ||
    any(full_refit_qa$max_abs_fitted_difference > 1e-10)) {
  stop("Full-sample refit did not reproduce the frozen Model C fits.", call. = FALSE)
}

log_line("All six full-sample refits reproduced the frozen models exactly.")

# Fixed full-sample panel definitions. LOCO fits are always compared on the
# same full-sample support and grid, so differences are not caused by moving
# percentile limits after deletion.
panel_definitions <- list()
panel_index <- 1L
for (phenotype in phenotype_order) {
  model <- models[[phenotype]]
  for (j in seq_len(nrow(v3_factors))) {
    factor_row <- v3_factors[j, ]
    variable <- factor_row$variable
    climate_value <- factor_row$climate
    x_obs <- model$model[[variable]]
    q <- as.numeric(quantile(x_obs, c(0.025, 0.975), type = 7, names = FALSE))
    q10_90 <- as.numeric(quantile(x_obs, c(0.10, 0.90), type = 7, names = FALSE))
    q25_75 <- as.numeric(quantile(x_obs, c(0.25, 0.75), type = 7, names = FALSE))
    curve_grid <- seq(q[1], q[2], length.out = 401L)
    search_grid <- seq(q[1], q[2], length.out = 1001L)
    full_curve <- v3_term_curve(model, variable, curve_grid, with_se = TRUE)
    full_search <- v3_term_curve(model, variable, search_grid, with_se = FALSE)
    full_stat <- v3_smooth_statistics(model, variable)
    reference_x <- median(x_obs, na.rm = TRUE)
    reference_eta <- approx(full_curve$x_scaled, full_curve$eta, xout = reference_x)$y
    map_row <- archived_curves %>%
      filter(
        .data$phenotype == .env$phenotype,
        .data$climate == .env$climate_value
      ) %>%
      slice(1)
    if (nrow(map_row) != 1L) stop("Missing physical-axis mapping.", call. = FALSE)
    panel_key <- paste(phenotype, factor_row$climate, sep = "__")
    panel_definitions[[panel_key]] <- list(
      phenotype = phenotype,
      climate = factor_row$climate,
      variable = variable,
      raw = factor_row$raw,
      colour = factor_row$colour,
      curve_grid = curve_grid,
      search_grid = search_grid,
      p2_5 = q[1],
      p97_5 = q[2],
      p10 = q10_90[1],
      p90 = q10_90[2],
      p25 = q25_75[1],
      p75 = q25_75[2],
      support_width = diff(q),
      reference_x = reference_x,
      reference_eta = reference_eta,
      map_intercept = map_row$map_intercept,
      map_slope = map_row$map_slope,
      reference_zone = map_row$reference_zone,
      full_curve = full_curve,
      full_search = full_search,
      full_stat = full_stat,
      full_candidates = full_candidates %>%
        filter(
          .data$phenotype == .env$phenotype,
          .data$climate == .env$climate_value
        )
    )
    panel_index <- panel_index + 1L
  }
}

extract_model_fit_summary <- function(fit) {
  s <- tryCatch(summary(fit), error = function(e) NULL)
  tibble(
    model_r_squared = if (is.null(s)) NA_real_ else as.numeric(s$r.sq),
    deviance_explained = if (is.null(s)) NA_real_ else as.numeric(s$dev.expl)
  )
}

run_one_country <- function(phenotype, country) {
  model <- models[[phenotype]]
  dat <- prepared[[phenotype]]
  keep_rows <- as.character(dat$NAME) != country
  n_removed <- sum(!keep_rows)
  result <- fit_model_c(model, model$model[keep_rows, , drop = FALSE])
  status <- tibble(
    phenotype = phenotype,
    country_left_out = country,
    n_full = nrow(model$model),
    n_removed = n_removed,
    fraction_removed = n_removed / nrow(model$model),
    n_remaining = sum(keep_rows),
    converged = result$converged,
    elapsed_seconds = result$elapsed_seconds,
    warning_text = result$warning_text,
    error_text = result$error_text
  )
  if (is.null(result$fit) || !result$converged) {
    return(list(
      status = status, model_fit = tibble(), curves = tibble(),
      metrics = tibble(), detected = tibble(), matches = tibble()
    ))
  }

  fit <- result$fit
  model_fit <- bind_cols(
    status %>% select(phenotype, country_left_out),
    extract_model_fit_summary(fit)
  )
  curve_rows <- list()
  metric_rows <- list()
  detected_rows <- list()
  match_rows <- list()

  for (climate in climate_order) {
    def <- panel_definitions[[paste(phenotype, climate, sep = "__")]]
    variable <- def$variable
    loco_curve <- v3_term_curve(fit, variable, def$curve_grid, with_se = FALSE)
    loco_reference_eta <- approx(
      loco_curve$x_scaled, loco_curve$eta, xout = def$reference_x
    )$y
    aligned_eta <- loco_curve$eta - loco_reference_eta + def$reference_eta
    full_eta <- def$full_curve$eta
    delta_raw <- loco_curve$eta - full_eta
    delta_aligned <- aligned_eta - full_eta
    central_10_90 <- def$curve_grid >= def$p10 & def$curve_grid <= def$p90
    full_p90_p10 <- approx(def$curve_grid, full_eta, xout = def$p90)$y -
      approx(def$curve_grid, full_eta, xout = def$p10)$y
    loco_p90_p10 <- approx(def$curve_grid, aligned_eta, xout = def$p90)$y -
      approx(def$curve_grid, aligned_eta, xout = def$p10)$y
    full_iqr <- approx(def$curve_grid, full_eta, xout = def$p75)$y -
      approx(def$curve_grid, full_eta, xout = def$p25)$y
    loco_iqr <- approx(def$curve_grid, aligned_eta, xout = def$p75)$y -
      approx(def$curve_grid, aligned_eta, xout = def$p25)$y
    loco_stat <- v3_smooth_statistics(fit, variable)
    full_sig <- is.finite(def$full_stat$p_value) && def$full_stat$p_value < 0.05
    loco_sig <- is.finite(loco_stat$p_value) && loco_stat$p_value < 0.05

    curve_rows[[climate]] <- tibble(
      phenotype = phenotype,
      climate = climate,
      variable = variable,
      country_left_out = country,
      n_removed = n_removed,
      grid_id = seq_along(def$curve_grid),
      x_scaled = def$curve_grid,
      x_physical_equivalent = def$map_intercept + def$map_slope * def$curve_grid,
      full_eta = full_eta,
      full_OR = exp(full_eta),
      loco_eta_termcentred = loco_curve$eta,
      loco_OR_termcentred = exp(loco_curve$eta),
      loco_eta_aligned_at_full_P50 = aligned_eta,
      loco_OR_aligned_at_full_P50 = exp(aligned_eta),
      delta_eta_termcentred = delta_raw,
      delta_eta_shape_aligned = delta_aligned,
      alignment_reference_x_scaled = def$reference_x,
      alignment_reference_OR_full = exp(def$reference_eta)
    )

    metric_rows[[climate]] <- tibble(
      phenotype = phenotype,
      climate = climate,
      variable = variable,
      country_left_out = country,
      n_removed = n_removed,
      fraction_removed = n_removed / nrow(model$model),
      full_edf = def$full_stat$edf,
      loco_edf = loco_stat$edf,
      delta_edf = loco_stat$edf - def$full_stat$edf,
      full_smooth_p_value = def$full_stat$p_value,
      loco_smooth_p_value = loco_stat$p_value,
      full_p_lt_0_05 = full_sig,
      loco_p_lt_0_05 = loco_sig,
      p_classification_preserved = identical(full_sig, loco_sig),
      full_edf_gt_1_05 = is.finite(def$full_stat$edf) && def$full_stat$edf > 1.05,
      loco_edf_gt_1_05 = is.finite(loco_stat$edf) && loco_stat$edf > 1.05,
      edf_gate_classification_preserved =
        (is.finite(def$full_stat$edf) && def$full_stat$edf > 1.05) ==
        (is.finite(loco_stat$edf) && loco_stat$edf > 1.05),
      max_abs_delta_logOR_termcentred_P2_5_P97_5 = max(abs(delta_raw)),
      median_abs_delta_logOR_termcentred_P2_5_P97_5 = median(abs(delta_raw)),
      max_abs_delta_logOR_shape_aligned_P2_5_P97_5 = max(abs(delta_aligned)),
      median_abs_delta_logOR_shape_aligned_P2_5_P97_5 = median(abs(delta_aligned)),
      max_abs_delta_logOR_shape_aligned_P10_P90 =
        max(abs(delta_aligned[central_10_90])),
      pearson_curve_correlation_shape_aligned = safe_curve_cor(
        full_eta, aligned_eta, method = "pearson"
      ),
      spearman_curve_correlation_shape_aligned = safe_curve_cor(
        full_eta, aligned_eta, method = "spearman"
      ),
      full_logOR_P90_minus_P10 = full_p90_p10,
      loco_logOR_P90_minus_P10 = loco_p90_p10,
      delta_logOR_P90_minus_P10 = loco_p90_p10 - full_p90_p10,
      P90_P10_direction_preserved =
        sign(loco_p90_p10) == sign(full_p90_p10) || abs(full_p90_p10) < 1e-8,
      full_logOR_P75_minus_P25 = full_iqr,
      loco_logOR_P75_minus_P25 = loco_iqr,
      delta_logOR_P75_minus_P25 = loco_iqr - full_iqr,
      P75_P25_direction_preserved =
        sign(loco_iqr) == sign(full_iqr) || abs(full_iqr) < 1e-8
    )

    search_curve <- v3_term_curve(fit, variable, def$search_grid, with_se = FALSE)
    detected <- v3_detect_turning_points(search_curve, loco_stat$edf)
    if (nrow(detected) > 0L) {
      detected_rows[[climate]] <- detected %>%
        mutate(
          phenotype = phenotype,
          climate = climate,
          variable = variable,
          country_left_out = country,
          n_removed = n_removed,
          x_physical_equivalent = def$map_intercept + def$map_slope * x_scaled
        )
    }
    if (nrow(def$full_candidates) > 0L) {
      matched <- v3_match_one_replicate(
        def$full_candidates, detected, def$support_width, tolerance_fraction = 0.10
      ) %>%
        rename(
          loco_x_scaled = bootstrap_x_scaled,
          loco_OR_termcentred = bootstrap_OR,
          loco_prominence_eta = bootstrap_prominence_eta
        ) %>%
        mutate(
          phenotype = phenotype,
          climate = climate,
          variable = variable,
          country_left_out = country,
          n_removed = n_removed,
          tolerance_scaled = 0.10 * def$support_width,
          loco_x_physical_equivalent =
            def$map_intercept + def$map_slope * loco_x_scaled
        )
      match_rows[[climate]] <- matched
    }
  }

  list(
    status = status,
    model_fit = model_fit,
    curves = bind_nonempty(curve_rows),
    metrics = bind_nonempty(metric_rows),
    detected = bind_nonempty(detected_rows),
    matches = bind_nonempty(match_rows)
  )
}

safe_country_job <- function(phenotype, country) {
  tryCatch(
    run_one_country(phenotype, country),
    error = function(e) {
      model <- models[[phenotype]]
      dat <- prepared[[phenotype]]
      n_removed <- sum(as.character(dat$NAME) == country)
      list(
        status = tibble(
          phenotype = phenotype, country_left_out = country,
          n_full = nrow(model$model), n_removed = n_removed,
          fraction_removed = n_removed / nrow(model$model),
          n_remaining = nrow(model$model) - n_removed,
          converged = FALSE, elapsed_seconds = NA_real_, warning_text = "",
          error_text = paste0("Unhandled job error: ", conditionMessage(e))
        ),
        model_fit = tibble(), curves = tibble(), metrics = tibble(),
        detected = tibble(), matches = tibble()
      )
    }
  )
}

all_results <- list()
for (phenotype in phenotype_order) {
  countries <- sort(unique(as.character(prepared[[phenotype]]$NAME)))
  checkpoint <- file.path(
    dirs$checkpoints,
    paste0("ModelC_LOCO_", v3_safe_stem(phenotype), "_checkpoint.rds")
  )
  if (file.exists(checkpoint)) {
    results <- readRDS(checkpoint)
    if (!setequal(names(results), countries)) {
      stop("Checkpoint country set does not match the analytical sample for ", phenotype, call. = FALSE)
    }
    results <- results[countries]
    log_line("Loaded completed checkpoint for ", phenotype, ": ", length(countries), " countries.")
  } else {
    log_line("LOCO refits for ", phenotype, ": ", length(countries), " countries.")
    results <- parallel::mclapply(
      countries,
      function(country) safe_country_job(phenotype, country),
      mc.cores = n_cores,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
    names(results) <- countries
    saveRDS(results, checkpoint, compress = "xz")
  }
  all_results[[phenotype]] <- results
  n_ok <- sum(map_lgl(results, ~ isTRUE(.x$status$converged[1])))
  log_line("Completed ", phenotype, ": ", n_ok, "/", length(countries), " converged.")
}

flat_results <- unlist(all_results, recursive = FALSE)
refit_status <- bind_rows(map(flat_results, "status"))
model_fit_results <- bind_nonempty(map(flat_results, "model_fit"))
loco_curves <- bind_nonempty(map(flat_results, "curves"))
loco_metrics <- bind_nonempty(map(flat_results, "metrics"))
loco_detected_turning_points <- bind_nonempty(map(flat_results, "detected"))

# Checkpoints contain model-scale curve and extremum results. Rebuild physical
# mappings and candidate matches here so all joins are explicit and immune to
# data-mask name collisions.
panel_mapping <- map_dfr(panel_definitions, function(def) {
  tibble(
    phenotype = def$phenotype,
    climate = def$climate,
    map_intercept = def$map_intercept,
    map_slope = def$map_slope,
    support_width = def$support_width
  )
})
loco_curves <- loco_curves %>%
  select(-x_physical_equivalent) %>%
  left_join(panel_mapping, by = c("phenotype", "climate")) %>%
  mutate(x_physical_equivalent = map_intercept + map_slope * x_scaled) %>%
  select(-map_intercept, -map_slope, -support_width)

if (nrow(loco_detected_turning_points) > 0L) {
  loco_detected_turning_points <- loco_detected_turning_points %>%
    select(-x_physical_equivalent) %>%
    left_join(panel_mapping, by = c("phenotype", "climate")) %>%
    mutate(x_physical_equivalent = map_intercept + map_slope * x_scaled) %>%
    select(-map_intercept, -map_slope, -support_width)
}

loco_turning_point_matches <- map_dfr(phenotype_order, function(phenotype_value) {
  successful_countries <- refit_status %>%
    filter(.data$phenotype == .env$phenotype_value, converged) %>%
    pull(country_left_out)
  map_dfr(climate_order, function(climate_value) {
    def <- panel_definitions[[paste(phenotype_value, climate_value, sep = "__")]]
    observed <- def$full_candidates
    if (nrow(observed) == 0L) return(tibble())
    map_dfr(successful_countries, function(country_value) {
      detected <- loco_detected_turning_points %>%
        filter(
          .data$phenotype == .env$phenotype_value,
          .data$climate == .env$climate_value,
          .data$country_left_out == .env$country_value
        )
      v3_match_one_replicate(
        observed, detected, def$support_width, tolerance_fraction = 0.10
      ) %>%
        rename(
          loco_x_scaled = bootstrap_x_scaled,
          loco_OR_termcentred = bootstrap_OR,
          loco_prominence_eta = bootstrap_prominence_eta
        ) %>%
        mutate(
          phenotype = phenotype_value,
          climate = climate_value,
          variable = def$variable,
          country_left_out = country_value,
          n_removed = sum(as.character(prepared[[phenotype_value]]$NAME) == country_value),
          tolerance_scaled = 0.10 * def$support_width,
          loco_x_physical_equivalent =
            def$map_intercept + def$map_slope * loco_x_scaled
        )
    })
  })
})

write_csv(refit_status, file.path(dirs$tables, "ModelC_LOCO_refit_status.csv"))
write_csv(model_fit_results, file.path(dirs$tables, "ModelC_LOCO_model_fit_statistics.csv"))
write_csv(loco_metrics, file.path(dirs$tables, "ModelC_LOCO_panel_curve_metrics.csv"))
write_csv(
  loco_detected_turning_points,
  file.path(dirs$source, "ModelC_LOCO_detected_turning_points.csv")
)
write_csv(
  loco_turning_point_matches,
  file.path(dirs$source, "ModelC_LOCO_full_sample_turning_point_matches.csv")
)
write_csv(
  loco_curves,
  file.path(dirs$source, "ModelC_LOCO_all_aligned_and_termcentred_curves.csv.gz")
)

successful_counts <- refit_status %>%
  group_by(phenotype) %>%
  summarise(
    n_attempted_refits = n(),
    n_successful_refits = sum(converged),
    n_failed_refits = sum(!converged),
    refit_success_fraction = mean(converged),
    .groups = "drop"
  )

panel_summary <- loco_metrics %>%
  group_by(phenotype, climate, variable) %>%
  summarise(
    n_successful_refits = n(),
    max_country_fraction_removed = max(fraction_removed),
    worst_case_max_abs_delta_logOR_shape_P2_5_P97_5 =
      max(max_abs_delta_logOR_shape_aligned_P2_5_P97_5),
    worst_case_country_shape_P2_5_P97_5 = country_left_out[
      which.max(max_abs_delta_logOR_shape_aligned_P2_5_P97_5)
    ],
    median_country_max_abs_delta_logOR_shape_P2_5_P97_5 =
      median(max_abs_delta_logOR_shape_aligned_P2_5_P97_5),
    worst_case_max_abs_delta_logOR_shape_P10_P90 =
      max(max_abs_delta_logOR_shape_aligned_P10_P90),
    worst_case_country_shape_P10_P90 = country_left_out[
      which.max(max_abs_delta_logOR_shape_aligned_P10_P90)
    ],
    minimum_pearson_curve_correlation = safe_min(
      pearson_curve_correlation_shape_aligned
    ),
    median_pearson_curve_correlation = safe_median(
      pearson_curve_correlation_shape_aligned
    ),
    minimum_spearman_curve_correlation = safe_min(
      spearman_curve_correlation_shape_aligned
    ),
    proportion_smooth_p_lt_0_05 = mean(loco_p_lt_0_05),
    proportion_p_classification_preserved = mean(p_classification_preserved),
    proportion_edf_gate_classification_preserved =
      mean(edf_gate_classification_preserved),
    proportion_P90_P10_direction_preserved = mean(P90_P10_direction_preserved),
    proportion_P75_P25_direction_preserved = mean(P75_P25_direction_preserved),
    full_edf = first(full_edf),
    full_smooth_p_value = first(full_smooth_p_value),
    full_logOR_P90_minus_P10 = first(full_logOR_P90_minus_P10),
    full_logOR_P75_minus_P25 = first(full_logOR_P75_minus_P25),
    .groups = "drop"
  ) %>%
  left_join(successful_counts, by = c("phenotype", "n_successful_refits"))

country_influence_ranking <- loco_metrics %>%
  group_by(phenotype, climate) %>%
  arrange(desc(max_abs_delta_logOR_shape_aligned_P2_5_P97_5), .by_group = TRUE) %>%
  mutate(influence_rank = row_number()) %>%
  ungroup()

if (nrow(loco_turning_point_matches) > 0L) {
  turning_point_stability <- loco_turning_point_matches %>%
    group_by(phenotype, climate, variable, observed_candidate_id) %>%
    summarise(
      n_successful_LOCO_refits = n(),
      n_matched_LOCO_refits = sum(matched),
      LOCO_match_fraction = mean(matched),
      location_median_scaled_LOCO = if (any(matched)) {
        median(loco_x_scaled[matched], na.rm = TRUE)
      } else NA_real_,
      location_lower_2_5pct_scaled_LOCO = if (any(matched)) {
        safe_quantile(loco_x_scaled[matched], 0.025)
      } else NA_real_,
      location_upper_97_5pct_scaled_LOCO = if (any(matched)) {
        safe_quantile(loco_x_scaled[matched], 0.975)
      } else NA_real_,
      location_min_scaled_LOCO = if (any(matched)) min(loco_x_scaled[matched]) else NA_real_,
      location_max_scaled_LOCO = if (any(matched)) max(loco_x_scaled[matched]) else NA_real_,
      maximum_match_distance_scaled = if (any(matched)) {
        max(match_distance_scaled[matched])
      } else NA_real_,
      .groups = "drop"
    ) %>%
    left_join(
      full_candidates %>%
        select(
          phenotype, climate, variable, observed_candidate_id, type, direction,
          x_scaled, x_physical_equivalent, observed_OR, observed_Lower_95CI,
          observed_Upper_95CI, support_attempted, bootstrap_retained_strict,
          map_intercept, map_slope, unit_display
        ),
      by = c("phenotype", "climate", "variable", "observed_candidate_id")
    ) %>%
    mutate(
      location_median_physical_LOCO =
        map_intercept + map_slope * location_median_scaled_LOCO,
      location_lower_2_5pct_physical_LOCO =
        map_intercept + map_slope * location_lower_2_5pct_scaled_LOCO,
      location_upper_97_5pct_physical_LOCO =
        map_intercept + map_slope * location_upper_97_5pct_scaled_LOCO,
      location_min_physical_LOCO = map_intercept + map_slope * location_min_scaled_LOCO,
      location_max_physical_LOCO = map_intercept + map_slope * location_max_scaled_LOCO
    )
} else {
  turning_point_stability <- tibble()
}

model_fit_summary <- model_fit_results %>%
  group_by(phenotype) %>%
  summarise(
    n_successful_refits = n(),
    LOCO_r_squared_min = min(model_r_squared, na.rm = TRUE),
    LOCO_r_squared_median = median(model_r_squared, na.rm = TRUE),
    LOCO_r_squared_max = max(model_r_squared, na.rm = TRUE),
    LOCO_deviance_explained_min = min(deviance_explained, na.rm = TRUE),
    LOCO_deviance_explained_median = median(deviance_explained, na.rm = TRUE),
    LOCO_deviance_explained_max = max(deviance_explained, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    map_dfr(phenotype_order, function(phenotype) {
      s <- summary(models[[phenotype]])
      tibble(
        phenotype = phenotype,
        full_r_squared = as.numeric(s$r.sq),
        full_deviance_explained = as.numeric(s$dev.expl)
      )
    }),
    by = "phenotype"
  )

write_csv(panel_summary, file.path(dirs$tables, "ModelC_LOCO_panel_summary.csv"))
write_csv(
  country_influence_ranking,
  file.path(dirs$tables, "ModelC_LOCO_country_influence_ranking.csv")
)
write_csv(
  turning_point_stability,
  file.path(dirs$tables, "ModelC_LOCO_turning_point_stability.csv")
)
write_csv(model_fit_summary, file.path(dirs$tables, "ModelC_LOCO_model_fit_summary.csv"))

full_curve_source <- map_dfr(panel_definitions, function(def) {
  def$full_curve %>%
    mutate(
      phenotype = def$phenotype,
      climate = def$climate,
      variable = def$variable,
      grid_id = row_number(),
      x_physical_equivalent = def$map_intercept + def$map_slope * x_scaled,
      reference_zone = def$reference_zone,
      p2_5_scaled = def$p2_5,
      p97_5_scaled = def$p97_5,
      alignment_reference_x_scaled = def$reference_x,
      alignment_reference_OR_full = exp(def$reference_eta)
    )
})

loco_envelopes <- loco_curves %>%
  group_by(phenotype, climate, variable, grid_id, x_scaled, x_physical_equivalent) %>%
  summarise(
    n_successful_refits = n(),
    loco_OR_aligned_min = min(loco_OR_aligned_at_full_P50),
    loco_OR_aligned_q2_5 = safe_quantile(loco_OR_aligned_at_full_P50, 0.025),
    loco_OR_aligned_median = median(loco_OR_aligned_at_full_P50),
    loco_OR_aligned_q97_5 = safe_quantile(loco_OR_aligned_at_full_P50, 0.975),
    loco_OR_aligned_max = max(loco_OR_aligned_at_full_P50),
    loco_OR_termcentred_q2_5 = safe_quantile(loco_OR_termcentred, 0.025),
    loco_OR_termcentred_median = median(loco_OR_termcentred),
    loco_OR_termcentred_q97_5 = safe_quantile(loco_OR_termcentred, 0.975),
    .groups = "drop"
  ) %>%
  left_join(
    full_curve_source %>%
      select(
        phenotype, climate, variable, grid_id, full_OR = OR,
        full_Lower_95CI = Lower_95CI, full_Upper_95CI = Upper_95CI,
        alignment_reference_x_scaled, alignment_reference_OR_full
      ),
    by = c("phenotype", "climate", "variable", "grid_id")
  )

write_csv(
  loco_envelopes,
  file.path(dirs$source, "ModelC_LOCO_country_deletion_curve_envelopes.csv")
)
write_csv(
  full_curve_source,
  file.path(dirs$source, "ModelC_LOCO_full_sample_reference_curves.csv")
)

theme_loco <- theme_classic(base_size = 8.3, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 9, hjust = 0),
    plot.subtitle = element_text(size = 7.4, colour = "grey25"),
    plot.caption = element_text(size = 7, colour = "grey25", hjust = 0),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    plot.margin = margin(5, 6, 5, 6)
  )

axis_label <- function(climate) {
  switch(
    climate,
    Temperature = "Temperature (°C)",
    Humidity = "Relative humidity (%)",
    Precipitation = "Precipitation (mm)",
    WetDays = "Wet days (d)"
  )
}

make_envelope_panel <- function(phenotype, climate) {
  z <- loco_envelopes %>%
    filter(
      .data$phenotype == .env$phenotype,
      .data$climate == .env$climate
    ) %>%
    arrange(x_scaled)
  colour <- climate_colours[[climate]]
  ggplot(z, aes(x = x_physical_equivalent)) +
    geom_ribbon(
      aes(ymin = loco_OR_aligned_q2_5, ymax = loco_OR_aligned_q97_5),
      fill = "grey68", alpha = 0.48, linewidth = 0
    ) +
    geom_ribbon(
      aes(ymin = full_Lower_95CI, ymax = full_Upper_95CI),
      fill = colour, alpha = 0.17, linewidth = 0
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey35", linewidth = 0.38) +
    geom_line(aes(y = loco_OR_aligned_median), colour = "grey38", linewidth = 0.45) +
    geom_line(aes(y = full_OR), colour = colour, linewidth = 0.82) +
    scale_x_continuous(breaks = breaks_pretty(n = 4)) +
    scale_y_continuous(breaks = breaks_pretty(n = 5), expand = expansion(mult = c(0.04, 0.08))) +
    labs(
      title = paste0(phenotype, " — ", climate),
      x = axis_label(climate), y = "Odds ratio"
    ) +
    theme_loco
}

envelope_panels <- map(phenotype_order, function(phenotype) {
  map(climate_order, ~ make_envelope_panel(phenotype, .x))
}) %>% flatten()

envelope_figure <- wrap_plots(envelope_panels, ncol = 4) +
  plot_annotation(
    title = "Country-deletion stability of Model C climate–AMR exposure–response functions",
    subtitle = str_wrap(paste0(
      "Coloured lines and ribbons are the frozen full-sample term-centred curves and pointwise 95% intervals. ",
      "Grey lines and bands are the median and central 95% range across leave-one-country-out refits, aligned to the full curve at the full-sample median exposure solely for shape comparison."
    ), width = 175),
    caption = paste0(
      "Curves are restricted to the phenotype–exposure-specific P2.5–P97.5 range. ",
      "The grey country-deletion envelope is descriptive and is not a confidence interval. ",
      "Physical-unit axes use the same phenotype-specific temperate-zone crosswalk as Figure 1."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, family = "Arial"),
      plot.subtitle = element_text(size = 9, family = "Arial", colour = "grey20"),
      plot.caption = element_text(size = 8, family = "Arial", colour = "grey25", hjust = 0)
    )
  )

save_figure_bundle <- function(plot, stem, width, height) {
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
  ragg::agg_png(
    file.path(dirs$figures, paste0(stem, "_preview.png")),
    width = width, height = height, units = "in", res = 180, background = "white"
  )
  print(plot)
  dev.off()
  ragg::agg_tiff(
    file.path(dirs$figures, paste0(stem, "_600dpi.tiff")),
    width = width, height = height, units = "in", res = 600,
    compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()
}

save_figure_bundle(
  envelope_figure,
  "Supplementary_Figure_ModelC_LOCO_curve_envelopes",
  width = 14.2, height = 19.2
)

heat_data <- panel_summary %>%
  mutate(
    phenotype = factor(phenotype, levels = rev(phenotype_order)),
    climate = factor(climate, levels = climate_order)
  )

heat_delta <- ggplot(
  heat_data,
  aes(x = climate, y = phenotype, fill = worst_case_max_abs_delta_logOR_shape_P2_5_P97_5)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", worst_case_max_abs_delta_logOR_shape_P2_5_P97_5)),
    size = 2.8
  ) +
  scale_fill_gradient(low = "#F7FBFF", high = "#B2182B", name = "Max |Δ log OR|") +
  labs(title = "a  Worst-case curve-shape change", x = NULL, y = NULL) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold"), panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

heat_central <- ggplot(
  heat_data,
  aes(x = climate, y = phenotype, fill = worst_case_max_abs_delta_logOR_shape_P10_P90)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", worst_case_max_abs_delta_logOR_shape_P10_P90)),
    size = 2.8
  ) +
  scale_fill_gradient(low = "#F7FBFF", high = "#B2182B", name = "Max |Δ log OR|") +
  labs(title = "b  Worst-case change within P10–P90", x = NULL, y = NULL) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold"), panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

heat_p <- ggplot(
  heat_data,
  aes(x = climate, y = phenotype, fill = proportion_smooth_p_lt_0_05)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = percent(proportion_smooth_p_lt_0_05, accuracy = 1)), size = 2.8) +
  scale_fill_gradient(low = "#F7F7F7", high = "#1B7837", limits = c(0, 1), name = "LOCO P<0.05") +
  labs(title = "c  Overall smooth-term P-value across refits", x = NULL, y = NULL) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold"), panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

summary_figure <- (heat_delta | heat_central | heat_p) +
  plot_annotation(
    title = "Summary of leave-one-country-out influence diagnostics",
    caption = paste0(
      "Curve differences use P50-aligned log-OR functions over the fixed full-sample support. Panel a uses P2.5–P97.5 and panel b uses P10–P90. ",
      "P values are the approximate overall smooth-term tests returned by mgcv after each country deletion."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, family = "Arial"),
      plot.caption = element_text(size = 8, family = "Arial", hjust = 0)
    )
  )
save_figure_bundle(
  summary_figure,
  "Supplementary_Figure_ModelC_LOCO_summary_heatmaps",
  width = 15.5, height = 6.2
)

if (nrow(turning_point_stability) > 0L) {
  turning_plot_data <- turning_point_stability %>%
    mutate(
      label_base = paste0(
        phenotype, " | ", climate, " | ", type, " | ",
        case_when(
          climate == "Temperature" ~ sprintf("%.1f °C", x_physical_equivalent),
          climate == "Humidity" ~ sprintf("%.1f%%", x_physical_equivalent),
          climate == "Precipitation" ~ sprintf("%.0f mm", x_physical_equivalent),
          TRUE ~ sprintf("%.0f d", x_physical_equivalent)
        )
      ),
      label = if_else(
        duplicated(label_base) | duplicated(label_base, fromLast = TRUE),
        paste0(label_base, " [", observed_candidate_id, "]"),
        label_base
      ),
      label = factor(label, levels = rev(label))
    )
  turning_plot <- ggplot(
    turning_plot_data,
    aes(x = LOCO_match_fraction, y = label, colour = climate)
  ) +
    geom_segment(aes(x = 0, xend = LOCO_match_fraction, yend = label), colour = "grey82") +
    geom_point(aes(shape = bootstrap_retained_strict), size = 2.7, stroke = 0.7) +
    geom_vline(xintercept = 0.90, linetype = "dashed", colour = "grey40") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.02)) +
    scale_colour_manual(values = climate_colours, name = "Climatic exposure") +
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      name = "Full-sample bootstrap support ≥90%"
    ) +
    labs(
      title = "Full-sample turning points reproduced after individual country deletion",
      subtitle = "Matching required the same extremum direction within 10% of the fixed P2.5–P97.5 search range.",
      x = "Proportion of successful LOCO refits with a matched turning point", y = NULL,
      caption = "This deletion fraction is a deterministic influence diagnostic, not a bootstrap probability or P value."
    ) +
    theme_classic(base_size = 9, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 7.2),
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = 8)
    )
  save_figure_bundle(
    turning_plot,
    "Supplementary_Figure_ModelC_LOCO_turning_point_retention",
    width = 11.5, height = 9.7
  )
}

# Targeted CR-Ab temperature diagnostic requested by the reviewer context.
crab_rank <- country_influence_ranking %>%
  filter(phenotype == "CR-Ab", climate == "Temperature") %>%
  arrange(influence_rank)
top_crab_countries <- crab_rank %>% slice_head(n = 5) %>% pull(country_left_out)
crab_curves <- loco_curves %>%
  filter(phenotype == "CR-Ab", climate == "Temperature")
crab_envelope <- loco_envelopes %>%
  filter(phenotype == "CR-Ab", climate == "Temperature")
highlight_colours <- setNames(
  c("#0072B2", "#009E73", "#E69F00", "#6A3D9A", "#5C4033"),
  top_crab_countries
)
crab_turning_point <- turning_point_stability %>%
  filter(phenotype == "CR-Ab", climate == "Temperature", type == "GMax") %>%
  slice(1)

crab_curve_plot <- ggplot() +
  geom_ribbon(
    data = crab_envelope,
    aes(
      x = x_physical_equivalent,
      ymin = loco_OR_aligned_q2_5, ymax = loco_OR_aligned_q97_5
    ),
    fill = "grey72", alpha = 0.5
  ) +
  geom_line(
    data = crab_curves %>% filter(!country_left_out %in% top_crab_countries),
    aes(x = x_physical_equivalent, y = loco_OR_aligned_at_full_P50, group = country_left_out),
    colour = "grey68", alpha = 0.30, linewidth = 0.28
  ) +
  geom_line(
    data = crab_curves %>% filter(country_left_out %in% top_crab_countries),
    aes(
      x = x_physical_equivalent, y = loco_OR_aligned_at_full_P50,
      colour = country_left_out, group = country_left_out
    ),
    linewidth = 0.72
  ) +
  geom_line(
    data = crab_envelope,
    aes(x = x_physical_equivalent, y = full_OR),
    colour = climate_colours[["Temperature"]], linewidth = 1.05
  ) +
  geom_vline(
    xintercept = crab_turning_point$x_physical_equivalent,
    linetype = "dotdash", colour = "grey28", linewidth = 0.55
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey35") +
  scale_colour_manual(values = highlight_colours, name = "Country omitted") +
  scale_x_continuous(breaks = breaks_pretty(n = 5)) +
  scale_y_continuous(breaks = breaks_pretty(n = 5)) +
  labs(
    title = "a  CR-Ab temperature curve under country deletion",
    subtitle = paste0(
      "Full-sample GMax ", sprintf("%.1f °C", crab_turning_point$x_physical_equivalent),
      "; matched in ", crab_turning_point$n_matched_LOCO_refits, "/",
      crab_turning_point$n_successful_LOCO_refits, " refits (",
      percent(crab_turning_point$LOCO_match_fraction, accuracy = 0.1), ")"
    ),
    x = "Temperature (°C)", y = "Odds ratio"
  ) +
  theme_loco +
  theme(legend.position = "bottom")

crab_bar_data <- crab_rank %>%
  slice_head(n = 15) %>%
  mutate(country_left_out = factor(
    country_left_out,
    levels = rev(country_left_out)
  ))
crab_rank_plot <- ggplot(
  crab_bar_data,
  aes(
    x = max_abs_delta_logOR_shape_aligned_P2_5_P97_5,
    y = country_left_out
  )
) +
  geom_col(fill = climate_colours[["Temperature"]], width = 0.72) +
  geom_text(
    aes(
      label = paste0(sprintf("%.2f", max_abs_delta_logOR_shape_aligned_P2_5_P97_5),
                     "  (n=", n_removed, ")")
    ),
    hjust = -0.05, size = 2.8
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28)), breaks = breaks_pretty(n = 5)) +
  labs(
    title = "b  Most influential country deletions",
    x = "Maximum absolute change in fitted log OR", y = NULL
  ) +
  theme_loco

crab_figure <- (crab_curve_plot | crab_rank_plot) +
  plot_annotation(
    title = "Targeted country-deletion assessment of the CR-Ab temperature association",
    caption = paste0(
      "All LOCO curves are aligned to the frozen full-sample curve at the full-sample median temperature solely for shape comparison. ",
      "n is the number of country-year observations removed. The red curve is the frozen Model C estimate."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, family = "Arial"),
      plot.caption = element_text(size = 8, family = "Arial", hjust = 0)
    )
  )
save_figure_bundle(
  crab_figure,
  "Supplementary_Figure_ModelC_LOCO_CRAb_temperature_targeted",
  width = 13.0, height = 6.8
)

# Compact workbook for manuscript and review use. The full curve matrix remains
# in compressed CSV because it exceeds a practical spreadsheet size.
workbook_path <- file.path(dirs$source, "Source_Data_ModelC_LOCO_sensitivity.xlsx")
wb <- createWorkbook()
sheet_data <- list(
  Refit_status = refit_status,
  Panel_summary = panel_summary,
  Country_ranking = country_influence_ranking,
  Turning_point_stability = turning_point_stability,
  Model_fit_summary = model_fit_summary,
  Curve_envelopes = loco_envelopes,
  Full_reference_curves = full_curve_source,
  Full_refit_QA = full_refit_qa,
  Input_alignment_QA = data_validation
)
for (sheet in names(sheet_data)) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, sheet_data[[sheet]])
  freezePane(wb, sheet, firstRow = TRUE)
  addFilter(wb, sheet, rows = 1, cols = seq_len(max(1, ncol(sheet_data[[sheet]]))))
  setColWidths(wb, sheet, cols = seq_len(max(1, ncol(sheet_data[[sheet]]))), widths = "auto")
}
saveWorkbook(wb, workbook_path, overwrite = TRUE)

contract_text <- c(
  "# Model C leave-one-country-out sensitivity analysis contract",
  "",
  "## Estimand and scope",
  "This analysis asks whether the fitted climate smooths from the six primary phenotype-specific Model C GAMMs are materially driven by any single country. It is a country-deletion influence analysis, not a held-out prediction or generalisation assessment.",
  "",
  "## Frozen elements",
  "- response: the same logit-transformed resistance outcome used by Model C",
  "- rows: the same lag-complete country-year analytical sample before deletion",
  "- predictors: the same precomputed within-climate-zone standardized and lagged climatic inputs and PLS components",
  "- model formula and basis dimensions: unchanged",
  "- fitting: Gaussian identity GAMM fitted by bam, REML and select=TRUE",
  "- primary models and downstream Results 2–3: unchanged",
  "",
  "## Country-deletion procedure",
  "For every phenotype, all rows from one country were removed and the complete model was re-estimated. The full-sample phenotype-exposure P2.5–P97.5 range was held fixed for every comparison.",
  "",
  "## Curve comparison",
  "Both raw term-centred changes and shape changes are reported. Shape comparisons align each LOCO smooth to the frozen smooth at the full-sample median exposure because mgcv imposes a separate centering constraint in every refit. This alignment is diagnostic only and does not alter the term-centred primary Figure 1.",
  "",
  "## Turning-point comparison",
  "Turning points were redetected with the same fixed geometry and EDF rules used for Figure 1 and matched to full-sample extrema by direction within 10% of the fixed P2.5–P97.5 search width. The LOCO match fraction is a deterministic deletion diagnostic, not a P value or bootstrap probability.",
  "",
  "## Interpretation limits",
  "The country-deletion envelope is descriptive and must not be called a confidence interval. LOCO influence analysis does not replace clustered uncertainty estimation, observation-level leverage diagnostics, residual checks or external predictive validation."
)
writeLines(contract_text, file.path(dirs$contract, "ModelC_LOCO_analysis_contract.md"))

methods_text <- paste0(
  "We assessed country-level influence using a leave-one-country-out refitting analysis. For each AMR phenotype, all country-year observations from one country were sequentially excluded and the complete Model C GAMM was re-estimated using the same preprocessed inputs, model formula, basis dimensions, Gaussian identity family, REML smoothing-parameter estimation and selection penalties as in the primary analysis. Climatic smooths from each restricted fit were evaluated over the corresponding full-sample phenotype–exposure-specific 2.5th–97.5th percentile range. Because penalized smooths are centred separately in each fit, curves were aligned to the complete-sample curve at the full-sample median exposure for shape comparisons; raw term-centred differences were retained separately. We summarized the maximum absolute change in fitted log odds, curve correlation, changes in effective degrees of freedom and overall smooth-term P values, and the preservation of full-sample turning points. Turning points were redetected using the same derivative-based procedure and matched by extremum direction within 10% of the fixed search range. Country-deletion envelopes and turning-point match fractions were treated as descriptive influence diagnostics rather than confidence intervals or probabilities."
)
writeLines(methods_text, file.path(dirs$report, "Recommended_Methods_ModelC_LOCO_sensitivity.txt"))

crab_summary <- panel_summary %>%
  filter(phenotype == "CR-Ab", climate == "Temperature")
all_panel_max <- max(panel_summary$worst_case_max_abs_delta_logOR_shape_P2_5_P97_5)
all_panel_central_max <- max(panel_summary$worst_case_max_abs_delta_logOR_shape_P10_P90)
full_significant_panels <- panel_summary %>% filter(full_smooth_p_value < 0.05)
n_sig_retained_90 <- sum(full_significant_panels$proportion_smooth_p_lt_0_05 >= 0.90)
n_turning_retained_90 <- sum(turning_point_stability$LOCO_match_fraction >= 0.90)
result_text <- c(
  "# Model C leave-one-country-out sensitivity analysis: result summary",
  "",
  paste0(
    "All ", sum(successful_counts$n_attempted_refits), " planned country-deletion models were attempted; ",
    sum(successful_counts$n_successful_refits), " converged successfully."
  ),
  paste0(
    "Across the 24 phenotype–climate panels, the largest worst-case P50-aligned curve change over the full-sample P2.5–P97.5 support was ",
    sprintf("%.3f", all_panel_max), " on the log-OR scale; the corresponding maximum within P10–P90 was ",
    sprintf("%.3f", all_panel_central_max), "."
  ),
  paste0(
    n_sig_retained_90, "/", nrow(full_significant_panels),
    " full-sample significant smooths remained P<0.05 in at least 90% of country-deletion refits, and ",
    n_turning_retained_90, "/", nrow(turning_point_stability),
    " full-sample turning points were matched in at least 90% of refits."
  ),
  paste0(
    "For CR-Ab temperature, the most influential deletion was ",
    crab_summary$worst_case_country_shape_P2_5_P97_5,
    " over P2.5–P97.5 (maximum aligned |delta log OR| ",
    sprintf("%.3f", crab_summary$worst_case_max_abs_delta_logOR_shape_P2_5_P97_5),
    "), whereas the worst-case change within P10–P90 was ",
    sprintf("%.3f", crab_summary$worst_case_max_abs_delta_logOR_shape_P10_P90),
    ". Its ", sprintf("%.1f °C", crab_turning_point$x_physical_equivalent),
    " maximum was matched in ", crab_turning_point$n_matched_LOCO_refits, "/",
    crab_turning_point$n_successful_LOCO_refits, " refits, with a descriptive central location range of ",
    sprintf("%.1f–%.1f °C", crab_turning_point$location_lower_2_5pct_physical_LOCO,
            crab_turning_point$location_upper_97_5pct_physical_LOCO), "."
  ),
  "",
  "These descriptive results must be interpreted together with the detailed panel table and the targeted CR-Ab figure. A high full-model R-squared is not validated merely by LOCO stability; the model-fit table only shows whether that statistic depends strongly on one country."
)
writeLines(result_text, file.path(dirs$report, "ModelC_LOCO_results_summary.md"))

paper_assessment <- c(
  "# Applicability of the published malaria-study hold-out strategy",
  "",
  "The malaria study uses hold-out analyses for two related but distinct purposes. Its published modelling supplement reports ordinary five-fold and sequential entire-survey hold-out validation for predictive performance, and compares learned fixed-effect coefficients across restricted fits. Its peer-review history additionally reports leave-one-country, leave-one-year and leave-one-month influence analyses in which the statistical model was re-estimated and the distributions of temperature coefficients and P values were compared with the full-sample estimates.",
  "",
  "The second use is directly transferable to the climate–AMR reviewer question because countries are the natural clustering unit and the concern is whether one country drives a fitted relationship. A direct coefficient comparison is not sufficient for Model C, however, because each climatic association is a penalized smooth rather than a single linear or quadratic coefficient. The present analysis therefore compares the complete smooth function, its EDF and overall test, and the locations of internal turning points.",
  "",
  "This LOCO refitting analysis is intentionally not labelled cross-validation because it does not score predictions for the omitted country. It diagnoses country-level influence on the fitted associations and complements, but does not replace, observation-level leverage/Cook-distance diagnostics, residual assessment or country-cluster uncertainty analyses."
)
writeLines(
  paper_assessment,
  file.path(dirs$report, "Scientific_assessment_of_malaria_paper_LOCO_applicability.md")
)

qa_summary <- tibble(
  check = c(
    "Six frozen models aligned to prepared data",
    "Full-sample refits exactly reproduce frozen models",
    "All planned country deletions attempted",
    "All country-deletion refits converged",
    "All 24 phenotype-climate panels summarized",
    "All full-sample turning points represented in LOCO summary",
    "All curve values finite",
    "All required figure PDF files created"
  ),
  passed = c(
    all(data_validation$passed),
    all(full_refit_qa$converged) &&
      all(full_refit_qa$max_abs_coefficient_difference <= 1e-10),
    nrow(refit_status) == sum(map_int(prepared, ~ n_distinct(.x$NAME))),
    all(refit_status$converged),
    nrow(panel_summary) == 24L,
    nrow(turning_point_stability) == nrow(full_candidates),
    all(is.finite(loco_curves$loco_eta_aligned_at_full_P50)),
    all(file.exists(file.path(
      dirs$figures,
      c(
        "Supplementary_Figure_ModelC_LOCO_curve_envelopes.pdf",
        "Supplementary_Figure_ModelC_LOCO_summary_heatmaps.pdf",
        "Supplementary_Figure_ModelC_LOCO_CRAb_temperature_targeted.pdf"
      )
    )))
  ),
  details = c(
    paste0(sum(data_validation$passed), "/", nrow(data_validation), " checks passed"),
    paste0("max coefficient delta=", max(full_refit_qa$max_abs_coefficient_difference)),
    paste0(nrow(refit_status), " attempted"),
    paste0(sum(refit_status$converged), "/", nrow(refit_status), " converged"),
    paste0(nrow(panel_summary), " panels"),
    paste0(nrow(turning_point_stability), "/", nrow(full_candidates), " candidates"),
    paste0(nrow(loco_curves), " finite curve rows"),
    "PDF, SVG, preview PNG and 600-dpi TIFF exported"
  )
)
write_csv(qa_summary, file.path(dirs$diagnostics, "ModelC_LOCO_QA_summary.csv"))

capture.output(
  sessionInfo(),
  file = file.path(dirs$logs, "ModelC_LOCO_sessionInfo.txt")
)

manifest_files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[!dir.exists(manifest_files)]
manifest <- tibble(
  relative_path = sub(paste0("^", output_root, "/?"), "", manifest_files),
  size_bytes = file.info(manifest_files)$size,
  md5 = unname(tools::md5sum(manifest_files))
) %>%
  arrange(relative_path)
write_csv(manifest, file.path(output_root, "MANIFEST.csv"))

log_line(
  "Completed. Successful refits: ", sum(refit_status$converged), "/", nrow(refit_status),
  "; panels: ", nrow(panel_summary), "; turning points: ", nrow(turning_point_stability), "."
)
cat("OUTPUT_ROOT=", output_root, "\n", sep = "")
