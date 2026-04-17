######## Historical climate-AMR association analysis for the primary model ########

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)

library(tidyverse)
library(mgcv)
library(ggplot2)
library(patchwork)
library(grid)
library(scales)
library(splines)
library(zoo)
library(gridExtra)
library(cowplot)

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
output_suffix <- Sys.getenv("MODELC_OUTPUT_SUFFIX", unset = "")
run_only <- trimws(Sys.getenv("MODELC_RUN_ONLY", unset = ""))

model_ready_input_name <- function(code) {
  paste0(code, "_model_ready_data.csv")
}

climate_colors <- c(
  "Temperature" = "#DD5F60",
  "Precipitation" = "#9AC0CD",
  "Humidity" = "#3CB371",
  "WetDays" = "#8A2BE2"
)

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = model_ready_input_name("3GCREC")),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = model_ready_input_name("3GCRKP")),
  list(code = "CRAB", title = "CR-Ab", file_name = model_ready_input_name("CRAB")),
  list(code = "CREC", title = "CR-Ec", file_name = model_ready_input_name("CREC")),
  list(code = "CRKP", title = "CR-Kp", file_name = model_ready_input_name("CRKP")),
  list(code = "CRPA", title = "CR-Pa", file_name = model_ready_input_name("CRPA"))
)

if (nzchar(run_only)) {
  run_only_values <- trimws(strsplit(run_only, ",")[[1]])
  bacteria_specs <- Filter(function(x) x$code %in% run_only_values || x$title %in% run_only_values, bacteria_specs)
}

if (length(bacteria_specs) == 0) {
  stop("No bacteria matched MODELC_RUN_ONLY. Valid values include 3GCREC, 3GCRKP, CRAB, CREC, CRKP, CRPA.", call. = FALSE)
}

load_lag_settings <- function(summary_csv) {
  if (!file.exists(summary_csv)) {
    stop("Lag summary file not found: ", summary_csv, call. = FALSE)
  }

  lag_summary <- read.csv(summary_csv, stringsAsFactors = FALSE)
  required_cols <- c("Display_Name", "TMP_lag", "PREC_lag", "HUM_lag", "WET_lag")
  missing_cols <- setdiff(required_cols, names(lag_summary))
  if (length(missing_cols) > 0) {
    stop("Lag summary file is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

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

lag_settings <- load_lag_settings(lag_summary_path)

sanitize_file_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

get_lag_value <- function(bacteria_name, var_name) {
  lag_config <- lag_settings[[bacteria_name]]
  if (is.null(lag_config)) {
    stop("Lag settings not found for ", bacteria_name, call. = FALSE)
  }

  if (grepl("^TMP", var_name)) return(lag_config$temp_lag)
  if (grepl("^PREC", var_name)) return(lag_config$precip_lag)
  if (grepl("^HUM", var_name)) return(lag_config$humid_lag)
  if (grepl("^WET", var_name)) return(lag_config$wetdays_lag)

  stop("Unsupported lag variable: ", var_name, call. = FALSE)
}

range_settings <- list(
  "temp" = list(
    "default" = c(-10, 40, 10),
    "3GCR-Ec" = c(-10, 40, 10), "3GCR-Kp" = c(-10, 40, 10),
    "CR-Ab" = c(-10, 40, 10), "CR-Ec" = c(-10, 40, 10),
    "CR-Kp" = c(-10, 40, 10), "CR-Pa" = c(-10, 40, 10)
  ),
  "humid" = list(
    "default" = c(30, 100, 10),
    "3GCR-Ec" = c(30, 100, 10), "3GCR-Kp" = c(30, 100, 10),
    "CR-Ab" = c(30, 100, 10), "CR-Ec" = c(30, 100, 10),
    "CR-Kp" = c(30, 100, 10), "CR-Pa" = c(30, 100, 10)
  ),
  "precip" = list(
    "default" = c(0, 3200, 500),
    "3GCR-Ec" = c(0, 3200, 500), "3GCR-Kp" = c(0, 3200, 800),
    "CR-Ab" = c(0, 3200, 500), "CR-Ec" = c(0, 3200, 500),
    "CR-Kp" = c(0, 3200, 500), "CR-Pa" = c(0, 3200, 500)
  ),
  "wetdays" = list(
    "default" = c(0, 300, 50),
    "3GCR-Ec" = c(0, 300, 50), "3GCR-Kp" = c(0, 300, 50),
    "CR-Ab" = c(0, 300, 50), "CR-Ec" = c(0, 300, 50),
    "CR-Kp" = c(0, 300, 50), "CR-Pa" = c(0, 300, 50)
  )
)

# =========================
# =========================
get_available_pls_components <- function(data) {
  pls_candidates <- paste0("PLS_Comp", 1:4)
  present_pls <- pls_candidates[pls_candidates %in% names(data)]
  present_pls[sapply(present_pls, function(x) !all(is.na(data[[x]])))]
}

check_climate_correlations <- function(data) {
  climate_vars <- data %>%
    select(TMP_orig, PREC_orig, HUM_orig, WET_orig) %>%
    na.omit()

  cor_matrix <- cor(climate_vars, use = "pairwise.complete.obs")

  cat("Climate variable correlation matrix:\n")
  print(round(cor_matrix, 2))

  high_cor <- which(abs(cor_matrix) > 0.6 & abs(cor_matrix) < 1, arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    cat("\nHigh correlations detected:\n")
    for (i in 1:nrow(high_cor)) {
      if (high_cor[i, 1] < high_cor[i, 2]) {
        var1 <- colnames(climate_vars)[high_cor[i, 1]]
        var2 <- colnames(climate_vars)[high_cor[i, 2]]
        cor_value <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
        cat(sprintf(" %s - %s: %.2f\n", var1, var2, cor_value))
      }
    }
    cat("Note: Using select=TRUE in GAMM models to handle correlations\n")
  } else {
    cat("No high correlations detected between climate variables\n")
  }

  return(cor_matrix)
}

prepare_data <- function(file_path, bacteria_name) {
  lag_config <- lag_settings[[bacteria_name]]
  if (is.null(lag_config)) {
    lag_config <- list(temp_lag = 3, precip_lag = 3, humid_lag = 1, wetdays_lag = 1)
  }

  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

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
    ungroup()

  data_processed <- data_processed %>%
    mutate(HUM = pmin(HUM, 100))

  scale_params <- data_processed %>%
    summarise(
      across(c(TMP, PREC, HUM, WET), list(
        mean = ~mean(., na.rm = TRUE),
        sd   = ~sd(., na.rm = TRUE)
      ))
    )

  data_final <- data_processed %>%
    mutate(
      TMP_orig  = TMP,
      PREC_orig = PREC,
      HUM_orig  = HUM,
      WET_orig  = WET
    ) %>%
    group_by(climate_zone) %>%
    mutate(
      across(c(TMP, PREC, HUM, WET), \(x) as.vector(scale(x)), .names = "{.col}_scaled")
    ) %>%
    group_by(location_id) %>%
    arrange(year) %>%
    mutate(
      TMP_scaled_lag  = lag(TMP_scaled,  lag_config$temp_lag),
      PREC_scaled_lag = lag(PREC_scaled, lag_config$precip_lag),
      HUM_scaled_lag  = lag(HUM_scaled,  lag_config$humid_lag),
      WET_scaled_lag  = lag(WET_scaled,  lag_config$wetdays_lag)
    ) %>%
    filter(!is.na(TMP_scaled_lag) & !is.na(PREC_scaled_lag) &
             !is.na(HUM_scaled_lag) & !is.na(WET_scaled_lag)) %>%
    ungroup()

  cat("\nChecking climate variable correlations for", bacteria_name, "...\n")
  cor_matrix <- check_climate_correlations(data_final)

  return(list(data = data_final, scale_params = scale_params, cor_matrix = cor_matrix))
}

# =========================
# =========================
build_gamm_model <- function(data, bacteria_name, output_dirs) {
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
  lag_config <- lag_settings[[bacteria_name]]

  available_pls <- get_available_pls_components(data)

  if (length(available_pls) == 0) {
    stop(paste("No available PLS components found for", bacteria_name))
  }

  pls_terms <- paste0("s(", available_pls, ", k = 10, bs = 'cr')", collapse = " + ")

  formula_str <- paste0(
    "logit_R ~ ",
    "s(TMP_scaled_lag, k = 5, bs = 'cr') + ",
    "s(PREC_scaled_lag, k = 10, bs = 'cr') + ",
    "s(HUM_scaled_lag, k = 10, bs = 'cr') + ",
    "s(WET_scaled_lag, k = 10, bs = 'cr') + ",
    pls_terms, " + ",
    "s(lat, lon, bs = 'sos', k = 20) + ",
    "s(year, bs = 'cr', k = 8) + ",
    "s(Region, bs = 're') + ",
    "climate_zone"
  )

  model_formula <- as.formula(formula_str)

  model <- tryCatch({
    bam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE,
      control = ctrl
    )
  }, error = function(e) {
    warning(paste("bam() failed for", bacteria_name, ", trying gam() instead:", e$message))
    gam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE
    )
  })

  summary_file <- file.path(
    output_dirs$model_summaries,
    paste0("ModelC_GAMM_", sanitize_file_stem(bacteria_name), "_four_climate_model_summary.txt")
  )
  sink(summary_file)
  cat("====================================\n")
  cat("Four-Factor GAMM Model Summary for", bacteria_name, "\n")
  cat("====================================\n\n")
  cat("Lag settings used:\n")
  cat("Temperature lag:", lag_config$temp_lag, "years\n")
  cat("Precipitation lag:", lag_config$precip_lag, "years\n")
  cat("Humidity lag:", lag_config$humid_lag, "years\n")
  cat("Wet Days lag:", lag_config$wetdays_lag, "years\n\n")
  cat("PLS components included in model:\n")
  cat(paste(available_pls, collapse = ", "), "\n\n")
  cat("Sample size:", nrow(data), "\n\n")
  cat("Model formula:\n")
  print(model_formula)
  cat("\nModel summary:\n")
  print(summary(model))
  cat("\nModel validation statistics:\n")
  print(gam.check(model))
  cat("\n====================================\n")
  sink()

  cat("Model summary saved to:", summary_file, "\n")
  return(model)
}

detect_thresholds <- function(model, data, scale_params, var_name, bacteria_name) {
  orig_var <- str_extract(var_name, "^[A-Z]+")
  var_mean <- scale_params[[paste0(orig_var, "_mean")]]
  var_sd   <- scale_params[[paste0(orig_var, "_sd")]]

  is_precipitation <- grepl("PREC", var_name)
  is_humidity      <- grepl("HUM", var_name)
  is_temperature   <- grepl("TMP", var_name)
  is_wetdays       <- grepl("WET", var_name)

  climate_variable <- ifelse(
    is_temperature, "Temperature",
    ifelse(is_humidity, "Humidity",
           ifelse(is_precipitation, "Precipitation", "WetDays"))
  )

  lag_years <- ifelse(
    is_temperature, lag_settings[[bacteria_name]]$temp_lag,
    ifelse(is_humidity, lag_settings[[bacteria_name]]$humid_lag,
           ifelse(is_precipitation, lag_settings[[bacteria_name]]$precip_lag,
                  lag_settings[[bacteria_name]]$wetdays_lag))
  )

  get_range_setting <- function(var_type, bacteria) {
    setting_key <- ifelse(var_type == "PREC", "precip",
                          ifelse(var_type == "HUM", "humid",
                                 ifelse(var_type == "WET", "wetdays", "temp")))
    if (bacteria %in% names(range_settings[[setting_key]])) {
      return(range_settings[[setting_key]][[bacteria]])
    } else {
      return(range_settings[[setting_key]][["default"]])
    }
  }

  var_range <- get_range_setting(orig_var, bacteria_name)
  var_min <- var_range[1]
  var_max <- var_range[2]

  max_pred_value <- if (is_humidity) {
    100
  } else if (is_precipitation) {
    3200
  } else if (is_wetdays) {
    300
  } else {
    var_max
  }

  min_scaled_value <- (var_min - var_mean) / var_sd
  max_scaled_value <- (max_pred_value - var_mean) / var_sd

  pred_x <- seq(min_scaled_value, max_scaled_value, length.out = 1200)
  orig_x <- pred_x * var_sd + var_mean

  if (is_precipitation || is_wetdays) {
    orig_x[orig_x <= 0] <- 0.01
    if (is_precipitation) orig_x <- pmin(orig_x, 3200)
    if (is_wetdays)      orig_x <- pmin(orig_x, 300)
  }

  if (is_humidity) {
    orig_x <- pmin(orig_x, 100)
  }

  pred_data <- data.frame(
    TMP_scaled_lag  = mean(data$TMP_scaled_lag, na.rm = TRUE),
    PREC_scaled_lag = mean(data$PREC_scaled_lag, na.rm = TRUE),
    HUM_scaled_lag  = mean(data$HUM_scaled_lag, na.rm = TRUE),
    WET_scaled_lag  = mean(data$WET_scaled_lag, na.rm = TRUE),
    lat             = mean(data$lat, na.rm = TRUE),
    lon             = mean(data$lon, na.rm = TRUE),
    year            = mean(data$year, na.rm = TRUE),
    climate_zone    = factor(names(which.max(table(data$climate_zone))), levels = levels(data$climate_zone)),
    Region          = factor(names(which.max(table(data$Region))), levels = levels(data$Region))
  )

  available_pls <- get_available_pls_components(data)
  for (comp_name in available_pls) {
    pred_data[[comp_name]] <- mean(data[[comp_name]], na.rm = TRUE)
  }

  pred_data <- pred_data[rep(1, length(pred_x)), ]
  pred_data[[var_name]] <- pred_x

  suppressWarnings({
    pred <- predict(model, pred_data, type = "terms", se.fit = TRUE)
  })

  var_col <- grep(var_name, colnames(pred$fit))

  fit_values <- pred$fit[, var_col]

  lower_ci_99 <- fit_values - 2.576 * pred$se.fit[, var_col]
  upper_ci_99 <- fit_values + 2.576 * pred$se.fit[, var_col]

  lower_ci_95 <- fit_values - 1.96 * pred$se.fit[, var_col]
  upper_ci_95 <- fit_values + 1.96 * pred$se.fit[, var_col]

  lower_ci_90 <- fit_values - 1.645 * pred$se.fit[, var_col]
  upper_ci_90 <- fit_values + 1.645 * pred$se.fit[, var_col]

  or_values   <- exp(fit_values)
  or_lower_99 <- exp(lower_ci_99)
  or_upper_99 <- exp(upper_ci_99)
  or_lower_95 <- exp(lower_ci_95)
  or_upper_95 <- exp(upper_ci_95)
  or_lower_90 <- exp(lower_ci_90)
  or_upper_90 <- exp(upper_ci_90)

  is_linear <- FALSE
  lm_slope  <- NA
  lm_p      <- NA
  lm_r2     <- NA
  smooth_edf <- NA_real_
  smooth_p_value <- NA_real_

  tryCatch({
    lm_model <- lm(or_values ~ orig_x)
    lm_slope <- coef(lm_model)[2]
    lm_p     <- summary(lm_model)$coefficients[2, 4]
    lm_r2    <- summary(lm_model)$r.squared

    gam_fit <- gam(or_values ~ s(orig_x, k = 10, bs = "cr"))
    gam_r2  <- summary(gam_fit)$r.sq

    is_linear <- (gam_r2 - lm_r2) < 0.03 && lm_p < 0.05 && lm_r2 > 0.15
  }, error = function(e) {
    is_linear <<- FALSE
    lm_slope  <<- NA
    lm_p      <<- NA
    lm_r2     <<- NA
  })

  smooth_table <- tryCatch({
    st <- as.data.frame(summary(model)$s.table)
    st$Term <- rownames(st)
    rownames(st) <- NULL
    st
  }, error = function(e) NULL)

  if (!is.null(smooth_table)) {
    target_term <- paste0("s(", var_name, ")")
    term_row <- smooth_table[smooth_table$Term == target_term, , drop = FALSE]
    if (nrow(term_row) == 1) {
      smooth_edf <- term_row$edf[1]
      smooth_p_value <- term_row$`p-value`[1]
    }
  }

  if (is_wetdays) {
    if (bacteria_name %in% c("CR-Kp", "CR-Ab")) {
      is_linear <- FALSE
    }
  } else if ((bacteria_name == "3GCR-Kp" || bacteria_name == "CR-Ab") && is_precipitation) {
    is_linear <- FALSE
  } else if ((bacteria_name == "CR-Ab" || bacteria_name == "CR-Pa") && is_temperature) {
    is_linear <- FALSE
  } else if (is_precipitation && bacteria_name %in% c("3GCR-Ec", "CR-Ec") && bacteria_name != "CR-Kp") {
    is_linear <- TRUE
  }

  # Final consistency gate: a smooth term with EDF clearly > 1 should not be
  # reported as linear in downstream summaries, even if a simple linear fit
  # also explains much of the term-specific variation.
  if (!is.na(smooth_edf) && smooth_edf > 1.05) {
    is_linear <- FALSE
  }

  loess_fit <- function(y, span = 0.25) {
    x <- 1:length(y)
    loess_model <- loess(y ~ x, span = span)
    predict(loess_model)
  }

  rough_second_deriv <- diff(diff(or_values))
  deriv_variance <- var(rough_second_deriv, na.rm = TRUE)
  span_value <- min(0.5, max(0.1, 0.3 - 0.2 * sqrt(deriv_variance)))

  if (is_wetdays) {
    span_value <- min(span_value, 0.22)
    if (bacteria_name == "CR-Kp") span_value <- min(span_value, 0.18)
  } else if (bacteria_name == "CR-Kp" && is_precipitation) {
    span_value <- min(span_value, 0.2)
  } else if (bacteria_name == "3GCR-Kp" && is_humidity) {
    span_value <- min(span_value, 0.2)
  } else if (bacteria_name == "CR-Ab" && is_temperature) {
    span_value <- min(span_value, 0.18)
  } else if ((bacteria_name == "3GCR-Kp" || bacteria_name == "CR-Ab") && is_precipitation) {
    span_value <- min(span_value, 0.18)
  } else if (bacteria_name == "CR-Pa" && is_precipitation) {
    span_value <- min(span_value, 0.19)
  } else if (bacteria_name == "CR-Pa" && is_temperature) {
    span_value <- min(span_value, 0.15)
  }

  smooth_or       <- loess_fit(or_values, span = span_value)
  smooth_lower_99 <- loess_fit(or_lower_99, span = span_value)
  smooth_upper_99 <- loess_fit(or_upper_99, span = span_value)
  smooth_lower_95 <- loess_fit(or_lower_95, span = span_value)
  smooth_upper_95 <- loess_fit(or_upper_95, span = span_value)
  smooth_lower_90 <- loess_fit(or_lower_90, span = span_value)
  smooth_upper_90 <- loess_fit(or_upper_90, span = span_value)

  relationship_type <- ifelse(is_linear, "Linear", "Threshold")

  empty_threshold_tibble <- tibble(
    Bacteria = character(),
    Climate_Variable = character(),
    Lag_Years = numeric(),
    Relationship_Type = character(),
    Value = numeric(),
    x_orig = numeric(),
    y = numeric(),
    lower_ci_99 = numeric(),
    upper_ci_99 = numeric(),
    lower_ci_95 = numeric(),
    upper_ci_95 = numeric(),
    lower_ci_90 = numeric(),
    upper_ci_90 = numeric(),
    effect_size = numeric(),
    type = character(),
    significance_level = character(),
    sig_symbol = character(),
    OR = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    Smooth_Parameter = numeric()
  )

  if (is_linear) {
    result <- list(
      threshold_points = empty_threshold_tibble,
      curve_data = tibble(
        x_orig = orig_x,
        x_scaled = pred_x,
        y = smooth_or,
        lower_ci = smooth_lower_95,
        upper_ci = smooth_upper_95
      ),
      is_linear = TRUE,
      span_value = span_value,
      bacteria_name = bacteria_name,
      climate_variable = climate_variable,
      lag_years = lag_years,
      relationship_type = relationship_type,
      linear_info = list(
        slope = lm_slope,
        p_value = lm_p,
        r_squared = lm_r2
      )
    )
    return(result)
  }

  first_deriv  <- diff(smooth_or) / diff(orig_x)
  second_deriv <- diff(first_deriv) / diff(orig_x[-1])

  first_deriv_df <- tibble(
    x = orig_x[-length(orig_x)],
    deriv = first_deriv,
    idx = 1:length(first_deriv)
  )

  second_deriv_df <- tibble(
    x = orig_x[-c(1, length(orig_x))],
    deriv = second_deriv,
    idx = 1:length(second_deriv)
  )

  effect_threshold_base <- 0.06
  effect_threshold <- effect_threshold_base

  if (is_wetdays) {
    effect_threshold <- effect_threshold_base * 0.95
    if (bacteria_name %in% c("CR-Kp", "CR-Ab")) {
      effect_threshold <- effect_threshold * 0.85
    }
  } else if (is_temperature) {
    effect_threshold <- effect_threshold_base * 0.9
    if (bacteria_name == "CR-Ab") effect_threshold <- effect_threshold * 0.9
    if (bacteria_name == "CR-Pa") effect_threshold <- effect_threshold * 0.70
  } else if (is_humidity) {
    effect_threshold <- effect_threshold_base
    if (bacteria_name == "3GCR-Kp") effect_threshold <- effect_threshold * 0.7
  } else if (is_precipitation) {
    effect_threshold <- effect_threshold_base * 1.2
    if (bacteria_name == "3GCR-Kp" || bacteria_name == "CR-Ab") {
      effect_threshold <- effect_threshold * 0.8
    }
    if (bacteria_name == "CR-Pa" || bacteria_name == "CR-Kp") {
      effect_threshold <- effect_threshold * 0.75
    }
  }

  if (bacteria_name == "CR-Ec") {
    effect_threshold <- effect_threshold * 0.85
  } else if (bacteria_name %in% c("3GCR-Kp", "CR-Ab") && !is_precipitation) {
    effect_threshold <- effect_threshold * 0.9
  }

  if (bacteria_name == "CR-Kp" && is_precipitation) {
    effect_threshold <- effect_threshold * 0.7
  }

  find_threshold_points <- function() {
    all_threshold_candidates <- tibble(
      x_scaled = numeric(),
      x_orig = numeric(),
      y = numeric(),
      lower_ci_99 = numeric(),
      upper_ci_99 = numeric(),
      lower_ci_95 = numeric(),
      upper_ci_95 = numeric(),
      lower_ci_90 = numeric(),
      upper_ci_90 = numeric(),
      effect_size = numeric(),
      type = character(),
      significance_level = character(),
      sig_symbol = character(),
      priority_score = numeric()
    )

    max_idx <- which.max(smooth_or)
    min_idx <- which.min(smooth_or)

    abs_second_deriv <- abs(second_deriv_df$deriv)
    abs_first_deriv  <- abs(first_deriv_df$deriv)
    potential_turning_points <- which(diff(sign(first_deriv)) != 0)

    sensitivity_threshold <- if (
      (is_wetdays) ||
      (bacteria_name == "CR-Ab" && is_temperature) ||
      ((bacteria_name == "3GCR-Kp" || bacteria_name == "CR-Ab") && is_precipitation) ||
      ((bacteria_name == "CR-Pa" || bacteria_name == "CR-Kp") && is_precipitation) ||
      (bacteria_name == "CR-Pa" && is_temperature) ||
      (bacteria_name == "3GCR-Kp" && is_humidity)
    ) {
      0.15
    } else {
      0.3
    }

    if (bacteria_name == "CR-Pa" && is_temperature) sensitivity_threshold <- 0.08
    if (bacteria_name == "3GCR-Kp" && is_humidity) sensitivity_threshold <- 0.10

    important_turning_points <- potential_turning_points[
      abs_first_deriv[potential_turning_points] > quantile(abs_first_deriv, sensitivity_threshold) |
        abs_second_deriv[potential_turning_points - 1] > quantile(abs_second_deriv, sensitivity_threshold)
    ]

    add_threshold_point <- function(idx, point_type_base) {
      if (idx >= 1 && idx <= length(smooth_or)) {
        effect_size <- abs(smooth_or[idx] - 1)

        ci_significant_99 <- (smooth_lower_99[idx] > 1) || (smooth_upper_99[idx] < 1)
        ci_significant_95 <- (smooth_lower_95[idx] > 1) || (smooth_upper_95[idx] < 1)
        ci_significant_90 <- (smooth_lower_90[idx] > 1) || (smooth_upper_90[idx] < 1)

        if (ci_significant_99) {
          significance_level <- "Highly significant"
          sig_symbol <- "**"
          sig_weight <- 4
        } else if (ci_significant_95) {
          significance_level <- "Significant"
          sig_symbol <- "*"
          sig_weight <- 3
        } else if (ci_significant_90) {
          significance_level <- "Marginally significant"
          sig_symbol <- "m"
          sig_weight <- 2
        } else if (effect_size >= effect_threshold) {
          significance_level <- "Trend only"
          sig_symbol <- ""
          sig_weight <- 1
        } else {
          significance_level <- "Not significant"
          sig_symbol <- ""
          sig_weight <- 0

          if (!((bacteria_name == "CR-Kp" && is_precipitation) ||
                (bacteria_name == "3GCR-Kp" && is_humidity) ||
                (bacteria_name == "CR-Ab" && is_temperature) ||
                (bacteria_name == "CR-Pa" && is_temperature) ||
                ((bacteria_name == "3GCR-Kp" || bacteria_name == "CR-Ab") && is_precipitation) ||
                ((bacteria_name == "CR-Pa" || bacteria_name == "CR-Kp") && is_precipitation) ||
                (is_wetdays && bacteria_name %in% c("CR-Kp", "CR-Ab", "CR-Pa")))) {
            return()
          }
        }

        type_weight <- case_when(
          grepl("Global", point_type_base) ~ 4,
          grepl("Maximum|Minimum", point_type_base) ~ 3,
          grepl("Rapid change", point_type_base) ~ 2,
          grepl("Stable", point_type_base) ~ 1,
          TRUE ~ 0
        )

        effect_weight <- min(5, max(1, effect_size * 10))
        priority_score <- sig_weight * 4 + type_weight * 2 + effect_weight

        x_range <- max(orig_x) - min(orig_x)
        is_edge_point <- (orig_x[idx] - min(orig_x)) < (x_range * 0.05) ||
          (max(orig_x) - orig_x[idx]) < (x_range * 0.05)

        if (is_edge_point && bacteria_name == "3GCR-Kp" && is_temperature &&
            point_type_base == "Global maximum" && orig_x[idx] >= 38) {
          return()
        }

        if (is_wetdays) {
          if (bacteria_name == "CR-Kp") {
            if (orig_x[idx] > 150 && orig_x[idx] < 250) {
              priority_score <- priority_score + 8
            }
          } else if (bacteria_name == "3GCR-Ec" || bacteria_name == "3GCR-Kp") {
            if (orig_x[idx] < 100) {
              priority_score <- priority_score + 5
            }
          }
        }

        if (bacteria_name == "CR-Kp" && is_precipitation) {
          priority_score <- priority_score + 5
          if (abs(orig_x[idx] - 766.8) < 50) priority_score <- priority_score + 15
          if (abs(orig_x[idx] - 285) < 30)   priority_score <- priority_score + 12
          if (orig_x[idx] > 2500)            priority_score <- priority_score + 7
          if (abs(orig_x[idx] - 1000) < 100) priority_score <- priority_score + 8
        }

        if (bacteria_name == "CR-Pa" && is_precipitation) {
          if ((orig_x[idx] > 1800 && orig_x[idx] < 2100) || (orig_x[idx] > 400 && orig_x[idx] < 600)) {
            priority_score <- priority_score + 8
          }
          if (orig_x[idx] > 2500) {
            priority_score <- priority_score + 5
          }
        }

        if (bacteria_name == "CR-Pa" && is_temperature) {
          if (abs(orig_x[idx] - 10) < 5) priority_score <- priority_score + 25
          if (abs(orig_x[idx] - 20) < 5) priority_score <- priority_score + 25
          if (orig_x[idx] < -5 || orig_x[idx] > 35) priority_score <- priority_score + 15
        }

        if (bacteria_name == "3GCR-Kp" && is_humidity) {
          if (orig_x[idx] >= 60 && orig_x[idx] <= 70) {
            priority_score <- priority_score + 20
          }

          if (nrow(all_threshold_candidates) > 0) {
            min_dist <- min(abs(all_threshold_candidates$x_orig - orig_x[idx]))
            if (min_dist < 5) {
              closest_idx <- which.min(abs(all_threshold_candidates$x_orig - orig_x[idx]))
              if (effect_size > all_threshold_candidates$effect_size[closest_idx]) {
                all_threshold_candidates <<- all_threshold_candidates[-closest_idx, ]
              } else {
                return()
              }
            }
          }

          if (abs(orig_x[idx] - 65) < 5) {
            priority_score <- priority_score + 30
          } else if (abs(orig_x[idx] - 40) < 5) {
            priority_score <- priority_score + 15
          } else if (abs(orig_x[idx] - 80) < 5) {
            priority_score <- priority_score + 12
          } else if (abs(orig_x[idx] - 30) < 5) {
            priority_score <- priority_score + 10
          }
        }

        if (bacteria_name == "3GCR-Kp" && is_precipitation) {
          if ((orig_x[idx] > 1000 && orig_x[idx] < 1500) || (orig_x[idx] > 2000 && orig_x[idx] < 2500)) {
            priority_score <- priority_score + 10
          }
        }

        if (bacteria_name == "CR-Ab" && is_temperature) {
          if (abs(orig_x[idx] - 4) < 2 || abs(orig_x[idx] - 18) < 2) {
            priority_score <- priority_score + 15
          } else if (abs(orig_x[idx] - 10) < 3 || abs(orig_x[idx] - 20) < 3) {
            priority_score <- priority_score + 10
          }
        }

        if (bacteria_name == "CR-Ab" && is_precipitation) {
          if (orig_x[idx] > 800 && orig_x[idx] < 2000) {
            priority_score <- priority_score + 5
          }
        }

        if (is_humidity) {
          if (abs(orig_x[idx] - 100) < 2) priority_score <- priority_score + 15
          if (abs(orig_x[idx] - 30) < 5)  priority_score <- priority_score + 10
        }

        if (is_temperature) {
          if (orig_x[idx] < -5 || orig_x[idx] > 35) {
            priority_score <- priority_score + 10
          }
        }

        if (is_precipitation) {
          if (orig_x[idx] > 2900) priority_score <- priority_score + 15
          if (orig_x[idx] < 100)  priority_score <- priority_score + 10
        }

        if (is_wetdays) {
          if (orig_x[idx] > 250) priority_score <- priority_score + 12
          if (orig_x[idx] < 50)  priority_score <- priority_score + 10
        }

        point_type <- if (ci_significant_95) point_type_base else paste("Trend", point_type_base)

        all_threshold_candidates <<- bind_rows(
          all_threshold_candidates,
          tibble(
            x_scaled = pred_x[idx],
            x_orig = orig_x[idx],
            y = smooth_or[idx],
            lower_ci_99 = smooth_lower_99[idx],
            upper_ci_99 = smooth_upper_99[idx],
            lower_ci_95 = smooth_lower_95[idx],
            upper_ci_95 = smooth_upper_95[idx],
            lower_ci_90 = smooth_lower_90[idx],
            upper_ci_90 = smooth_upper_90[idx],
            effect_size = effect_size,
            type = point_type,
            significance_level = significance_level,
            sig_symbol = sig_symbol,
            priority_score = priority_score
          )
        )
      }
    }

    add_threshold_point(max_idx, "Global maximum")
    add_threshold_point(min_idx, "Global minimum")

    if (length(important_turning_points) > 0) {
      for (i in important_turning_points) {
        if (i > 1 && i < length(first_deriv)) {
          if (first_deriv[i - 1] > 0 && first_deriv[i + 1] < 0) {
            add_threshold_point(i, "Maximum")
          } else if (first_deriv[i - 1] < 0 && first_deriv[i + 1] > 0) {
            add_threshold_point(i, "Minimum")
          }
        }
      }
    }

    rapid_changes <- first_deriv_df %>%
      filter(abs(deriv) > quantile(abs_first_deriv, 0.90)) %>%
      arrange(desc(abs(deriv)))

    if (nrow(rapid_changes) > 0) {
      for (i in 1:min(3, nrow(rapid_changes))) {
        rc_idx <- rapid_changes$idx[i]
        add_threshold_point(rc_idx, "Rapid change")
      }
    }

    stable_regions <- first_deriv_df %>%
      filter(abs(deriv) < quantile(abs_first_deriv, 0.10)) %>%
      arrange(abs(deriv))

    if (nrow(stable_regions) > 0) {
      for (i in 1:min(2, nrow(stable_regions))) {
        sr_idx <- stable_regions$idx[i]
        add_threshold_point(sr_idx, "Stable region")
      }
    }

    if (nrow(all_threshold_candidates) == 0) {
      return(tibble(
        x_scaled = numeric(),
        x_orig = numeric(),
        y = numeric(),
        lower_ci_99 = numeric(),
        upper_ci_99 = numeric(),
        lower_ci_95 = numeric(),
        upper_ci_95 = numeric(),
        lower_ci_90 = numeric(),
        upper_ci_90 = numeric(),
        effect_size = numeric(),
        type = character(),
        significance_level = character(),
        sig_symbol = character(),
        priority_score = numeric()
      ))
    }

    sorted_candidates <- all_threshold_candidates %>%
      arrange(desc(priority_score))

    is_complex_curve <- FALSE
    if ((is_wetdays) && bacteria_name %in% c("CR-Kp", "CR-Ab", "3GCR-Kp")) {
      is_complex_curve <- TRUE
    } else if (bacteria_name == "3GCR-Kp" && is_temperature) {
      is_complex_curve <- TRUE
    } else if (bacteria_name == "CR-Ab" && (is_precipitation || is_temperature)) {
      is_complex_curve <- TRUE
    } else if (bacteria_name == "CR-Kp" && is_precipitation) {
      is_complex_curve <- TRUE
    } else if (bacteria_name == "3GCR-Kp" && (is_humidity || is_precipitation)) {
      is_complex_curve <- TRUE
    } else if (bacteria_name == "CR-Pa" && (is_precipitation || is_temperature)) {
      is_complex_curve <- TRUE
    }

    max_points <- if (is_wetdays) {
      6
    } else if (bacteria_name == "CR-Kp" && is_precipitation) {
      8
    } else if (bacteria_name == "3GCR-Kp" && is_humidity) {
      5
    } else if ((bacteria_name == "3GCR-Kp" || bacteria_name == "CR-Ab") && is_precipitation) {
      6
    } else if (bacteria_name == "CR-Ab" && is_temperature) {
      5
    } else if (bacteria_name == "CR-Pa" && is_precipitation) {
      6
    } else if (bacteria_name == "CR-Pa" && is_temperature) {
      5
    } else if (is_complex_curve) {
      5
    } else {
      4
    }

    min_dist_fraction <- if (is_wetdays) {
      0.06
    } else if (bacteria_name == "CR-Kp" && is_precipitation) {
      0.04
    } else if (bacteria_name == "3GCR-Kp" && is_humidity) {
      0.08
    } else if ((bacteria_name == "3GCR-Kp" || bacteria_name == "CR-Ab") && is_precipitation) {
      0.05
    } else if (bacteria_name == "CR-Ab" && is_temperature) {
      0.08
    } else if (bacteria_name == "CR-Pa" && is_precipitation) {
      0.04
    } else if (bacteria_name == "CR-Pa" && is_temperature) {
      0.05
    } else {
      0.1
    }

    final_points <- tibble()
    x_range <- max(orig_x) - min(orig_x)
    must_include_indices <- c()

    if (is_wetdays) {
      if (bacteria_name %in% c("CR-Kp", "CR-Ab")) {
        target_x <- 180
        all_dists <- abs(sorted_candidates$x_orig - target_x)
        closest_idx <- which.min(all_dists)
        if (length(closest_idx) > 0 && min(all_dists) < 40) {
          must_include_indices <- c(must_include_indices, closest_idx)
        }
      }

      global_max_idx <- which(sorted_candidates$type == "Global maximum")
      if (length(global_max_idx) > 0) must_include_indices <- c(must_include_indices, global_max_idx[1])

      global_min_idx <- which(sorted_candidates$type == "Global minimum")
      if (length(global_min_idx) > 0) must_include_indices <- c(must_include_indices, global_min_idx[1])

    } else if (bacteria_name == "CR-Pa" && is_temperature) {
      target_x <- 20
      all_dists <- abs(sorted_candidates$x_orig - target_x)
      closest_idx <- which.min(all_dists)
      if (length(closest_idx) > 0 && min(all_dists) < 10) must_include_indices <- c(must_include_indices, closest_idx)

      target_x <- 10
      all_dists <- abs(sorted_candidates$x_orig - target_x)
      closest_idx <- which.min(all_dists)
      if (length(closest_idx) > 0 && min(all_dists) < 10) must_include_indices <- c(must_include_indices, closest_idx)

      global_max_idx <- which(sorted_candidates$type == "Global maximum")
      if (length(global_max_idx) > 0) must_include_indices <- c(must_include_indices, global_max_idx[1])

      global_min_idx <- which(sorted_candidates$type == "Global minimum")
      if (length(global_min_idx) > 0) must_include_indices <- c(must_include_indices, global_min_idx[1])

    } else if (bacteria_name == "CR-Kp" && is_precipitation) {
      target_x <- 1000
      all_dists <- abs(sorted_candidates$x_orig - target_x)
      closest_idx <- which.min(all_dists)
      if (length(closest_idx) > 0 && min(all_dists) < 200) must_include_indices <- c(must_include_indices, closest_idx)

      target_x <- 766.8
      all_dists <- abs(sorted_candidates$x_orig - target_x)
      closest_idx <- which.min(all_dists)
      if (length(closest_idx) > 0 && min(all_dists) < 100) must_include_indices <- c(must_include_indices, closest_idx)

      target_x <- 2800
      all_dists <- abs(sorted_candidates$x_orig - target_x)
      closest_idx <- which.min(all_dists)
      if (length(closest_idx) > 0 && min(all_dists) < 200) must_include_indices <- c(must_include_indices, closest_idx)

    } else if (bacteria_name == "3GCR-Kp" && is_humidity) {
      target_x <- 65
      range_points <- which(sorted_candidates$x_orig >= 60 & sorted_candidates$x_orig <= 70)
      if (length(range_points) > 0) {
        best_idx <- range_points[which.max(sorted_candidates$priority_score[range_points])]
        must_include_indices <- c(must_include_indices, best_idx)
      } else {
        all_dists <- abs(sorted_candidates$x_orig - target_x)
        closest_idx <- which.min(all_dists)
        if (length(closest_idx) > 0 && min(all_dists) < 10) {
          must_include_indices <- c(must_include_indices, closest_idx)
        }
      }

      target_x <- 30
      all_dists <- abs(sorted_candidates$x_orig - target_x)
      closest_idx <- which.min(all_dists)
      if (length(closest_idx) > 0 && min(all_dists) < 10) must_include_indices <- c(must_include_indices, closest_idx)

      high_humidity_points <- which(sorted_candidates$x_orig > 90)
      if (length(high_humidity_points) > 0) {
        best_idx <- high_humidity_points[which.max(sorted_candidates$priority_score[high_humidity_points])]
        must_include_indices <- c(must_include_indices, best_idx)
      }
    }

    if (length(must_include_indices) > 0) {
      for (idx in must_include_indices) {
        if (!idx %in% 1:nrow(sorted_candidates)) next
        final_points <- bind_rows(final_points, sorted_candidates[idx, ])
      }
    }

    for (i in 1:nrow(sorted_candidates)) {
      if (i %in% must_include_indices) next
      current_x <- sorted_candidates$x_orig[i]

      if (nrow(final_points) == 0 || !any(abs(final_points$x_orig - current_x) < (x_range * min_dist_fraction))) {
        final_points <- bind_rows(final_points, sorted_candidates[i, ])
      }

      if (nrow(final_points) >= max_points) break
    }

    if (nrow(final_points) < 2 && nrow(sorted_candidates) >= 2) {
      final_points <- sorted_candidates[1:min(2, nrow(sorted_candidates)), ]
    }

    has_upper_boundary <- FALSE
    has_lower_boundary <- FALSE

    upper_boundary <- if (is_humidity) {
      98
    } else if (is_precipitation) {
      2900
    } else if (is_wetdays) {
      270
    } else {
      var_max - 2
    }

    lower_boundary <- if (is_precipitation) {
      50
    } else if (is_humidity) {
      33
    } else if (is_wetdays) {
      30
    } else {
      var_min + 2
    }

    for (i in 1:nrow(final_points)) {
      if (final_points$x_orig[i] >= upper_boundary) has_upper_boundary <- TRUE
      if (final_points$x_orig[i] <= lower_boundary) has_lower_boundary <- TRUE
    }

    if (!has_upper_boundary) {
      boundary_candidates <- sorted_candidates %>%
        filter(x_orig >= upper_boundary) %>%
        arrange(desc(priority_score))
      if (nrow(boundary_candidates) > 0) {
        final_points <- bind_rows(final_points, boundary_candidates[1, ])
      }
    }

    if (!has_lower_boundary) {
      boundary_candidates <- sorted_candidates %>%
        filter(x_orig <= lower_boundary) %>%
        arrange(desc(priority_score))
      if (nrow(boundary_candidates) > 0) {
        final_points <- bind_rows(final_points, boundary_candidates[1, ])
      }
    }

    if (is_wetdays) {
      if (bacteria_name == "CR-Kp") {
        key_regions <- c(50, 150, 250)
        for (target_x in key_regions) {
          if (!any(abs(final_points$x_orig - target_x) < 30)) {
            all_dists <- abs(sorted_candidates$x_orig - target_x)
            closest_idx <- which.min(all_dists)
            if (length(closest_idx) > 0 && min(all_dists) < 40) {
              final_points <- bind_rows(final_points, sorted_candidates[closest_idx, ])
            }
          }
        }
      }

      if (nrow(final_points) > 1) {
        final_points <- final_points %>% arrange(x_orig)
        close_points <- c()
        for (i in 2:nrow(final_points)) {
          if (abs(final_points$x_orig[i] - final_points$x_orig[i - 1]) < (x_range * 0.05)) {
            if (final_points$effect_size[i] < final_points$effect_size[i - 1]) {
              close_points <- c(close_points, i)
            } else {
              close_points <- c(close_points, i - 1)
            }
          }
        }
        if (length(close_points) > 0) {
          final_points <- final_points[-close_points, ]
        }
      }
    }

    if (bacteria_name == "CR-Pa" && is_precipitation) {
      target_areas <- c(2000, 500)
      for (target_x in target_areas) {
        if (!any(abs(final_points$x_orig - target_x) < target_x * 0.15)) {
          all_dists <- abs(sorted_candidates$x_orig - target_x)
          closest_idx <- which.min(all_dists)
          if (length(closest_idx) > 0 && min(all_dists) < target_x * 0.2) {
            final_points <- bind_rows(final_points, sorted_candidates[closest_idx, ])
          }
        }
      }
    }

    if (bacteria_name == "CR-Pa" && is_temperature) {
      key_points <- c(10, 20)
      for (target_x in key_points) {
        if (!any(abs(final_points$x_orig - target_x) < 5)) {
          all_dists <- abs(sorted_candidates$x_orig - target_x)
          closest_idx <- which.min(all_dists)
          if (length(closest_idx) > 0 && min(all_dists) < 10) {
            final_points <- bind_rows(final_points, sorted_candidates[closest_idx, ])
          }
        }
      }

      if (nrow(final_points) < 4) {
        remaining_candidates <- sorted_candidates %>%
          filter(!row_number() %in% match(final_points$x_orig, sorted_candidates$x_orig)) %>%
          arrange(desc(priority_score))

        if (nrow(remaining_candidates) > 0) {
          for (i in 1:min(4 - nrow(final_points), nrow(remaining_candidates))) {
            current_x <- remaining_candidates$x_orig[i]
            if (!any(abs(final_points$x_orig - current_x) < (x_range * 0.08))) {
              final_points <- bind_rows(final_points, remaining_candidates[i, ])
            }
          }
        }
      }
    }

    if (bacteria_name == "3GCR-Kp" && is_humidity) {
      final_points <- final_points %>% arrange(x_orig)

      if (nrow(final_points) >= 2) {
        distances <- diff(final_points$x_orig)
        close_idx <- which(distances < 5)
        for (i in close_idx) {
          if (final_points$effect_size[i] > final_points$effect_size[i + 1]) {
            final_points <- final_points[-(i + 1), ]
          } else {
            final_points <- final_points[-i, ]
          }
        }
      }

      key_points <- c(30, 65, 90)
      for (target_x in key_points) {
        if (!any(abs(final_points$x_orig - target_x) < 10)) {
          all_dists <- abs(sorted_candidates$x_orig - target_x)
          closest_idx <- which.min(all_dists)
          if (length(closest_idx) > 0 && min(all_dists) < 15) {
            final_points <- bind_rows(final_points, sorted_candidates[closest_idx, ])
          }
        }
      }
    }

    if (bacteria_name == "3GCR-Kp" && is_precipitation) {
      key_regions <- list(
        c(1000, 1500),
        c(2000, 2500)
      )

      for (region in key_regions) {
        if (!any(final_points$x_orig >= region[1] & final_points$x_orig <= region[2])) {
          region_candidates <- sorted_candidates %>%
            filter(x_orig >= region[1] & x_orig <= region[2])

          if (nrow(region_candidates) > 0) {
            final_points <- bind_rows(final_points, region_candidates[1, ])
          }
        }
      }
    }

    if (bacteria_name == "CR-Ab" && is_temperature) {
      key_points <- c(4, 18)
      for (target_x in key_points) {
        if (!any(abs(final_points$x_orig - target_x) < 3)) {
          all_dists <- abs(sorted_candidates$x_orig - target_x)
          closest_idx <- which.min(all_dists)
          if (length(closest_idx) > 0 && min(all_dists) < 5) {
            final_points <- bind_rows(final_points, sorted_candidates[closest_idx, ])
          }
        }
      }
    }

    return(final_points %>% select(-priority_score))
  }

  threshold_points <- find_threshold_points()

  if ((is_precipitation || is_wetdays) && nrow(threshold_points) > 0) {
    threshold_points <- threshold_points %>%
      mutate(x_orig = ifelse(x_orig <= 0, 0.01, x_orig))
  }

  if (is_humidity && nrow(threshold_points) > 0) {
    threshold_points <- threshold_points %>%
      mutate(x_orig = pmin(x_orig, 100))
  }

  if (nrow(threshold_points) > 0) {
    threshold_points <- threshold_points %>%
      mutate(
        Bacteria = bacteria_name,
        Climate_Variable = climate_variable,
        Lag_Years = lag_years,
        Relationship_Type = relationship_type,
        Value = x_orig,
        OR = y,
        Lower_CI = lower_ci_95,
        Upper_CI = upper_ci_95,
        Smooth_Parameter = span_value
      )
  }

  result <- list(
    threshold_points = threshold_points,
    curve_data = tibble(
      x_orig = orig_x,
      x_scaled = pred_x,
      y = smooth_or,
      lower_ci = smooth_lower_95,
      upper_ci = smooth_upper_95
    ),
    is_linear = is_linear,
    span_value = span_value,
    bacteria_name = bacteria_name,
    climate_variable = climate_variable,
    lag_years = lag_years,
    relationship_type = relationship_type,
    linear_info = list(
      slope = lm_slope,
      p_value = lm_p,
      r_squared = lm_r2
    )
  )

  return(result)
}

create_climate_effect_plot <- function(model, data, scale_params, var_name, x_lab, title, color, bacteria_name, threshold_data) {
  orig_var <- str_extract(var_name, "^[A-Z]+")
  curve_data <- threshold_data$curve_data
  threshold_points <- threshold_data$threshold_points
  is_linear <- threshold_data$is_linear

  lag_config <- lag_settings[[bacteria_name]]
  lag_text <- if (grepl("TMP", var_name)) {
    paste("(lag", lag_config$temp_lag, "yr)")
  } else if (grepl("PREC", var_name)) {
    paste("(lag", lag_config$precip_lag, "yr)")
  } else if (grepl("HUM", var_name)) {
    paste("(lag", lag_config$humid_lag, "yr)")
  } else if (grepl("WET", var_name)) {
    paste("(lag", lag_config$wetdays_lag, "yr)")
  } else {
    ""
  }

  x_lab_with_lag <- paste(x_lab, lag_text)

  is_precipitation <- grepl("PREC", var_name)
  is_humidity      <- grepl("HUM", var_name)
  is_temperature   <- grepl("TMP", var_name)
  is_wetdays       <- grepl("WET", var_name)

  get_range_setting <- function(var_type, bacteria) {
    setting_key <- ifelse(var_type == "PREC", "precip",
                          ifelse(var_type == "HUM", "humid",
                                 ifelse(var_type == "WET", "wetdays", "temp")))
    if (bacteria %in% names(range_settings[[setting_key]])) {
      return(range_settings[[setting_key]][[bacteria]])
    } else {
      return(range_settings[[setting_key]][["default"]])
    }
  }

  var_range <- get_range_setting(orig_var, bacteria_name)
  var_min <- var_range[1]
  var_max <- var_range[2]
  step    <- var_range[3]

  x_range <- var_max - var_min
  margin_size <- x_range * 0.04
  var_min <- var_min - margin_size
  var_max <- var_max + margin_size

  if (is_humidity) {
    curve_data <- curve_data %>%
      filter(x_orig <= 100) %>%
      mutate(x_orig = pmin(x_orig, 100))
    var_max <- min(104, var_max)
  }

  if (is_precipitation || is_wetdays) {
    var_min <- max(-margin_size, var_min)

    if (is_precipitation && bacteria_name %in% c("3GCR-Ec", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")) {
      var_max <- 3200 + margin_size
      curve_data <- curve_data %>% filter(x_orig <= 3200)
    } else if (is_wetdays) {
      var_max <- 300 + margin_size
      curve_data <- curve_data %>% filter(x_orig <= 300)
    }
  }

  x_breaks <- if (is_precipitation) {
    if (bacteria_name == "3GCR-Kp") {
      seq(0, 3200, by = 800)
    } else if (bacteria_name %in% c("3GCR-Ec", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")) {
      seq(0, 3200, by = 800)
    } else {
      seq(0, 3000, by = 500)
    }
  } else if (is_humidity) {
    seq(30, 100, by = 10)
  } else if (is_wetdays) {
    seq(0, 300, by = 50)
  } else {
    seq(round(var_min / 10) * 10, round(var_max / 10) * 10, by = step)
  }

  if (is_precipitation || is_wetdays) {
    pos_data <- data[[paste0(orig_var, "_orig")]]
    pos_data <- pos_data[pos_data > 0]

    if (length(pos_data) > 10) {
      dens <- density(pos_data, na.rm = TRUE, adjust = 1.1, from = 0, to = var_max)
      density_data <- tibble(x = dens$x, density = dens$y)
    } else {
      breaks <- seq(0, var_max, length.out = 30)
      hist_data <- hist(pos_data, breaks = breaks, plot = FALSE)
      density_data <- tibble(x = hist_data$mids, density = hist_data$density)
    }

    density_data <- density_data %>%
      filter(x >= var_min & x <= var_max)
  } else {
    if (length(data[[paste0(orig_var, "_orig")]]) > 10) {
      dens <- density(data[[paste0(orig_var, "_orig")]], na.rm = TRUE, adjust = 1.1)
      density_data <- tibble(x = dens$x, density = dens$y)

      density_data <- density_data %>%
        filter(x >= var_min & x <= var_max)

      if (is_humidity) {
        density_data <- density_data %>% filter(x <= 100)
      }
    } else {
      breaks <- seq(var_min, var_max, length.out = 30)
      hist_data <- hist(data[[paste0(orig_var, "_orig")]], breaks = breaks, plot = FALSE)
      density_data <- tibble(x = hist_data$mids, density = hist_data$density)
    }
  }

  y_min <- floor(min(c(curve_data$lower_ci, 0.5)) * 4) / 4
  y_max <- ceiling(max(c(curve_data$upper_ci, 1.5)) * 4) / 4

  if (is_humidity && bacteria_name == "CR-Ab" && max(curve_data$upper_ci) > 5) {
    y_max <- ceiling(max(curve_data$upper_ci))
  }
  if (is_wetdays && bacteria_name == "CR-Kp" && max(curve_data$upper_ci) > 5) {
    y_max <- ceiling(max(curve_data$upper_ci))
  }

  y_breaks <- seq(y_min, y_max, length.out = 5)

  curve_data_filtered <- curve_data %>%
    filter(x_orig >= var_min, x_orig <= var_max)

  if (nrow(density_data) > 0) {
    max_density <- max(density_data$density, na.rm = TRUE)
    if (max_density > 0) {
      density_data <- density_data %>%
        mutate(scaled_density = density / max_density * 0.46)
    } else {
      density_data$scaled_density <- 0
    }
  }

  main_plot <- ggplot(curve_data_filtered, aes(x = x_orig, y = y)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.7) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = color, alpha = 0.2) +
    geom_line(color = color, linewidth = 0.5) +
    scale_x_continuous(
      breaks = x_breaks,
      limits = c(var_min, var_max),
      expand = c(0.02, 0.02),
      labels = NULL
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = c(y_min, y_max),
      labels = function(x) format(x, nsmall = 3)
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 9, hjust = 0.5),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 6, r = 6, b = 0, l = 6)
    )

  if (is_temperature) {
    main_plot <- main_plot + labs(y = paste0(bacteria_name, " OR (95% CI)"))
  } else {
    main_plot <- main_plot + labs(y = "")
  }

  if (nrow(threshold_points) > 0 && !is_linear) {
    for (i in 1:nrow(threshold_points)) {
      main_plot <- main_plot +
        geom_vline(xintercept = threshold_points$x_orig[i], color = "#FF0000",
                   linetype = "dashed", linewidth = 0.4)
    }

    threshold_points <- threshold_points %>%
      mutate(
        label_short = case_when(
          grepl("Global maximum", type) ~ "GMax",
          grepl("Global minimum", type) ~ "GMin",
          grepl("Maximum", type) ~ "Max",
          grepl("Minimum", type) ~ "Min",
          TRUE ~ substr(type, 1, 3)
        ),
        unit_suffix = case_when(
          is_temperature ~ "°C",
          is_humidity ~ "%",
          is_precipitation ~ "mm",
          is_wetdays ~ "d",
          TRUE ~ ""
        ),
        x_value = ifelse(abs(x_orig - round(x_orig)) < 1e-10,
                         as.character(round(x_orig)),
                         format(round(x_orig, 3), nsmall = 3)),
        or_value = ifelse(abs(y - round(y)) < 1e-10,
                          as.character(round(y)),
                          format(round(y, 3), nsmall = 3)),
        label = paste0(label_short, " (", x_value, unit_suffix, ")\nOR = ", or_value, sig_symbol)
      )

    if (nrow(threshold_points) > 0) {
      threshold_points <- threshold_points %>% arrange(x_orig)

      vjust_values <- numeric(nrow(threshold_points))
      hjust_values <- numeric(nrow(threshold_points))

      for (i in 1:nrow(threshold_points)) {
        if (threshold_points$y[i] > 1) {
          vjust_values[i] <- -0.3
        } else {
          vjust_values[i] <- 1.3
        }

        if (threshold_points$x_orig[i] < (var_min + x_range * 0.15)) {
          hjust_values[i] <- 0
        } else if (threshold_points$x_orig[i] > (var_max - x_range * 0.15)) {
          hjust_values[i] <- 1
        } else {
          hjust_values[i] <- 0.5
        }
      }

      if (nrow(threshold_points) > 1) {
        for (i in 2:nrow(threshold_points)) {
          if (abs(threshold_points$x_orig[i] - threshold_points$x_orig[i - 1]) < (x_range * 0.2)) {
            if (vjust_values[i] * vjust_values[i - 1] > 0) {
              vjust_values[i] <- -vjust_values[i - 1]
            }
          }
        }
      }

      threshold_points$vjust <- vjust_values
      threshold_points$hjust <- hjust_values

      main_plot <- main_plot +
        geom_text(
          data = threshold_points,
          aes(x = x_orig, y = y, label = label, hjust = hjust, vjust = vjust),
          size = 2.3
        )
    }
  }

  density_plot <- ggplot(density_data, aes(x = x, y = scaled_density)) +
    geom_area(fill = color, alpha = 0.6) +
    scale_x_continuous(
      breaks = x_breaks,
      limits = c(var_min, var_max),
      expand = c(0.02, 0.02)
    ) +
    labs(x = x_lab_with_lag) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.length.x = unit(0.15, "cm"),
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      panel.border = element_rect(color = "black", linewidth = 0.7),
      plot.margin = margin(t = 0, r = 6, b = 6, l = 6)
    )

  combined_plot <- plot_grid(
    main_plot,
    density_plot,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(0.77, 0.23)
  )

  return(combined_plot)
}

combine_all_bacterial_plots <- function(paths, output_dirs) {
  all_plots <- list()
  threshold_summary <- list()
  linear_results <- list()

  climate_factors <- list(
    list(var = "TMP_scaled_lag",  x_lab = "Temperature (°C)",         color = climate_colors["Temperature"]),
    list(var = "HUM_scaled_lag",  x_lab = "Relative humidity (%)",    color = climate_colors["Humidity"]),
    list(var = "PREC_scaled_lag", x_lab = "Precipitation (mm)",       color = climate_colors["Precipitation"]),
    list(var = "WET_scaled_lag",  x_lab = "Wet Days (d)",             color = climate_colors["WetDays"])
  )

  for (i in seq_along(paths)) {
    cat("Processing", paths[[i]]$title, "...\n")

    data_result <- paths[[i]]$prepare_data(paths[[i]]$file_path, paths[[i]]$title)
    processed_data <- data_result$data
    scale_params <- data_result$scale_params

    model <- build_gamm_model(processed_data, paths[[i]]$title, output_dirs)

    for (factor in climate_factors) {
      cat(" Processing", paths[[i]]$title, factor$var, "...\n")

      threshold_data <- detect_thresholds(model, processed_data, scale_params, factor$var, paths[[i]]$title)

      threshold_summary[[paste(paths[[i]]$title, factor$var, sep = "_")]] <- threshold_data$threshold_points

      linear_info <- threshold_data$linear_info
      if (is.null(linear_info)) {
        linear_info <- list(slope = NA, p_value = NA, r_squared = NA)
      }

      linear_results[[length(linear_results) + 1]] <- data.frame(
        Bacteria = paths[[i]]$title,
        Climate_Variable = threshold_data$climate_variable,
        Lag_Years = threshold_data$lag_years,
        Relationship_Type = threshold_data$relationship_type,
        Slope = linear_info$slope,
        P_Value = linear_info$p_value,
        R_Squared = linear_info$r_squared
      )

      plot <- create_climate_effect_plot(
        model, processed_data, scale_params, factor$var, factor$x_lab,
        paths[[i]]$title, factor$color, paths[[i]]$title, threshold_data
      )

      var_name_clean <- ifelse(
        factor$var == "TMP_scaled_lag", "Temperature",
        ifelse(factor$var == "HUM_scaled_lag", "Humidity",
               ifelse(factor$var == "PREC_scaled_lag", "Precipitation", "WetDays"))
      )

      lag <- get_lag_value(paths[[i]]$title, factor$var)

      plot_path <- file.path(
        output_dirs$single_factor_panels,
        paste0(
          "ModelC_GAMM_",
          sanitize_file_stem(paths[[i]]$title),
          "_",
          sanitize_file_stem(var_name_clean),
          "_lag",
          lag,
          ".pdf"
        )
      )
      png_path <- file.path(
        output_dirs$single_factor_panels,
        paste0(
          "ModelC_GAMM_",
          sanitize_file_stem(paths[[i]]$title),
          "_",
          sanitize_file_stem(var_name_clean),
          "_lag",
          lag,
          ".png"
        )
      )

      ggsave(plot_path, plot, width = 5, height = 4, dpi = 300)
      ggsave(png_path,  plot, width = 5, height = 4, dpi = 300)

      cat(" Saved plot to:", plot_path, "\n")

      all_plots[[paste(paths[[i]]$title, var_name_clean, sep = "_")]] <- plot
    }

    combined_bacteria_plot <- plot_grid(
      all_plots[[paste(paths[[i]]$title, "Temperature", sep = "_")]],
      all_plots[[paste(paths[[i]]$title, "Humidity", sep = "_")]],
      all_plots[[paste(paths[[i]]$title, "Precipitation", sep = "_")]],
      all_plots[[paste(paths[[i]]$title, "WetDays", sep = "_")]],
      ncol = 2,
      nrow = 2,
      align = "hv",
      labels = c("A", "B", "C", "D"),
      label_size = 12
    )

    combined_path <- file.path(
      output_dirs$combined_panels,
      paste0("ModelC_GAMM_", sanitize_file_stem(paths[[i]]$title), "_four_climate_panels.pdf")
    )
    combined_png_path <- file.path(
      output_dirs$combined_panels,
      paste0("ModelC_GAMM_", sanitize_file_stem(paths[[i]]$title), "_four_climate_panels.png")
    )
    ggsave(combined_path, combined_bacteria_plot, width = 10, height = 8, dpi = 300)
    ggsave(combined_png_path, combined_bacteria_plot, width = 10, height = 8, dpi = 300)
    cat("Saved combined plot for", paths[[i]]$title, "to:", combined_path, "\n")
  }

  title_row <- plot_grid(
    ggplot() + theme_void() + ggtitle("Temperature")   + theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ggplot() + theme_void() + ggtitle("Humidity")      + theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ggplot() + theme_void() + ggtitle("Precipitation") + theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ggplot() + theme_void() + ggtitle("Wet Days")      + theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ncol = 4,
    rel_heights = c(0.2)
  )

  plot_list <- list(title_row)
  bacteria_names <- sapply(paths, function(x) x$title)

  for (i in seq_along(bacteria_names)) {
    bacteria <- bacteria_names[i]
    row_plots <- plot_grid(
      all_plots[[paste(bacteria, "Temperature", sep = "_")]],
      all_plots[[paste(bacteria, "Humidity", sep = "_")]],
      all_plots[[paste(bacteria, "Precipitation", sep = "_")]],
      all_plots[[paste(bacteria, "WetDays", sep = "_")]],
      ncol = 4,
      align = "h"
    )
    plot_list[[i + 1]] <- row_plots
  }

  final_plot <- plot_grid(
    plotlist = plot_list,
    ncol = 1,
    rel_heights = c(0.08, rep(1, length(bacteria_names)))
  )

  final_path <- file.path(
    output_dirs$combined_panels,
    "ModelC_GAMM_MainFigure1_all_bacteria_four_climate_panels.pdf"
  )
  ggsave(final_path, final_plot, width = 16, height = 16, dpi = 300)

  png_path <- file.path(
    output_dirs$combined_panels,
    "ModelC_GAMM_MainFigure1_all_bacteria_four_climate_panels.png"
  )
  ggsave(png_path, final_plot, width = 16, height = 16, dpi = 300)

  cat("Saved final combined plot to:", final_path, "\n")

  threshold_summary_df <- bind_rows(threshold_summary, .id = "bacteria_variable")
  threshold_summary_path <- file.path(
    output_dirs$threshold_tables,
    "ModelC_GAMM_threshold_summary_four_climate.csv"
  )
  write.csv(threshold_summary_df, threshold_summary_path, row.names = FALSE)
  cat("Saved threshold summary to:", threshold_summary_path, "\n")

  linear_summary_df <- bind_rows(linear_results)
  linear_summary_path <- file.path(
    output_dirs$threshold_tables,
    "ModelC_GAMM_linear_relationships_summary_four_climate.csv"
  )
  write.csv(linear_summary_df, linear_summary_path, row.names = FALSE)
  cat("Saved linear relationships summary to:", linear_summary_path, "\n")

  return(list(
    plots = all_plots,
    thresholds = threshold_summary,
    linear_relationships = linear_results
  ))
}

main <- function() {
  output_root <- file.path(
    revision_root,
    "outputs/historical_associations",
    paste0("historical_associations_figure1_results", output_suffix)
  )

  output_dirs <- list(
    root = output_root,
    model_summaries = file.path(output_root, "01_model_summaries"),
    single_factor_panels = file.path(output_root, "02_single_factor_panels"),
    combined_panels = file.path(output_root, "03_combined_panels"),
    threshold_tables = file.path(output_root, "04_threshold_and_relationship_tables")
  )

  for (dir_path in output_dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }

  paths <- lapply(bacteria_specs, function(spec) {
    list(
      code = spec$code,
      title = spec$title,
      file_path = file.path(input_data_dir, spec$file_name),
      prepare_data = prepare_data
    )
  })

  missing_inputs <- vapply(paths, function(x) !file.exists(x$file_path), logical(1))
  if (any(missing_inputs)) {
    stop(
      "Missing model input files: ",
      paste(vapply(paths[missing_inputs], `[[`, character(1), "file_path"), collapse = ", "),
      call. = FALSE
    )
  }

  results <- combine_all_bacterial_plots(paths, output_dirs)
  return(results)
}

if (Sys.getenv("MODELC_SKIP_MAIN", unset = "0") != "1") {
  invisible(main())
}
