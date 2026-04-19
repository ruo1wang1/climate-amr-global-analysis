cat("\n================================================================\n")
cat("FOUR-SPECIFICATION COVARIATE ROBUSTNESS -- historical climate inputs\n")
cat("================================================================\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(writexl)
})

# ============================================================================
# 1. SETTINGS
# ============================================================================
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

output_root <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "10_CovariateRobustness_4Specs_historical"
)

source_data_root <- file.path(
  revision_root,
  "data/source_data",
  "covariate_adjustment_robustness"
)

dir.create(file.path(output_root, "01_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "02_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "03_workbook"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "04_metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_data_root, "01_csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_data_root, "02_workbook"), recursive = TRUE, showWarnings = FALSE)

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = "3GCREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = "3GCRKP_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRAB", title = "CR-Ab",   file_name = "CRAB_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CREC", title = "CR-Ec",   file_name = "CREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRKP", title = "CR-Kp",   file_name = "CRKP_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRPA", title = "CR-Pa",   file_name = "CRPA_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv")
)

bacteria_levels <- vapply(bacteria_specs, `[[`, character(1), "title")
climate_levels <- c("Temperature", "Humidity", "Precipitation", "Wet days")

all_14_covariates <- c(
  "CPI", "OPE", "HEG", "HD", "IC", "UP",
  "log_PD", "log_IT", "log_AMC", "log_PMP",
  "SR", "SSR", "TCL", "WASH"
)

theory_6 <- c("log_AMC", "CPI", "HEG", "UP", "log_PMP", "WASH")

top6_per_bact <- list(
  "3GCR-Ec" = c("CPI", "OPE", "UP", "log_PMP", "HEG", "HD"),
  "3GCR-Kp" = c("CPI", "UP", "HEG", "log_PMP", "OPE", "HD"),
  "CR-Ab"   = c("CPI", "HEG", "log_PMP", "UP", "OPE", "WASH"),
  "CR-Ec"   = c("CPI", "OPE", "UP", "HEG", "log_PMP", "WASH"),
  "CR-Kp"   = c("CPI", "OPE", "UP", "HEG", "log_PMP", "WASH"),
  "CR-Pa"   = c("CPI", "UP", "log_PMP", "HEG", "OPE", "WASH")
)

spec_colors <- c(
  "PLS (main)" = "#4A708B",
  "PCA" = "#2E8B57",
  "Theory-6" = "#CD853F",
  "Top6 from PLS1" = "#8B668B"
)

spec_display <- c(
  "PLS (main)" = "PLS (Model C)",
  "PCA" = "PCA",
  "Theory-6" = "Theory-selected covariates",
  "Top6 from PLS1" = "PLS-informed top covariates"
)

spec_colors_display <- setNames(unname(spec_colors), unname(spec_display[names(spec_colors)]))

variable_display <- c(
  "Temperature" = "Temperature (\u00B0C)",
  "Humidity" = "Relative Humidity (%)",
  "Precipitation" = "Precipitation (mm)",
  "Wet days" = "Wet days (d)"
)

range_settings <- list(
  temp = c(-10, 40, 10),
  humid = c(30, 100, 10),
  precip = c(0, 3200, 500),
  wetdays = c(0, 300, 50)
)

ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)

# ============================================================================
# 2. HELPERS
# ============================================================================
load_lag_settings <- function(summary_csv) {
  if (!file.exists(summary_csv)) {
    stop("Lag summary file not found: ", summary_csv, call. = FALSE)
  }
  lag_df <- read.csv(summary_csv, stringsAsFactors = FALSE)
  required_cols <- c("Display_Name", "TMP_lag", "PREC_lag", "HUM_lag", "WET_lag")
  missing_cols <- setdiff(required_cols, names(lag_df))
  if (length(missing_cols) > 0) {
    stop("Lag summary missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  out <- setNames(vector("list", nrow(lag_df)), lag_df$Display_Name)
  for (i in seq_len(nrow(lag_df))) {
    out[[lag_df$Display_Name[i]]] <- list(
      temp_lag = as.integer(lag_df$TMP_lag[i]),
      precip_lag = as.integer(lag_df$PREC_lag[i]),
      humid_lag = as.integer(lag_df$HUM_lag[i]),
      wetdays_lag = as.integer(lag_df$WET_lag[i])
    )
  }
  out
}

lag_settings <- load_lag_settings(lag_summary_path)

sig_code <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.10 ~ "\u2020",
    TRUE ~ "ns"
  )
}

safe_first <- function(df, col) {
  if (nrow(df) == 0) return(NA)
  df[[col]][1]
}

get_available_pls_components <- function(data) {
  pls_candidates <- paste0("PLS_Comp", 1:4)
  present <- pls_candidates[pls_candidates %in% names(data)]
  present[sapply(present, function(x) !all(is.na(data[[x]])))]
}

identify_available_covariates <- function(data) {
  present <- intersect(all_14_covariates, names(data))
  keep <- present[sapply(present, function(v) {
    x <- data[[v]]
    all(is.na(x)) == FALSE && sd(x, na.rm = TRUE) > 0
  })]
  keep
}

impute_scale_vec <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  sx <- sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

get_range_setting <- function(var_name) {
  if (grepl("^TMP", var_name)) return(range_settings$temp)
  if (grepl("^HUM", var_name)) return(range_settings$humid)
  if (grepl("^PREC", var_name)) return(range_settings$precip)
  if (grepl("^WET", var_name)) return(range_settings$wetdays)
  range_settings$temp
}

get_scale_key <- function(var_name) {
  if (grepl("^TMP", var_name)) return("TMP")
  if (grepl("^HUM", var_name)) return("HUM")
  if (grepl("^PREC", var_name)) return("PREC")
  if (grepl("^WET", var_name)) return("WET")
  "TMP"
}

compute_vif_table <- function(data, covariates_sc, bacteria_name, spec_label) {
  covariates_sc <- covariates_sc[covariates_sc %in% names(data)]
  if (length(covariates_sc) == 0) {
    return(tibble())
  }
  if (length(covariates_sc) == 1) {
    return(tibble(
      Bacteria = bacteria_name,
      Spec = spec_label,
      Covariate_sc = covariates_sc,
      VIF = 1
    ))
  }

  out <- lapply(covariates_sc, function(v) {
    others <- setdiff(covariates_sc, v)
    fit <- lm(reformulate(others, response = v), data = data)
    r2 <- summary(fit)$r.squared
    vif_val <- ifelse(is.na(r2) || r2 >= 1, NA_real_, 1 / (1 - r2))
    tibble(Bacteria = bacteria_name, Spec = spec_label, Covariate_sc = v, VIF = vif_val)
  })
  bind_rows(out)
}

extract_climate_terms <- function(model, spec_label) {
  st <- as.data.frame(summary(model)$s.table)
  st$term <- rownames(st)
  st %>%
    filter(grepl("TMP_scaled_lag|PREC_scaled_lag|HUM_scaled_lag|WET_scaled_lag", term)) %>%
    mutate(
      Spec = spec_label,
      Variable = case_when(
        grepl("TMP", term) ~ "Temperature",
        grepl("PREC", term) ~ "Precipitation",
        grepl("HUM", term) ~ "Humidity",
        grepl("WET", term) ~ "Wet days"
      )
    ) %>%
    rename(EDF = edf, F_stat = F, p_value = `p-value`) %>%
    select(Spec, Variable, EDF, F_stat, p_value, term)
}

get_curve <- function(model, data, var_name, spec_label, n = 300) {
  base <- data[rep(1, n), , drop = FALSE]
  for (col in names(data)) {
    if (is.numeric(data[[col]])) {
      base[[col]] <- mean(data[[col]], na.rm = TRUE)
    } else if (is.factor(data[[col]])) {
      base[[col]] <- factor(names(which.max(table(data[[col]]))), levels = levels(data[[col]]))
    }
  }

  scale_params <- attr(data, "scale_params")
  raw_key <- get_scale_key(var_name)
  raw_mean <- scale_params[[paste0(raw_key, "_mean")]]
  raw_sd <- scale_params[[paste0(raw_key, "_sd")]]
  var_range <- get_range_setting(var_name)
  x_orig <- seq(var_range[1], var_range[2], length.out = n)
  x_scaled <- (x_orig - raw_mean) / raw_sd
  base[[var_name]] <- x_scaled

  tryCatch({
    pred <- predict(model, base, type = "terms", se.fit = TRUE)
    idx <- grep(var_name, colnames(pred$fit), fixed = TRUE)
    if (length(idx) == 0) return(NULL)

    fit <- pred$fit[, idx[1]]
    se <- pred$se.fit[, idx[1]]
    tibble(
      x = x_scaled,
      x_orig = x_orig,
      fit = fit,
      se = se,
      or = exp(fit),
      lo = exp(fit - 1.96 * se),
      hi = exp(fit + 1.96 * se),
      Spec = spec_label
    )
  }, error = function(e) NULL)
}

compute_or_range <- function(model, data, var_name) {
  curve <- get_curve(model, data, var_name, "tmp")
  if (is.null(curve) || nrow(curve) == 0) return(NA_real_)

  q10 <- quantile(curve$x_orig, 0.10, na.rm = TRUE)
  q90 <- quantile(curve$x_orig, 0.90, na.rm = TRUE)
  sub <- curve %>% filter(x_orig >= q10, x_orig <= q90)
  if (nrow(sub) < 10) return(NA_real_)
  round(max(sub$or, na.rm = TRUE) / min(sub$or, na.rm = TRUE), 3)
}

compute_shape_similarity <- function(curve_main, curve_alt, bacteria_name, variable, alt_spec) {
  if (is.null(curve_main) || is.null(curve_alt) || nrow(curve_main) < 20 || nrow(curve_alt) < 20) {
    return(NULL)
  }

  xmin <- max(min(curve_main$x_orig, na.rm = TRUE), min(curve_alt$x_orig, na.rm = TRUE))
  xmax <- min(max(curve_main$x_orig, na.rm = TRUE), max(curve_alt$x_orig, na.rm = TRUE))
  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin) return(NULL)

  grid <- seq(xmin, xmax, length.out = 250)
  fit_main <- approx(curve_main$x_orig, curve_main$fit, xout = grid, rule = 2)$y
  fit_alt <- approx(curve_alt$x_orig, curve_alt$fit, xout = grid, rule = 2)$y
  rho <- suppressWarnings(cor(fit_main, fit_alt, use = "complete.obs"))
  slope_main <- coef(lm(fit_main ~ grid))[2]
  slope_alt <- coef(lm(fit_alt ~ grid))[2]

  tibble(
    Bacteria = bacteria_name,
    Variable = variable,
    Alt_Spec = alt_spec,
    rho = round(rho, 3),
    Main_direction = ifelse(slope_main > 0, "+", ifelse(slope_main < 0, "-", "~0")),
    Alt_direction = ifelse(slope_alt > 0, "+", ifelse(slope_alt < 0, "-", "~0")),
    Direction_match = Main_direction == Alt_direction
  )
}

fit_model <- function(formula_str, data, label) {
  cat("    Fitting:", label, "...")
  m <- tryCatch(
    bam(as.formula(formula_str), data = data, family = gaussian(), method = "REML", select = TRUE, control = ctrl),
    error = function(e) {
      warning(paste(label, "bam failed, using gam:", conditionMessage(e)))
      gam(as.formula(formula_str), data = data, family = gaussian(), method = "REML", select = TRUE)
    }
  )
  cat(" R2=", round(summary(m)$r.sq, 4), "\n")
  m
}

# ============================================================================
# 3. DATA PREPARATION
# ============================================================================
prepare_data_robust <- function(file_path, bacteria_name) {
  lag_config <- lag_settings[[bacteria_name]]
  if (is.null(lag_config)) {
    stop("Lag settings not found for ", bacteria_name, call. = FALSE)
  }
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path, call. = FALSE)
  }

  data <- read.csv(file_path, stringsAsFactors = FALSE) %>%
    mutate(
      year = as.numeric(as.character(year)),
      Region = factor(Region),
      NAME = factor(NAME),
      climate_zone = case_when(
        abs(lat) > 66.5 ~ "Polar Zone",
        abs(lat) > 23.5 ~ "Temperate Zone",
        TRUE ~ "Tropical Zone"
      ),
      climate_zone = factor(climate_zone, levels = c("Polar Zone", "Temperate Zone", "Tropical Zone")),
      HUM = pmin(HUM, 100)
    ) %>%
    group_by(NAME) %>%
    mutate(location_id = cur_group_id()) %>%
    ungroup()

  scale_params <- data %>%
    summarise(
      across(c(TMP, PREC, HUM, WET), list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)))
    )

  data <- data %>%
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
    ungroup() %>%
    group_by(location_id) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      TMP_scaled_lag  = lag(TMP_scaled, lag_config$temp_lag),
      PREC_scaled_lag = lag(PREC_scaled, lag_config$precip_lag),
      HUM_scaled_lag  = lag(HUM_scaled, lag_config$humid_lag),
      WET_scaled_lag  = lag(WET_scaled, lag_config$wetdays_lag)
    ) %>%
    ungroup() %>%
    filter(
      !is.na(TMP_scaled_lag),
      !is.na(PREC_scaled_lag),
      !is.na(HUM_scaled_lag),
      !is.na(WET_scaled_lag),
      !is.na(logit_R)
    )

  avail_covs <- identify_available_covariates(data)
  for (v in avail_covs) {
    data[[paste0(v, "_sc")]] <- impute_scale_vec(data[[v]])
  }

  cov_matrix <- data %>%
    select(all_of(paste0(avail_covs, "_sc"))) %>%
    as.data.frame()
  pca_result <- prcomp(cov_matrix, center = FALSE, scale. = FALSE)
  cumvar <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
  n_pc <- min(which(cumvar >= 0.80))
  for (j in seq_len(n_pc)) {
    data[[paste0("PC", j)]] <- pca_result$x[, j]
  }

  top6_vec <- intersect(top6_per_bact[[bacteria_name]], avail_covs)
  top6_df <- tibble(
    Covariate = top6_vec,
    Covariate_sc = paste0(top6_vec, "_sc"),
    Source = "original covariate set definition"
  )

  attr(data, "avail_covs") <- avail_covs
  attr(data, "n_pc") <- n_pc
  attr(data, "pca_cumvar") <- cumvar[n_pc]
  attr(data, "top6_df") <- top6_df
  attr(data, "scale_params") <- scale_params
  data
}

# ============================================================================
# 4. FOUR SPECIFICATIONS
# ============================================================================
climate_block <- paste(
  "s(TMP_scaled_lag, k = 5, bs = 'cr')",
  "s(PREC_scaled_lag, k = 10, bs = 'cr')",
  "s(HUM_scaled_lag, k = 10, bs = 'cr')",
  "s(WET_scaled_lag, k = 10, bs = 'cr')",
  sep = " + "
)

spatial_block <- "s(lat, lon, bs = 'sos', k = 20) + s(year, bs = 'cr', k = 8) + s(Region, bs = 're') + climate_zone"

build_specA <- function(data, bact) {
  pls <- get_available_pls_components(data)
  if (length(pls) == 0) {
    stop("No PLS components available for ", bact, call. = FALSE)
  }
  pls_terms <- paste0("s(", pls, ", k = 10, bs = 'cr')", collapse = " + ")
  fit_model(
    paste("logit_R ~", climate_block, "+", pls_terms, "+", spatial_block),
    data,
    paste("Spec A (PLS) -", bact)
  )
}

build_specB <- function(data, bact) {
  n_pc <- attr(data, "n_pc")
  pc_terms <- paste0("s(PC", seq_len(n_pc), ", k = 5, bs = 'cr')", collapse = " + ")
  fit_model(
    paste("logit_R ~", climate_block, "+", pc_terms, "+", spatial_block),
    data,
    paste("Spec B (PCA) -", bact)
  )
}

build_specC <- function(data, bact) {
  avail <- intersect(theory_6, attr(data, "avail_covs"))
  lin_terms <- paste(paste0(avail, "_sc"), collapse = " + ")
  fit_model(
    paste("logit_R ~", climate_block, "+", lin_terms, "+", spatial_block),
    data,
    paste("Spec C (Theory-6) -", bact)
  )
}

build_specD <- function(data, bact) {
  top6_df <- attr(data, "top6_df")
  avail <- top6_df$Covariate
  lin_terms <- paste(paste0(avail, "_sc"), collapse = " + ")
  fit_model(
    paste("logit_R ~", climate_block, "+", lin_terms, "+", spatial_block),
    data,
    paste("Spec D (Top6 from PLS1) -", bact)
  )
}

# ============================================================================
# 5. RUN ONE BACTERIUM
# ============================================================================
run_bacterium <- function(file_path, bact) {
  cat("\n", strrep("=", 60), "\n", bact, "\n", strrep("=", 60), "\n", sep = "")
  data <- prepare_data_robust(file_path, bact)

  mA <- build_specA(data, bact)
  mB <- build_specB(data, bact)
  mC <- build_specC(data, bact)
  mD <- build_specD(data, bact)

  models <- list(
    "PLS (main)" = mA,
    "PCA" = mB,
    "Theory-6" = mC,
    "Top6 from PLS1" = mD
  )

  vars <- c("TMP_scaled_lag", "PREC_scaled_lag", "HUM_scaled_lag", "WET_scaled_lag")
  vlabs <- c("Temperature", "Precipitation", "Humidity", "Wet days")

  all_stats <- bind_rows(lapply(names(models), function(spec_label) {
    extract_climate_terms(models[[spec_label]], spec_label)
  }))

  or_ranges <- expand.grid(Spec = names(models), Var = vars, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      OR_range = compute_or_range(models[[Spec]], data, Var),
      Variable = vlabs[match(Var, vars)]
    ) %>%
    ungroup() %>%
    select(Spec, Variable, OR_range)

  curve_store <- list()
  shape_store <- list()
  for (i in seq_along(vars)) {
    cv <- bind_rows(lapply(names(models), function(spec_label) {
      out <- get_curve(models[[spec_label]], data, vars[i], spec_label)
      if (is.null(out)) return(NULL)
      out %>% mutate(
        Bacteria = bact,
        Variable = vlabs[i],
        Variable_label = recode(vlabs[i], !!!variable_display)
      )
    }))
    curve_store[[vlabs[i]]] <- cv

    main_curve <- cv %>% filter(Spec == "PLS (main)")
    for (alt in c("PCA", "Theory-6", "Top6 from PLS1")) {
      alt_curve <- cv %>% filter(Spec == alt)
      sm <- compute_shape_similarity(main_curve, alt_curve, bact, vlabs[i], alt)
      if (!is.null(sm)) {
        shape_store[[length(shape_store) + 1L]] <- sm
      }
    }
  }

  top6_df <- attr(data, "top6_df") %>%
    mutate(Bacteria = bact, Rank = row_number())

  theory_avail <- intersect(theory_6, attr(data, "avail_covs"))
  vif_df <- bind_rows(
    compute_vif_table(data, paste0(theory_avail, "_sc"), bact, "Theory-6"),
    compute_vif_table(data, paste0(top6_df$Covariate, "_sc"), bact, "Top6 from PLS1")
  )

  model_info <- tibble(
    Bacteria = bact,
    N = nrow(data),
    Countries = dplyr::n_distinct(data$NAME),
    PLS_components_used = paste(get_available_pls_components(data), collapse = ", "),
    PCA_components_used = attr(data, "n_pc"),
    PCA_cumulative_variance = round(100 * attr(data, "pca_cumvar"), 1),
    Theory6_used = paste(theory_avail, collapse = ", "),
    Top6_used = paste(top6_df$Covariate, collapse = ", "),
    TMP_lag = lag_settings[[bact]]$temp_lag,
    PREC_lag = lag_settings[[bact]]$precip_lag,
    HUM_lag = lag_settings[[bact]]$humid_lag,
    WET_lag = lag_settings[[bact]]$wetdays_lag,
    A_Rsq = round(summary(mA)$r.sq, 4),
    A_Dev = round(100 * summary(mA)$dev.expl, 2),
    B_Rsq = round(summary(mB)$r.sq, 4),
    B_Dev = round(100 * summary(mB)$dev.expl, 2),
    C_Rsq = round(summary(mC)$r.sq, 4),
    C_Dev = round(100 * summary(mC)$dev.expl, 2),
    D_Rsq = round(summary(mD)$r.sq, 4),
    D_Dev = round(100 * summary(mD)$dev.expl, 2)
  )

  comp <- all_stats %>%
    left_join(or_ranges, by = c("Spec", "Variable")) %>%
    mutate(
      Bacteria = bact,
      Sig = sig_code(p_value),
      Retained = ifelse(EDF >= 0.5, "Yes", "No"),
      EDF = round(EDF, 2),
      F_stat = round(F_stat, 3),
      OR_range = round(OR_range, 3)
    ) %>%
    select(Bacteria, Spec, Variable, term, EDF, F_stat, p_value, Sig, Retained, OR_range)

  table_b <- comp %>%
    transmute(
      Bacteria,
      Variable,
      Spec,
      Cell = paste0("EDF=", EDF, " ", Sig, " [", OR_range, "]")
    ) %>%
    pivot_wider(names_from = Spec, values_from = Cell)

  table_c <- comp %>%
    transmute(
      Bacteria,
      Variable,
      Spec,
      Cell = paste0(Retained, " (", Sig, ")")
    ) %>%
    pivot_wider(names_from = Spec, values_from = Cell)

  list(
    comparison = comp,
    model_info = model_info,
    curves = curve_store,
    top6 = top6_df,
    vif = vif_df,
    shape = bind_rows(shape_store),
    table_b = table_b,
    table_c = table_c
  )
}

# ============================================================================
# 6. VISUALIZATION
# ============================================================================
plot_overlay <- function(curve_list, output_path) {
  df <- bind_rows(lapply(curve_list, function(x) bind_rows(x)))
  if (nrow(df) == 0) return(NULL)

  df$Bacteria <- factor(df$Bacteria, levels = bacteria_levels)
  df$Variable_label <- factor(df$Variable_label, levels = variable_display[climate_levels])
  df$Spec_display <- factor(recode(df$Spec, !!!spec_display), levels = unname(spec_display))

  p <- ggplot(df, aes(x = x_orig, y = or, color = Spec_display)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.6, alpha = 0.9) +
    facet_grid(rows = vars(Bacteria), cols = vars(Variable_label), scales = "free_x") +
    scale_color_manual(values = spec_colors_display) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = "Odds Ratio",
      color = "Covariate specification"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.3),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
      strip.text.x = element_text(face = "bold", size = 9),
      strip.text.y = element_text(face = "bold", size = 9, angle = 270),
      strip.placement = "outside",
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 8, colour = "black"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(file.path(output_path, "covariate_adjustment_curve_overlay.pdf"), p, width = 14, height = 16, dpi = 300, bg = "white")
  ggsave(file.path(output_path, "covariate_adjustment_curve_overlay.png"), p, width = 14, height = 16, dpi = 300, bg = "white")
  p
}

plot_edf_heat <- function(full_comp, output_path) {
  heat <- full_comp %>%
    mutate(
      Bacteria = factor(Bacteria, levels = rev(bacteria_levels)),
      Variable = factor(Variable, levels = climate_levels),
      Spec_display = factor(recode(Spec, !!!spec_display), levels = unname(spec_display)),
      Label = paste0(EDF, Sig),
      Fill_val = pmin(EDF, 8)
    )

  p <- ggplot(heat, aes(x = Variable, y = Bacteria, fill = Fill_val)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = Label), size = 2.5) +
    facet_wrap(~Spec_display, nrow = 1) +
    scale_fill_gradient2(low = "#F4B4B4", mid = "#FFFACD", high = "#A8D5BA", midpoint = 1.5, limits = c(0, 8), name = "EDF") +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(file.path(output_path, "covariate_adjustment_edf_heatmap.pdf"), p, width = 16, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(output_path, "covariate_adjustment_edf_heatmap.png"), p, width = 16, height = 6, dpi = 300, bg = "white")
  p
}

plot_or_range <- function(full_comp, output_path) {
  or_data <- full_comp %>%
    filter(!is.na(OR_range)) %>%
    mutate(
      Bacteria = factor(Bacteria, levels = bacteria_levels),
      Variable = factor(Variable, levels = climate_levels),
      Spec_display = factor(recode(Spec, !!!spec_display), levels = unname(spec_display))
    )

  p <- ggplot(or_data, aes(x = Bacteria, y = OR_range, fill = Spec_display)) +
    geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
    facet_wrap(~Variable, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = spec_colors_display) +
    labs(
      title = NULL,
      subtitle = NULL,
      y = "OR range (max/min)",
      x = NULL,
      fill = "Covariate specification"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(file.path(output_path, "covariate_adjustment_or_range.pdf"), p, width = 14, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(output_path, "covariate_adjustment_or_range.png"), p, width = 14, height = 5.5, dpi = 300, bg = "white")
  p
}

plot_shape_heat <- function(shape_df, output_path) {
  if (nrow(shape_df) == 0) return(NULL)
  heat <- shape_df %>%
    mutate(
      Bacteria = factor(Bacteria, levels = bacteria_levels),
      Variable = factor(Variable, levels = climate_levels),
      Alt_Spec_display = factor(recode(Alt_Spec, !!!spec_display), levels = unname(spec_display[c("PCA", "Theory-6", "Top6 from PLS1")])),
      Label = paste0("rho=", sprintf("%.2f", rho))
    )

  p <- ggplot(heat, aes(x = Variable, y = Bacteria, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = Label), size = 2.6) +
    facet_wrap(~Alt_Spec_display, nrow = 1) +
    scale_fill_gradient2(low = "#F4B4B4", mid = "#FFFACD", high = "#A8D5BA", midpoint = 0.8, limits = c(-1, 1), name = "rho") +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  ggsave(file.path(output_path, "covariate_adjustment_shape_concordance.pdf"), p, width = 14, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(output_path, "covariate_adjustment_shape_concordance.png"), p, width = 14, height = 6, dpi = 300, bg = "white")
  p
}

create_docx_tables <- function(table_a, table_b, table_c, output_path) {
  if (!requireNamespace("flextable", quietly = TRUE) || !requireNamespace("officer", quietly = TRUE)) {
    return(invisible(NULL))
  }
  library(flextable)
  library(officer)

  ft_a <- flextable(table_a) %>% theme_booktabs() %>% autofit()
  ft_b <- flextable(table_b) %>% theme_booktabs() %>% autofit()
  ft_c <- flextable(table_c) %>% theme_booktabs() %>% autofit()

  doc <- read_docx() %>%
    body_add_par("covariate-adjustment robustness summary. Sensitivity of climate-AMR associations to covariate specification", style = "heading 1") %>%
    body_add_par("Table A. Model-level goodness-of-fit across four covariate adjustment specifications", style = "heading 2") %>%
    body_add_flextable(ft_a) %>%
    body_add_par("Table B. Climate smooth term statistics across specifications", style = "heading 2") %>%
    body_add_flextable(ft_b) %>%
    body_add_par("Table C. Smooth term retention matrix", style = "heading 2") %>%
    body_add_flextable(ft_c)

  print(doc, target = file.path(output_path, "Covariate_Robustness_Tables.docx"))
}

# ============================================================================
# 7. MAIN EXECUTION
# ============================================================================
cat("[START]\n")

all_comp <- list()
all_info <- list()
all_curves <- list()
all_top6 <- list()
all_vif <- list()
all_shape <- list()
all_table_b <- list()
all_table_c <- list()

for (spec in bacteria_specs) {
  file_path <- file.path(input_data_dir, spec$file_name)
  if (!file.exists(file_path)) {
    warning("Skipping missing file: ", file_path)
    next
  }
  res <- run_bacterium(file_path, spec$title)
  all_comp[[spec$title]] <- res$comparison
  all_info[[spec$title]] <- res$model_info
  all_curves[[spec$title]] <- res$curves
  all_top6[[spec$title]] <- res$top6
  all_vif[[spec$title]] <- res$vif
  all_shape[[spec$title]] <- res$shape
  all_table_b[[spec$title]] <- res$table_b
  all_table_c[[spec$title]] <- res$table_c
}

full_comp <- bind_rows(all_comp)
full_info <- bind_rows(all_info)
full_top6 <- bind_rows(all_top6)
full_vif <- bind_rows(all_vif)
full_shape <- bind_rows(all_shape)
full_curve_data <- bind_rows(lapply(all_curves, function(x) bind_rows(x)))
full_table_b <- bind_rows(all_table_b)
full_table_c <- bind_rows(all_table_c)

vif_max <- full_vif %>%
  group_by(Bacteria, Spec) %>%
  summarise(Max_VIF = round(max(VIF, na.rm = TRUE), 2), .groups = "drop") %>%
  mutate(
    Spec = case_when(
      Spec == "Theory-6" ~ "Theory6_Max_VIF",
      Spec == "Top6 from PLS1" ~ "Top6_Max_VIF",
      TRUE ~ Spec
    )
  ) %>%
  pivot_wider(names_from = Spec, values_from = Max_VIF)

table_a <- full_info %>%
  left_join(vif_max, by = "Bacteria") %>%
  transmute(
    `AMR strain` = Bacteria,
    N,
    Countries,
    `PLS (Model C) R2` = A_Rsq,
    `PLS (Model C) Dev%` = A_Dev,
    `PCA R2` = B_Rsq,
    `PCA Dev%` = B_Dev,
    `Theory-selected covariates R2` = C_Rsq,
    `Theory-selected covariates Dev%` = C_Dev,
    `PLS-informed top covariates R2` = D_Rsq,
    `PLS-informed top covariates Dev%` = D_Dev,
    `Theory-selected covariates Max VIF` = Theory6_Max_VIF,
    `PLS-informed top covariates Max VIF` = Top6_Max_VIF
  )

table_b_export <- full_table_b %>%
  rename(`AMR strain` = Bacteria, `Climate variable` = Variable) %>%
  rename(
    `PLS (Model C)` = `PLS (main)`,
    `Theory-selected covariates` = `Theory-6`,
    `PLS-informed top covariates` = `Top6 from PLS1`
  )

table_c_export <- full_table_c %>%
  rename(`AMR strain` = Bacteria, `Climate variable` = Variable) %>%
  rename(
    `PLS (Model C)` = `PLS (main)`,
    `Theory-selected covariates` = `Theory-6`,
    `PLS-informed top covariates` = `Top6 from PLS1`
  )

write.csv(full_comp, file.path(output_root, "01_tables", "covariate_robustness_4specs_comparison.csv"), row.names = FALSE)
write.csv(full_info, file.path(output_root, "01_tables", "covariate_robustness_4specs_model_info.csv"), row.names = FALSE)
write.csv(full_top6, file.path(output_root, "01_tables", "covariate_robustness_top6_manifest.csv"), row.names = FALSE)
write.csv(full_vif, file.path(output_root, "01_tables", "covariate_robustness_linear_vif.csv"), row.names = FALSE)
write.csv(full_shape, file.path(output_root, "01_tables", "covariate_robustness_shape_similarity.csv"), row.names = FALSE)
write.csv(full_curve_data, file.path(output_root, "01_tables", "covariate_robustness_curve_plot_data.csv"), row.names = FALSE)
write.csv(table_a, file.path(output_root, "01_tables", "covariate_adjustment_model_level_model_level.csv"), row.names = FALSE)
write.csv(table_b_export, file.path(output_root, "01_tables", "covariate_adjustment_climate_smooth_statistics_climate_smooth_statistics.csv"), row.names = FALSE)
write.csv(table_c_export, file.path(output_root, "01_tables", "covariate_adjustment_smooth_term_retention_smooth_term_retention.csv"), row.names = FALSE)

cat("\n--- Generating figures ---\n")
p_overlay <- plot_overlay(all_curves, file.path(output_root, "02_figures"))
p_edf <- plot_edf_heat(full_comp, file.path(output_root, "02_figures"))
p_or <- plot_or_range(full_comp, file.path(output_root, "02_figures"))
p_shape <- plot_shape_heat(full_shape, file.path(output_root, "02_figures"))

if (!is.null(p_edf) && !is.null(p_or)) {
  p_comb <- plot_grid(p_edf, p_or, ncol = 1, rel_heights = c(1, 0.9), labels = c("A", "B"), label_size = 16, label_fontface = "bold")
  ggsave(file.path(output_root, "02_figures", "covariate_adjustment_combined_summary.pdf"), p_comb, width = 16, height = 11, dpi = 300, bg = "white")
  ggsave(file.path(output_root, "02_figures", "covariate_adjustment_combined_summary.png"), p_comb, width = 16, height = 11, dpi = 300, bg = "white")
}

wb_path <- file.path(output_root, "03_workbook", "covariate_adjustment_robustness.xlsx")
write_xlsx(
  list(
    Table_A_Model_Level = table_a,
    Table_B_Smooth_Statistics = table_b_export,
    Table_C_Retention = table_c_export,
    Model_Info = full_info,
    Climate_Term_Comparison = full_comp,
    Shape_Similarity = full_shape,
    Linear_VIF = full_vif,
    Top6_Manifest = full_top6,
    Curve_Plot_Data = full_curve_data
  ),
  wb_path
)

create_docx_tables(table_a, table_b_export, table_c_export, file.path(output_root, "03_workbook"))

write.csv(full_info, file.path(source_data_root, "01_csv", "CovariateRobustness_Model_Info.csv"), row.names = FALSE)
write.csv(full_comp, file.path(source_data_root, "01_csv", "CovariateRobustness_Climate_Term_Comparison.csv"), row.names = FALSE)
write.csv(full_shape, file.path(source_data_root, "01_csv", "CovariateRobustness_Shape_Similarity.csv"), row.names = FALSE)
write.csv(full_vif, file.path(source_data_root, "01_csv", "CovariateRobustness_Linear_VIF.csv"), row.names = FALSE)
write.csv(full_top6, file.path(source_data_root, "01_csv", "CovariateRobustness_Top6_Manifest.csv"), row.names = FALSE)
write.csv(full_curve_data, file.path(source_data_root, "01_csv", "CovariateRobustness_Curve_Plot_Data.csv"), row.names = FALSE)
write.csv(table_a, file.path(source_data_root, "01_csv", "covariate_adjustment_model_level_model_level.csv"), row.names = FALSE)
write.csv(table_b_export, file.path(source_data_root, "01_csv", "covariate_adjustment_climate_smooth_statistics_climate_smooth_statistics.csv"), row.names = FALSE)
write.csv(table_c_export, file.path(source_data_root, "01_csv", "covariate_adjustment_smooth_term_retention_smooth_term_retention.csv"), row.names = FALSE)

write_xlsx(
  list(
    Table_A_Model_Level = table_a,
    Table_B_Smooth_Statistics = table_b_export,
    Table_C_Retention = table_c_export,
    Model_Info = full_info,
    Climate_Term_Comparison = full_comp,
    Shape_Similarity = full_shape,
    Linear_VIF = full_vif,
    Top6_Manifest = full_top6,
    Curve_Plot_Data = full_curve_data
  ),
  file.path(source_data_root, "02_workbook", "SourceData_covariate_adjustment_robustness.xlsx")
)

sink(file.path(output_root, "04_metadata", "SUMMARY.txt"))
cat("================================================================\n")
cat("FOUR-SPECIFICATION COVARIATE ROBUSTNESS SUMMARY\n")
cat("================================================================\n\n")
cat("Model-level summary:\n")
print(as.data.frame(full_info))
cat("\nRetention summary (EDF >= 0.5):\n")
print(as.data.frame(full_comp %>% count(Spec, Retained)))
cat("\nSignificance summary (p < 0.05):\n")
print(as.data.frame(full_comp %>% mutate(Significant = p_value < 0.05) %>% count(Spec, Significant)))
cat("\nVIF maxima by specification:\n")
print(as.data.frame(full_vif %>% group_by(Spec) %>% summarise(Max_VIF = round(max(VIF, na.rm = TRUE), 3), .groups = "drop")))
  cat("\nMean shape concordance versus main model:\n")
print(as.data.frame(full_shape %>% group_by(Alt_Spec) %>% summarise(mean_rho = round(mean(rho, na.rm = TRUE), 3), .groups = "drop")))
sink()

cat("\n================================================================\n")
cat("DONE\n")
cat(output_root, "\n")
cat("================================================================\n")
