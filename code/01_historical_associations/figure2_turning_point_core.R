#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

options(stringsAsFactors = FALSE, warn = 1)

v3_specs <- tribble(
  ~code,     ~phenotype, ~input_file,                                                            ~model_file,
  "3GCREC", "3GCR-Ec", "3GCREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv", "ModelC_3GCR_Ec.rds",
  "3GCRKP", "3GCR-Kp", "3GCRKP_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv", "ModelC_3GCR_Kp.rds",
  "CRAB",   "CR-Ab",   "CRAB_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv",   "ModelC_CR_Ab.rds",
  "CREC",   "CR-Ec",   "CREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv",   "ModelC_CR_Ec.rds",
  "CRKP",   "CR-Kp",   "CRKP_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv",   "ModelC_CR_Kp.rds",
  "CRPA",   "CR-Pa",   "CRPA_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv",   "ModelC_CR_Pa.rds"
)

v3_factors <- tribble(
  ~variable,          ~raw,   ~climate,        ~title,                 ~unit, ~colour,
  "TMP_scaled_lag",  "TMP",  "Temperature",   "Temperature",         "degC", "#DD5F60",
  "HUM_scaled_lag",  "HUM",  "Humidity",      "Relative humidity",   "%",    "#3CB371",
  "PREC_scaled_lag", "PREC", "Precipitation", "Precipitation",       "mm",   "#9AC0CD",
  "WET_scaled_lag",  "WET",  "WetDays",       "Wet days",            "d",    "#8A2BE2"
)

v3_safe_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

v3_mode_value <- function(x) {
  tab <- table(x, useNA = "no")
  if (length(tab) == 0L) return(NA)
  names(tab)[which.max(tab)]
}

v3_load_lags <- function(path) {
  x <- read.csv(path, check.names = FALSE)
  required <- c("Display_Name", "TMP_lag", "PREC_lag", "HUM_lag", "WET_lag")
  if (length(setdiff(required, names(x))) > 0L) {
    stop("Lag summary is missing required columns.", call. = FALSE)
  }
  setNames(lapply(seq_len(nrow(x)), function(i) {
    c(
      TMP = as.integer(x$TMP_lag[i]),
      PREC = as.integer(x$PREC_lag[i]),
      HUM = as.integer(x$HUM_lag[i]),
      WET = as.integer(x$WET_lag[i])
    )
  }), x$Display_Name)
}

v3_prepare_data <- function(path, phenotype, lag_settings) {
  if (!file.exists(path)) stop("Input not found: ", path, call. = FALSE)
  lag_cfg <- lag_settings[[phenotype]]
  if (is.null(lag_cfg)) stop("Missing lag settings for ", phenotype, call. = FALSE)

  dat <- read.csv(path, check.names = FALSE) %>%
    mutate(
      year = as.numeric(as.character(year)),
      Region = factor(Region),
      climate_zone = case_when(
        abs(lat) > 66.5 ~ "Polar Zone",
        abs(lat) > 23.5 ~ "Temperate Zone",
        TRUE ~ "Tropical Zone"
      ),
      climate_zone = factor(climate_zone),
      HUM = pmin(HUM, 100)
    ) %>%
    group_by(NAME) %>%
    mutate(location_id = cur_group_id()) %>%
    ungroup() %>%
    mutate(
      TMP_orig = TMP,
      PREC_orig = PREC,
      HUM_orig = HUM,
      WET_orig = WET
    ) %>%
    group_by(climate_zone) %>%
    mutate(
      across(c(TMP, PREC, HUM, WET), ~ as.vector(scale(.x)), .names = "{.col}_scaled")
    ) %>%
    group_by(location_id) %>%
    # Preserve the preprocessing order used by the fitted model: arrange by year while
    # grouped, without .by_group=TRUE.
    arrange(year) %>%
    mutate(
      TMP_scaled_lag = dplyr::lag(TMP_scaled, lag_cfg[["TMP"]]),
      PREC_scaled_lag = dplyr::lag(PREC_scaled, lag_cfg[["PREC"]]),
      HUM_scaled_lag = dplyr::lag(HUM_scaled, lag_cfg[["HUM"]]),
      WET_scaled_lag = dplyr::lag(WET_scaled, lag_cfg[["WET"]]),
      TMP_orig_lag = dplyr::lag(TMP_orig, lag_cfg[["TMP"]]),
      PREC_orig_lag = dplyr::lag(PREC_orig, lag_cfg[["PREC"]]),
      HUM_orig_lag = dplyr::lag(HUM_orig, lag_cfg[["HUM"]]),
      WET_orig_lag = dplyr::lag(WET_orig, lag_cfg[["WET"]])
    ) %>%
    filter(
      !is.na(TMP_scaled_lag), !is.na(PREC_scaled_lag),
      !is.na(HUM_scaled_lag), !is.na(WET_scaled_lag)
    ) %>%
    ungroup()

  dat
}

v3_validate_model_data <- function(model, prepared, phenotype, tolerance = 1e-10) {
  mf <- model$model
  if (nrow(mf) != nrow(prepared)) {
    return(tibble(
      phenotype = phenotype,
      check = "row_count",
      passed = FALSE,
      max_abs_difference = NA_real_,
      details = paste0("model=", nrow(mf), "; prepared=", nrow(prepared))
    ))
  }

  common <- intersect(names(mf), names(prepared))
  rows <- lapply(common, function(nm) {
    a <- mf[[nm]]
    b <- prepared[[nm]]
    if (is.factor(a) || is.character(a)) {
      same <- identical(as.character(a), as.character(b))
      tibble(
        phenotype = phenotype, check = nm, passed = same,
        max_abs_difference = NA_real_, details = ifelse(same, "exact", "factor/character mismatch")
      )
    } else {
      delta <- suppressWarnings(max(abs(as.numeric(a) - as.numeric(b)), na.rm = TRUE))
      if (!is.finite(delta)) delta <- NA_real_
      same <- isTRUE(delta <= tolerance)
      tibble(
        phenotype = phenotype, check = nm, passed = same,
        max_abs_difference = delta,
        details = ifelse(same, "within tolerance", "numeric mismatch")
      )
    }
  })
  bind_rows(rows)
}

v3_prediction_template <- function(model, variable, x_grid) {
  mf <- model$model
  response <- all.vars(formula(model))[1]
  predictor_names <- setdiff(names(mf), response)
  base <- mf[1, predictor_names, drop = FALSE]

  for (nm in predictor_names) {
    x <- mf[[nm]]
    if (is.factor(x)) {
      base[[nm]] <- factor(v3_mode_value(x), levels = levels(x))
    } else if (is.character(x)) {
      base[[nm]] <- v3_mode_value(x)
    } else if (is.numeric(x) || is.integer(x)) {
      base[[nm]] <- mean(x, na.rm = TRUE)
    }
  }

  out <- base[rep(1, length(x_grid)), , drop = FALSE]
  out[[variable]] <- x_grid
  out
}

v3_target_term_index <- function(prediction, variable) {
  target <- paste0("s(", variable, ")")
  idx <- match(target, colnames(prediction))
  if (is.na(idx)) {
    idx <- grep(paste0("^s\\(", variable, "\\)"), colnames(prediction))
  }
  if (length(idx) != 1L || is.na(idx)) {
    stop("Could not uniquely identify smooth term for ", variable, call. = FALSE)
  }
  idx
}

v3_smooth_statistics <- function(model, variable) {
  st <- as.data.frame(summary(model)$s.table)
  st$term <- rownames(st)
  rownames(st) <- NULL
  target <- paste0("s(", variable, ")")
  z <- st[st$term == target, , drop = FALSE]
  if (nrow(z) != 1L) {
    return(tibble(edf = NA_real_, ref_df = NA_real_, statistic = NA_real_, p_value = NA_real_))
  }
  statistic_col <- intersect(c("F", "Chi.sq"), names(z))[1]
  p_col <- grep("p-value", names(z), fixed = TRUE, value = TRUE)[1]
  ref_col <- intersect(c("Ref.df", "ref.df"), names(z))[1]
  tibble(
    edf = as.numeric(z$edf[1]),
    ref_df = if (!is.na(ref_col)) as.numeric(z[[ref_col]][1]) else NA_real_,
    statistic = if (!is.na(statistic_col)) as.numeric(z[[statistic_col]][1]) else NA_real_,
    p_value = if (!is.na(p_col)) as.numeric(z[[p_col]][1]) else NA_real_
  )
}

v3_term_curve <- function(model, variable, x_grid, with_se = TRUE) {
  nd <- v3_prediction_template(model, variable, x_grid)
  if (with_se) {
    pr <- suppressWarnings(predict(model, nd, type = "terms", se.fit = TRUE, unconditional = FALSE))
    idx <- v3_target_term_index(pr$fit, variable)
    eta <- as.numeric(pr$fit[, idx])
    se <- as.numeric(pr$se.fit[, idx])
  } else {
    pr <- suppressWarnings(predict(model, nd, type = "terms", se.fit = FALSE))
    idx <- v3_target_term_index(pr, variable)
    eta <- as.numeric(pr[, idx])
    se <- rep(NA_real_, length(eta))
  }

  tibble(
    x_scaled = x_grid,
    eta = eta,
    se = se,
    OR = exp(eta),
    Lower_95CI = exp(eta - 1.96 * se),
    Upper_95CI = exp(eta + 1.96 * se)
  )
}

v3_fill_zero_signs <- function(signs) {
  if (!any(signs != 0L)) return(signs)
  out <- signs
  # Forward and backward passes bridge the near-zero derivative neighborhood
  # without imposing a direction at a completely flat curve.
  for (i in seq_along(out)) {
    if (out[i] == 0L && i > 1L && out[i - 1L] != 0L) out[i] <- out[i - 1L]
  }
  for (i in rev(seq_along(out))) {
    if (out[i] == 0L && i < length(out) && out[i + 1L] != 0L) out[i] <- out[i + 1L]
  }
  out
}

v3_empty_candidates <- function() {
  tibble(
    x_scaled = numeric(), eta = numeric(), OR = numeric(),
    direction = character(), prominence_eta = numeric(),
    type = character(), point_id = character()
  )
}

v3_detect_turning_points <- function(curve, edf, amplitude_epsilon = 1e-4, edf_cutoff = 1.05) {
  curve <- curve %>% arrange(x_scaled)
  x <- curve$x_scaled
  eta <- curve$eta
  if (length(x) < 21L || any(!is.finite(x)) || any(!is.finite(eta))) {
    return(v3_empty_candidates())
  }

  eta_span <- diff(range(eta))
  if (!is.finite(edf) || edf <= edf_cutoff || eta_span < amplitude_epsilon) {
    return(v3_empty_candidates())
  }

  derivative <- c(
    (eta[2] - eta[1]) / (x[2] - x[1]),
    (eta[3:length(eta)] - eta[1:(length(eta) - 2L)]) /
      (x[3:length(x)] - x[1:(length(x) - 2L)]),
    (eta[length(eta)] - eta[length(eta) - 1L]) /
      (x[length(x)] - x[length(x) - 1L])
  )

  derivative_tolerance <- max(1e-10, max(abs(derivative), na.rm = TRUE) * 1e-5)
  signs <- sign(derivative)
  signs[abs(derivative) <= derivative_tolerance] <- 0L
  signs <- v3_fill_zero_signs(signs)

  change_idx <- which(signs[-length(signs)] * signs[-1L] < 0L)
  if (length(change_idx) == 0L) return(v3_empty_candidates())

  rows <- lapply(change_idx, function(i) {
    d1 <- derivative[i]
    d2 <- derivative[i + 1L]
    denom <- abs(d1) + abs(d2)
    frac <- if (is.finite(denom) && denom > 0) abs(d1) / denom else 0.5
    x0 <- x[i] + frac * (x[i + 1L] - x[i])
    eta0 <- approx(x, eta, xout = x0, ties = "ordered")$y
    direction <- if (signs[i] > 0L && signs[i + 1L] < 0L) "Max" else "Min"
    tibble(x_scaled = x0, eta = eta0, OR = exp(eta0), direction = direction)
  })
  roots <- bind_rows(rows) %>% arrange(x_scaled)

  # Local prominence is calculated against the lowest/highest fitted value on
  # each side up to the adjacent turning point or support boundary.
  roots$prominence_eta <- vapply(seq_len(nrow(roots)), function(i) {
    left_bound <- if (i == 1L) min(x) else roots$x_scaled[i - 1L]
    right_bound <- if (i == nrow(roots)) max(x) else roots$x_scaled[i + 1L]
    left_eta <- eta[x >= left_bound & x <= roots$x_scaled[i]]
    right_eta <- eta[x >= roots$x_scaled[i] & x <= right_bound]
    if (roots$direction[i] == "Max") {
      roots$eta[i] - max(min(left_eta), min(right_eta))
    } else {
      min(max(left_eta), max(right_eta)) - roots$eta[i]
    }
  }, numeric(1))

  roots <- roots %>%
    filter(is.finite(prominence_eta), prominence_eta >= amplitude_epsilon) %>%
    mutate(
      type = direction,
      point_id = paste0(direction, "_", row_number())
    )
  if (nrow(roots) == 0L) return(roots)

  global_max_i <- which.max(eta)
  global_min_i <- which.min(eta)
  grid_step <- median(diff(x))
  roots <- roots %>%
    mutate(
      type = case_when(
        direction == "Max" & abs(x_scaled - x[global_max_i]) <= 2 * grid_step ~ "GMax",
        direction == "Min" & abs(x_scaled - x[global_min_i]) <= 2 * grid_step ~ "GMin",
        TRUE ~ direction
      ),
      point_id = paste0(type, "_", row_number())
    )
  roots
}

v3_observed_analysis <- function(model, prepared, phenotype, factor_row, n_curve = 801L, n_search = 1001L) {
  variable <- factor_row$variable
  x_obs <- prepared[[variable]]
  x_obs <- x_obs[is.finite(x_obs)]
  q <- as.numeric(quantile(x_obs, c(0.025, 0.975), type = 7, names = FALSE))
  observed_range <- range(x_obs)
  lag_years <- unique(v3_load_lags_cache[[phenotype]][[factor_row$raw]])

  curve <- v3_term_curve(model, variable, seq(observed_range[1], observed_range[2], length.out = n_curve), with_se = TRUE) %>%
    mutate(
      phenotype = phenotype,
      climate = factor_row$climate,
      variable = variable,
      lag_years = lag_years,
      inside_P2_5_P97_5 = x_scaled >= q[1] & x_scaled <= q[2]
    )
  search_curve <- v3_term_curve(model, variable, seq(q[1], q[2], length.out = n_search), with_se = FALSE)
  stat <- v3_smooth_statistics(model, variable)
  candidates <- v3_detect_turning_points(search_curve, stat$edf) %>%
    mutate(
      phenotype = phenotype,
      climate = factor_row$climate,
      variable = variable,
      lag_years = lag_years,
      observed_candidate_id = paste0(v3_safe_stem(phenotype), "__", factor_row$climate, "__", point_id),
      support_width_scaled = diff(q),
      p2_5_scaled = q[1],
      p97_5_scaled = q[2]
    )

  list(
    curve = curve,
    search_curve = search_curve,
    candidates = candidates,
    support = tibble(
      phenotype = phenotype,
      climate = factor_row$climate,
      variable = variable,
      lag_years = lag_years,
      observed_min_scaled = observed_range[1],
      p2_5_scaled = q[1],
      p97_5_scaled = q[2],
      observed_max_scaled = observed_range[2],
      support_width_scaled = diff(q),
      n_model = nrow(prepared),
      n_countries = n_distinct(prepared$NAME),
      edf = stat$edf,
      ref_df = stat$ref_df,
      statistic = stat$statistic,
      p_value = stat$p_value
    )
  )
}

v3_refit_once <- function(model, seed) {
  set.seed(seed)
  dat <- model$model
  response <- all.vars(formula(model))[1]
  sigma <- sqrt(model$sig2)
  if (!is.finite(sigma) || sigma <= 0) stop("Invalid Gaussian residual SD.", call. = FALSE)
  dat[[response]] <- fitted(model) + rnorm(nrow(dat), mean = 0, sd = sigma)

  warnings <- character()
  started <- proc.time()[[3]]
  fit <- tryCatch(
    withCallingHandlers(
      bam(
        formula(model), data = dat, family = gaussian(), method = "REML",
        select = TRUE, discrete = FALSE, use.chol = FALSE,
        control = gam.control(nthreads = 1, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
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
    error_text = if (converged) "" else "Fit returned but convergence/finite-coefficient check failed"
  )
}

v3_match_one_replicate <- function(observed, boot, support_width, tolerance_fraction) {
  if (nrow(observed) == 0L) return(tibble())
  ans <- observed %>%
    transmute(
      observed_candidate_id,
      observed_direction = direction,
      observed_x_scaled = x_scaled,
      matched = FALSE,
      bootstrap_x_scaled = NA_real_,
      bootstrap_OR = NA_real_,
      bootstrap_prominence_eta = NA_real_,
      match_distance_scaled = NA_real_
    )
  if (nrow(boot) == 0L) return(ans)

  tolerance <- tolerance_fraction * support_width

  # Exact maximum-cardinality bipartite matching with minimum total distance
  # as the tie-breaker. Candidate counts are small, so exhaustive recursion is
  # transparent and avoids the under-matching that a nearest-pair greedy rule
  # can produce.
  best <- new.env(parent = emptyenv())
  best$n_match <- -1L
  best$cost <- Inf
  best$pairs <- matrix(integer(), ncol = 2L)

  recurse <- function(obs_position, used_boot, pairs, cost) {
    if (obs_position > nrow(observed)) {
      n_match <- nrow(pairs)
      if (n_match > best$n_match || (n_match == best$n_match && cost < best$cost - 1e-15)) {
        best$n_match <- n_match
        best$cost <- cost
        best$pairs <- pairs
      }
      return(invisible(NULL))
    }

    # Option 1: leave this observed candidate unmatched.
    recurse(obs_position + 1L, used_boot, pairs, cost)

    # Option 2: match it to each eligible unused bootstrap candidate.
    eligible <- which(
      !(seq_len(nrow(boot)) %in% used_boot) &
        boot$direction == observed$direction[obs_position] &
        abs(boot$x_scaled - observed$x_scaled[obs_position]) <= tolerance
    )
    for (bi in eligible) {
      recurse(
        obs_position + 1L,
        c(used_boot, bi),
        rbind(pairs, c(obs_position, bi)),
        cost + abs(boot$x_scaled[bi] - observed$x_scaled[obs_position])
      )
    }
    invisible(NULL)
  }
  recurse(1L, integer(), matrix(integer(), ncol = 2L), 0)

  if (nrow(best$pairs) > 0L) {
    for (i in seq_len(nrow(best$pairs))) {
      oi <- best$pairs[i, 1]
      bi <- best$pairs[i, 2]
      ans$matched[oi] <- TRUE
      ans$bootstrap_x_scaled[oi] <- boot$x_scaled[bi]
      ans$bootstrap_OR[oi] <- boot$OR[bi]
      ans$bootstrap_prominence_eta[oi] <- boot$prominence_eta[bi]
      ans$match_distance_scaled[oi] <- abs(boot$x_scaled[bi] - observed$x_scaled[oi])
    }
  }
  ans
}

# Populated by the driver after reading the packaged lag summary. It is kept as
# a dedicated cache so every function uses the same pre-specified lag values.
v3_load_lags_cache <- NULL
