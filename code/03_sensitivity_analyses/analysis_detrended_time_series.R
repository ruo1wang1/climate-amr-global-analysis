###############################################################################
# Analysis of climate-AMR associations using detrended time series
###############################################################################

cat("\n============================================================\n")
cat("ANALYSIS OF CLIMATE-AMR ASSOCIATIONS USING DETRENDED TIME SERIES\n")
cat("============================================================\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(scales)
  library(writexl)
})

# 1. Paths and global settings
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

default_output_root <- file.path(
  revision_root,
  "outputs",
  "analysis_detrended_time_series"
)
output_root <- Sys.getenv("DETREND_OUTPUT_ROOT", unset = default_output_root)

default_source_data_root <- file.path(
  revision_root,
  "data/source_data",
  "detrended_time_series"
)
source_data_root <- Sys.getenv("DETREND_SOURCE_DATA_ROOT", unset = default_source_data_root)

model_ready_input_name <- function(code) {
  paste0(code, "_model_ready_data.csv")
}

dir.create(file.path(output_root, "01_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "02_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "03_workbook"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "04_metadata"), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(source_data_root, "01_csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_data_root, "02_workbook"), recursive = TRUE, showWarnings = FALSE)

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = model_ready_input_name("3GCREC")),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = model_ready_input_name("3GCRKP")),
  list(code = "CRAB", title = "CR-Ab", file_name = model_ready_input_name("CRAB")),
  list(code = "CREC", title = "CR-Ec", file_name = model_ready_input_name("CREC")),
  list(code = "CRKP", title = "CR-Kp", file_name = model_ready_input_name("CRKP")),
  list(code = "CRPA", title = "CR-Pa", file_name = model_ready_input_name("CRPA"))
)

bacteria_levels <- vapply(bacteria_specs, `[[`, character(1), "title")
climate_levels <- c("Temperature", "Humidity", "Precipitation", "Wet Days")
min_obs_country <- 5L

climate_colors <- c(
  "Temperature" = "#C0392B",
  "Humidity" = "#27AE60",
  "Precipitation" = "#2980B9",
  "Wet Days" = "#8E44AD"
)

robustness_pal <- c(
  "Retained" = "#2166AC",
  "Direction changed" = "#E08214",
  "Attenuated" = "#92C5DE",
  "Lost" = "#D6604D",
  "Gained" = "#4DAF4A",
  "Concordant ns" = "#B8B8B8"
)

fig_width_in <- 7.2
fig_height_ds1 <- 8.5
fig_height_ds2 <- 9.0
fig_height_ds3 <- 5.5

phys_limits <- list(
  TMP = c(-10, 40),
  HUM = c(30, 100),
  PREC = c(0, 3200),
  WET = c(0, 300)
)

climate_specs <- list(
  list(var = "TMP_scaled_lag", orig = "TMP", label = "Temperature"),
  list(var = "HUM_scaled_lag", orig = "HUM", label = "Humidity"),
  list(var = "PREC_scaled_lag", orig = "PREC", label = "Precipitation"),
  list(var = "WET_scaled_lag", orig = "WET", label = "Wet Days")
)

x_axis_labels <- c(
  "Temperature" = "Temperature (\u00B0C)",
  "Humidity" = "Relative humidity (%)",
  "Precipitation" = "Precipitation (mm)",
  "Wet Days" = "Wet days (days)"
)

# 2. Helpers
load_lag_settings <- function(summary_csv) {
  if (!file.exists(summary_csv)) {
    stop("Lag summary file not found: ", summary_csv, call. = FALSE)
  }

  lag_df <- read.csv(summary_csv, stringsAsFactors = FALSE)
  required_cols <- c("Display_Name", "TMP_lag", "PREC_lag", "HUM_lag", "WET_lag")
  missing_cols <- setdiff(required_cols, names(lag_df))
  if (length(missing_cols) > 0) {
    stop("Lag summary file missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  lag_settings <- setNames(vector("list", nrow(lag_df)), lag_df$Display_Name)
  for (i in seq_len(nrow(lag_df))) {
    lag_settings[[lag_df$Display_Name[i]]] <- list(
      temp_lag = as.integer(lag_df$TMP_lag[i]),
      precip_lag = as.integer(lag_df$PREC_lag[i]),
      humid_lag = as.integer(lag_df$HUM_lag[i]),
      wetdays_lag = as.integer(lag_df$WET_lag[i])
    )
  }
  lag_settings
}

lag_settings <- load_lag_settings(lag_summary_path)

get_available_pls_components <- function(data) {
  pls_candidates <- paste0("PLS_Comp", 1:4)
  present <- pls_candidates[pls_candidates %in% names(data)]
  present[sapply(present, function(x) !all(is.na(data[[x]])))]
}

safe_linear_detrend <- function(y, t) {
  keep <- is.finite(y) & is.finite(t)
  if (sum(keep) < 3 || length(unique(t[keep])) < 2) {
    return(y)
  }

  fit <- lm(y[keep] ~ t[keep])
  slope <- coef(fit)[2]
  centered_t <- t - mean(t[keep], na.rm = TRUE)
  y - slope * centered_t
}

read_base_data <- function(file_path) {
  read.csv(file_path, stringsAsFactors = FALSE) %>%
    mutate(
      year = as.numeric(as.character(year)),
      Region = factor(Region),
      HUM = pmin(HUM, 100),
      climate_zone = case_when(
        abs(lat) > 66.5 ~ "Polar Zone",
        abs(lat) > 23.5 ~ "Temperate Zone",
        TRUE ~ "Tropical Zone"
      ),
      climate_zone = factor(climate_zone, levels = c("Polar Zone", "Temperate Zone", "Tropical Zone"))
    ) %>%
    group_by(NAME) %>%
    mutate(location_id = cur_group_id()) %>%
    ungroup()
}

compute_scale_params <- function(data) {
  data %>%
    summarise(
      across(
        c(TMP, PREC, HUM, WET),
        list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))
      )
    )
}

apply_common_detrend <- function(data) {
  detrended <- data
  year_mean <- mean(detrended$year, na.rm = TRUE)

  for (var in c("logit_R", "TMP", "PREC", "HUM", "WET")) {
    fit <- lm(reformulate(c("year", "factor(location_id)"), response = var), data = detrended)
    slope <- coef(fit)[["year"]]
    detrended[[var]] <- detrended[[var]] - slope * (detrended$year - year_mean)
  }

  detrended
}

apply_country_detrend <- function(data, min_obs = min_obs_country) {
  data %>%
    group_by(location_id) %>%
    mutate(n_obs_country = n()) %>%
    ungroup() %>%
    filter(n_obs_country >= min_obs) %>%
    group_by(location_id) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      logit_R = safe_linear_detrend(logit_R, year),
      TMP = safe_linear_detrend(TMP, year),
      PREC = safe_linear_detrend(PREC, year),
      HUM = safe_linear_detrend(HUM, year),
      WET = safe_linear_detrend(WET, year)
    ) %>%
    ungroup() %>%
    mutate(HUM = pmin(HUM, 100))
}

prepare_model_data <- function(data, lag_config, use_year_factor = FALSE) {
  out <- data %>%
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
      TMP_scaled_lag = lag(TMP_scaled, lag_config$temp_lag),
      PREC_scaled_lag = lag(PREC_scaled, lag_config$precip_lag),
      HUM_scaled_lag = lag(HUM_scaled, lag_config$humid_lag),
      WET_scaled_lag = lag(WET_scaled, lag_config$wetdays_lag)
    ) %>%
    ungroup() %>%
    filter(
      !is.na(TMP_scaled_lag),
      !is.na(PREC_scaled_lag),
      !is.na(HUM_scaled_lag),
      !is.na(WET_scaled_lag),
      !is.na(logit_R)
    )

  if (use_year_factor) {
    out <- out %>% mutate(year_factor = factor(year))
  }
  out
}

build_formula <- function(data, use_year_factor = FALSE) {
  pls_terms <- get_available_pls_components(data)
  pls_rhs <- paste0("s(", pls_terms, ", k = 10, bs = 'cr')")

  time_term <- if (use_year_factor) {
    "year_factor"
  } else {
    "s(year, bs = 'cr', k = 8)"
  }

  rhs_terms <- c(
    "s(TMP_scaled_lag, k = 5, bs = 'cr')",
    "s(PREC_scaled_lag, k = 10, bs = 'cr')",
    "s(HUM_scaled_lag, k = 10, bs = 'cr')",
    "s(WET_scaled_lag, k = 10, bs = 'cr')",
    pls_rhs,
    "s(lat, lon, bs = 'sos', k = 20)",
    time_term,
    "s(Region, bs = 're')",
    "climate_zone"
  )

  as.formula(paste("logit_R ~", paste(rhs_terms, collapse = " + ")))
}

fit_gamm <- function(data, use_year_factor = FALSE) {
  form <- build_formula(data, use_year_factor = use_year_factor)
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)

  tryCatch(
    bam(
      form,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE,
      control = ctrl
    ),
    error = function(e) {
      message("bam() failed, falling back to gam(): ", conditionMessage(e))
      gam(
        form,
        data = data,
        family = gaussian(),
        method = "REML",
        select = TRUE
      )
    }
  )
}

extract_stats <- function(model, bacteria, model_label) {
  s_tbl <- as.data.frame(summary(model)$s.table)
  s_tbl$term <- rownames(s_tbl)

  s_tbl %>%
    as_tibble() %>%
    filter(grepl("TMP_scaled_lag|PREC_scaled_lag|HUM_scaled_lag|WET_scaled_lag", term)) %>%
    mutate(
      Bacteria = bacteria,
      Model = model_label,
      Climate = case_when(
        grepl("TMP", term) ~ "Temperature",
        grepl("HUM", term) ~ "Humidity",
        grepl("PREC", term) ~ "Precipitation",
        grepl("WET", term) ~ "Wet Days",
        TRUE ~ NA_character_
      ),
      Climate = factor(Climate, levels = climate_levels)
    ) %>%
    rename(EDF = edf, F_stat = F, p_value = `p-value`) %>%
    select(Bacteria, Model, Climate, EDF, F_stat, p_value)
}

predict_smoothed <- function(model, data, var_name, climate_orig, scale_params, n = 500) {
  gm <- scale_params[[paste0(climate_orig, "_mean")]]
  gs <- scale_params[[paste0(climate_orig, "_sd")]]
  limits <- phys_limits[[climate_orig]]
  if (climate_orig == "HUM") limits[2] <- 100

  x_lo <- (limits[1] - gm) / gs
  x_hi <- (limits[2] - gm) / gs
  xs <- seq(x_lo, x_hi, length.out = n)

  av <- names(model$var.summary)
  pred_data <- as.data.frame(setNames(lapply(av, function(v) {
    if (v == var_name) return(xs)

    val <- data[[v]]
    if (is.factor(val)) {
      return(factor(rep(names(which.max(table(val))), n), levels = levels(val)))
    }
    rep(mean(val, na.rm = TRUE), n)
  }), av))

  pr <- suppressWarnings(predict(model, pred_data, type = "terms", se.fit = TRUE))
  var_col <- grep(var_name, colnames(pr$fit), fixed = TRUE)
  if (length(var_col) == 0) {
    return(tibble(x = xs, or = rep(1, n), lo = rep(1, n), hi = rep(1, n)))
  }

  fit_vals <- pr$fit[, var_col[1]]
  se_vals <- pr$se.fit[, var_col[1]]

  tibble(
    x = xs,
    x_orig = xs * gs + gm,
    or = exp(fit_vals),
    lo = exp(fit_vals - 1.96 * se_vals),
    hi = exp(fit_vals + 1.96 * se_vals)
  )
}

compute_shape <- function(model_primary, model_detrended, data_primary, data_detrended, scale_params) {
  out <- list()

  for (ci in climate_specs) {
    c1 <- predict_smoothed(model_primary, data_primary, ci$var, ci$orig, scale_params)
    c2 <- predict_smoothed(model_detrended, data_detrended, ci$var, ci$orig, scale_params)
    z_lo <- max(min(c1$x), min(c2$x))
    z_hi <- min(max(c1$x), max(c2$x))
    zs <- seq(z_lo, z_hi, length.out = 200)
    o1 <- approx(c1$x, c1$or, xout = zs)$y
    o2 <- approx(c2$x, c2$or, xout = zs)$y

    p10 <- max(1L, round(length(zs) * 0.10))
    p90 <- min(length(zs), round(length(zs) * 0.90))
    out[[ci$label]] <- list(
      rho = suppressWarnings(cor(o1, o2, method = "spearman", use = "complete.obs")),
      dir_match = sign(o1[p90] - o1[p10]) == sign(o2[p90] - o2[p10])
    )
  }

  out
}

compute_common_diagnostics <- function(data, bacteria) {
  raw <- data
  year_mean <- mean(raw$year, na.rm = TRUE)
  vars <- c("logit_R", "TMP", "PREC", "HUM", "WET")
  out <- vector("list", length(vars))

  for (i in seq_along(vars)) {
    var <- vars[i]
    fit_before <- lm(reformulate(c("year", "factor(location_id)"), response = var), data = raw)
    slope_before <- summary(fit_before)$coefficients["year", ]

    raw_after <- raw
    raw_after[[var]] <- raw_after[[var]] - slope_before["Estimate"] * (raw_after$year - year_mean)
    fit_after <- lm(reformulate(c("year", "factor(location_id)"), response = var), data = raw_after)
    slope_after <- summary(fit_after)$coefficients["year", ]

    out[[i]] <- tibble(
      Bacteria = bacteria,
      Variable = var,
      Slope_Before = slope_before["Estimate"],
      t_Before = slope_before["t value"],
      p_Before = slope_before["Pr(>|t|)"],
      Slope_After = slope_after["Estimate"],
      t_After = slope_after["t value"],
      p_After = slope_after["Pr(>|t|)"]
    )
  }

  bind_rows(out)
}

# 3. Figures
build_robustness_lookup <- function(stats_primary_common, shape_all) {
  shape_df <- bind_rows(lapply(names(shape_all), function(b) {
    shape <- shape_all[[b]]
    if (is.null(shape)) return(NULL)
    tibble(
      Bacteria = b,
      Climate = factor(names(shape), levels = climate_levels),
      rho = vapply(shape, `[[`, numeric(1), "rho"),
      dir_match = vapply(shape, `[[`, logical(1), "dir_match")
    )
  }))

  primary_ref <- stats_primary_common %>%
    filter(Model == "Primary") %>%
    select(Bacteria, Climate, p_primary = p_value)

  stats_primary_common %>%
    filter(Model == "Common-DT") %>%
    left_join(primary_ref, by = c("Bacteria", "Climate")) %>%
    left_join(shape_df, by = c("Bacteria", "Climate")) %>%
    mutate(
      robustness = case_when(
        p_value < 0.05 & p_primary < 0.05 & (is.na(dir_match) | dir_match) ~ "Retained",
        p_value < 0.05 & p_primary < 0.05 & !dir_match ~ "Direction changed",
        p_value >= 0.10 & p_primary < 0.05 ~ "Lost",
        p_value >= 0.05 & p_value < 0.10 & p_primary < 0.05 ~ "Attenuated",
        p_value < 0.05 & p_primary >= 0.05 ~ "Gained",
        TRUE ~ "Concordant ns"
      ),
      robustness = factor(
        robustness,
        levels = c("Retained", "Direction changed", "Attenuated", "Lost", "Gained", "Concordant ns")
      )
    )
}

theme_nature <- function(base_size = 7) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      text = element_text(family = "Helvetica", colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "pt"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.4),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.key.size = unit(8, "pt"),
      legend.background = element_blank(),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "grey40", hjust = 0, margin = margin(b = 3)),
      plot.margin = margin(3, 4, 2, 3)
    )
}

build_fig_ds1 <- function(all_results, shape_all, stats_primary_common) {
  rob_lkp <- build_robustness_lookup(stats_primary_common, shape_all)

  rob_fill_alpha <- c(
    "Retained" = 0.22,
    "Direction changed" = 0.22,
    "Attenuated" = 0.15,
    "Lost" = 0.15,
    "Gained" = 0.15,
    "Concordant ns" = 0.10
  )

  panels <- list()

  for (b in bacteria_levels) {
    res <- all_results[[b]]
    if (is.null(res)) next

    for (ci in climate_specs) {
      ri <- rob_lkp %>%
        filter(Bacteria == b, as.character(Climate) == ci$label)

      rho_val <- if (nrow(ri) > 0 && !is.na(ri$rho[1])) ri$rho[1] else NA_real_
      rob_cat <- if (nrow(ri) > 0) as.character(ri$robustness[1]) else "Concordant ns"
      p_prim <- if (nrow(ri) > 0) ri$p_primary[1] else 1.0
      p_dt <- if (nrow(ri) > 0) ri$p_value[1] else 1.0

      clim_col <- climate_colors[[ci$label]]
      rob_col <- robustness_pal[[rob_cat]]
      rib_alpha <- rob_fill_alpha[[rob_cat]]

      c1 <- predict_smoothed(
        res$primary_model, res$primary_data, ci$var, ci$orig, res$scale_params
      ) %>% mutate(Model = "Primary (Model C)")

      c2 <- predict_smoothed(
        res$common_model, res$common_data, ci$var, ci$orig, res$scale_params
      ) %>% mutate(Model = "Common detrended")

      z_lo <- max(min(c1$x), min(c2$x))
      z_hi <- min(max(c1$x), max(c2$x))

      curves <- bind_rows(c1, c2) %>%
        filter(x >= z_lo, x <= z_hi) %>%
        mutate(Model = factor(Model, levels = c("Primary (Model C)", "Common detrended")))

      show_y <- identical(ci$orig, "TMP")
      show_x <- identical(b, tail(bacteria_levels, 1))

      x_limits <- phys_limits[[ci$orig]]

      p <- ggplot(curves, aes(x = x_orig, y = or)) +
        geom_hline(yintercept = 1, linetype = "longdash", colour = "grey60", linewidth = 0.25) +
        geom_ribbon(
          data = filter(curves, Model == "Primary (Model C)"),
          aes(ymin = lo, ymax = hi),
          fill = rob_col, alpha = rib_alpha
        ) +
        geom_ribbon(
          data = filter(curves, Model == "Common detrended"),
          aes(ymin = lo, ymax = hi),
          fill = "grey50", alpha = 0.08
        ) +
        geom_line(
          data = filter(curves, Model == "Primary (Model C)"),
          colour = clim_col, linewidth = 0.65
        ) +
        geom_line(
          data = filter(curves, Model == "Common detrended"),
          colour = "grey20", linewidth = 0.45, linetype = "dashed"
        ) +
        scale_x_continuous(
          limits = x_limits,
          breaks = pretty(x_limits, n = 4),
          expand = expansion(mult = c(0.01, 0.01))
        ) +
        theme_nature(base_size = 6.5) +
        theme(plot.margin = margin(1.5, 2, 1.5, 2))

      if (!is.na(rho_val)) {
        p <- p + annotate(
          "text", x = Inf, y = Inf, label = sprintf("rho=%.2f", rho_val),
          hjust = 1.05, vjust = 1.45, size = 1.8, colour = rob_col, fontface = "italic"
        )
      }

      if (p_prim < 0.05) {
        sig_lab <- if (p_dt < 0.05) "*" else "o"
        sig_col <- if (p_dt < 0.05) clim_col else "grey55"
        p <- p + annotate(
          "text", x = -Inf, y = Inf, label = sig_lab,
          hjust = -0.20, vjust = 1.30, size = 2.2, colour = sig_col
        )
      }

      y_lab <- if (show_y) paste0(b, "\nOR (95% CI)") else NULL
      x_lab <- x_axis_labels[[ci$label]]
      p <- p + labs(x = x_lab, y = y_lab)

      if (!show_y) {
        p <- p + theme(
          axis.text.y = element_text(colour = "transparent"),
          axis.ticks.y = element_line(colour = "transparent"),
          axis.title.y = element_blank()
        )
      }
      if (!show_x) {
        p <- p + theme(
          axis.text.x = element_text(colour = "transparent"),
          axis.ticks.x = element_line(colour = "transparent"),
          axis.title.x = element_text(colour = "transparent")
        )
      }

      panels[[paste0(b, "_", ci$orig)]] <- p
    }
  }

  make_header <- function(label) {
    col <- climate_colors[[label]]
    ggplot() +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = col) +
      annotate("text", x = 0, y = 0, label = label, colour = "white", fontface = "bold", size = 2.5) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  headers <- lapply(climate_specs, function(ci) make_header(ci$label))

  legend_grob <- get_legend(
    ggplot(
      tibble(
        x = c(1, 2, 3, 4),
        y = c(1, 1, 1, 1),
        lt = factor(
          c("Primary (Model C)", "Primary (Model C)", "Common detrended", "Common detrended"),
          levels = c("Primary (Model C)", "Common detrended")
        )
      ),
      aes(x, y, linetype = lt, group = lt)
    ) +
      geom_line(linewidth = 0.7, colour = "black") +
      scale_linetype_manual(
        values = c("Primary (Model C)" = "solid", "Common detrended" = "dashed"),
        name = NULL
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 6.5, family = "Helvetica"),
        legend.key = element_blank()
      )
  )

  footnote_grob <- ggdraw() +
    draw_label(
      "* both models significant (P < 0.05); o primary model only; rho, Spearman correlation on the shared exposure grid; tint, robustness class.",
      size = 5.5, colour = "grey35", x = 0.01, hjust = 0, fontfamily = "Helvetica"
    )

  header_row <- plot_grid(plotlist = headers, ncol = 4, align = "h")
  panel_order <- unlist(lapply(bacteria_levels, function(b) {
    paste0(b, "_", vapply(climate_specs, `[[`, character(1), "orig"))
  }))
  body_grid <- plot_grid(plotlist = panels[panel_order], ncol = 4, align = "hv", axis = "tblr")

  final_plot <- plot_grid(
    header_row, body_grid, legend_grob, footnote_grob,
    ncol = 1, rel_heights = c(0.03, 1, 0.03, 0.028)
  )

  ggsave(file.path(output_root, "02_figures", "detrended_time_series_curves.pdf"), final_plot, width = fig_width_in, height = fig_height_ds1)
  ggsave(file.path(output_root, "02_figures", "detrended_time_series_curves.png"), final_plot, width = fig_width_in, height = fig_height_ds1, dpi = 600)
}

build_fig_ds2 <- function(stats_primary_common, shape_all) {
  bacteria_rev <- rev(bacteria_levels)
  rob_lkp <- build_robustness_lookup(stats_primary_common, shape_all)

  sig_label <- function(p) {
    case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.10 ~ ".",
      TRUE ~ ""
    )
  }

  df_a <- stats_primary_common %>%
    mutate(
      sig = sig_label(p_value),
      cell_lab = paste0(sprintf("%.1f", EDF), sig),
      Bacteria = factor(Bacteria, levels = bacteria_levels),
      Model = factor(Model, levels = c("Primary", "Common-DT"), labels = c("Primary model", "Common detrended"))
    )

  f_max <- max(df_a$F_stat, na.rm = TRUE)

  p_a <- ggplot(df_a, aes(x = Climate, y = Bacteria)) +
    geom_tile(aes(fill = F_stat), colour = "white", linewidth = 0.5) +
    geom_text(aes(label = cell_lab), size = 2.3, colour = "grey10") +
    facet_wrap(~Model, ncol = 2) +
    scale_fill_gradientn(
      colours = c("#F7F7F7", "#FEE08B", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
      values = rescale(c(0, 1, 3, 6, 12, f_max)),
      limits = c(0, f_max),
      name = "F-statistic"
    ) +
    scale_y_discrete(limits = bacteria_rev) +
    scale_x_discrete(labels = c("Temperature" = "Temp.", "Humidity" = "Humid.", "Precipitation" = "Precip.", "Wet Days" = "Wet\nDays")) +
    labs(
      x = NULL,
      y = NULL,
      title = "(a) F-statistic and EDF: Primary vs Common detrended",
      subtitle = "Cell labels show EDF and significance"
    ) +
    theme_nature(base_size = 7) +
    theme(
      strip.text = element_text(face = "bold", size = 7.5),
      axis.text.x = element_text(size = 6.5)
    )

  dc <- rob_lkp %>%
    mutate(
      rho_lab = if_else(!is.na(rho), sprintf("rho=%.2f", rho), "NA"),
      Bacteria = factor(Bacteria, levels = bacteria_levels)
    )

  text_colour_map <- c(
    "Retained" = "white",
    "Direction changed" = "white",
    "Attenuated" = "grey15",
    "Lost" = "white",
    "Gained" = "white",
    "Concordant ns" = "grey25"
  )

  p_b <- ggplot(dc, aes(x = Climate, y = Bacteria)) +
    geom_tile(aes(fill = robustness), colour = "white", linewidth = 0.5) +
    geom_text(aes(label = rho_lab, colour = robustness), size = 2.1, show.legend = FALSE) +
    scale_fill_manual(values = robustness_pal, name = "Robustness", drop = FALSE) +
    scale_colour_manual(values = text_colour_map) +
    scale_y_discrete(limits = bacteria_rev) +
    scale_x_discrete(labels = c("Temperature" = "Temp.", "Humidity" = "Humid.", "Precipitation" = "Precip.", "Wet Days" = "Wet\nDays")) +
    labs(
      x = NULL,
      y = NULL,
      title = "(b) Robustness of detrended dose-response shapes"
    ) +
    theme_nature(base_size = 7) +
    theme(axis.text.x = element_text(size = 6.5))

  rob_counts <- dc %>%
    count(robustness) %>%
    mutate(
      pct = n / sum(n) * 100,
      bar_lab = sprintf("%d (%.0f%%)", n, pct),
      robustness = factor(robustness, levels = rev(levels(rob_lkp$robustness)))
    )

  p_c <- ggplot(rob_counts, aes(x = robustness, y = n, fill = robustness)) +
    geom_col(width = 0.62, show.legend = FALSE) +
    geom_text(aes(label = bar_lab), hjust = -0.10, size = 2.3, colour = "grey20", family = "Helvetica") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = robustness_pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(
      x = NULL,
      y = "Number of associations",
      title = "(c) Overall robustness summary",
      subtitle = sprintf("Total = %d pathogen x climate associations", sum(rob_counts$n))
    ) +
    theme_nature(base_size = 7) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "black", linewidth = 0.3),
      axis.text.y = element_text(size = 6.5),
      plot.subtitle = element_text(size = 6, colour = "grey40")
    )

  bottom_row <- plot_grid(p_b, p_c, ncol = 2, align = "h", rel_widths = c(1.15, 0.85))
  final_plot <- plot_grid(p_a, bottom_row, ncol = 1, rel_heights = c(1, 0.68))

  ggsave(file.path(output_root, "02_figures", "detrended_time_series_summary.pdf"), final_plot, width = fig_width_in, height = fig_height_ds2)
  ggsave(file.path(output_root, "02_figures", "detrended_time_series_summary.png"), final_plot, width = fig_width_in, height = fig_height_ds2, dpi = 600)

  write.csv(
    dc %>% select(Bacteria, Climate, EDF, F_stat, p_value, rho, dir_match, robustness),
    file.path(output_root, "01_tables", "detrended_time_series_robustness.csv"),
    row.names = FALSE
  )
}

build_fig_ds3 <- function(diag_all) {
  df <- diag_all %>%
    mutate(
      Variable = factor(
        Variable,
        levels = c("logit_R", "TMP", "PREC", "HUM", "WET"),
        labels = c("logit(R)", "Temp.", "Precip.", "Humid.", "Wet Days")
      ),
      Bacteria = factor(Bacteria, levels = bacteria_levels)
    ) %>%
    pivot_longer(
      cols = c(t_Before, t_After),
      names_to = "Stage",
      values_to = "t_stat"
    ) %>%
    mutate(
      Stage = factor(
        if_else(Stage == "t_Before", "Before detrending", "After detrending"),
        levels = c("Before detrending", "After detrending")
      )
    )

  p <- ggplot(df, aes(x = Variable, y = t_stat, fill = Stage)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.30) +
    geom_hline(yintercept = 1.96, colour = "grey60", linewidth = 0.25, linetype = "dotted") +
    geom_hline(yintercept = -1.96, colour = "grey60", linewidth = 0.25, linetype = "dotted") +
    geom_col(position = position_dodge(0.72), width = 0.65, alpha = 0.88) +
    facet_wrap(~Bacteria, ncol = 3) +
    scale_fill_manual(
      values = c("Before detrending" = "#4393C3", "After detrending" = "#D6604D"),
      name = NULL
    ) +
    labs(
      x = NULL,
      y = "Year-trend t-statistic (pooled fixed-effects regression)",
      title = "Detrending effectiveness: year-trend t-statistic before vs after common detrending",
      subtitle = "Dotted lines show +/-1.96; near-zero t values after detrending indicate successful trend removal."
    ) +
    theme_nature(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
      strip.text = element_text(face = "bold", size = 7),
      legend.position = "bottom",
      legend.margin = margin(t = 2),
      plot.subtitle = element_text(size = 6, colour = "grey40")
    )

  ggsave(file.path(output_root, "02_figures", "common_detrending_diagnostics.pdf"), p, width = fig_width_in, height = fig_height_ds3)
  ggsave(file.path(output_root, "02_figures", "common_detrending_diagnostics.png"), p, width = fig_width_in, height = fig_height_ds3, dpi = 600)
}

# 4. Main execution
cat("[START]\n")

all_results <- list()
all_stats <- list()
all_shapes <- list()
all_diags <- list()
model_manifest_rows <- list()

for (spec in bacteria_specs) {
  bacteria <- spec$title
  file_path <- file.path(input_data_dir, spec$file_name)
  if (!file.exists(file_path)) {
    warning("Missing file for ", bacteria, ": ", file_path)
    next
  }

  cat("\n", strrep("=", 60), "\n ", bacteria, "\n", strrep("=", 60), "\n", sep = "")

  lag_config <- lag_settings[[bacteria]]
  raw_data <- read_base_data(file_path)
  scale_params <- compute_scale_params(raw_data)

  cat("  [A] Primary primary model ...\n")
  primary_data <- prepare_model_data(raw_data, lag_config, use_year_factor = FALSE)
  primary_model <- fit_gamm(primary_data, use_year_factor = FALSE)
  primary_stats <- extract_stats(primary_model, bacteria, "Primary")

  cat("  [B] Common linear detrending ...\n")
  common_raw <- apply_common_detrend(raw_data)
  common_data <- prepare_model_data(common_raw, lag_config, use_year_factor = TRUE)
  common_model <- fit_gamm(common_data, use_year_factor = TRUE)
  common_stats <- extract_stats(common_model, bacteria, "Common-DT")

  cat("  [C] Country-specific linear detrending ...\n")
  country_raw <- apply_country_detrend(raw_data, min_obs = min_obs_country)
  country_data <- prepare_model_data(country_raw, lag_config, use_year_factor = TRUE)
  country_model <- fit_gamm(country_data, use_year_factor = TRUE)
  country_stats <- extract_stats(country_model, bacteria, "Country-DT")

  cat("  [D] Shape consistency ...\n")
  shape <- compute_shape(primary_model, common_model, primary_data, common_data, scale_params)
  all_shapes[[bacteria]] <- shape

  cat("  [E] Common-trend diagnostics ...\n")
  diagnostics <- compute_common_diagnostics(raw_data, bacteria)
  all_diags[[bacteria]] <- diagnostics

  model_manifest_rows[[length(model_manifest_rows) + 1L]] <- tibble(
    Bacteria = bacteria,
    input_file = file_path,
    n_primary = nrow(primary_data),
    n_common_dt = nrow(common_data),
    n_country_dt = nrow(country_data),
    TMP_lag = lag_config$temp_lag,
    PREC_lag = lag_config$precip_lag,
    HUM_lag = lag_config$humid_lag,
    WET_lag = lag_config$wetdays_lag,
    PLS_components_used = paste(get_available_pls_components(primary_data), collapse = ", "),
    primary_formula = paste(deparse(formula(primary_model)), collapse = " "),
    common_dt_formula = paste(deparse(formula(common_model)), collapse = " "),
    country_dt_formula = paste(deparse(formula(country_model)), collapse = " ")
  )

  stats_combined <- bind_rows(primary_stats, common_stats, country_stats)
  all_stats[[bacteria]] <- stats_combined

  all_results[[bacteria]] <- list(
    lag_config = lag_config,
    scale_params = scale_params,
    primary_data = primary_data,
    common_data = common_data,
    country_data = country_data,
    primary_model = primary_model,
    common_model = common_model,
    country_model = country_model
  )
}

stats_all <- bind_rows(all_stats) %>%
  mutate(
    Bacteria = factor(Bacteria, levels = bacteria_levels),
    Model = factor(Model, levels = c("Primary", "Common-DT", "Country-DT"))
  ) %>%
  arrange(Bacteria, Climate, Model)

diagnostics_all <- bind_rows(all_diags) %>%
  mutate(Bacteria = factor(Bacteria, levels = bacteria_levels))

shape_summary <- bind_rows(lapply(names(all_shapes), function(b) {
  shape <- all_shapes[[b]]
  if (is.null(shape)) return(NULL)
  tibble(
    Bacteria = b,
    Climate = factor(names(shape), levels = climate_levels),
    rho_primary_vs_common = sapply(shape, `[[`, "rho"),
    direction_match_primary_vs_common = sapply(shape, `[[`, "dir_match")
  )
})) %>%
  mutate(Bacteria = factor(Bacteria, levels = bacteria_levels))

model_manifest <- bind_rows(model_manifest_rows) %>%
  mutate(Bacteria = factor(Bacteria, levels = bacteria_levels)) %>%
  arrange(Bacteria)

write.csv(stats_all, file.path(output_root, "01_tables", "detrending_sensitivity_summary.csv"), row.names = FALSE)
write.csv(diagnostics_all, file.path(output_root, "01_tables", "detrending_diagnostics_common_linear.csv"), row.names = FALSE)
write.csv(shape_summary, file.path(output_root, "01_tables", "detrending_shape_consistency_primary_vs_common.csv"), row.names = FALSE)
write.csv(model_manifest, file.path(output_root, "04_metadata", "detrending_model_manifest.csv"), row.names = FALSE)

build_fig_ds1(all_results, all_shapes, stats_all %>% filter(Model %in% c("Primary", "Common-DT")))
build_fig_ds2(stats_all %>% filter(Model %in% c("Primary", "Common-DT")), all_shapes)
build_fig_ds3(diagnostics_all)

robustness_table <- read.csv(
  file.path(output_root, "01_tables", "detrended_time_series_robustness.csv"),
  stringsAsFactors = FALSE
)

wb_path <- file.path(output_root, "03_workbook", "detrended_time_series.xlsx")
write_xlsx(
  list(
    Summary = stats_all,
    Shape_Consistency = shape_summary,
    Common_Diagnostics = diagnostics_all,
    Robustness = robustness_table,
    Model_Manifest = model_manifest
  ),
  wb_path
)

write.csv(stats_all, file.path(source_data_root, "01_csv", "Detrending_Sensitivity_summary.csv"), row.names = FALSE)
write.csv(shape_summary, file.path(source_data_root, "01_csv", "Detrending_Shape_Consistency.csv"), row.names = FALSE)
write.csv(diagnostics_all, file.path(source_data_root, "01_csv", "Detrending_Common_Diagnostics.csv"), row.names = FALSE)
write.csv(robustness_table, file.path(source_data_root, "01_csv", "Detrending_Robustness_Table.csv"), row.names = FALSE)
write.csv(model_manifest, file.path(source_data_root, "01_csv", "Detrending_Model_Manifest.csv"), row.names = FALSE)
write_xlsx(
  list(
    Summary = stats_all,
    Shape_Consistency = shape_summary,
    Common_Diagnostics = diagnostics_all,
    Robustness = robustness_table,
    Model_Manifest = model_manifest
  ),
  file.path(source_data_root, "02_workbook", "SourceData_detrended_time_series.xlsx")
)

method_note <- c(
  "analysis of climate-AMR associations using detrended time series",
  "",
  "Primary model:",
  "- Full main Model C aligned to the current historical analysis pipeline.",
  "- Uses model_ready_inputs and public lag-selection source data for the primary historical model.",
  "- Dynamically includes all available PLS components (3 or 4 depending on phenotype).",
  "",
  "Common detrending:",
  "- Removes a shared linear year trend from logit_R and climate variables using year + country fixed effects.",
  "- Refits the same model structure but replaces s(year) with year_factor to avoid reintroducing the common trend.",
  "",
  "Country-specific detrending:",
  "- Removes within-country linear trends from logit_R and climate variables.",
  "- Uses the same lag settings and same broad adjustment structure as the detrended sensitivity model.",
  "",
  "Shape diagnostics and figures:",
  "- Response-shape comparisons are computed from raw GAMM term predictions.",
  "- No loess or post-estimation curve smoothing is used in the robustness metrics.",
  "- The common-detrending diagnostics figure reports pooled fixed-effects year-trend t-statistics before and after common detrending.",
  "",
  "Current optimal lags:",
  paste(
    apply(model_manifest[, c("Bacteria", "TMP_lag", "PREC_lag", "HUM_lag", "WET_lag")], 1, function(x) {
      sprintf("%s: TMP=%s, PREC=%s, HUM=%s, WET=%s", x[1], x[2], x[3], x[4], x[5])
    }),
    collapse = "\n"
  )
)
writeLines(method_note, file.path(output_root, "04_metadata", "README_detrending_sensitivity.txt"))

cat("\n", strrep("=", 60), "\n", sep = "")
cat("DONE\n")
cat(output_root, "\n")
cat(strrep("=", 60), "\n")
