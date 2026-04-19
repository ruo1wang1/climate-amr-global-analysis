suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(openxlsx)
  library(flextable)
  library(officer)
  library(cowplot)
  library(patchwork)
  library(scales)
})

set.seed(20260408)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
input_data_dir <- file.path(revision_root, "outputs/historical_associations", "model_ready_inputs")
helper_script_path <- file.path(
  revision_root,
  "code",
  "01_historical_associations",
  "analysis_historical_associations_main_model.R"
)

output_suffix <- Sys.getenv("SPECMODELFX_OUTPUT_SUFFIX", unset = "")
output_root_override <- trimws(Sys.getenv("SPECMODELFX_OUTPUT_ROOT", unset = ""))
source_data_root_override <- trimws(Sys.getenv("SPECMODELFX_SOURCE_DATA_ROOT", unset = ""))
run_only <- trimws(Sys.getenv("SPECMODELFX_RUN_ONLY", unset = ""))

default_output_root <- file.path(
  revision_root,
  "outputs",
  paste0("climate_specification_checks", output_suffix)
)
default_source_data_root <- file.path(
  revision_root,
  "data/source_data",
  paste0("climate_variable_specification", output_suffix)
)

output_root <- if (nzchar(output_root_override)) output_root_override else default_output_root
source_data_root <- if (nzchar(source_data_root_override)) source_data_root_override else default_source_data_root

model_ready_input_name <- function(code) {
  paste0(code, "_model_ready_data.csv")
}

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = model_ready_input_name("3GCREC"), pls_components = 3),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = model_ready_input_name("3GCRKP"), pls_components = 4),
  list(code = "CRAB", title = "CR-Ab", file_name = model_ready_input_name("CRAB"), pls_components = 4),
  list(code = "CREC", title = "CR-Ec", file_name = model_ready_input_name("CREC"), pls_components = 3),
  list(code = "CRKP", title = "CR-Kp", file_name = model_ready_input_name("CRKP"), pls_components = 3),
  list(code = "CRPA", title = "CR-Pa", file_name = model_ready_input_name("CRPA"), pls_components = 4)
)

if (nzchar(run_only)) {
  run_only_values <- trimws(strsplit(run_only, ",")[[1]])
  bacteria_specs <- Filter(
    function(x) x$code %in% run_only_values || x$title %in% run_only_values,
    bacteria_specs
  )
}

if (length(bacteria_specs) == 0) {
  stop("No bacteria matched SPECMODELFX_RUN_ONLY.", call. = FALSE)
}

bacteria_order <- vapply(bacteria_specs, `[[`, character(1), "title")
variable_table_order <- c("Humidity", "Precipitation", "Temperature", "Wet Days")
variable_heatmap_order <- c("HUM", "PREC", "TMP", "WET")
model_order <- c("Model A", "Model B", "Model C")
lag_name_map <- c(TMP = "temp_lag", PREC = "precip_lag", HUM = "humid_lag", WET = "wetdays_lag")

model_specs <- list(
  A = list(
    key = "A",
    short_label = "Model A",
    display_label = "Model A (TMP, PREC, HUM)",
    climate_vars = c("TMP", "PREC", "HUM"),
    climate_labels = c(TMP = "Temperature", PREC = "Precipitation", HUM = "Humidity"),
    lag_summary_path = file.path(
      revision_root,
      "data",
      "source_data",
      "lag_selection",
      "historical_model_a",
      "01_csv",
      "historical_lag_summary_model_a.csv"
    ),
    figure_id = "A",
    figure_basename = "model_a_climate_response",
    panel_b_factor_order = c("Temperature", "Humidity", "Precipitation"),
    climate_factors = list(
      list(var = "TMP_scaled_lag", original = "TMP", climate = "Temperature", x_lab = "Temperature (°C)", color = "#DD5F60"),
      list(var = "HUM_scaled_lag", original = "HUM", climate = "Humidity", x_lab = "Relative humidity (%)", color = "#3CB371"),
      list(var = "PREC_scaled_lag", original = "PREC", climate = "Precipitation", x_lab = "Precipitation (mm)", color = "#9AC0CD")
    )
  ),
  B = list(
    key = "B",
    short_label = "Model B",
    display_label = "Model B (TMP, PREC, WET)",
    climate_vars = c("TMP", "PREC", "WET"),
    climate_labels = c(TMP = "Temperature", PREC = "Precipitation", WET = "Wet Days"),
    lag_summary_path = file.path(
      revision_root,
      "data",
      "source_data",
      "lag_selection",
      "historical_model_b",
      "01_csv",
      "historical_lag_summary_model_b.csv"
    ),
    figure_id = "B",
    figure_basename = "model_b_climate_response",
    panel_b_factor_order = c("Temperature", "Wet Days", "Precipitation"),
    climate_factors = list(
      list(var = "TMP_scaled_lag", original = "TMP", climate = "Temperature", x_lab = "Temperature (°C)", color = "#DD5F60"),
      list(var = "WET_scaled_lag", original = "WET", climate = "Wet Days", x_lab = "Wet Days (d)", color = "#8A2BE2"),
      list(var = "PREC_scaled_lag", original = "PREC", climate = "Precipitation", x_lab = "Precipitation (mm)", color = "#9AC0CD")
    )
  ),
  C = list(
    key = "C",
    short_label = "Model C",
    display_label = "Model C (TMP, PREC, HUM, WET)",
    climate_vars = c("TMP", "PREC", "HUM", "WET"),
    climate_labels = c(TMP = "Temperature", PREC = "Precipitation", HUM = "Humidity", WET = "Wet Days"),
    lag_summary_path = file.path(
      revision_root,
      "data",
      "source_data",
      "lag_selection",
      "historical_model_c",
      "01_csv",
      "historical_lag_summary_model_c.csv"
    )
  )
)

debug_print <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s][%s] %s\n", timestamp, level, message))
}

ensure_dirs <- function(paths) {
  for (dir_path in paths) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

sanitize_file_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) x <- "file"
  x
}

sanitize_sheet_name <- function(x) {
  x <- gsub("[\\\\/:*?\\[\\]]", "_", x)
  x <- substr(x, 1, 31)
  if (!nzchar(x)) x <- "Sheet"
  x
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(0)
  }
  mean(x)
}

format_p_value_display <- function(p_value) {
  ifelse(
    is.na(p_value),
    "",
    ifelse(
      p_value < 0.001,
      "<0.001",
      formatC(p_value, digits = 5, format = "fg", flag = "#")
    )
  )
}

compute_significance <- function(p_value) {
  case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ ".",
    TRUE ~ ""
  )
}

get_available_pls_components <- function(data) {
  pls_candidates <- paste0("PLS_Comp", 1:4)
  present_pls <- pls_candidates[pls_candidates %in% names(data)]
  present_pls[vapply(present_pls, function(x) !all(is.na(data[[x]])), logical(1))]
}

check_factor_levels <- function(data) {
  factor_vars <- names(data)[vapply(data, is.factor, logical(1))]
  problems <- character(0)
  for (var in factor_vars) {
    if (length(levels(data[[var]])) < 2) {
      problems <- c(problems, var)
    }
  }
  problems
}

compute_vif_values <- function(data, climate_vars) {
  vars <- paste0(climate_vars, "_scaled_lag")
  x <- data[, vars, drop = FALSE]
  x <- x[complete.cases(x), , drop = FALSE]
  out <- setNames(rep(NA_real_, length(vars)), vars)

  if (nrow(x) < 10) {
    return(out)
  }

  zero_var <- vapply(x, function(col) sd(col) == 0, logical(1))
  if (any(zero_var)) {
    x <- x[, !zero_var, drop = FALSE]
  }

  if (ncol(x) < 2) {
    return(out)
  }

  for (var in names(x)) {
    others <- setdiff(names(x), var)
    if (length(others) == 0) {
      out[var] <- 1
      next
    }
    fit <- try(lm(reformulate(others, response = var), data = x), silent = TRUE)
    if (inherits(fit, "try-error")) next
    r2 <- summary(fit)$r.squared
    if (is.na(r2)) next
    out[var] <- if (r2 >= 0.999999) Inf else 1 / (1 - r2)
  }

  out
}

load_modelc_helpers <- function(helper_path) {
  if (!file.exists(helper_path)) {
    stop("Helper script not found: ", helper_path, call. = FALSE)
  }

  helper_env <- new.env(parent = globalenv())
  old_skip <- Sys.getenv("MODELC_SKIP_MAIN", unset = NA_character_)

  on.exit({
    if (is.na(old_skip)) {
      Sys.unsetenv("MODELC_SKIP_MAIN")
    } else {
      Sys.setenv(MODELC_SKIP_MAIN = old_skip)
    }
  }, add = TRUE)

  Sys.setenv(MODELC_SKIP_MAIN = "1")
  source(helper_path, local = helper_env)
  helper_env
}

helper_env <- load_modelc_helpers(helper_script_path)

prepare_base_data <- function(file_path) {
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
      climate_zone = factor(climate_zone, levels = c("Tropical Zone", "Temperate Zone", "Polar Zone"))
    ) %>%
    group_by(NAME) %>%
    mutate(location_id = cur_group_id()) %>%
    ungroup() %>%
    mutate(HUM = pmin(HUM, 100))

  scale_params <- data_processed %>%
    summarise(
      across(
        c(TMP, PREC, HUM, WET),
        list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))
      )
    )

  base_data <- data_processed %>%
    mutate(
      TMP_orig = TMP,
      PREC_orig = PREC,
      HUM_orig = HUM,
      WET_orig = WET
    ) %>%
    group_by(climate_zone) %>%
    mutate(
      across(c(TMP, PREC, HUM, WET), \(x) as.numeric(scale(x)), .names = "{.col}_scaled")
    ) %>%
    ungroup()

  list(base_data = base_data, scale_params = scale_params)
}

read_best_lag_settings <- function(model_spec) {
  if (!file.exists(model_spec$lag_summary_path)) {
    stop("Lag summary not found for ", model_spec$short_label, ": ", model_spec$lag_summary_path, call. = FALSE)
  }

  lag_df <- read.csv(model_spec$lag_summary_path, stringsAsFactors = FALSE) %>%
    mutate(
      is_best_flag = if ("is_best" %in% names(.)) {
        is_best %in% c(TRUE, "TRUE", 1, "1")
      } else {
        FALSE
      }
    ) %>%
    filter(is_best_flag | rank == 1) %>%
    group_by(Display_Name) %>%
    arrange(AIC, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()

  lag_settings <- list()
  for (i in seq_len(nrow(lag_df))) {
    row <- lag_df[i, ]
    lag_settings[[row$Display_Name]] <- list(
      temp_lag = if (!is.na(row$TMP_lag)) as.integer(row$TMP_lag) else NA_integer_,
      precip_lag = if (!is.na(row$PREC_lag)) as.integer(row$PREC_lag) else NA_integer_,
      humid_lag = if (!is.na(row$HUM_lag)) as.integer(row$HUM_lag) else NA_integer_,
      wetdays_lag = if (!is.na(row$WET_lag)) as.integer(row$WET_lag) else NA_integer_
    )
  }

  lag_settings
}

lag_settings_by_model <- lapply(model_specs, read_best_lag_settings)

create_final_dataset <- function(base_data, lag_config, climate_vars) {
  data_lagged <- base_data %>%
    group_by(location_id) %>%
    arrange(year, .by_group = TRUE)

  for (var_name in names(lag_name_map)) {
    lag_field <- lag_name_map[[var_name]]
    lag_value <- lag_config[[lag_field]]
    if (is.null(lag_value) || is.na(lag_value) || lag_value < 1) {
      lag_value <- 1L
    }
    data_lagged <- data_lagged %>%
      mutate(!!paste0(var_name, "_scaled_lag") := dplyr::lag(.data[[paste0(var_name, "_scaled")]], n = lag_value))
  }

  required_cols <- paste0(climate_vars, "_scaled_lag")
  data_lagged <- data_lagged %>%
    filter(if_all(all_of(required_cols), ~ !is.na(.x))) %>%
    ungroup() %>%
    droplevels()

  optional_vars <- setdiff(names(lag_name_map), climate_vars)
  for (var_name in optional_vars) {
    lagged_col <- paste0(var_name, "_scaled_lag")
    replacement_value <- safe_mean(data_lagged[[lagged_col]])
    data_lagged[[lagged_col]][is.na(data_lagged[[lagged_col]])] <- replacement_value
  }

  data_lagged
}

build_final_formula <- function(data, climate_vars) {
  factor_problems <- check_factor_levels(data)
  available_pls <- get_available_pls_components(data)

  if (length(available_pls) == 0) {
    stop("No available PLS components found in final dataset.", call. = FALSE)
  }

  climate_terms <- vapply(
    climate_vars,
    function(var_name) {
      k_value <- if (var_name == "TMP") 5 else 10
      sprintf("s(%s_scaled_lag, k = %d, bs = 'cr')", var_name, k_value)
    },
    character(1)
  )

  pls_terms <- vapply(
    available_pls,
    function(var_name) sprintf("s(%s, k = 10, bs = 'cr')", var_name),
    character(1)
  )

  formula_parts <- c(climate_terms, pls_terms)

  if (nrow(data) > 30) {
    formula_parts <- c(formula_parts, "s(lat, lon, bs = 'sos', k = 20)")
  }

  year_levels <- length(unique(data$year))
  if (year_levels > 2) {
    year_k <- min(8, year_levels - 1)
    formula_parts <- c(formula_parts, sprintf("s(year, bs = 'cr', k = %d)", year_k))
  } else if (year_levels == 2) {
    formula_parts <- c(formula_parts, "year")
  }

  if ("Region" %in% names(data) && !("Region" %in% factor_problems)) {
    formula_parts <- c(formula_parts, "s(Region, bs = 're')")
  }

  if ("climate_zone" %in% names(data) && !("climate_zone" %in% factor_problems)) {
    formula_parts <- c(formula_parts, "climate_zone")
  }

  list(
    formula_str = paste("logit_R ~", paste(formula_parts, collapse = " + ")),
    available_pls = available_pls
  )
}

fit_final_model <- function(data, climate_vars, bacteria_name) {
  model_setup <- build_final_formula(data, climate_vars)
  model_formula <- as.formula(model_setup$formula_str)
  ctrl <- gam.control(nthreads = 1, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)

  model_obj <- tryCatch({
    bam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE,
      control = ctrl
    )
  }, error = function(e) {
    warning("bam() failed for ", bacteria_name, "; trying gam(): ", e$message)
    gam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE
    )
  })

  list(
    model = model_obj,
    formula = model_setup$formula_str,
    available_pls = model_setup$available_pls
  )
}

extract_smooth_stats <- function(model, bacteria_name, model_spec, lag_config) {
  s_table <- as.data.frame(summary(model)$s.table)
  s_table$Term <- rownames(s_table)
  rownames(s_table) <- NULL

  climate_map <- setNames(unname(model_spec$climate_labels[model_spec$climate_vars]), paste0("s(", model_spec$climate_vars, "_scaled_lag)"))

  s_table %>%
    filter(Term %in% names(climate_map)) %>%
    mutate(
      Bacterial_Species = bacteria_name,
      Model = model_spec$short_label,
      Climate_Variable = unname(climate_map[Term]),
      Lag_Years = vapply(
        names(climate_map)[match(Term, names(climate_map))],
        function(term_name) {
          var_name <- sub("^s\\((.*)_scaled_lag\\)$", "\\1", term_name)
          lag_field <- lag_name_map[[var_name]]
          as.integer(lag_config[[lag_field]])
        },
        integer(1)
      ),
      Effective_DF = edf,
      Ref_df = Ref.df,
      F_value = `F`,
      p_value = `p-value`,
      Significance = compute_significance(p_value)
    ) %>%
    select(
      Bacterial_Species,
      Model,
      Climate_Variable,
      Lag_Years,
      Term,
      Effective_DF,
      Ref_df,
      F_value,
      p_value,
      Significance
    ) %>%
    mutate(
      Climate_Variable = factor(Climate_Variable, levels = variable_table_order),
      Model = factor(Model, levels = model_order),
      Bacterial_Species = factor(Bacterial_Species, levels = bacteria_order)
    ) %>%
    arrange(Bacterial_Species, Model, Climate_Variable) %>%
    mutate(
      Climate_Variable = as.character(Climate_Variable),
      Model = as.character(Model),
      Bacterial_Species = as.character(Bacterial_Species)
    )
}

build_density_source_data <- function(data, orig_var, bacteria_name) {
  is_precipitation <- orig_var == "PREC"
  is_humidity <- orig_var == "HUM"
  is_wetdays <- orig_var == "WET"

  get_range_setting <- function(var_type, bacteria) {
    setting_key <- ifelse(var_type == "PREC", "precip", ifelse(var_type == "HUM", "humid", ifelse(var_type == "WET", "wetdays", "temp")))
    if (bacteria %in% names(helper_env$range_settings[[setting_key]])) {
      helper_env$range_settings[[setting_key]][[bacteria]]
    } else {
      helper_env$range_settings[[setting_key]][["default"]]
    }
  }

  var_range <- get_range_setting(orig_var, bacteria_name)
  var_min <- var_range[1]
  var_max <- var_range[2]
  x_range <- var_max - var_min
  margin_size <- x_range * 0.04
  var_min <- var_min - margin_size
  var_max <- var_max + margin_size

  values <- data[[paste0(orig_var, "_orig")]]

  if (is_humidity) {
    values <- pmin(values, 100)
    var_max <- min(104, var_max)
  }

  if (is_precipitation || is_wetdays) {
    var_min <- max(-margin_size, var_min)
  }

  if (is_precipitation || is_wetdays) {
    pos_data <- values[values > 0]
    if (length(pos_data) > 10) {
      dens <- density(pos_data, na.rm = TRUE, adjust = 1.1, from = 0, to = var_max)
      density_data <- tibble(x = dens$x, density = dens$y)
    } else {
      breaks <- seq(0, var_max, length.out = 30)
      hist_data <- hist(pos_data, breaks = breaks, plot = FALSE)
      density_data <- tibble(x = hist_data$mids, density = hist_data$density)
    }
  } else {
    if (length(values) > 10) {
      dens <- density(values, na.rm = TRUE, adjust = 1.1)
      density_data <- tibble(x = dens$x, density = dens$y)
    } else {
      breaks <- seq(var_min, var_max, length.out = 30)
      hist_data <- hist(values, breaks = breaks, plot = FALSE)
      density_data <- tibble(x = hist_data$mids, density = hist_data$density)
    }
  }

  density_data <- density_data %>%
    filter(x >= var_min & x <= var_max)

  if (is_humidity) {
    density_data <- density_data %>% filter(x <= 100)
  }

  max_density <- max(density_data$density, na.rm = TRUE)
  density_data %>%
    mutate(
      scaled_density = if (is.finite(max_density) && max_density > 0) density / max_density * 0.46 else 0
    )
}

build_threshold_labels <- function(threshold_points, climate_variable) {
  if (nrow(threshold_points) == 0) {
    return(threshold_points)
  }

  unit_suffix <- if (climate_variable == "Temperature") {
    "°C"
  } else if (climate_variable == "Humidity") {
    "%"
  } else if (climate_variable == "Precipitation") {
    "mm"
  } else {
    "d"
  }

  threshold_points %>%
    mutate(
      label_short = case_when(
        grepl("Global maximum", type) ~ "GMax",
        grepl("Global minimum", type) ~ "GMin",
        grepl("Maximum", type) ~ "Max",
        grepl("Minimum", type) ~ "Min",
        TRUE ~ substr(type, 1, 3)
      ),
      x_value = ifelse(abs(x_orig - round(x_orig)) < 1e-10, as.character(round(x_orig)), format(round(x_orig, 3), nsmall = 3)),
      or_value = ifelse(abs(y - round(y)) < 1e-10, as.character(round(y)), format(round(y, 3), nsmall = 3)),
      Figure_Label = paste0(label_short, " (", x_value, unit_suffix, ")\nOR = ", or_value, sig_symbol)
    )
}

prepare_panel_b_plot_data <- function(smooth_df, factor_order) {
  smooth_df %>%
    transmute(
      Bacteria = factor(Bacterial_Species, levels = bacteria_order),
      Factor = factor(Climate_Variable, levels = rev(factor_order)),
      Effective_DF,
      F_value,
      Significance,
      label_text = sprintf("F=%.2f%s\nEDF=%.2f", F_value, Significance, Effective_DF)
    ) %>%
    arrange(Bacteria, Factor)
}

make_panel_b_plot <- function(plot_data, factor_order, panel_label = "B") {
  factor_colors <- c(
    "Temperature" = "#DD5F60",
    "Humidity" = "#3CB371",
    "Precipitation" = "#9AC0CD",
    "Wet Days" = "#8A2BE2"
  )
  factor_shapes <- c(
    "Temperature" = 16,
    "Humidity" = 17,
    "Precipitation" = 18,
    "Wet Days" = 15
  )

  max_f <- max(plot_data$F_value, na.rm = TRUE)
  x_limit <- max(15, ceiling(max_f / 5) * 5 + 5)
  x_breaks <- seq(0, x_limit, by = 5)
  label_nudge <- max(1.6, x_limit * 0.055)
  tick_color <- "gray60"

  base_plot <- ggplot(plot_data, aes(x = F_value, y = Factor)) +
    theme_bw(base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray96", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 11, color = "gray25"),
      axis.text.y = element_text(size = 12, face = "bold", color = "gray25"),
      strip.text = element_text(size = 13, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      legend.position = "right",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      panel.spacing = unit(0.25, "lines"),
      legend.box = "vertical",
      legend.key.size = unit(1.15, "lines"),
      legend.key = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", linewidth = 0.8),
      axis.ticks.length = unit(0.20, "cm"),
      axis.ticks = element_line(color = tick_color),
      axis.ticks.x = element_line(color = tick_color),
      axis.ticks.y = element_line(color = tick_color),
      plot.margin = margin(t = 8, r = 14, b = 8, l = 8)
    ) +
    geom_point(
      aes(size = pmin(Effective_DF, 6), shape = Factor, color = Factor),
      alpha = 0.9
    ) +
    geom_text(
      aes(label = label_text),
      nudge_x = label_nudge,
      hjust = 0,
      size = 3.9,
      family = "serif",
      lineheight = 1.05
    ) +
    facet_wrap(~Bacteria, nrow = 1, strip.position = "top") +
    scale_color_manual(values = factor_colors[factor_order], name = "Climate Factor") +
    scale_shape_manual(values = factor_shapes[factor_order], name = "Climate Factor") +
    scale_size_continuous(
      name = "EDF Value",
      breaks = c(0, 2, 4, 6),
      labels = c("0", "2", "4", "6"),
      range = c(2.4, 6.6),
      limits = c(0, max(6, max(plot_data$Effective_DF, na.rm = TRUE)))
    ) +
    scale_x_continuous(
      limits = c(0, x_limit),
      breaks = x_breaks,
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    labs(x = "F value") +
    guides(
      shape = guide_legend(order = 1, override.aes = list(size = 5, alpha = 1, color = factor_colors[factor_order])),
      size = guide_legend(order = 2),
      color = "none"
    )

  ggdraw() +
    draw_label(panel_label, x = 0.018, y = 0.92, hjust = 0, vjust = 1, fontfamily = "serif", size = 34) +
    draw_plot(base_plot, x = 0.08, y = 0, width = 0.92, height = 1)
}

save_plot_pair <- function(plot_obj, pdf_path, png_path, width, height) {
  ggsave(pdf_path, plot = plot_obj, width = width, height = height, dpi = 300, device = grDevices::cairo_pdf)
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = 300)
}

save_model_summary <- function(model_fit, bacteria_name, model_spec, lag_config, output_dir) {
  summary_path <- file.path(
    output_dir,
    paste0(model_spec$key, "_", sanitize_file_stem(bacteria_name), "_model_summary.txt")
  )

  sink(summary_path)
  cat("====================================\n")
  cat("Model Specification Sensitivity Final Model Summary\n")
  cat("====================================\n\n")
  cat("Bacteria:", bacteria_name, "\n")
  cat("Model:", model_spec$display_label, "\n\n")
  cat("Lag settings used:\n")
  cat("Temperature lag:", lag_config$temp_lag, "years\n")
  cat("Precipitation lag:", lag_config$precip_lag, "years\n")
  cat("Humidity lag:", lag_config$humid_lag, "years\n")
  cat("Wet days lag:", lag_config$wetdays_lag, "years\n\n")
  cat("PLS components included:\n")
  cat(paste(model_fit$available_pls, collapse = ", "), "\n\n")
  cat("Model formula:\n")
  cat(model_fit$formula, "\n\n")
  cat("Model summary:\n")
  print(summary(model_fit$model))
  cat("\nModel validation statistics:\n")
  print(gam.check(model_fit$model))
  cat("\n====================================\n")
  sink()

  summary_path
}

fit_all_final_models <- function(metadata_dir) {
  model_summary_dir <- file.path(metadata_dir, "model_summaries")
  ensure_dirs(model_summary_dir)

  fitted_results <- list()
  smooth_rows <- list()
  vif_rows <- list()
  manifest_rows <- list()

  for (bacteria in bacteria_specs) {
    debug_print(paste("Preparing final-model datasets for", bacteria$title))
    base_prepared <- prepare_base_data(file.path(input_data_dir, bacteria$file_name))
    base_data <- base_prepared$base_data
    scale_params <- base_prepared$scale_params

    fitted_results[[bacteria$title]] <- list()

    for (model_key in names(model_specs)) {
      model_spec <- model_specs[[model_key]]
      lag_config <- lag_settings_by_model[[model_key]][[bacteria$title]]
      if (is.null(lag_config)) {
        stop("Missing lag settings for ", bacteria$title, " / ", model_spec$short_label, call. = FALSE)
      }

      debug_print(paste("  Fitting", model_spec$short_label, "for", bacteria$title))
      final_data <- create_final_dataset(base_data, lag_config, model_spec$climate_vars)
      model_fit <- fit_final_model(final_data, model_spec$climate_vars, bacteria$title)
      summary_path <- save_model_summary(model_fit, bacteria$title, model_spec, lag_config, model_summary_dir)

      fitted_results[[bacteria$title]][[model_key]] <- list(
        data = final_data,
        scale_params = scale_params,
        model = model_fit$model,
        formula = model_fit$formula,
        available_pls = model_fit$available_pls,
        lag_config = lag_config,
        summary_path = summary_path,
        model_spec = model_spec
      )

      smooth_rows[[length(smooth_rows) + 1]] <- extract_smooth_stats(model_fit$model, bacteria$title, model_spec, lag_config)

      vif_values <- compute_vif_values(final_data, model_spec$climate_vars)
      vif_named <- c(
        TMP = unname(vif_values["TMP_scaled_lag"]),
        PREC = unname(vif_values["PREC_scaled_lag"]),
        HUM = unname(vif_values["HUM_scaled_lag"]),
        WET = unname(vif_values["WET_scaled_lag"])
      )
      for (var_name in variable_heatmap_order) {
        vif_rows[[length(vif_rows) + 1]] <- tibble(
          Bacterial_Species = bacteria$title,
          Model = model_spec$short_label,
          Climate_Code = var_name,
          Climate_Variable = case_when(
            var_name == "TMP" ~ "Temperature",
            var_name == "PREC" ~ "Precipitation",
            var_name == "HUM" ~ "Humidity",
            TRUE ~ "Wet Days"
          ),
          VIF = vif_named[[var_name]]
        )
      }

      manifest_rows[[length(manifest_rows) + 1]] <- tibble(
        Bacteria = bacteria$title,
        Model = model_spec$short_label,
        Display_Model = model_spec$display_label,
        Input_File = file.path(input_data_dir, bacteria$file_name),
        Sample_Size = nrow(final_data),
        Temp_Lag = lag_config$temp_lag,
        Prec_Lag = lag_config$precip_lag,
        Hum_Lag = lag_config$humid_lag,
        Wet_Lag = lag_config$wetdays_lag,
        Available_PLS = paste(model_fit$available_pls, collapse = ", "),
        Formula = model_fit$formula,
        AIC = AIC(model_fit$model),
        Explained_Deviance_pct = summary(model_fit$model)$dev.expl * 100,
        Max_VIF = suppressWarnings(max(vif_values, na.rm = TRUE)),
        Summary_Path = summary_path
      )
    }
  }

  list(
    fitted_results = fitted_results,
    smooth_df = bind_rows(smooth_rows),
    vif_df = bind_rows(vif_rows),
    manifest_df = bind_rows(manifest_rows)
  )
}

build_climate_specification_summary_table <- function(smooth_df) {
  smooth_df %>%
    mutate(
      Bacterial_Species = factor(Bacterial_Species, levels = bacteria_order),
      Model = factor(Model, levels = model_order),
      Climate_Variable = factor(Climate_Variable, levels = variable_table_order)
    ) %>%
    arrange(Bacterial_Species, Model, Climate_Variable) %>%
    mutate(
      Bacterial_Species = as.character(Bacterial_Species),
      Model = as.character(Model),
      Climate_Variable = as.character(Climate_Variable)
    )
}

write_climate_specification_summary_outputs <- function(summary_table, vif_df, manifest_df, output_root) {
  tables_dir <- file.path(output_root, "01_tables")
  workbook_dir <- file.path(output_root, "02_workbook")
  doc_dir <- file.path(output_root, "03_docx")
  metadata_dir <- file.path(output_root, "06_metadata")
  ensure_dirs(c(output_root, tables_dir, workbook_dir, doc_dir, metadata_dir))

  csv_path <- file.path(tables_dir, "climate_specification_summary.csv")
  xlsx_path <- file.path(workbook_dir, "climate_specification_summary.xlsx")
  docx_path <- file.path(doc_dir, "climate_specification_summary.docx")
  manifest_path <- file.path(metadata_dir, "climate_specification_model_manifest.csv")

  write.csv(summary_table, csv_path, row.names = FALSE)
  write.csv(manifest_df, manifest_path, row.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
  addWorksheet(wb, "VIF_Values")
  addWorksheet(wb, "Model_Manifest")
  writeData(wb, "Summary", summary_table)
  writeData(wb, "VIF_Values", vif_df)
  writeData(wb, "Model_Manifest", manifest_df)
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  docx_table <- summary_table %>%
    transmute(
      `Bacterial Species` = Bacterial_Species,
      Model,
      `Climate Variable` = Climate_Variable,
      `Effective DF` = sprintf("%.3f", Effective_DF),
      `F-value` = sprintf("%.3f", F_value),
      `p-value` = format_p_value_display(p_value),
      Significance
    )

  ft <- flextable(docx_table)
  ft <- theme_booktabs(ft)
  ft <- autofit(ft)
  ft <- add_header_lines(
    ft,
    values = "climate specification comparison summary. Climate Variable Effects Across Three Models"
  )
  ft <- add_footer_lines(
    ft,
    values = c(
      "Significance codes: '***' p<0.001, '**' p<0.01, '*' p<0.05, '.' p<0.1",
      "Note: Smooth-term statistics were extracted from refitted GAMMs using the selected lag combinations from model-structure comparison."
    )
  )
  save_as_docx(ft, path = docx_path)

  invisible(list(csv = csv_path, xlsx = xlsx_path, docx = docx_path, manifest = manifest_path))
}

build_vif_heatmap <- function(vif_df) {
  heatmap_df <- expand.grid(
    Bacterial_Species = bacteria_order,
    Model = rev(model_order),
    Climate_Code = variable_heatmap_order,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(vif_df, by = c("Bacterial_Species", "Model", "Climate_Code")) %>%
    mutate(
      Climate_Code = factor(Climate_Code, levels = variable_heatmap_order),
      Model = factor(Model, levels = rev(model_order)),
      Bacterial_Species = factor(Bacterial_Species, levels = bacteria_order)
    )

  max_vif <- max(heatmap_df$VIF, na.rm = TRUE)
  max_vif <- ifelse(is.finite(max_vif), max(5, max_vif), 5)

  ggplot(heatmap_df, aes(x = Climate_Code, y = Model, fill = VIF)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(
      aes(label = ifelse(is.na(VIF), "", sprintf("%.2f", VIF))),
      color = "black",
      size = 3.2
    ) +
    facet_wrap(~Bacterial_Species, ncol = 3) +
    scale_fill_gradientn(
      colours = c("#176c2f", "#7ba63f", "#c6d400", "#fff431"),
      values = rescale(c(1, 2, 3, max_vif)),
      limits = c(1, max_vif),
      na.value = "white",
      name = "VIF Value"
    ) +
    labs(
      x = "Climate Variable",
      y = "Model Structure"
    ) +
    theme_bw(base_family = "serif") +
    theme(
      axis.text = element_text(size = 11, color = "gray25"),
      axis.title = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray55", linewidth = 0.8),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

write_vif_outputs <- function(vif_df, output_root, source_data_root) {
  figures_dir <- file.path(output_root, "04_figures")
  metadata_dir <- file.path(output_root, "06_metadata")
  ensure_dirs(c(figures_dir, metadata_dir, source_data_root))

  pdf_path <- file.path(figures_dir, "climate_collinearity_heatmap.pdf")
  png_path <- file.path(figures_dir, "climate_collinearity_heatmap.png")
  csv_path <- file.path(metadata_dir, "climate_collinearity_heatmap_values.csv")
  source_csv_dir <- file.path(source_data_root, "climate_collinearity_heatmap", "01_csv")
  source_wb_dir <- file.path(source_data_root, "climate_collinearity_heatmap", "02_workbook")
  ensure_dirs(c(source_csv_dir, source_wb_dir))

  plot_obj <- build_vif_heatmap(vif_df)
  save_plot_pair(plot_obj, pdf_path, png_path, width = 11.5, height = 7.8)
  write.csv(vif_df, csv_path, row.names = FALSE)
  write.csv(vif_df, file.path(source_csv_dir, "climate_collinearity_heatmap_values.csv"), row.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "VIF_Heatmap")
  writeData(wb, "VIF_Heatmap", vif_df)
  saveWorkbook(
    wb,
    file.path(source_wb_dir, "source_data_climate_collinearity_heatmap.xlsx"),
    overwrite = TRUE
  )

  invisible(list(pdf = pdf_path, png = png_path, csv = csv_path))
}

build_model_curve_figure <- function(model_key, fitted_results, output_root, source_data_root) {
  model_spec <- model_specs[[model_key]]
  helper_env$lag_settings <- lag_settings_by_model[[model_key]]

  figures_dir <- file.path(output_root, "04_figures")
  metadata_dir <- file.path(output_root, "06_metadata")
  model_metadata_dir <- file.path(metadata_dir, model_spec$figure_basename)
  source_model_root <- file.path(source_data_root, model_spec$figure_basename)
  source_csv_dir <- file.path(source_model_root, "01_csv")
  source_wb_dir <- file.path(source_model_root, "02_workbook")
  ensure_dirs(c(figures_dir, model_metadata_dir, source_csv_dir, source_wb_dir))

  row_plots <- list()
  curves_rows <- list()
  density_rows <- list()
  threshold_rows <- list()
  linear_rows <- list()

  panel_b_input <- list()

  for (bacteria in bacteria_specs) {
    fit_obj <- fitted_results[[bacteria$title]][[model_key]]
    factor_plots <- list()
    panel_b_input[[length(panel_b_input) + 1]] <- extract_smooth_stats(
      fit_obj$model,
      bacteria$title,
      model_spec,
      fit_obj$lag_config
    )

    for (factor in model_spec$climate_factors) {
      debug_print(paste("Building", model_spec$short_label, factor$climate, "curve for", bacteria$title))
      threshold_data <- helper_env$detect_thresholds(
        fit_obj$model,
        fit_obj$data,
        fit_obj$scale_params,
        factor$var,
        bacteria$title
      )

      factor_plot <- helper_env$create_climate_effect_plot(
        fit_obj$model,
        fit_obj$data,
        fit_obj$scale_params,
        factor$var,
        factor$x_lab,
        bacteria$title,
        factor$color,
        bacteria$title,
        threshold_data
      )
      factor_plots[[factor$climate]] <- factor_plot

      curves_rows[[length(curves_rows) + 1]] <- threshold_data$curve_data %>%
        transmute(
          Bacteria = bacteria$title,
          Model = model_spec$short_label,
          Climate_Variable = factor$climate,
          Lag_Years = threshold_data$lag_years,
          x_value = x_orig,
          x_scaled,
          OR = y,
          Lower_95CI = lower_ci,
          Upper_95CI = upper_ci,
          Relationship_Type = threshold_data$relationship_type
        )

      density_rows[[length(density_rows) + 1]] <- build_density_source_data(
        fit_obj$data,
        factor$original,
        bacteria$title
      ) %>%
        transmute(
          Bacteria = bacteria$title,
          Model = model_spec$short_label,
          Climate_Variable = factor$climate,
          x_value = x,
          Density = density,
          Scaled_Density = scaled_density
        )

      threshold_rows[[length(threshold_rows) + 1]] <- build_threshold_labels(
        threshold_data$threshold_points,
        factor$climate
      ) %>%
        mutate(
          Bacteria = bacteria$title,
          Model = model_spec$short_label,
          Climate_Variable = factor$climate
        )

      linear_rows[[length(linear_rows) + 1]] <- tibble(
        Bacteria = bacteria$title,
        Model = model_spec$short_label,
        Climate_Variable = factor$climate,
        Lag_Years = threshold_data$lag_years,
        Relationship_Type = threshold_data$relationship_type,
        Slope = threshold_data$linear_info$slope,
        P_Value = threshold_data$linear_info$p_value,
        R_Squared = threshold_data$linear_info$r_squared
      )
    }

    row_plots[[length(row_plots) + 1]] <- plot_grid(
      plotlist = factor_plots,
      ncol = length(model_spec$climate_factors),
      align = "h",
      axis = "tb"
    )
  }

  panel_a_grid <- plot_grid(plotlist = row_plots, ncol = 1, align = "v")
  panel_a_plot <- ggdraw() +
    draw_label("A", x = 0.005, y = 0.995, hjust = 0, vjust = 1, fontfamily = "serif", size = 34) +
    draw_plot(panel_a_grid, x = 0.03, y = 0, width = 0.97, height = 1)

  panel_b_df <- bind_rows(panel_b_input) %>%
    filter(Model == model_spec$short_label)
  panel_b_plot_data <- prepare_panel_b_plot_data(panel_b_df, model_spec$panel_b_factor_order)
  panel_b_plot <- make_panel_b_plot(panel_b_plot_data, model_spec$panel_b_factor_order, panel_label = "B")

  final_plot <- plot_grid(
    panel_a_plot,
    panel_b_plot,
    ncol = 1,
    rel_heights = c(6.8, 1.4)
  )

  pdf_path <- file.path(figures_dir, paste0(model_spec$figure_basename, ".pdf"))
  png_path <- file.path(figures_dir, paste0(model_spec$figure_basename, ".png"))
  save_plot_pair(final_plot, pdf_path, png_path, width = 15, height = 18.8)

  curves_df <- bind_rows(curves_rows)
  density_df <- bind_rows(density_rows)
  thresholds_df <- bind_rows(threshold_rows)
  linear_df <- bind_rows(linear_rows)

  write.csv(curves_df, file.path(model_metadata_dir, paste0(model_spec$figure_basename, "_curves.csv")), row.names = FALSE)
  write.csv(density_df, file.path(model_metadata_dir, paste0(model_spec$figure_basename, "_density.csv")), row.names = FALSE)
  write.csv(thresholds_df, file.path(model_metadata_dir, paste0(model_spec$figure_basename, "_thresholds.csv")), row.names = FALSE)
  write.csv(panel_b_df, file.path(model_metadata_dir, paste0(model_spec$figure_basename, "_panelB_stats.csv")), row.names = FALSE)
  write.csv(linear_df, file.path(model_metadata_dir, paste0(model_spec$figure_basename, "_linear_relationships.csv")), row.names = FALSE)

  write.csv(curves_df, file.path(source_csv_dir, paste0(model_spec$figure_basename, "_curves.csv")), row.names = FALSE)
  write.csv(density_df, file.path(source_csv_dir, paste0(model_spec$figure_basename, "_density.csv")), row.names = FALSE)
  write.csv(thresholds_df, file.path(source_csv_dir, paste0(model_spec$figure_basename, "_thresholds.csv")), row.names = FALSE)
  write.csv(panel_b_df, file.path(source_csv_dir, paste0(model_spec$figure_basename, "_panelB_stats.csv")), row.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "Curves")
  addWorksheet(wb, "Density")
  addWorksheet(wb, "Thresholds")
  addWorksheet(wb, "PanelB")
  writeData(wb, "Curves", curves_df)
  writeData(wb, "Density", density_df)
  writeData(wb, "Thresholds", thresholds_df)
  writeData(wb, "PanelB", panel_b_df)
  saveWorkbook(
    wb,
    file.path(source_wb_dir, paste0("source_data_", model_spec$figure_basename, ".xlsx")),
    overwrite = TRUE
  )

  invisible(
    list(
      pdf = pdf_path,
      png = png_path,
      curves = curves_df,
      density = density_df,
      thresholds = thresholds_df,
      panel_b = panel_b_df,
      linear = linear_df
    )
  )
}

main <- function() {
  tables_dir <- file.path(output_root, "01_tables")
  workbook_dir <- file.path(output_root, "02_workbook")
  doc_dir <- file.path(output_root, "03_docx")
  figures_dir <- file.path(output_root, "04_figures")
  source_dir <- file.path(output_root, "05_source_data")
  metadata_dir <- file.path(output_root, "06_metadata")
  ensure_dirs(c(output_root, tables_dir, workbook_dir, doc_dir, figures_dir, source_dir, metadata_dir, source_data_root))

  debug_print("Fitting GAMMs for Models A, B, and C using the selected lag settings")
  all_results <- fit_all_final_models(metadata_dir)

  summary_table <- build_climate_specification_summary_table(all_results$smooth_df)
  debug_print("Writing climate specification comparison summary outputs")
  summary_paths <- write_climate_specification_summary_outputs(summary_table, all_results$vif_df, all_results$manifest_df, output_root)

  debug_print("Writing climate collinearity figure outputs")
  vif_paths <- write_vif_outputs(all_results$vif_df, output_root, source_data_root)

  debug_print("Writing model A climate-response figure outputs")
  model_a_paths <- build_model_curve_figure("A", all_results$fitted_results, output_root, source_data_root)

  debug_print("Writing model B climate-response figure outputs")
  model_b_paths <- build_model_curve_figure("B", all_results$fitted_results, output_root, source_data_root)

  manifest_summary <- tibble(
    Output = c(
      "summary_csv",
      "summary_xlsx",
      "summary_docx",
      "climate_collinearity_pdf",
      "climate_collinearity_png",
      "model_a_pdf",
      "model_a_png",
      "model_b_pdf",
      "model_b_png"
    ),
    Path = c(
      summary_paths$csv,
      summary_paths$xlsx,
      summary_paths$docx,
      vif_paths$pdf,
      vif_paths$png,
      model_a_paths$pdf,
      model_a_paths$png,
      model_b_paths$pdf,
      model_b_paths$png
    )
  )
  write.csv(manifest_summary, file.path(metadata_dir, "climate_specification_output_manifest.csv"), row.names = FALSE)

  debug_print("Completed climate specification summary and figure generation")
  invisible(
    list(
      summary = summary_paths,
      climate_collinearity = vif_paths,
      model_a = model_a_paths,
      model_b = model_b_paths
    )
  )
}

invisible(main())
