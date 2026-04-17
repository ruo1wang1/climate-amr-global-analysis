######## Figure 1 Figure 1 source data export: primary model ########

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(openxlsx)
})

Sys.setenv(MODELC_SKIP_MAIN = "1")

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "analysis_historical_associations_main_model.R"), local = FALSE)

fig1_output_suffix <- Sys.getenv("MODELC_FIG1SRC_OUTPUT_SUFFIX", unset = "")
fig1_run_only <- trimws(Sys.getenv("MODELC_FIG1SRC_RUN_ONLY", unset = ""))

fit_primary_model_for_source_data <- function(data, bacteria_name) {
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
  available_pls <- get_available_pls_components(data)

  if (length(available_pls) == 0) {
    stop("No available PLS components found for ", bacteria_name, call. = FALSE)
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
    warning("bam() failed for ", bacteria_name, "; trying gam(): ", e$message)
    gam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE
    )
  })

  list(model = model, formula = formula_str, available_pls = available_pls)
}

build_density_source_data <- function(data, orig_var, bacteria_name) {
  is_precipitation <- orig_var == "PREC"
  is_humidity <- orig_var == "HUM"
  is_wetdays <- orig_var == "WET"

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

  if (is_precipitation) {
    pos_data <- values[values > 0]
    if (length(pos_data) > 10) {
      dens <- density(pos_data, na.rm = TRUE, adjust = 1.1, from = 0, to = var_max)
      density_data <- tibble(x = dens$x, density = dens$y)
    } else {
      breaks <- seq(0, var_max, length.out = 30)
      hist_data <- hist(pos_data, breaks = breaks, plot = FALSE)
      density_data <- tibble(x = hist_data$mids, density = hist_data$density)
    }
  } else if (is_wetdays) {
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
  density_data <- density_data %>%
    mutate(
      scaled_density = if (is.finite(max_density) && max_density > 0) density / max_density * 0.46 else 0
    )
}

extract_panel_b_stats <- function(model, bacteria_name) {
  s_table <- as.data.frame(summary(model)$s.table)
  s_table$Term <- rownames(s_table)
  rownames(s_table) <- NULL

  climate_map <- c(
    "s(TMP_scaled_lag)" = "Temperature",
    "s(HUM_scaled_lag)" = "Humidity",
    "s(PREC_scaled_lag)" = "Precipitation",
    "s(WET_scaled_lag)" = "WetDays"
  )

  s_table %>%
    filter(Term %in% names(climate_map)) %>%
    mutate(
      Bacteria = bacteria_name,
      Climate_Variable = unname(climate_map[Term]),
      Lag_Years = vapply(
        Climate_Variable,
        function(x) {
          if (x == "Temperature") return(lag_settings[[bacteria_name]]$temp_lag)
          if (x == "Humidity") return(lag_settings[[bacteria_name]]$humid_lag)
          if (x == "Precipitation") return(lag_settings[[bacteria_name]]$precip_lag)
          lag_settings[[bacteria_name]]$wetdays_lag
        },
        numeric(1)
      ),
      F_statistic = `F`,
      EDF = edf,
      Ref_df = Ref.df,
      P_Value = `p-value`,
      Significance = case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01 ~ "**",
        P_Value < 0.05 ~ "*",
        P_Value < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    ) %>%
    select(Bacteria, Climate_Variable, Lag_Years, Term, EDF, Ref_df, F_statistic, P_Value, Significance)
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
      x_value = ifelse(abs(x_orig - round(x_orig)) < 1e-10,
                       as.character(round(x_orig)),
                       format(round(x_orig, 3), nsmall = 3)),
      or_value = ifelse(abs(y - round(y)) < 1e-10,
                        as.character(round(y)),
                        format(round(y, 3), nsmall = 3)),
      Figure_Label = paste0(label_short, " (", x_value, unit_suffix, ")\nOR = ", or_value, sig_symbol)
    )
}

export_fig1_source_data <- function(specs_to_run, output_root) {
  curves_all <- list()
  densities_all <- list()
  thresholds_all <- list()
  panel_b_all <- list()
  metadata_rows <- list()

  climate_factors <- list(
    list(var = "TMP_scaled_lag", original = "TMP", climate = "Temperature"),
    list(var = "HUM_scaled_lag", original = "HUM", climate = "Humidity"),
    list(var = "PREC_scaled_lag", original = "PREC", climate = "Precipitation"),
    list(var = "WET_scaled_lag", original = "WET", climate = "WetDays")
  )

  for (spec in specs_to_run) {
    cat("Exporting Figure 1 source data for", spec$title, "...\n")
    prepared <- prepare_data(file.path(input_data_dir, spec$file_name), spec$title)
    data_ready <- prepared$data
    scale_params <- prepared$scale_params
    model_obj <- fit_primary_model_for_source_data(data_ready, spec$title)

    panel_b_all[[length(panel_b_all) + 1]] <- extract_panel_b_stats(model_obj$model, spec$title)
    metadata_rows[[length(metadata_rows) + 1]] <- data.frame(
      Bacteria = spec$title,
      Input_File = file.path(input_data_dir, spec$file_name),
      n_samples = nrow(data_ready),
      PLS_Components = paste(model_obj$available_pls, collapse = ", "),
      Model_Formula = model_obj$formula,
      stringsAsFactors = FALSE
    )

    for (factor in climate_factors) {
      threshold_data <- detect_thresholds(model_obj$model, data_ready, scale_params, factor$var, spec$title)
      lag_years <- get_lag_value(spec$title, factor$var)

      curves_all[[length(curves_all) + 1]] <- threshold_data$curve_data %>%
        transmute(
          Bacteria = spec$title,
          Climate_Variable = factor$climate,
          Lag_Years = lag_years,
          x_value = x_orig,
          x_scaled,
          OR = y,
          Lower_95CI = lower_ci,
          Upper_95CI = upper_ci,
          Relationship_Type = threshold_data$relationship_type
        )

      densities_all[[length(densities_all) + 1]] <- build_density_source_data(data_ready, factor$original, spec$title) %>%
        transmute(
          Bacteria = spec$title,
          Climate_Variable = factor$climate,
          Lag_Years = lag_years,
          x_value = x,
          Density = density,
          Scaled_Density = scaled_density
        )

      thresholds_all[[length(thresholds_all) + 1]] <- build_threshold_labels(
        threshold_data$threshold_points,
        factor$climate
      ) %>%
        mutate(
          Bacteria = spec$title,
          Climate_Variable = factor$climate,
          Lag_Years = lag_years
        )
    }
  }

  curves_df <- bind_rows(curves_all)
  density_df <- bind_rows(densities_all)
  thresholds_df <- bind_rows(thresholds_all)
  panel_b_df <- bind_rows(panel_b_all)
  metadata_df <- bind_rows(metadata_rows)

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  csv_dir <- file.path(output_root, "01_csv")
  workbook_dir <- file.path(output_root, "02_workbook")
  for (dir_path in c(csv_dir, workbook_dir)) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }

  write.csv(curves_df, file.path(csv_dir, "Figure1_ModelC_Fig1A_curves.csv"), row.names = FALSE)
  write.csv(density_df, file.path(csv_dir, "Figure1_ModelC_Fig1A_density.csv"), row.names = FALSE)
  write.csv(thresholds_df, file.path(csv_dir, "Figure1_ModelC_Fig1A_thresholds.csv"), row.names = FALSE)
  write.csv(panel_b_df, file.path(csv_dir, "Figure1_ModelC_Fig1B_summary.csv"), row.names = FALSE)
  write.csv(metadata_df, file.path(csv_dir, "Figure1_ModelC_metadata.csv"), row.names = FALSE)

  workbook_path <- file.path(workbook_dir, "SourceData_figure1_historical_associations.xlsx")
  wb <- createWorkbook()
  addWorksheet(wb, "Fig1A_curves")
  addWorksheet(wb, "Fig1A_density")
  addWorksheet(wb, "Fig1A_thresholds")
  addWorksheet(wb, "Fig1B_summary")
  addWorksheet(wb, "Metadata")
  addWorksheet(wb, "README")

  writeData(wb, "Fig1A_curves", curves_df)
  writeData(wb, "Fig1A_density", density_df)
  writeData(wb, "Fig1A_thresholds", thresholds_df)
  writeData(wb, "Fig1B_summary", panel_b_df)
  writeData(wb, "Metadata", metadata_df)
  writeData(
    wb,
    "README",
    data.frame(
      Item = c(
        "Figure",
        "Panels covered",
        "Model",
        "Input data",
        "Lag source",
        "How to use"
      ),
      Value = c(
        "Figure 1",
        "Fig1A response curves, density panels, and threshold annotations; Fig1B summary statistics",
        "primary model with four climate factors and dynamic PLS components",
        input_data_dir,
        lag_summary_path,
        "This workbook is intended as figure-specific source data for journal submission. The plotted curve and density values correspond to the data actually used to generate the revised Figure 1."
      ),
      stringsAsFactors = FALSE
    )
  )
  saveWorkbook(wb, workbook_path, overwrite = TRUE)

  readme_lines <- c(
    "# Figure 1 Source Data",
    "",
    "This folder contains a publication-style source data workbook for primary model Figure 1.",
    "",
    "Sheets in the workbook:",
    "- `Fig1A_curves`",
    "- `Fig1A_density`",
    "- `Fig1A_thresholds`",
    "- `Fig1B_summary`",
    "- `Metadata`",
    "- `README`"
  )
  writeLines(readme_lines, file.path(output_root, "README.md"))

  list(
    workbook = workbook_path,
    curves_csv = file.path(csv_dir, "Figure1_ModelC_Fig1A_curves.csv"),
    density_csv = file.path(csv_dir, "Figure1_ModelC_Fig1A_density.csv"),
    thresholds_csv = file.path(csv_dir, "Figure1_ModelC_Fig1A_thresholds.csv"),
    panel_b_csv = file.path(csv_dir, "Figure1_ModelC_Fig1B_summary.csv")
  )
}

main <- function() {
  specs_to_run <- bacteria_specs
  if (nzchar(fig1_run_only)) {
    run_only_values <- trimws(strsplit(fig1_run_only, ",")[[1]])
    specs_to_run <- Filter(function(x) x$code %in% run_only_values || x$title %in% run_only_values, specs_to_run)
  }

  if (length(specs_to_run) == 0) {
    stop("No bacteria matched MODELC_FIG1SRC_RUN_ONLY.", call. = FALSE)
  }

  output_root <- file.path(
    revision_root,
    "data/source_data",
    paste0("figure1_historical_associations", fig1_output_suffix)
  )

  output_paths <- export_fig1_source_data(specs_to_run, output_root)
  cat("Figure 1 source data saved to:\n")
  print(output_paths)

  invisible(output_paths)
}

invisible(main())
