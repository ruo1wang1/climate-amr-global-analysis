######## Analysis of climate-AMR associations under alternative climatic-variable specifications ########

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(openxlsx)
  library(flextable)
  library(officer)
})

set.seed(20260407)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
input_data_dir <- file.path(revision_root, "outputs/historical_associations", "model_ready_inputs")
existing_modelc_root <- file.path(
  revision_root,
  "data",
  "source_data",
  "lag_selection",
  "historical_model_c",
  "01_csv"
)

output_suffix <- Sys.getenv("SPECMODEL_OUTPUT_SUFFIX", unset = "")
output_root_override <- trimws(Sys.getenv("SPECMODEL_OUTPUT_ROOT", unset = ""))
run_only <- trimws(Sys.getenv("SPECMODEL_RUN_ONLY", unset = ""))
run_models <- trimws(Sys.getenv("SPECMODEL_RUN_MODELS", unset = "A,B,C"))
reuse_existing_modelc <- Sys.getenv("SPECMODEL_REUSE_MODELC", unset = "0") != "0"
max_lag_to_run <- suppressWarnings(as.integer(Sys.getenv("SPECMODEL_MAX_LAG", unset = "3")))
cv_scope <- trimws(Sys.getenv("SPECMODEL_CV_SCOPE", unset = "best_only"))

if (is.na(max_lag_to_run) || max_lag_to_run < 1 || max_lag_to_run > 3) {
  stop("SPECMODEL_MAX_LAG must be an integer between 1 and 3.", call. = FALSE)
}

if (!cv_scope %in% c("best_only", "all", "none")) {
  stop("SPECMODEL_CV_SCOPE must be one of: best_only, all, none.", call. = FALSE)
}

default_output_base <- file.path(
  revision_root,
  "outputs"
)
output_base <- if (nzchar(output_root_override)) output_root_override else default_output_base

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
  stop("No bacteria matched SPECMODEL_RUN_ONLY.", call. = FALSE)
}

model_specs <- list(
  A = list(
    key = "A",
    short_label = "Model A",
    display_label = "Model A (TMP, PREC, HUM)",
    climate_vars = c("TMP", "PREC", "HUM"),
    omitted_vars = c("WET"),
    output_root = file.path(output_base, "historical_associations", paste0("model_a_lag_selection", output_suffix))
  ),
  B = list(
    key = "B",
    short_label = "Model B",
    display_label = "Model B (TMP, PREC, WET)",
    climate_vars = c("TMP", "PREC", "WET"),
    omitted_vars = c("HUM"),
    output_root = file.path(output_base, "historical_associations", paste0("model_b_lag_selection", output_suffix))
  ),
  C = list(
    key = "C",
    short_label = "Model C",
    display_label = "Model C (TMP, PREC, HUM, WET)",
    climate_vars = c("TMP", "PREC", "HUM", "WET"),
    omitted_vars = character(0),
    output_root = file.path(output_base, "historical_associations", "model_c_lag_selection")
  )
)

run_model_values <- trimws(strsplit(run_models, ",")[[1]])
run_model_values <- run_model_values[nzchar(run_model_values)]
model_specs <- model_specs[names(model_specs) %in% run_model_values]

if (length(model_specs) == 0) {
  stop("No model matched SPECMODEL_RUN_MODELS. Use a comma-separated subset of A,B,C.", call. = FALSE)
}

model_structure_output_root <- file.path(
  output_base,
  paste0("climate_variable_specification", output_suffix)
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

sanitize_sheet_name <- function(x) {
  x <- gsub("[\\\\/:*?\\[\\]]", "_", x)
  x <- substr(x, 1, 31)
  if (!nzchar(x)) x <- "Sheet"
  x
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

build_safe_formula <- function(data, pls_components, climate_vars) {
  factor_problems <- check_factor_levels(data)

  climate_terms <- paste(
    vapply(
      climate_vars,
      function(var_name) sprintf("s(%s_scaled_lag, k = 10, bs = 'cr')", var_name),
      character(1)
    ),
    collapse = " + "
  )

  available_pls <- paste0("PLS_Comp", seq_len(pls_components))
  available_pls <- available_pls[available_pls %in% names(data)]
  available_pls <- available_pls[vapply(available_pls, function(x) !all(is.na(data[[x]])), logical(1))]

  pls_terms <- if (length(available_pls) > 0) {
    paste(
      vapply(available_pls, function(x) sprintf("s(%s, k = 10, bs = 'cr')", x), character(1)),
      collapse = " + "
    )
  } else {
    NULL
  }

  formula_parts <- c(climate_terms, pls_terms)

  if (nrow(data) > 30) {
    spatial_k <- max(10, min(20, floor(nrow(data) / 10)))
    formula_parts <- c(formula_parts, sprintf("s(lat, lon, bs = 'sos', k = %d)", spatial_k))
  }

  year_levels <- length(unique(data$year))
  if (year_levels > 2) {
    year_k <- min(8, year_levels - 1)
    formula_parts <- c(formula_parts, sprintf("s(year, bs = 'cr', k = %d)", year_k))
  } else if (year_levels == 2) {
    formula_parts <- c(formula_parts, "year")
  }

  if (!("Region" %in% factor_problems) && "Region" %in% names(data)) {
    formula_parts <- c(formula_parts, "s(Region, bs = 're')")
  }

  if (!("climate_zone" %in% factor_problems) && "climate_zone" %in% names(data)) {
    formula_parts <- c(formula_parts, "climate_zone")
  }

  paste("logit_R ~", paste(formula_parts, collapse = " + "))
}

align_test_factor_levels <- function(train_data, test_data) {
  factor_cols <- names(train_data)[vapply(train_data, is.factor, logical(1))]
  factor_cols <- intersect(factor_cols, names(test_data))

  for (col in factor_cols) {
    train_levels <- levels(train_data[[col]])
    if (length(train_levels) == 0) next

    train_values <- as.character(train_data[[col]])
    reference_level <- names(sort(table(train_values), decreasing = TRUE))[1]
    test_values <- as.character(test_data[[col]])
    unseen <- is.na(test_values) | !(test_values %in% train_levels)
    test_values[unseen] <- reference_level
    test_data[[col]] <- factor(test_values, levels = train_levels)
  }

  test_data
}

safe_rmse <- function(obs, pred) {
  keep <- is.finite(obs) & is.finite(pred)
  if (sum(keep) < 5) {
    return(NA_real_)
  }
  sqrt(mean((obs[keep] - pred[keep])^2))
}

compute_best_flags <- function(values, mode = c("min", "max"), digits = 3) {
  mode <- match.arg(mode)
  rounded_values <- round(values, digits)
  out <- rep("", length(values))
  valid <- !is.na(rounded_values)

  if (!any(valid)) {
    return(out)
  }

  target <- if (mode == "min") min(rounded_values[valid]) else max(rounded_values[valid])
  out[valid & rounded_values == target] <- "\u2713"
  out
}

compute_max_vif <- function(data, climate_vars) {
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

prepare_base_data <- function(file_path) {
  data <- read.csv(file_path)

  data %>%
    mutate(
      year = as.numeric(as.character(year)),
      Region = if ("Region" %in% names(.)) factor(Region) else factor(1),
      climate_zone = factor(
        case_when(
          abs(lat) > 66.5 ~ "Polar",
          abs(lat) > 23.5 ~ "Temperate",
          TRUE ~ "Tropical"
        ),
        levels = c("Tropical", "Temperate", "Polar")
      ),
      hemisphere = factor(if_else(lat >= 0, "North", "South"))
    ) %>%
    filter(!is.na(year), !is.na(lat), !is.na(lon), !is.na(logit_R)) %>%
    group_by(NAME) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      location_id = cur_group_id(),
      across(c(TMP, PREC, HUM, WET), ~ as.numeric(scale(.x)), .names = "{.col}_scaled")
    ) %>%
    ungroup()
}

create_lagged_dataset <- function(base_data, lag_combo) {
  lagged_data <- base_data %>%
    group_by(location_id) %>%
    arrange(year, .by_group = TRUE)

  for (var_name in names(lag_combo)) {
    lag_value <- as.integer(lag_combo[[var_name]])
    scaled_name <- paste0(var_name, "_scaled")
    lagged_name <- paste0(var_name, "_scaled_lag")
    lagged_data <- lagged_data %>%
      mutate(!!lagged_name := dplyr::lag(.data[[scaled_name]], n = lag_value))
  }

  required_cols <- paste0(names(lag_combo), "_scaled_lag")

  lagged_data %>%
    filter(if_all(all_of(required_cols), ~ !is.na(.x))) %>%
    droplevels() %>%
    ungroup()
}

fit_gamm_model <- function(data, formula_str) {
  model_formula <- as.formula(formula_str)

  tryCatch(
    bam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "fREML",
      discrete = TRUE,
      select = TRUE,
      nthreads = 1
    ),
    error = function(e) {
      warning("bam() failed, trying gam(): ", e$message)
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

make_seed_from_key <- function(key) {
  chars <- utf8ToInt(enc2utf8(key))
  if (length(chars) == 0) {
    return(20260407L)
  }
  seed <- sum(chars * seq_along(chars)) + 20260407L
  as.integer(seed %% .Machine$integer.max)
}

assign_country_folds <- function(countries, k, cv_key) {
  countries <- sort(unique(as.character(countries)))
  if (length(countries) == 0 || k < 2) {
    return(setNames(integer(0), character(0)))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(make_seed_from_key(cv_key))
  shuffled_countries <- sample(countries, length(countries), replace = FALSE)
  setNames(rep(seq_len(k), length.out = length(shuffled_countries)), shuffled_countries)
}

perform_cv <- function(data, pls_components, climate_vars, cv_key, k = 5) {
  n <- nrow(data)

  if (n < 50) {
    k <- min(k, max(2, floor(n / 10)))
  }

  countries <- unique(data$NAME)
  if (length(countries) < k) {
    k <- max(2, min(length(countries), k))
  }
  if (k < 2) {
    return(NA_real_)
  }

  country_folds <- assign_country_folds(countries, k = k, cv_key = cv_key)

  data$fold <- NA_integer_
  for (country in countries) {
    data$fold[data$NAME == country] <- country_folds[country]
  }

  cv_rmse <- numeric(k)
  cv_rmse[] <- NA_real_

  for (fold_id in seq_len(k)) {
    train_data <- data %>% filter(fold != !!fold_id) %>% droplevels()
    test_data <- data %>% filter(fold == !!fold_id) %>% droplevels()

    if (nrow(train_data) < 20 || nrow(test_data) < 5) next

    test_data <- align_test_factor_levels(train_data, test_data) %>% droplevels()
    formula_str <- build_safe_formula(train_data, pls_components, climate_vars)

    fold_model <- tryCatch(
      suppressWarnings(fit_gamm_model(train_data, formula_str)),
      error = function(e) NULL
    )

    if (is.null(fold_model)) next

    pred <- tryCatch(
      predict(fold_model, newdata = test_data),
      error = function(e) rep(NA_real_, nrow(test_data))
    )

    cv_rmse[fold_id] <- safe_rmse(test_data$logit_R, pred)
  }

  if (all(is.na(cv_rmse))) {
    return(NA_real_)
  }

  mean(cv_rmse, na.rm = TRUE)
}

generate_lag_grid <- function(climate_vars, max_lag = 3) {
  grid_values <- rep(list(seq_len(max_lag)), length(climate_vars))
  names(grid_values) <- paste0(climate_vars, "_lag")
  as_tibble(expand.grid(grid_values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
}

run_single_model_spec <- function(spec) {
  debug_print(paste("Running lag search for", spec$display_label))

  output_root <- spec$output_root
  tables_dir <- file.path(output_root, "01_tables")
  workbook_dir <- file.path(output_root, "02_workbook")
  metadata_dir <- file.path(output_root, "03_metadata")
  ensure_dirs(c(output_root, tables_dir, workbook_dir, metadata_dir))

  lag_grid <- generate_lag_grid(spec$climate_vars, max_lag = max_lag_to_run)

  all_metrics <- list()
  summary_rows <- list()

  for (bacteria in bacteria_specs) {
    debug_print(paste("  Processing", bacteria$title, "for", spec$short_label))
    file_path <- file.path(input_data_dir, bacteria$file_name)
    if (!file.exists(file_path)) {
      warning("Missing input file: ", file_path)
      next
    }

    base_data <- prepare_base_data(file_path)
    metrics_rows <- list()

    for (combo_idx in seq_len(nrow(lag_grid))) {
      lag_row <- lag_grid[combo_idx, , drop = FALSE]
      lag_combo <- as.list(lag_row)
      names(lag_combo) <- gsub("_lag$", "", names(lag_combo))

      current_data <- create_lagged_dataset(base_data, lag_combo)
      if (nrow(current_data) < max(40, bacteria$pls_components * 12)) {
        next
      }

      formula_str <- build_safe_formula(current_data, bacteria$pls_components, spec$climate_vars)
      model_fit <- tryCatch(
        fit_gamm_model(current_data, formula_str),
        error = function(e) NULL
      )

      if (is.null(model_fit)) next

      model_summary <- summary(model_fit)
      vif_values <- compute_max_vif(current_data, spec$climate_vars)
      cv_rmse <- if (cv_scope == "all") {
        perform_cv(current_data, bacteria$pls_components, spec$climate_vars, cv_key = bacteria$code)
      } else {
        NA_real_
      }

      max_vif <- suppressWarnings(max(vif_values, na.rm = TRUE))
      if (!is.finite(max_vif)) {
        max_vif <- NA_real_
      }

      metrics_rows[[length(metrics_rows) + 1]] <- tibble(
        Model_Key = spec$key,
        Model_Label = spec$display_label,
        Bacteria = bacteria$code,
        Display_Name = bacteria$title,
        combination_id = combo_idx,
        TMP_lag = if ("TMP" %in% spec$climate_vars) as.integer(lag_combo$TMP) else NA_integer_,
        PREC_lag = if ("PREC" %in% spec$climate_vars) as.integer(lag_combo$PREC) else NA_integer_,
        HUM_lag = if ("HUM" %in% spec$climate_vars) as.integer(lag_combo$HUM) else NA_integer_,
        WET_lag = if ("WET" %in% spec$climate_vars) as.integer(lag_combo$WET) else NA_integer_,
        AIC = AIC(model_fit),
        BIC = BIC(model_fit),
        Dev_explained = model_summary$dev.expl * 100,
        GCV = model_fit$gcv.ubre,
        CV_RMSE = cv_rmse,
        Max_VIF = max_vif,
        n_samples = nrow(current_data),
        n_countries = n_distinct(current_data$NAME),
        Model_Formula = formula_str
      )
    }

    bacteria_metrics <- bind_rows(metrics_rows)
    if (nrow(bacteria_metrics) == 0) {
      warning("No valid lag combinations were fitted for ", bacteria$title, " in ", spec$short_label)
      next
    }

    best_aic <- min(bacteria_metrics$AIC, na.rm = TRUE)
    best_model_idx <- which.min(bacteria_metrics$AIC)
    bacteria_metrics <- bacteria_metrics %>%
      mutate(
        AIC_diff = AIC - best_aic,
        rank = rank(AIC, ties.method = "min"),
        is_best = row_number() == best_model_idx
      ) %>%
      arrange(rank, AIC, CV_RMSE)

    best_row <- bacteria_metrics %>%
      filter(is_best) %>%
      slice(1)

    if (cv_scope == "best_only") {
      best_lag_combo <- list()
      for (var_name in spec$climate_vars) {
        best_lag_combo[[var_name]] <- best_row[[paste0(var_name, "_lag")]][[1]]
      }
      best_data <- create_lagged_dataset(base_data, best_lag_combo)
      best_cv_rmse <- perform_cv(best_data, bacteria$pls_components, spec$climate_vars, cv_key = bacteria$code)
      bacteria_metrics <- bacteria_metrics %>%
        mutate(CV_RMSE = ifelse(is_best, best_cv_rmse, CV_RMSE))
      best_row$CV_RMSE <- best_cv_rmse
    }

    all_metrics[[length(all_metrics) + 1]] <- bacteria_metrics
    summary_rows[[length(summary_rows) + 1]] <- best_row
  }

  all_metrics_df <- bind_rows(all_metrics)
  summary_df <- bind_rows(summary_rows) %>%
    arrange(factor(Display_Name, levels = vapply(bacteria_specs, `[[`, character(1), "title")))

  metrics_csv_path <- file.path(tables_dir, paste0("all_bacteria_", spec$key, "_lag_metrics.csv"))
  summary_csv_path <- file.path(tables_dir, paste0("bacteria_", spec$key, "_lag_summary.csv"))
  workbook_path <- file.path(workbook_dir, paste0("LagSelection_", spec$key, "_historical.xlsx"))
  manifest_path <- file.path(metadata_dir, paste0("LagSelection_", spec$key, "_manifest.csv"))

  write.csv(all_metrics_df, metrics_csv_path, row.names = FALSE)
  write.csv(summary_df, summary_csv_path, row.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "Optimal_Summary")
  addWorksheet(wb, "All_Metrics")
  writeData(wb, "Optimal_Summary", summary_df)
  writeData(wb, "All_Metrics", all_metrics_df)
  saveWorkbook(wb, workbook_path, overwrite = TRUE)

  bacteria_lookup <- tibble(
    Display_Name = vapply(bacteria_specs, `[[`, character(1), "title"),
    Input_File = file.path(input_data_dir, vapply(bacteria_specs, `[[`, character(1), "file_name"))
  )

  manifest_df <- summary_df %>%
    left_join(bacteria_lookup, by = "Display_Name") %>%
    transmute(
      Model = Model_Label,
      Bacteria = Display_Name,
      Input_File,
      n_samples,
      n_countries,
      Model_Formula
    )
  write.csv(manifest_df, manifest_path, row.names = FALSE)

  list(
    summary = summary_df,
    all_metrics = all_metrics_df,
    paths = list(
      summary_csv = summary_csv_path,
      metrics_csv = metrics_csv_path,
      workbook = workbook_path,
      manifest = manifest_path
    )
  )
}

import_existing_modelc_results <- function() {
  summary_path <- file.path(existing_modelc_root, "historical_lag_summary_model_c.csv")
  metrics_path <- file.path(existing_modelc_root, "historical_lag_metrics_model_c.csv")

  if (!file.exists(summary_path) || !file.exists(metrics_path)) {
    stop("Existing Model C lag-selection outputs were not found under: ", existing_modelc_root, call. = FALSE)
  }

  summary_df <- read.csv(summary_path, stringsAsFactors = FALSE) %>%
    mutate(
      Model_Key = "C",
      Model_Label = "Model C (TMP, PREC, HUM, WET)"
    ) %>%
    rename(
      Bacteria = Bacteria,
      Display_Name = Display_Name
    ) %>%
    select(
      Model_Key,
      Model_Label,
      Bacteria,
      Display_Name,
      TMP_lag,
      PREC_lag,
      HUM_lag,
      WET_lag,
      AIC,
      BIC,
      Dev_explained,
      GCV,
      CV_RMSE,
      Max_VIF,
      rank,
      n_samples
    ) %>%
    mutate(
      n_countries = NA_integer_,
      Model_Formula = NA_character_,
      AIC_diff = 0,
      is_best = TRUE
    ) %>%
    filter(Display_Name %in% vapply(bacteria_specs, `[[`, character(1), "title"))

  metrics_df <- read.csv(metrics_path, stringsAsFactors = FALSE) %>%
    mutate(
      Model_Key = "C",
      Model_Label = "Model C (TMP, PREC, HUM, WET)",
      Display_Name = recode(
        bacteria_name,
        "3GCREC" = "3GCR-Ec",
        "3GCRKP" = "3GCR-Kp",
        "CREC" = "CR-Ec",
        "CRKP" = "CR-Kp",
        "CRAB" = "CR-Ab",
        "CRPA" = "CR-Pa"
      ),
      Bacteria = bacteria_name,
      n_countries = NA_integer_,
      Model_Formula = NA_character_
    ) %>%
    filter(Display_Name %in% vapply(bacteria_specs, `[[`, character(1), "title")) %>%
    select(
      Model_Key,
      Model_Label,
      Bacteria,
      Display_Name,
      combination_id,
      TMP_lag,
      PREC_lag,
      HUM_lag,
      WET_lag,
      AIC,
      BIC,
      Dev_explained,
      GCV,
      CV_RMSE,
      Max_VIF,
      n_samples,
      n_countries,
      Model_Formula,
      AIC_diff,
      is_best,
      rank
    )

  list(
    summary = summary_df,
    all_metrics = metrics_df,
    paths = list(
      summary_csv = summary_path,
      metrics_csv = metrics_path
    )
  )
}

format_model_structure_summary <- function(summary_df) {
  summary_df %>%
    group_by(Display_Name) %>%
    mutate(
      Best_AIC = compute_best_flags(AIC, mode = "min", digits = 3),
      Best_ED = compute_best_flags(Dev_explained, mode = "max", digits = 3),
      Best_RMSE = compute_best_flags(CV_RMSE, mode = "min", digits = 3)
    ) %>%
    ungroup() %>%
    arrange(
      factor(Display_Name, levels = vapply(bacteria_specs, `[[`, character(1), "title")),
      factor(Model_Key, levels = c("A", "B", "C"))
    ) %>%
    transmute(
      AMR_Strain = Display_Name,
      Sample_Size = n_samples,
      Model = Model_Label,
      TMP_Lag = ifelse(is.na(TMP_lag), "", as.character(TMP_lag)),
      PREC_Lag = ifelse(is.na(PREC_lag), "", as.character(PREC_lag)),
      HUM_Lag = ifelse(is.na(HUM_lag), "", as.character(HUM_lag)),
      WET_Lag = ifelse(is.na(WET_lag), "", as.character(WET_lag)),
      AIC = round(AIC, 3),
      Explained_Deviance_pct = round(Dev_explained, 3),
      CV_RMSE = round(CV_RMSE, 3),
      Best_AIC,
      Best_ED,
      Best_RMSE
    )
}

write_model_structure_outputs <- function(summary_table, model_results, output_root) {
  tables_dir <- file.path(output_root, "01_tables")
  workbook_dir <- file.path(output_root, "02_workbook")
  doc_dir <- file.path(output_root, "03_docx")
  metadata_dir <- file.path(output_root, "04_metadata")
  ensure_dirs(c(output_root, tables_dir, workbook_dir, doc_dir, metadata_dir))

  csv_path <- file.path(tables_dir, "model_structure_comparison_ThreeClimateModels_OptimalLag_Comparison.csv")
  xlsx_path <- file.path(workbook_dir, "model_structure_comparison_ThreeClimateModels_OptimalLag_Comparison.xlsx")
  docx_path <- file.path(doc_dir, "model_structure_comparison_ThreeClimateModels_OptimalLag_Comparison.docx")
  manifest_path <- file.path(metadata_dir, "model_structure_comparison_run_manifest.csv")

  write.csv(summary_table, csv_path, row.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_table)

  manifest_rows <- list()
  for (model_key in names(model_results)) {
    result_obj <- model_results[[model_key]]
    summary_sheet <- sanitize_sheet_name(paste0("Summary_", model_key))
    metrics_sheet <- sanitize_sheet_name(paste0("AllMetrics_", model_key))
    addWorksheet(wb, summary_sheet)
    writeData(wb, summary_sheet, result_obj$summary)
    if (!is.null(result_obj$all_metrics)) {
      addWorksheet(wb, metrics_sheet)
      writeData(wb, metrics_sheet, result_obj$all_metrics)
    }

    manifest_rows[[length(manifest_rows) + 1]] <- tibble(
      Model_Key = model_key,
      Model_Label = unique(result_obj$summary$Model_Label),
      Summary_Path = unname(result_obj$paths$summary_csv),
      Metrics_Path = unname(result_obj$paths$metrics_csv),
      Imported = model_key == "C" && reuse_existing_modelc
    )
  }
  manifest_df <- bind_rows(manifest_rows)
  write.csv(manifest_df, manifest_path, row.names = FALSE)
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  docx_table <- summary_table %>%
    mutate(
      across(c(TMP_Lag, PREC_Lag, HUM_Lag, WET_Lag), ~ ifelse(is.na(.x), "", as.character(.x))),
      `Explained Deviance (%)` = sprintf("%.3f", Explained_Deviance_pct)
    ) %>%
    select(
      `AMR Strain` = AMR_Strain,
      `Sample Size` = Sample_Size,
      Model,
      `TMP Lag` = TMP_Lag,
      `PREC Lag` = PREC_Lag,
      `HUM Lag` = HUM_Lag,
      `WET Lag` = WET_Lag,
      AIC,
      `Explained Deviance (%)`,
      `CV-RMSE` = CV_RMSE,
      `Best AIC` = Best_AIC,
      `Best ED` = Best_ED,
      `Best RMSE` = Best_RMSE
    )

  ft <- flextable(docx_table)
  ft <- theme_booktabs(ft)
  ft <- autofit(ft)
  ft <- add_header_lines(
    ft,
    values = "model-structure comparison summary. Comparison of Three Climate Models' Optimal Lag Structures and Performance Metrics for AMR Resistance Prediction"
  )
  ft <- add_footer_lines(
    ft,
    values = "Note: Model A includes temperature, precipitation, and relative humidity; Model B includes temperature, precipitation, and wet days; Model C includes all four climate factors. Best AIC, Best ED, and Best RMSE are indicated within each AMR phenotype across the three model specifications."
  )
  save_as_docx(ft, path = docx_path)

  list(
    csv = csv_path,
    xlsx = xlsx_path,
    docx = docx_path,
    manifest = manifest_path
  )
}

main <- function() {
  debug_print("Starting model-structure comparison analysis")

  model_results <- list()

  for (model_key in names(model_specs)) {
    if (model_key == "C" && reuse_existing_modelc) {
      debug_print("Importing existing Model C lag-selection outputs")
      model_results[[model_key]] <- import_existing_modelc_results()
    } else {
      model_results[[model_key]] <- run_single_model_spec(model_specs[[model_key]])
    }
  }

  summary_table <- bind_rows(lapply(model_results, `[[`, "summary"))
  model_structure_table <- format_model_structure_summary(summary_table)
  model_structure_paths <- write_model_structure_outputs(model_structure_table, model_results, model_structure_output_root)

  debug_print("Model-specification sensitivity analysis completed")
  print(model_structure_paths)

  invisible(
    list(
      model_results = model_results,
      summary_table = model_structure_table,
      summary_paths = model_structure_paths
    )
  )
}

main()
