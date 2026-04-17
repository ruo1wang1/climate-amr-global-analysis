###############################################################################
# Projection lag-selection analysis for the simplified future-projection model
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(grid)
  library(gridExtra)
  library(writexl)
  library(officer)
  library(flextable)
})

set.seed(20260330)

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
output_base <- file.path(
  revision_root,
  "outputs",
  "ModelC_Full",
  "projection_preparation",
  "02_simplified_model_lag_selection_S18_S19_S20"
)

dir.create(file.path(output_base, "01_full_results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "02_toprank_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "03_pdf_figure_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "03_docx_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "04_workbook"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "05_model_summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_base, "06_metadata"), recursive = TRUE, showWarnings = FALSE)

bacteria_specs <- list(
  list(code = "3GCREC", title = "3GCR-Ec", file_name = "3GCREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "3GCRKP", title = "3GCR-Kp", file_name = "3GCRKP_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRAB", title = "CR-Ab", file_name = "CRAB_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CREC", title = "CR-Ec", file_name = "CREC_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRKP", title = "CR-Kp", file_name = "CRKP_data_with_PLS_components_3_with_MERRA2_and_ERA5_climate.csv"),
  list(code = "CRPA", title = "CR-Pa", file_name = "CRPA_data_with_PLS_components_4_with_MERRA2_and_ERA5_climate.csv")
)

load_lag_summary <- function(summary_csv) {
  df <- read.csv(summary_csv, stringsAsFactors = FALSE)
  df %>%
    select(Bacteria, Display_Name, TMP_lag, PREC_lag, HUM_lag, WET_lag, AIC, Dev_explained, n_samples)
}

prepare_base_data <- function(file_path) {
  data <- read.csv(file_path)

  data %>%
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
    ungroup() %>%
    mutate(HUM = pmin(HUM, 100)) %>%
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
    ungroup()
}

apply_lags <- function(data, lag_config) {
  data %>%
    group_by(location_id) %>%
    arrange(year) %>%
    mutate(
      TMP_scaled_lag = lag(TMP_scaled, lag_config$temp_lag),
      PREC_scaled_lag = lag(PREC_scaled, lag_config$prec_lag),
      HUM_scaled_lag = lag(HUM_scaled, lag_config$hum_lag),
      WET_scaled_lag = lag(WET_scaled, lag_config$wet_lag)
    ) %>%
    filter(
      !is.na(TMP_scaled_lag),
      !is.na(PREC_scaled_lag),
      !is.na(HUM_scaled_lag),
      !is.na(WET_scaled_lag)
    ) %>%
    ungroup()
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

calculate_model_performance <- function(model, data, bacteria_name, lag_config, combination_id) {
  tibble(
    Bacteria = bacteria_name,
    Lag_Combination = sprintf("%d-%d-%d-%d", lag_config$temp_lag, lag_config$prec_lag, lag_config$hum_lag, lag_config$wet_lag),
    combination_id = combination_id,
    AIC = AIC(model),
    BIC = BIC(model),
    Deviance_Explained_pct = summary(model)$dev.expl * 100,
    GCV = model$gcv.ubre,
    sample_size = nrow(data),
    temp_lag = lag_config$temp_lag,
    prec_lag = lag_config$prec_lag,
    hum_lag = lag_config$hum_lag,
    wet_lag = lag_config$wet_lag
  )
}

align_newdata_factor_levels <- function(train_data, test_data) {
  factor_cols <- names(train_data)[vapply(train_data, is.factor, logical(1))]
  factor_cols <- intersect(factor_cols, names(test_data))

  for (col in factor_cols) {
    train_levels <- levels(train_data[[col]])
    if (length(train_levels) == 0) {
      next
    }

    train_values <- as.character(train_data[[col]])
    reference_level <- names(sort(table(train_values), decreasing = TRUE))[1]
    test_values <- as.character(test_data[[col]])
    unseen_mask <- is.na(test_values) | !(test_values %in% train_levels)
    test_values[unseen_mask] <- reference_level
    test_data[[col]] <- factor(test_values, levels = train_levels)
  }

  test_data
}

perform_cv <- function(data, lag_config, k = 5) {
  set.seed(123)
  countries <- unique(data$NAME)
  folds <- split(countries, sample(rep(seq_len(k), length.out = length(countries))))
  rmse_values <- rep(NA_real_, k)

  for (i in seq_len(k)) {
    test_data <- data %>% filter(NAME %in% folds[[i]])
    train_data <- data %>% filter(!(NAME %in% folds[[i]]))

    if (nrow(test_data) < 10 || nrow(train_data) < 10) {
      next
    }

    train_model <- tryCatch(build_simplified_model(train_data), error = function(e) NULL)
    if (is.null(train_model)) {
      next
    }

    test_data_aligned <- align_newdata_factor_levels(train_data, test_data)
    pred <- tryCatch(
      predict(train_model, newdata = test_data_aligned),
      error = function(e) rep(NA_real_, nrow(test_data_aligned))
    )
    rmse_values[i] <- sqrt(mean((test_data$logit_R - pred)^2, na.rm = TRUE))
  }

  tibble(
    Lag_Combination = sprintf("%d-%d-%d-%d", lag_config$temp_lag, lag_config$prec_lag, lag_config$hum_lag, lag_config$wet_lag),
    CV_RMSE = mean(rmse_values, na.rm = TRUE),
    CV_RMSE_sd = sd(rmse_values, na.rm = TRUE)
  )
}

create_lag_combinations <- function(reference_lag, search_range = 1, max_lag = 3) {
  temp_range <- max(1, reference_lag$TMP_lag - search_range):min(max_lag, reference_lag$TMP_lag + search_range)
  prec_range <- max(1, reference_lag$PREC_lag - search_range):min(max_lag, reference_lag$PREC_lag + search_range)
  hum_range <- max(1, reference_lag$HUM_lag - search_range):min(max_lag, reference_lag$HUM_lag + search_range)
  wet_range <- max(1, reference_lag$WET_lag - search_range):min(max_lag, reference_lag$WET_lag + search_range)

  expand.grid(
    temp_lag = temp_range,
    prec_lag = prec_range,
    hum_lag = hum_range,
    wet_lag = wet_range
  ) %>%
    as_tibble() %>%
    mutate(
      combination_id = paste0("T", temp_lag, "_P", prec_lag, "_H", hum_lag, "_W", wet_lag),
      Lag_Combination = sprintf("%d-%d-%d-%d", temp_lag, prec_lag, hum_lag, wet_lag),
      is_original_optimal = temp_lag == reference_lag$TMP_lag &
        prec_lag == reference_lag$PREC_lag &
        hum_lag == reference_lag$HUM_lag &
        wet_lag == reference_lag$WET_lag
    )
}

significance_adjustment <- function(p_value) {
  case_when(
    is.na(p_value) ~ 0.1,
    p_value < 0.01 ~ 1.0,
    p_value < 0.05 ~ 0.7,
    p_value < 0.10 ~ 0.3,
    TRUE ~ 0.1
  )
}

bound_and_normalize_weights <- function(weights, min_w = 0.05, max_w = 0.7) {
  if (sum(weights, na.rm = TRUE) <= 0) {
    return(rep(0.25, length(weights)))
  }

  w <- weights / sum(weights, na.rm = TRUE)
  names(w) <- names(weights)

  for (iter in 1:20) {
    lower_mask <- w < min_w
    upper_mask <- w > max_w

    if (!any(lower_mask | upper_mask)) {
      break
    }

    w[lower_mask] <- min_w
    w[upper_mask] <- max_w

    fixed_mask <- lower_mask | upper_mask
    free_mask <- !fixed_mask
    remainder <- 1 - sum(w[fixed_mask], na.rm = TRUE)

    if (remainder <= 0 || !any(free_mask)) {
      w <- w / sum(w, na.rm = TRUE)
      break
    }

    free_weights <- weights[free_mask]
    if (sum(free_weights, na.rm = TRUE) <= 0) {
      w[free_mask] <- remainder / sum(free_mask)
    } else {
      w[free_mask] <- remainder * free_weights / sum(free_weights, na.rm = TRUE)
    }
  }

  w / sum(w, na.rm = TRUE)
}

derive_climate_weights <- function(model, bacteria_name) {
  s_table <- as.data.frame(summary(model)$s.table)
  s_table$term <- rownames(s_table)
  rownames(s_table) <- NULL

  climate_stats <- s_table %>%
    filter(term %in% c("s(TMP_scaled_lag)", "s(PREC_scaled_lag)", "s(HUM_scaled_lag)", "s(WET_scaled_lag)")) %>%
    transmute(
      Bacteria = bacteria_name,
      Climate_Factor = recode(
        term,
        "s(TMP_scaled_lag)" = "Temperature",
        "s(PREC_scaled_lag)" = "Precipitation",
        "s(HUM_scaled_lag)" = "Humidity",
        "s(WET_scaled_lag)" = "Wet Days"
      ),
      EDF = edf,
      F_stat = F,
      P_value = `p-value`,
      Sig_Adjustment = significance_adjustment(`p-value`),
      Raw_Weight = pmax(F, 0) * pmax(edf, 0) * Sig_Adjustment
    )

  ordered_factors <- c("Temperature", "Precipitation", "Humidity", "Wet Days")
  raw_weights <- climate_stats$Raw_Weight
  names(raw_weights) <- climate_stats$Climate_Factor
  raw_weights <- raw_weights[ordered_factors]
  final_weights <- bound_and_normalize_weights(raw_weights)

  climate_stats %>%
    mutate(
      Final_Weight = as.numeric(final_weights[Climate_Factor])
    ) %>%
    arrange(match(Climate_Factor, ordered_factors))
}

recommend_lag_update <- function(orig_rank, delta_aic, delta_dev) {
  if (is.na(orig_rank)) {
    return("Best simplified lag selected; review original-lag rank manually")
  }

  if (orig_rank <= 5) {
    if (abs(delta_aic) < 2) {
      return("Best simplified lag selected; original lag also competitive")
    }
    if (delta_aic <= 5) {
      return("Best simplified lag selected; original lag still top-5")
    }
    return("Best simplified lag selected; original lag still top-5 but less optimal")
  }

  if (delta_aic > 10) {
    return("Best simplified lag selected (clear improvement vs original)")
  }
  "Best simplified lag selected (moderate improvement vs original)"
}

format_numeric <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}

render_table_pdf <- function(df, pdf_path, title, footnote = NULL, width = 12, row_height = 0.32) {
  table_theme <- ttheme_minimal(
    base_size = 10,
    core = list(
      fg_params = list(fontfamily = "serif", cex = 0.9),
      bg_params = list(fill = rep(c("#F7F9FC", "#FFFFFF"), length.out = nrow(df)), col = NA)
    ),
    colhead = list(
      fg_params = list(fontfamily = "serif", fontface = "bold", col = "white", cex = 0.95),
      bg_params = list(fill = "#4C76C2", col = "#1E2A38")
    )
  )

  tg <- tableGrob(df, rows = NULL, theme = table_theme)
  height <- max(4.5, nrow(df) * row_height + ifelse(is.null(footnote), 1.5, 2.2))

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE, family = "serif")
  grid.newpage()
  grid.text(title, x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontsize = 15, fontface = "bold"))
  grid.draw(editGrob(tg, vp = viewport(x = 0.5, y = 0.52, width = 0.96, height = 0.82)))
  if (!is.null(footnote)) {
    grid.text(footnote, x = 0.02, y = 0.04, just = c("left", "bottom"), gp = gpar(fontsize = 10, fontface = "italic"))
  }
  dev.off()
}

make_ft_theme <- function(ft, blue_header = TRUE) {
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- fontsize(ft, size = 10.5, part = "body")
  ft <- fontsize(ft, size = 11, part = "header")
  ft <- align(ft, align = "center", part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- padding(ft, padding.top = 3, padding.bottom = 3, padding.left = 3, padding.right = 3)

  if (blue_header) {
    ft <- bg(ft, bg = "#4C76C2", part = "header")
    ft <- color(ft, color = "white", part = "header")
  }

  ft <- bold(ft, part = "header")
  ft <- border_remove(ft)

  outer_border <- fp_border(color = "black", width = 1.2)
  inner_border <- fp_border(color = "#333333", width = 0.6)

  ft <- border_outer(ft, border = outer_border, part = "all")
  ft <- hline_top(ft, border = outer_border, part = "all")
  ft <- hline_bottom(ft, border = outer_border, part = "all")
  ft <- hline(ft, border = inner_border, part = "header")
  ft <- hline(ft, border = inner_border, part = "body")
  ft <- vline(ft, border = inner_border, part = "all")
  ft
}

save_table_docx <- function(ft, title, note, path) {
  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(doc, "", style = "Normal")
  doc <- body_add_flextable(doc, value = paginate(ft, init = TRUE, hdr_ftr = TRUE))
  if (!is.null(note)) {
    doc <- body_add_par(doc, "", style = "Normal")
    doc <- body_add_par(doc, note, style = "Normal")
  }
  print(doc, target = path)
}

build_s18_flextable <- function(df) {
  ft <- flextable(df)
  ft <- make_ft_theme(ft)
  ft <- autofit(ft)
  ft <- width(ft, j = 1, width = 1.0)
  ft <- width(ft, j = 2, width = 1.0)
  ft <- width(ft, j = c(3, 4), width = 0.65)
  ft <- width(ft, j = 5:9, width = 0.9)
  ft <- width(ft, j = 10, width = 0.75)
  ft
}

build_s19_flextable <- function(df) {
  display_df <- df %>%
    transmute(
      `AMR Pathogens` = AMR_Strain,
      TMP = Original_TMP,
      PREC = Original_PREC,
      HUM = Original_HUM,
      WET = Original_WET,
      TMP_best = Best_TMP,
      PREC_best = Best_PREC,
      HUM_best = Best_HUM,
      WET_best = Best_WET,
      `Orig. Rank` = Original_Rank,
      `ΔAIC` = format_numeric(Delta_AIC),
      `ΔDev%` = format_numeric(Delta_Dev_pct),
      `Original AIC` = format_numeric(Original_AIC),
      `Best AIC` = format_numeric(Best_AIC),
      `Orig. Dev%` = format_numeric(Original_Dev_pct),
      `Best Dev%` = format_numeric(Best_Dev_pct),
      Recommendation = Recommendation
    )

  ft <- flextable(display_df, col_keys = names(display_df))
  ft <- add_header_row(
    ft,
    values = c("", "ORIGINAL OPTIMAL LAG", "BEST SIMPLIFIED MODEL LAG", "PERFORMANCE COMPARISON"),
    colwidths = c(1, 4, 4, 8)
  )
  ft <- set_header_labels(
    ft,
    `AMR Pathogens` = "AMR Pathogens",
    TMP = "TMP",
    PREC = "PREC",
    HUM = "HUM",
    WET = "WET",
    TMP_best = "TMP",
    PREC_best = "PREC",
    HUM_best = "HUM",
    WET_best = "WET",
    `Orig. Rank` = "Orig.\nRank",
    `ΔAIC` = "ΔAIC",
    `ΔDev%` = "ΔDev%",
    `Original AIC` = "Original\nAIC",
    `Best AIC` = "Best AIC",
    `Orig. Dev%` = "Orig.\nDev%",
    `Best Dev%` = "Best\nDev%",
    Recommendation = "Recommendation"
  )
  ft <- make_ft_theme(ft)
  ft <- height(ft, i = 1:2, height = 0.45, part = "header")
  ft <- bold(ft, j = c(1, 10, 17), bold = TRUE, part = "body")

  rank_vals <- df$Original_Rank
  competitive_vals <- grepl("competitive|top-5", df$Recommendation, ignore.case = TRUE)
  clear_vals <- grepl("clear improvement", df$Recommendation, ignore.case = TRUE)

  ft <- bg(ft, i = which(rank_vals <= 5), j = 10, bg = "#DDEBD2", part = "body")
  ft <- bg(ft, i = which(rank_vals > 5 & rank_vals <= 10), j = 10, bg = "#FCECC2", part = "body")
  ft <- bg(ft, i = which(rank_vals > 10), j = 10, bg = "#F7C9A9", part = "body")

  ft <- bg(ft, i = which(competitive_vals), j = 17, bg = "#E7F2DA", part = "body")
  ft <- bg(ft, i = which(!competitive_vals & !clear_vals), j = 17, bg = "#FFF0C7", part = "body")
  ft <- bg(ft, i = which(clear_vals), j = 17, bg = "#F7D7B7", part = "body")

  ft <- color(ft, i = which(as.numeric(df$Delta_AIC) > 0), j = 11, color = "#E67E00", part = "body")
  ft <- align(ft, j = 17, align = "left", part = "body")
  ft <- autofit(ft)
  ft <- width(ft, j = 1, width = 1.15)
  ft <- width(ft, j = 2:9, width = 0.56)
  ft <- width(ft, j = 10:16, width = 0.82)
  ft <- width(ft, j = 17, width = 1.9)
  ft
}

build_s20_flextable <- function(df) {
  ft <- flextable(df)
  ft <- make_ft_theme(ft)
  ft <- autofit(ft)
  ft <- width(ft, j = 1, width = 1.0)
  ft <- width(ft, j = 2:5, width = 1.1)
  ft
}

lag_summary <- load_lag_summary(lag_summary_path)
all_results <- list()
best_models <- list()
weight_tables <- list()
s19_rows <- list()

for (spec in bacteria_specs) {
  message("Processing simplified lag search for ", spec$title)
  base_data <- prepare_base_data(file.path(input_data_dir, spec$file_name))
  reference_row <- lag_summary %>% filter(Display_Name == spec$title)
  lag_grid <- create_lag_combinations(reference_row)

  perf_list <- list()

  for (i in seq_len(nrow(lag_grid))) {
    lag_config <- list(
      temp_lag = lag_grid$temp_lag[i],
      prec_lag = lag_grid$prec_lag[i],
      hum_lag = lag_grid$hum_lag[i],
      wet_lag = lag_grid$wet_lag[i]
    )

    lagged_data <- apply_lags(base_data, lag_config)
    if (nrow(lagged_data) < 100) {
      next
    }

    model <- tryCatch(build_simplified_model(lagged_data), error = function(e) NULL)
    if (is.null(model)) {
      next
    }

    perf <- calculate_model_performance(model, lagged_data, spec$title, lag_config, lag_grid$combination_id[i]) %>%
      mutate(is_original_optimal = lag_grid$is_original_optimal[i])
    perf_list[[length(perf_list) + 1]] <- perf
  }

  perf_df <- bind_rows(perf_list) %>%
    arrange(AIC) %>%
    mutate(
      Rank = row_number(),
      Optimal = ifelse(Rank == 1, "Yes", "No")
    )

  if (nrow(perf_df) == 0) {
    next
  }

  top_ranked <- perf_df %>% slice_head(n = 10)
  eval_combos <- union(top_ranked$Lag_Combination, perf_df %>% filter(is_original_optimal) %>% pull(Lag_Combination))

  cv_rows <- list()
  for (lag_string in eval_combos) {
    split_lag <- as.integer(strsplit(lag_string, "-", fixed = TRUE)[[1]])
    lag_config <- list(temp_lag = split_lag[1], prec_lag = split_lag[2], hum_lag = split_lag[3], wet_lag = split_lag[4])
    lagged_data <- apply_lags(base_data, lag_config)
    cv_rows[[length(cv_rows) + 1]] <- perform_cv(lagged_data, lag_config, k = 5)
  }
  cv_df <- bind_rows(cv_rows)

  perf_df <- perf_df %>% left_join(cv_df, by = "Lag_Combination")
  top_ranked <- perf_df %>% slice_head(n = 10)

  best_row <- perf_df %>% slice(1)
  best_lag_vec <- as.integer(strsplit(best_row$Lag_Combination, "-", fixed = TRUE)[[1]])
  best_model_data <- apply_lags(
    base_data,
    list(temp_lag = best_lag_vec[1], prec_lag = best_lag_vec[2], hum_lag = best_lag_vec[3], wet_lag = best_lag_vec[4])
  )
  best_model <- build_simplified_model(best_model_data)
  best_models[[spec$title]] <- best_model

  model_summary_file <- file.path(output_base, "05_model_summaries", paste0("Simplified_ModelC_", gsub("-", "_", spec$title), "_best_model_summary.txt"))
  sink(model_summary_file)
  cat("Simplified Projection Model Summary:", spec$title, "\n")
  cat("Best lag combination:", best_row$Lag_Combination, "\n\n")
  print(summary(best_model))
  cat("\n")
  print(gam.check(best_model))
  sink()

  weight_tables[[spec$title]] <- derive_climate_weights(best_model, spec$title)

  original_row <- perf_df %>% filter(is_original_optimal)
  s19_rows[[length(s19_rows) + 1]] <- tibble(
    AMR_Strain = spec$title,
    Original_TMP = reference_row$TMP_lag,
    Original_PREC = reference_row$PREC_lag,
    Original_HUM = reference_row$HUM_lag,
    Original_WET = reference_row$WET_lag,
    Best_TMP = best_lag_vec[1],
    Best_PREC = best_lag_vec[2],
    Best_HUM = best_lag_vec[3],
    Best_WET = best_lag_vec[4],
    Original_Rank = original_row$Rank,
    Delta_AIC = original_row$AIC - best_row$AIC,
    Delta_Dev_pct = best_row$Deviance_Explained_pct - original_row$Deviance_Explained_pct,
    Original_AIC = original_row$AIC,
    Best_AIC = best_row$AIC,
    Original_Dev_pct = original_row$Deviance_Explained_pct,
    Best_Dev_pct = best_row$Deviance_Explained_pct,
    Recommendation = recommend_lag_update(original_row$Rank, original_row$AIC - best_row$AIC, best_row$Deviance_Explained_pct - original_row$Deviance_Explained_pct)
  )

  all_results[[spec$title]] <- perf_df
}

full_results_df <- bind_rows(all_results)
s18_df <- full_results_df %>%
  group_by(Bacteria) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  transmute(
    AMR_Strain = Bacteria,
    Lag_Combination = Lag_Combination,
    Optimal = Optimal,
    Rank = Rank,
    AIC = round(AIC, 3),
    BIC = round(BIC, 3),
    Deviance_Explained_pct = round(Deviance_Explained_pct, 3),
    GCV = round(GCV, 3),
    CV_RMSE = round(CV_RMSE, 3),
    sample_size = sample_size
  )

s19_df <- bind_rows(s19_rows) %>%
  mutate(
    Delta_AIC = round(Delta_AIC, 3),
    Delta_Dev_pct = round(Delta_Dev_pct, 3),
    Original_AIC = round(Original_AIC, 3),
    Best_AIC = round(Best_AIC, 3),
    Original_Dev_pct = round(Original_Dev_pct, 3),
    Best_Dev_pct = round(Best_Dev_pct, 3)
  )

projection_lag_settings_df <- s19_df %>%
  transmute(
    Bacteria = AMR_Strain,
    TMP_lag = Best_TMP,
    PREC_lag = Best_PREC,
    HUM_lag = Best_HUM,
    WET_lag = Best_WET,
    Source = "Best simplified projection GAMM lag",
    Original_Lag_Rank_In_Simplified_Search = Original_Rank,
    Delta_AIC_vs_Original = Delta_AIC,
    Delta_Dev_pct_vs_Original = Delta_Dev_pct
  )

weights_long <- bind_rows(weight_tables)
s20_df <- weights_long %>%
  select(Bacteria, Climate_Factor, Final_Weight) %>%
  mutate(Climate_Factor = factor(Climate_Factor, levels = c("Temperature", "Precipitation", "Humidity", "Wet Days"))) %>%
  pivot_wider(names_from = Climate_Factor, values_from = Final_Weight) %>%
  arrange(match(Bacteria, c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa"))) %>%
  mutate(across(-Bacteria, ~ round(.x, 4)))

write.csv(full_results_df, file.path(output_base, "01_full_results", "projection_lag_validation_all_simplified_lag_results.csv"), row.names = FALSE)
write.csv(s18_df, file.path(output_base, "02_toprank_tables", "projection_lag_validation_top10_simplified_lag_results.csv"), row.names = FALSE)
write.csv(s19_df, file.path(output_base, "02_toprank_tables", "projection_lag_comparison_simplified_vs_original_lag_comparison.csv"), row.names = FALSE)
write.csv(s20_df, file.path(output_base, "02_toprank_tables", "projection_simplified_climate_weights.csv"), row.names = FALSE)
write.csv(weights_long, file.path(output_base, "02_toprank_tables", "projection_simplified_climate_weight_details.csv"), row.names = FALSE)
write.csv(
  projection_lag_settings_df,
  file.path(output_base, "02_toprank_tables", "projection_simplified_final_lag_settings.csv"),
  row.names = FALSE
)

write_xlsx(
  list(
    S18_All_Results = full_results_df,
    S18_Top10 = s18_df,
    S19_Comparison = s19_df,
    Projection_Final_Lags = projection_lag_settings_df,
    S20_Weights = s20_df,
    S20_Derivation = weights_long
  ),
  file.path(output_base, "04_workbook", "source_data_projection_simplified_lag_selection.xlsx")
)

s18_display <- s18_df %>%
  mutate(
    AIC = format_numeric(AIC),
    BIC = format_numeric(BIC),
    Deviance_Explained_pct = format_numeric(Deviance_Explained_pct),
    GCV = format_numeric(GCV),
    CV_RMSE = ifelse(is.na(CV_RMSE), "NA", format_numeric(CV_RMSE)),
    sample_size = as.character(sample_size)
  )

s19_display <- s19_df %>%
  mutate(
    Original_Lag = sprintf("%d-%d-%d-%d", Original_TMP, Original_PREC, Original_HUM, Original_WET),
    Best_Lag = sprintf("%d-%d-%d-%d", Best_TMP, Best_PREC, Best_HUM, Best_WET)
  ) %>%
  transmute(
    AMR_Strain,
    Original_Lag,
    Best_Lag,
    Original_Rank,
    Delta_AIC = format_numeric(Delta_AIC),
    Delta_Dev_pct = format_numeric(Delta_Dev_pct),
    Original_AIC = format_numeric(Original_AIC),
    Best_AIC = format_numeric(Best_AIC),
    Original_Dev_pct = format_numeric(Original_Dev_pct),
    Best_Dev_pct = format_numeric(Best_Dev_pct),
    Recommendation
  )

s20_display <- s20_df %>%
  mutate(
    Temperature = format(round(Temperature, 4), nsmall = 4),
    Precipitation = format(round(Precipitation, 4), nsmall = 4),
    Humidity = format(round(Humidity, 4), nsmall = 4),
    `Wet Days` = format(round(`Wet Days`, 4), nsmall = 4)
  )

render_table_pdf(
  s18_display,
  file.path(output_base, "03_pdf_figure_tables", "projection_lag_validation_simplified_lag_results_figurestyle.pdf"),
  "projection lag validation summary. Validation Results of Climate Lag Combinations for Simplified Model C",
  "Note: Rows are ranked by AIC within each AMR phenotype. CV_RMSE is reported for the top-ranked combinations and the original main-model lag when evaluated in the simplified model."
)

render_table_pdf(
  s19_display,
  file.path(output_base, "03_pdf_figure_tables", "projection_lag_comparison_simplified_vs_original_lag_comparison_figurestyle.pdf"),
  "projection lag comparison summary. Climate Lag Comparison Between Full Model C and Simplified Projection Model",
  "Note: Delta_AIC is calculated as Original AIC minus Best Simplified AIC. Positive values indicate improved fit for the best simplified-model lag, which is the lag setting used in downstream projection modeling."
)

render_table_pdf(
  s20_display,
  file.path(output_base, "03_pdf_figure_tables", "projection_climate_weight_climate_factor_weights_figurestyle.pdf"),
  "projection climate-weight summary. Climatic Factor Weights for Different AMR Pathogens Based on Simplified GAMM Analysis",
  "Note: Final weights are normalized from F-statistic x EDF x significance-adjustment scores, then bounded to the biologically plausible interval 0.05-0.70 with renormalization."
)

s18_ft <- build_s18_flextable(s18_display)
s19_ft <- build_s19_flextable(s19_df)
s20_ft <- build_s20_flextable(s20_display)

save_table_docx(
  s18_ft,
  "projection lag validation summary. Validation Results of Climate Lag Combinations for Various Bacterial AMR Strains (Model C - Four Climatic Variables Simplified Model Lag)",
  "Note: Rows are ranked by AIC within each AMR phenotype. CV_RMSE is reported for the top-ranked combinations and for the original main-model lag when evaluated in the simplified model.",
  file.path(output_base, "03_docx_tables", "projection_lag_validation_simplified_lag_results.docx")
)

save_table_docx(
  s19_ft,
  "projection lag comparison summary. Climate Lag Combination Comparison and Validation in Simplified Models",
  "Note: This table compares optimal climate lag combinations and performance differences between the original GAMM model C and the simplified projection model C (with PLS components removed). Positive ΔAIC values indicate improved fit for the best simplified-model lag combination. The BEST SIMPLIFIED MODEL LAG column is the lag setting carried forward into downstream projection modeling.",
  file.path(output_base, "03_docx_tables", "projection_lag_comparison_simplified_vs_original_lag_comparison.docx")
)

save_table_docx(
  s20_ft,
  "projection climate-weight summary. Climatic Factor Weights for Different AMR Pathogens Based on GAMM Analysis",
  "Note: Values represent normalized climatic-factor weights derived from the best simplified GAMM for each AMR phenotype using the rule Raw weight = F × EDF × significance adjustment, followed by bounding to the interval 0.05-0.70 and renormalization.",
  file.path(output_base, "03_docx_tables", "projection_climate_weight_climate_factor_weights.docx")
)

combined_doc <- read_docx()
combined_doc <- body_add_par(combined_doc, "projection lag summaries", style = "heading 1")
combined_doc <- body_add_par(combined_doc, "", style = "Normal")
combined_doc <- body_add_par(combined_doc, "projection lag validation summary. Validation Results of Climate Lag Combinations for Various Bacterial AMR Strains (Model C - Four Climatic Variables Simplified Model Lag)", style = "heading 2")
combined_doc <- body_add_flextable(combined_doc, value = paginate(s18_ft, init = TRUE, hdr_ftr = TRUE))
combined_doc <- body_add_par(combined_doc, "Note: Rows are ranked by AIC within each AMR phenotype. CV_RMSE is reported for the top-ranked combinations and for the original main-model lag when evaluated in the simplified model.", style = "Normal")
combined_doc <- body_add_break(combined_doc)
combined_doc <- body_add_par(combined_doc, "projection lag comparison summary. Climate Lag Combination Comparison and Validation in Simplified Models", style = "heading 2")
combined_doc <- body_add_flextable(combined_doc, value = paginate(s19_ft, init = TRUE, hdr_ftr = TRUE))
combined_doc <- body_add_par(combined_doc, "Note: This table compares optimal climate lag combinations and performance differences between the original GAMM model C and the simplified projection model C (with PLS components removed). Positive ΔAIC values indicate improved fit for the best simplified-model lag combination. The BEST SIMPLIFIED MODEL LAG column is the lag setting carried forward into downstream projection modeling.", style = "Normal")
combined_doc <- body_add_break(combined_doc)
combined_doc <- body_add_par(combined_doc, "projection climate-weight summary. Climatic Factor Weights for Different AMR Pathogens Based on GAMM Analysis", style = "heading 2")
combined_doc <- body_add_flextable(combined_doc, value = paginate(s20_ft, init = TRUE, hdr_ftr = TRUE))
combined_doc <- body_add_par(combined_doc, "Note: Values represent normalized climatic-factor weights derived from the best simplified GAMM for each AMR phenotype using the rule Raw weight = F × EDF × significance adjustment, followed by bounding to the interval 0.05-0.70 and renormalization.", style = "Normal")
print(combined_doc, target = file.path(output_base, "03_docx_tables", "projection_lag_summaries.docx"))

methods_lines <- c(
  "# Future Projection Preparation Methods",
  "",
  "## Simplified projection model",
  "The projection-phase simplified GAMM retained the four climatic smooth terms and the same spatiotemporal structure as the full Model C, while removing all PLS components.",
  "",
  "## Lag re-selection",
  "For each AMR phenotype, candidate lag combinations were searched within a +/-1 year window around the full-model lag combination, bounded to 1-3 years for each climatic factor. Candidate models were ranked by AIC after fitting the simplified GAMM.",
  "The original full-model lag combination was explicitly retained in the evaluated set and its ranking within the simplified-model search space was recorded for comparison only.",
  "Downstream projection modeling uses the AIC-optimal BEST SIMPLIFIED MODEL LAG for each AMR phenotype rather than mechanically retaining the original full-model lag.",
  "",
  "## Climate factor weights",
  "Climatic factor weights were derived from the best simplified model using a pragmatic weighting rule aligned with the original methods documentation:",
  "- Raw weight_i = F_i x EDF_i x S_i",
  "- Significance adjustment S_i: p < 0.01 -> 1.0; p < 0.05 -> 0.7; p < 0.10 -> 0.3; p >= 0.10 -> 0.1",
  "- Raw weights were normalized to sum to 1",
  "- A lower bound of 0.05 and upper bound of 0.70 were imposed for biological plausibility, followed by renormalization"
)
writeLines(methods_lines, file.path(output_base, "06_metadata", "METHODS.md"))
writeLines(capture.output(sessionInfo()), file.path(output_base, "06_metadata", "sessionInfo.txt"))

cat("Saved simplified projection lag-selection outputs to:\n")
cat(output_base, "\n")
