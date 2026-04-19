######## Analysis of climate-AMR associations across alternative basis dimensions ########

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(openxlsx)
  library(flextable)
  library(officer)
})

Sys.setenv(MODELC_SKIP_MAIN = "1")

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "..", "01_historical_associations", "analysis_historical_associations_main_model.R"), local = FALSE)

basis_output_suffix <- Sys.getenv("MODELC_BASIS_OUTPUT_SUFFIX", unset = "")
basis_run_only <- trimws(Sys.getenv("MODELC_BASIS_RUN_ONLY", unset = ""))
basis_k_values <- c(4, 5, 8, 10, 12)

fit_k_sensitivity_model <- function(data, bacteria_name, basis_k) {
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
  available_pls <- get_available_pls_components(data)

  if (length(available_pls) == 0) {
    stop("No available PLS components found for ", bacteria_name, call. = FALSE)
  }

  pls_terms <- paste0("s(", available_pls, ", k = ", basis_k, ", bs = 'cr')", collapse = " + ")
  spatial_k <- min(20, basis_k * 2)
  year_k <- min(8, basis_k)

  formula_str <- paste0(
    "logit_R ~ ",
    "s(TMP_scaled_lag, k = ", basis_k, ", bs = 'cr') + ",
    "s(PREC_scaled_lag, k = ", basis_k, ", bs = 'cr') + ",
    "s(HUM_scaled_lag, k = ", basis_k, ", bs = 'cr') + ",
    "s(WET_scaled_lag, k = ", basis_k, ", bs = 'cr') + ",
    pls_terms, " + ",
    "s(lat, lon, bs = 'sos', k = ", spatial_k, ") + ",
    "s(year, bs = 'cr', k = ", year_k, ") + ",
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
    warning("bam() failed for ", bacteria_name, " at k=", basis_k, "; trying gam(): ", e$message)
    gam(
      model_formula,
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE
    )
  })

  list(
    model = model,
    formula = formula_str,
    available_pls = available_pls,
    spatial_k = spatial_k,
    year_k = year_k
  )
}

extract_basis_row <- function(model_fit, bacteria_name, basis_k, n_samples, available_pls) {
  model_summary <- summary(model_fit)
  data.frame(
    AMR_Strain = bacteria_name,
    k_Value = basis_k,
    R2 = round(model_summary$r.sq, 3),
    Explained_Deviance_pct = round(model_summary$dev.expl * 100, 3),
    AIC = round(AIC(model_fit), 3),
    BIC = round(BIC(model_fit), 3),
    Effective_df = round(sum(model_fit$edf), 3),
    n_samples = n_samples,
    PLS_Components = paste(available_pls, collapse = ", "),
    stringsAsFactors = FALSE
  )
}

build_basis_dimension_table <- function(specs_to_run) {
  results <- list()
  model_details <- list()

  for (spec in specs_to_run) {
    cat("Running basis-dimension sensitivity for", spec$title, "...\n")
    prepared <- prepare_data(file.path(input_data_dir, spec$file_name), spec$title)
    data_ready <- prepared$data

    for (basis_k in basis_k_values) {
      cat("  Testing k =", basis_k, "\n")
      fit_obj <- fit_k_sensitivity_model(data_ready, spec$title, basis_k)
      results[[length(results) + 1]] <- extract_basis_row(
        model_fit = fit_obj$model,
        bacteria_name = spec$title,
        basis_k = basis_k,
        n_samples = nrow(data_ready),
        available_pls = fit_obj$available_pls
      )
      model_details[[length(model_details) + 1]] <- data.frame(
        AMR_Strain = spec$title,
        k_Value = basis_k,
        Model_Formula = fit_obj$formula,
        Spatial_k = fit_obj$spatial_k,
        Year_k = fit_obj$year_k,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    summary_table = bind_rows(results),
    model_details = bind_rows(model_details)
  )
}

write_basis_dimension_outputs <- function(summary_table, model_details, output_root) {
  tables_dir <- file.path(output_root, "01_tables")
  workbook_dir <- file.path(output_root, "02_workbook")
  doc_dir <- file.path(output_root, "03_docx")
  metadata_dir <- file.path(output_root, "04_metadata")

  for (dir_path in c(output_root, tables_dir, workbook_dir, doc_dir, metadata_dir)) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }

  csv_path <- file.path(tables_dir, "basis_dimension_sensitivity_summary.csv")
  xlsx_path <- file.path(workbook_dir, "source_data_basis_dimension_sensitivity.xlsx")
  docx_path <- file.path(doc_dir, "basis_dimension_sensitivity_summary.docx")
  details_csv_path <- file.path(metadata_dir, "basis_dimension_sensitivity_model_formulas.csv")

  write.csv(summary_table, csv_path, row.names = FALSE)
  write.csv(model_details, details_csv_path, row.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "BasisDimension")
  addWorksheet(wb, "Model_Formulas")

  writeData(wb, "BasisDimension", summary_table)
  writeData(wb, "Model_Formulas", model_details)
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  table_for_docx <- summary_table %>%
    mutate(
      `Explained Deviance (%)` = sprintf("%.3f%%", Explained_Deviance_pct)
    ) %>%
    select(
      `AMR Strain` = AMR_Strain,
      `k Value` = k_Value,
      `R2` = R2,
      `Explained Deviance (%)`,
      AIC,
      BIC,
      `Effective df` = Effective_df
    )

  ft <- flextable(table_for_docx)
  ft <- theme_booktabs(ft)
  ft <- autofit(ft)
  ft <- add_header_lines(
    ft,
    values = "Basis-dimension sensitivity summary. Sensitivity analysis of GAMM model performance across alternative basis dimensions"
  )
  ft <- add_footer_lines(
    ft,
    values = "Note: Effective df represents the effective degrees of freedom used by the model. Results are based on the primary historical model inputs, phenotype-specific optimal lags, and dynamic PLS-component inclusion."
  )

  save_as_docx(ft, path = docx_path)

  list(
    csv = csv_path,
    xlsx = xlsx_path,
    docx = docx_path,
    model_formulas = details_csv_path
  )
}

main <- function() {
  specs_to_run <- bacteria_specs
  if (nzchar(basis_run_only)) {
    run_only_values <- trimws(strsplit(basis_run_only, ",")[[1]])
    specs_to_run <- Filter(function(x) x$code %in% run_only_values || x$title %in% run_only_values, specs_to_run)
  }

  if (length(specs_to_run) == 0) {
    stop("No bacteria matched MODELC_BASIS_RUN_ONLY.", call. = FALSE)
  }

  output_root <- file.path(
    revision_root,
    "outputs",
    "historical_associations",
    paste0("basis_dimension_sensitivity", basis_output_suffix)
  )

  basis_results <- build_basis_dimension_table(specs_to_run)
  output_paths <- write_basis_dimension_outputs(
    summary_table = basis_results$summary_table,
    model_details = basis_results$model_details,
    output_root = output_root
  )

  cat("Basis-dimension outputs saved to:\n")
  print(output_paths)

  invisible(list(
    summary_table = basis_results$summary_table,
    model_details = basis_results$model_details,
    output_paths = output_paths
  ))
}

invisible(main())
