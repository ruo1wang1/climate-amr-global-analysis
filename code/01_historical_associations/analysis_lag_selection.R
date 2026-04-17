
#================================================================= 
# Climate Lag Combination Analysis for Six Bacterial Species (Four Climate Factors)
# 
# Purpose: To provide scientific basis for lag combination selection
# Author: [Your Name]
# Date: 2025
#================================================================= 

# Install necessary packages (if not already installed)
required_packages <- c("tidyverse", "mgcv", "viridis", "patchwork", "gridExtra",
                      "grid", "scales", "extrafont", "cowplot", "kableExtra", "tidytext")
for(pkg in required_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
}

# Load necessary packages
library(tidyverse)
library(mgcv)
library(viridis)
library(patchwork)
library(gridExtra)
library(grid)
library(scales)
library(cowplot)
library(kableExtra)
library(tidytext)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
base_path <- Sys.getenv("CLIMATE_AMR_ANALYSIS_INPUT_DIR", unset = file.path(repo_root, "data", "analysis_inputs"))
code_base_dir <- file.path(repo_root, "outputs")

# Make grouped CV reproducible across reruns.
set.seed(20260329)

# Try loading extrafont, continue if failed
tryCatch({
  library(extrafont)
  # Uncomment the following line to import system fonts if first use
  # font_import()
  # loadfonts(device = "win") # Windows system
  # loadfonts(device = "pdf") # For exporting PDF
}, error = function(e) {
  warning("extrafont package failed to load, will use system default fonts. Error message: ", e$message)
})

output_suffix <- Sys.getenv("LAGSEL_OUTPUT_SUFFIX", unset = "")
output_dir <- file.path(
  code_base_dir,
  paste0("analysis_lag_selection", output_suffix)
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define bacteria configuration list - using new abbreviated naming
bacteria_config <- list(
  list(name = "3GCREC", pls_components = 3, display_name = "3GCR-Ec"),
  list(name = "3GCRKP", pls_components = 4, display_name = "3GCR-Kp"),
  list(name = "CREC", pls_components = 3, display_name = "CR-Ec"),
  list(name = "CRKP", pls_components = 3, display_name = "CR-Kp"),
  list(name = "CRAB", pls_components = 4, display_name = "CR-Ab"),
  list(name = "CRPA", pls_components = 4, display_name = "CR-Pa")
)

run_only <- Sys.getenv("LAGSEL_RUN_ONLY", unset = "")
if(nzchar(run_only)) {
  run_only_values <- trimws(strsplit(run_only, ",")[[1]])
  bacteria_config <- bacteria_config[vapply(bacteria_config, function(x) x$name %in% run_only_values, logical(1))]
}

max_lag_to_run <- suppressWarnings(as.integer(Sys.getenv("LAGSEL_MAX_LAG", unset = "3")))
if(is.na(max_lag_to_run) || max_lag_to_run < 1 || max_lag_to_run > 3) {
  stop("LAGSEL_MAX_LAG must be an integer between 1 and 3.")
}

# Set theme and colors
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_rect(colour = "gray80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
      panel.grid.minor = element_line(colour = "gray95", linewidth = 0.1),
      axis.title = element_text(size = rel(0.9), face = "bold"),
      strip.background = element_rect(fill = "gray95", colour = "gray80"),
      strip.text = element_text(face = "bold", size = rel(0.9)),
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0.5),
      legend.title = element_text(face = "bold", size = rel(0.9)),
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = rel(0.8)),
      plot.caption = element_text(size = rel(0.7), hjust = 0),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# Define richer color schemes
bacteria_colors <- c(
  "3GCREC" = "#1f77b4", "3GCRKP" = "#ff7f0e",
  "CREC" = "#2ca02c", "CRKP" = "#d62728",
  "CRAB" = "#9467bd", "CRPA" = "#8c564b"
)

# Lag year color schemes - different color schemes for each climate variable
temp_lag_colors <- c("1" = "#e41a1c", "2" = "#377eb8", "3" = "#4daf4a")
prec_lag_shapes <- c("1" = 16, "2" = 17, "3" = 15) # Circle, triangle, square
hum_lag_sizes <- c("1" = 2, "2" = 3.5, "3" = 5) # Small, medium, large
wet_lag_colors <- c("1" = "#984ea3", "2" = "#ff7f00", "3" = "#a65628") # New color scheme for wet days

# Simple lag year colors - for simple charts
lag_colors <- c("1" = "#3f007d", "2" = "#35978f", "3" = "#f7e269")

# AIC difference color gradient
aic_diff_colors <- c(
  "0-2" = "#4575b4",   # Deep blue - Very close to optimal
  "2-5" = "#74add1",   # Medium blue - Close to optimal
  "5-10" = "#abd9e9",  # Light blue - Acceptable
  "10-15" = "#fee090", # Yellow - Not ideal
  "15-20" = "#fdae61", # Orange - Poor
  "20+" = "#d73027"    # Red - Very poor
)

#------------------------------------------------------------------------
# Debug and utility functions
#------------------------------------------------------------------------
# Debug info function
debug_print <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s][%s] %s\n", timestamp, level, message))
}

# Check factor levels count
check_factor_levels <- function(data) {
  factor_vars <- sapply(data, is.factor)
  factor_names <- names(data)[factor_vars]
  problems <- list()
  for(var in factor_names) {
    levels_count <- length(levels(data[[var]]))
    if(levels_count < 2) {
      problems[[var]] <- levels_count
    }
  }
  return(problems)
}

# Dynamically build model formula with four climate variables
build_safe_formula <- function(data, pls_components) {
  # Check categorical variables
  factor_problems <- check_factor_levels(data)
  
  # Build base formula part with FOUR climate variables
  base_formula <- paste("logit_R ~ s(TMP_scaled_lag, k = 10, bs = 'cr') + ",
                       "s(PREC_scaled_lag, k = 10, bs = 'cr') + ",
                       "s(HUM_scaled_lag, k = 10, bs = 'cr') + ",
                       "s(WET_scaled_lag, k = 10, bs = 'cr')")
  
  # Add PLS components
  pls_terms <- paste0(sapply(1:pls_components, function(j) {
    sprintf("s(PLS_Comp%d, k = 10, bs = 'cr')", j)
  }), collapse = " + ")
  
  formula_parts <- c(base_formula, pls_terms)
  
  # Add spatial term if enough data points
  if(nrow(data) > 30) {
    formula_parts <- c(formula_parts, "s(lat, lon, bs = 'sos', k = min(20, nrow(data)/10))")
  }
  
  # Add year term
  year_levels <- length(unique(data$year))
  if(year_levels > 2) {
    k_year <- min(8, year_levels - 1)
    formula_parts <- c(formula_parts, sprintf("s(year, bs = 'cr', k = %d)", k_year))
  } else if(year_levels == 2) {
    formula_parts <- c(formula_parts, "year") # Linear effect
  }
  
  # Add safe categorical variables
  if(!("Region" %in% names(factor_problems)) && "Region" %in% names(data)) {
    formula_parts <- c(formula_parts, "s(Region, bs = 're')")
  }
  
  if(!("climate_zone" %in% names(factor_problems)) && "climate_zone" %in% names(data)) {
    formula_parts <- c(formula_parts, "climate_zone")
  }
  
  # Combine all parts
  full_formula <- paste(formula_parts, collapse = " + ")
  return(full_formula)
}

# Align factor levels between training and test sets so grouped CV can predict
# when a held-out country introduces unseen Region/climate levels.
align_test_factor_levels <- function(train_data, test_data) {
  factor_cols <- names(train_data)[vapply(train_data, is.factor, logical(1))]
  factor_cols <- intersect(factor_cols, names(test_data))
  
  for(col in factor_cols) {
    train_levels <- levels(train_data[[col]])
    if(length(train_levels) == 0) next
    
    train_values <- as.character(train_data[[col]])
    reference_level <- names(sort(table(train_values), decreasing = TRUE))[1]
    test_values <- as.character(test_data[[col]])
    unseen <- is.na(test_values) | !(test_values %in% train_levels)
    test_values[unseen] <- reference_level
    test_data[[col]] <- factor(test_values, levels = train_levels)
  }
  
  return(test_data)
}

safe_rmse <- function(obs, pred) {
  keep <- is.finite(obs) & is.finite(pred)
  if(sum(keep) < 5) {
    return(NA_real_)
  }
  sqrt(mean((obs[keep] - pred[keep])^2))
}

compute_max_vif <- function(data) {
  vars <- c("TMP_scaled_lag", "PREC_scaled_lag", "HUM_scaled_lag", "WET_scaled_lag")
  x <- data[, vars, drop = FALSE]
  
  # Remove non-finite rows for the VIF calculation only.
  x <- x[stats::complete.cases(x), , drop = FALSE]
  out <- setNames(rep(NA_real_, length(vars)), vars)
  
  if(nrow(x) < 10) {
    return(out)
  }
  
  zero_var <- vapply(x, function(col) stats::sd(col) == 0, logical(1))
  if(any(zero_var)) {
    x <- x[, !zero_var, drop = FALSE]
  }
  
  if(ncol(x) < 2) {
    return(out)
  }
  
  for(var in names(x)) {
    others <- setdiff(names(x), var)
    if(length(others) == 0) {
      out[var] <- 1
      next
    }
    fit <- try(lm(reformulate(others, response = var), data = x), silent = TRUE)
    if(inherits(fit, "try-error")) next
    r2 <- summary(fit)$r.squared
    if(is.na(r2)) next
    out[var] <- if(r2 >= 0.999999) Inf else 1 / (1 - r2)
  }
  
  return(out)
}

# Improved cross-validation function with robustness
perform_cv <- function(data, pls_components, k = 5) {
  n <- nrow(data)
  
  # If data too small, reduce folds
  if(n < 50) {
    k <- min(k, max(2, floor(n/10)))
    debug_print(sprintf("Small sample size (%d), adjusted cross-validation folds to %d", n, k), "INFO")
  }
  
  # Stratified sampling - grouped by country
  countries <- unique(data$NAME)
  country_folds <- sample(rep(1:k, length.out = length(countries)))
  names(country_folds) <- countries
  
  data$fold <- NA
  for(country in countries) {
    data$fold[data$NAME == country] <- country_folds[country]
  }
  
  cv_rmse <- numeric(k)
  
  for(i in 1:k) {
    test_idx <- which(data$fold == i)
    train_data <- droplevels(data[-test_idx, ])
    test_data <- data[test_idx, ]
    
    # Check if training and test data have enough samples
    if(nrow(train_data) < 20 || nrow(test_data) < 5) {
      debug_print(sprintf("Fold %d: Training data (%d) or test data (%d) sample size insufficient", 
                         i, nrow(train_data), nrow(test_data)), "WARNING")
      cv_rmse[i] <- NA
      next
    }
    
    # Check factor levels in training data
    factor_problems <- check_factor_levels(train_data)
    if(length(factor_problems) > 0) {
      debug_print(sprintf("Fold %d: Factor level issues: %s", 
                         i, paste(names(factor_problems), collapse=", ")), "WARNING")
    }
    
    # Build safe formula for each training set
    formula_str <- build_safe_formula(train_data, pls_components)
    
    cv_model <- try({
      bam(
        as.formula(formula_str),
        data = train_data,
        family = gaussian(),
        method = "fREML",
        discrete = TRUE,
        select = TRUE,  # Use select=TRUE to handle correlated predictors
        nthreads = 1    # Avoid OpenMP issues, use single thread
      )
    }, silent = TRUE)
    
    if(!inherits(cv_model, "try-error")) {
      aligned_test_data <- align_test_factor_levels(train_data, test_data)
      pred <- try(predict(cv_model, newdata = aligned_test_data), silent = TRUE)
      
      if(!inherits(pred, "try-error")) {
        fold_rmse <- safe_rmse(aligned_test_data$logit_R, pred)
        if(is.na(fold_rmse)) {
          cv_rmse[i] <- NA
          debug_print("Cross-validation prediction produced too few finite values", "WARNING")
        } else {
          cv_rmse[i] <- fold_rmse
        }
      } else {
        cv_rmse[i] <- NA
        debug_print("Cross-validation prediction failed", "WARNING")
      }
    } else {
      cv_rmse[i] <- NA
      debug_print("Cross-validation model fitting failed", "WARNING")
    }
  }
  
  result <- mean(cv_rmse, na.rm = TRUE)
  
  if(is.nan(result) || is.na(result) || length(cv_rmse) == 0 || all(is.na(cv_rmse))) {
    debug_print("All cross-validation folds failed, returning default large value", "WARNING")
    return(999) # Return a large value indicating validation failure
  }
  
  return(result)
}

#------------------------------------------------------------------------
# Main analysis functions
#------------------------------------------------------------------------
# Analyze all lag combinations for a single bacteria - FOUR Climate Factors
analyze_bacteria_lag_combinations_four_factors <- function(bacteria_name, pls_components, max_lag = 3) {
  debug_print(paste("Analyzing", bacteria_name, "lag combinations with four climate factors"))
  
  # Set data file path
  data_file <- file.path(
    base_path,
    paste0(
      bacteria_name,
      "_data_with_PLS_components_",
      pls_components,
      "_with_MERRA2_and_ERA5_climate.csv"
    )
  )
  
  # Load data
  if(!file.exists(data_file)) {
    debug_print(paste("Data file not found:", data_file), "ERROR")
    return(NULL)
  }
  
  data <- read.csv(data_file)
  debug_print("Data loaded successfully")
  
  # Process data
  data_processed <- data %>%
    mutate(
      year = as.numeric(as.character(year)),
      # Ensure categorical variables set correctly
      Region = if("Region" %in% names(.)) factor(Region) else factor(1),
      climate_zone = factor(case_when(
        abs(lat) > 66.5 ~ "Polar",
        abs(lat) > 23.5 ~ "Temperate",
        TRUE ~ "Tropical"
      ), levels = c("Tropical", "Temperate", "Polar")),
      hemisphere = factor(if_else(lat >= 0, "North", "South"))
    ) %>%
    # Remove missing values
    filter(!is.na(year), !is.na(lat), !is.na(lon), !is.na(logit_R)) %>%
    group_by(NAME) %>%
    arrange(year) %>%
    # Scale all FOUR climate variables within each country
    mutate(across(c(TMP, PREC, HUM, WET), list(scaled = ~scale(.)), .names = "{.col}_scaled")) %>%
    ungroup()
  
  # Check data
  debug_print(sprintf("Processed data contains %d rows from %d countries", 
                     nrow(data_processed), length(unique(data_processed$NAME))))
  
  # Check for missing WET data
  if(sum(is.na(data_processed$WET)) > 0) {
    debug_print(sprintf("WARNING: %d missing values in WET variable", sum(is.na(data_processed$WET))), "WARNING")
  }
  
  # Generate all possible lag combinations for four factors
  lag_combinations <- expand.grid(
    TMP_lag = 1:max_lag,
    PREC_lag = 1:max_lag,
    HUM_lag = 1:max_lag,
    WET_lag = 1:max_lag
  )
  
  # Initialize result storage
  results <- list()
  model_metrics <- data.frame()
  
  # Iterate through all lag combinations
  for(i in 1:nrow(lag_combinations)) {
    tmp_lag <- lag_combinations$TMP_lag[i]
    prec_lag <- lag_combinations$PREC_lag[i]
    hum_lag <- lag_combinations$HUM_lag[i]
    wet_lag <- lag_combinations$WET_lag[i]
    
    debug_print(sprintf("Testing combination %d/%d: TMP=%d, PREC=%d, HUM=%d, WET=%d", 
                       i, nrow(lag_combinations), tmp_lag, prec_lag, hum_lag, wet_lag))
    
    # Create lag variables for all four climate factors
    current_data <- data_processed %>%
      group_by(NAME) %>%
      mutate(
        TMP_scaled_lag = lag(TMP_scaled, tmp_lag),
        PREC_scaled_lag = lag(PREC_scaled, prec_lag),
        HUM_scaled_lag = lag(HUM_scaled, hum_lag),
        WET_scaled_lag = lag(WET_scaled, wet_lag)
      ) %>%
      ungroup() %>%
      filter(!is.na(TMP_scaled_lag), !is.na(PREC_scaled_lag), 
             !is.na(HUM_scaled_lag), !is.na(WET_scaled_lag))
    
    # Check data after lag
    if(nrow(current_data) < 30) {
      debug_print(sprintf("Combination %d insufficient data after lag (%d rows)", i, nrow(current_data)), "WARNING")
      next
    }
    
    # Check factor levels
    problems <- check_factor_levels(current_data)
    if(length(problems) > 0) {
      debug_print(sprintf("Combination %d factor issues: %s", i, paste(names(problems), collapse=", ")), "WARNING")
    }
    
    # Build safe model formula including WET variable
    formula_str <- build_safe_formula(current_data, pls_components)
    debug_print(sprintf("Using formula: %s", formula_str), "INFO")
    
    # Fit model with select=TRUE to handle correlated predictors
    model <- try({
      bam(
        as.formula(formula_str),
        data = current_data,
        family = gaussian(),
        method = "fREML",
        discrete = TRUE,
        select = TRUE,  # Important: helps with correlated predictors
        nthreads = 1    # Avoid OpenMP issues
      )
    }, silent = TRUE)
    
    # Store results
    if(!inherits(model, "try-error")) {
      # Calculate VIF to monitor multicollinearity
      vif_values <- try(compute_max_vif(current_data), silent = TRUE)
      
      if(inherits(vif_values, "try-error")) {
        vif_values <- rep(NA, 4)
        names(vif_values) <- c("TMP_scaled_lag", "PREC_scaled_lag", "HUM_scaled_lag", "WET_scaled_lag")
      }
      
      # Max VIF as indicator of multicollinearity severity
      finite_vif <- vif_values[is.finite(vif_values)]
      max_vif <- if(any(is.infinite(vif_values))) {
        Inf
      } else if(length(finite_vif) == 0) {
        NA_real_
      } else {
        max(finite_vif)
      }
      
      # Perform cross-validation
      cv_rmse <- perform_cv(current_data, pls_components)
      
      results[[i]] <- model
      model_metrics <- rbind(model_metrics, data.frame(
        bacteria_name = bacteria_name,
        combination_id = i,
        TMP_lag = tmp_lag,
        PREC_lag = prec_lag,
        HUM_lag = hum_lag,
        WET_lag = wet_lag,
        AIC = AIC(model),
        BIC = BIC(model),
        Dev_explained = summary(model)$dev.expl * 100,
        GCV = model$gcv.ubre,
        CV_RMSE = cv_rmse,
        Max_VIF = max_vif,
        lag_combo = paste(tmp_lag, prec_lag, hum_lag, wet_lag, sep=","),
        n_samples = nrow(current_data)
      ))
      
      debug_print(sprintf("Model %d fit successfully: AIC=%.2f, Dev=%.2f%%, CV-RMSE=%.4f, Max_VIF=%.2f", 
                         i, AIC(model), summary(model)$dev.expl * 100, cv_rmse, max_vif))
    } else {
      debug_print(sprintf("Model %d fitting failed: %s", i, as.character(model)), "WARNING")
    }
  }
  
  if(nrow(model_metrics) == 0) {
    debug_print("No successful model fits", "ERROR")
    return(NULL)
  }
  
  # Calculate relative AIC values
  min_aic <- min(model_metrics$AIC)
  model_metrics$AIC_diff <- model_metrics$AIC - min_aic
  
  # Mark best combination
  best_model_idx <- which.min(model_metrics$AIC)
  model_metrics$is_best <- (1:nrow(model_metrics) == best_model_idx)
  
  # Rank combinations
  model_metrics$rank <- rank(model_metrics$AIC)
  
  return(list(
    model_metrics = model_metrics,
    best_model = results[[best_model_idx]],
    best_combination = lag_combinations[best_model_idx, ],
    all_models = results
  ))
}

#------------------------------------------------------------------------
# Analyze and collect data for all bacteria with four climate factors
#------------------------------------------------------------------------
# Run analysis for all bacteria
analyze_all_bacteria_four_factors <- function(max_lag = 3) {
  all_results <- list()
  all_metrics <- data.frame()
  summary_table <- data.frame()
  
  for(i in 1:length(bacteria_config)) {
    bacteria_name <- bacteria_config[[i]]$name
    pls_components <- bacteria_config[[i]]$pls_components
    display_name <- bacteria_config[[i]]$display_name
    
    debug_print(paste("Analyzing", bacteria_name, "with four climate factors"))
    
    # Run analysis
    results <- analyze_bacteria_lag_combinations_four_factors(
      bacteria_name,
      pls_components,
      max_lag = max_lag
    )
    
    if(!is.null(results)) {
      all_results[[bacteria_name]] <- results
      all_metrics <- rbind(all_metrics, results$model_metrics)
      
      # Extract best combination
      best_idx <- which(results$model_metrics$is_best)
      best_metrics <- results$model_metrics[best_idx, ]
      
      # Add to summary table
      summary_table <- rbind(summary_table, data.frame(
        Bacteria = bacteria_name,
        Display_Name = display_name,
        TMP_lag = best_metrics$TMP_lag,
        PREC_lag = best_metrics$PREC_lag,
        HUM_lag = best_metrics$HUM_lag,
        WET_lag = best_metrics$WET_lag,  # Added WET lag
        AIC = best_metrics$AIC,
        AIC_diff = best_metrics$AIC_diff, # Always 0
        BIC = best_metrics$BIC,
        Dev_explained = best_metrics$Dev_explained,
        GCV = best_metrics$GCV,
        CV_RMSE = best_metrics$CV_RMSE,
        Max_VIF = best_metrics$Max_VIF,  # Added VIF
        rank = 1, # Best combination always rank 1
        n_samples = best_metrics$n_samples
      ))
    } else {
      debug_print(paste("Failed to analyze", bacteria_name), "WARNING")
    }
  }
  
  # Save all data
  saveRDS(all_results, file.path(output_dir, "historical_lag_selection_results_model_c.rds"))
  write.csv(all_metrics, file.path(output_dir, "historical_lag_metrics_model_c.csv"), row.names = FALSE)
  write.csv(summary_table, file.path(output_dir, "historical_lag_summary_model_c.csv"), row.names = FALSE)
  
  return(list(
    all_results = all_results,
    all_metrics = all_metrics,
    summary_table = summary_table
  ))
}

#------------------------------------------------------------------------
# Create high-quality visualization functions for four climate factors
#------------------------------------------------------------------------

# 1. Create faceted AIC heatmap - adapted for four factors
create_facet_aic_heatmap_four_factors <- function(all_metrics, summary_table) {
  # For four factors, we need to create separate visualizations
  # This is a challenge because we can't easily show 4D in one plot
  
  # Create heatmap for TMP and PREC, faceted by HUM and WET
  for(wet_lag in 1:3) {
    plot_data <- all_metrics %>% 
      filter(WET_lag == wet_lag) %>%
      left_join(summary_table %>% select(Bacteria, Display_Name), by = c("bacteria_name" = "Bacteria"))
    
    # Skip if no data for this WET lag
    if(nrow(plot_data) == 0) next
    
    # Create relative AIC values for each bacteria
    plot_data <- plot_data %>%
      group_by(bacteria_name) %>%
      mutate(
        AIC_relative = AIC - min(AIC),
        best_combo = ifelse(is_best, "Optimal", "Other")
      ) %>%
      ungroup()
    
    # Calculate min and max relative AIC for color scale
    min_rel_aic <- 0 # Best combination's relative AIC is always 0
    max_rel_aic <- max(min(20, max(plot_data$AIC_relative, na.rm = TRUE)), 5) # Limit max value
    
    # Create heatmap for this WET lag
    p <- ggplot(plot_data, aes(x = factor(TMP_lag), y = factor(PREC_lag), fill = AIC_relative)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_point(data = filter(plot_data, is_best), 
                 aes(x = factor(TMP_lag), y = factor(PREC_lag)), 
                 shape = 21, fill = NA, color = "black", size = 5, stroke = 1) +
      scale_fill_viridis_c(
        option = "plasma",
        direction = -1,
        limits = c(min_rel_aic, max_rel_aic),
        breaks = c(0, 2, 5, 10, 15, 20),
        labels = c("0 (Best)", "2", "5", "10", "15", "20+"),
        name = "AIC Difference"
      ) +
      facet_grid(Display_Name ~ HUM_lag, 
                 labeller = labeller(HUM_lag = function(x) paste("Humidity Lag:", x, "years"))) +
      labs(
        title = paste("AIC Evaluation Heatmap - Wet Days Lag:", wet_lag, "years"),
        subtitle = "Black circles indicate optimal combinations (minimum AIC)",
        x = "Temperature Lag (years)",
        y = "Precipitation Lag (years)"
      ) +
      theme_publication() +
      theme(
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "right"
      )
    
    # Save image
    ggsave(file.path(output_dir, paste0("facet_aic_heatmap_wet_lag_", wet_lag, ".pdf")), 
           p, width = 10, height = 12, dpi = 300)
  }
  
  # Also create a summary heatmap showing only the best combinations
  best_combs <- all_metrics %>% 
    filter(is_best) %>%
    left_join(summary_table %>% select(Bacteria, Display_Name), by = c("bacteria_name" = "Bacteria"))
  
  # Create a summary plot
  p_summary <- ggplot(best_combs, aes(x = Display_Name, y = "Best Combination")) +
    geom_tile(fill = "lightblue") +
    geom_text(aes(label = paste("T:", TMP_lag, "P:", PREC_lag, "H:", HUM_lag, "W:", WET_lag)), 
              fontface = "bold") +
    labs(
      title = "Optimal Lag Combinations (Four Climate Factors)",
      subtitle = "T=Temperature, P=Precipitation, H=Humidity, W=Wet Days",
      x = NULL, y = NULL
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  # Save summary image
  ggsave(file.path(output_dir, "best_combinations_summary.pdf"), 
         p_summary, width = 10, height = 4, dpi = 300)
  
  return(p_summary)
}

# 2. Multi-criteria comparison plot for four factors - simplified version
create_multicriteria_plot_four_factors <- function(all_metrics, summary_table) {
  # Prepare data
  plot_data <- all_metrics %>%
    left_join(summary_table %>% select(Bacteria, Display_Name), by = c("bacteria_name" = "Bacteria")) %>%
    group_by(bacteria_name) %>%
    mutate(
      AIC_min = min(AIC),
      is_within_2 = AIC <= (AIC_min + 2),
      # Create a combined factor for visualization
      TMP_PREC = paste("T", TMP_lag, "P", PREC_lag, sep=""),
      HUM_WET = paste("H", HUM_lag, "W", WET_lag, sep="")
    ) %>%
    ungroup()
  
  # Create simplified plot using combined factors
  p <- ggplot(plot_data, aes(x = Dev_explained, y = CV_RMSE)) +
    # Main data points with combined visual elements
    geom_point(
      aes(color = TMP_PREC, shape = HUM_WET),
      size = 3, alpha = 0.7
    ) +
    # Mark combinations with AIC difference < 2 with square outline
    geom_point(
      data = filter(plot_data, is_within_2 & !is_best),
      aes(x = Dev_explained, y = CV_RMSE),
      shape = 0, size = 5, color = "black", stroke = 1
    ) +
    # Highlight optimal combinations with circular outline
    geom_point(
      data = filter(plot_data, is_best),
      aes(x = Dev_explained, y = CV_RMSE),
      shape = 1, size = 5, color = "black", stroke = 1.5
    ) +
    # Facet by bacteria with free scales
    facet_wrap(~ Display_Name, scales = "free") +
    # Labels and theme
    labs(
      title = "Multi-Criteria Evaluation of Four Climate Factor Lag Combinations",
      subtitle = "Black circles mark optimal combinations (lowest AIC). Square outlines indicate AIC difference < 2",
      x = "Explained Deviance (%)",
      y = "Cross-Validation RMSE",
      color = "Temp & Prec Lags",
      shape = "Humidity & Wet Lags"
    ) +
    theme_publication() +
    theme(
      strip.background = element_rect(fill = "grey95"),
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin(6, 6, 6, 6)
    )
  
  # Save image
  ggsave(file.path(output_dir, "multicriteria_comparison_four_factors.pdf"), 
         p, width = 12, height = 8, dpi = 300)
  
  return(p)
}

# 3. Best lag combination visualization for four factors
create_best_lag_visualization_four_factors <- function(summary_table) {
  # Copy data and ensure order
  plot_data <- summary_table %>%
    arrange(Bacteria) %>%
    mutate(Display_Name = factor(Display_Name, levels = Display_Name))
  
  # Create best lag combination visualization with four climate factors
  p <- ggplot(plot_data) +
    # Temperature lag
    geom_segment(aes(x = 0, xend = TMP_lag, y = Display_Name, yend = Display_Name, color = Bacteria), 
                 size = 2.5, lineend = "round") +
    geom_point(aes(x = TMP_lag, y = Display_Name, color = Bacteria), size = 5, shape = 16) +
    # Precipitation lag
    geom_segment(aes(x = 3, xend = 3 + PREC_lag, y = Display_Name, yend = Display_Name, color = Bacteria), 
                 size = 2.5, lineend = "round") +
    geom_point(aes(x = 3 + PREC_lag, y = Display_Name, color = Bacteria), size = 5, shape = 17) +
    # Humidity lag
    geom_segment(aes(x = 6, xend = 6 + HUM_lag, y = Display_Name, yend = Display_Name, color = Bacteria), 
                 size = 2.5, lineend = "round") +
    geom_point(aes(x = 6 + HUM_lag, y = Display_Name, color = Bacteria), size = 5, shape = 15) +
    # Wet days lag (new)
    geom_segment(aes(x = 9, xend = 9 + WET_lag, y = Display_Name, yend = Display_Name, color = Bacteria), 
                 size = 2.5, lineend = "round") +
    geom_point(aes(x = 9 + WET_lag, y = Display_Name, color = Bacteria), size = 5, shape = 18) +
    # Add labels and reference lines
    scale_color_manual(values = bacteria_colors, guide = "none") +
    scale_x_continuous(
      breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
      labels = c("0", "1", "2", "0", "1", "2", "0", "1", "2", "0", "1", "2", "3"),
      sec.axis = sec_axis(~., breaks = c(1.5, 4.5, 7.5, 10.5), 
                         labels = c("Temperature", "Precipitation", "Humidity", "Wet Days"))
    ) +
    annotate("rect", xmin = 0, xmax = 3, ymin = 0.5, ymax = 6.5, fill = "gray97", alpha = 0.3) +
    annotate("rect", xmin = 3, xmax = 6, ymin = 0.5, ymax = 6.5, fill = "gray94", alpha = 0.3) +
    annotate("rect", xmin = 6, xmax = 9, ymin = 0.5, ymax = 6.5, fill = "gray97", alpha = 0.3) +
    annotate("rect", xmin = 9, xmax = 12, ymin = 0.5, ymax = 6.5, fill = "gray94", alpha = 0.3) +
    labs(
      title = "Optimal Climate Lag Periods for Six Bacterial Species (Four Climate Factors)",
      x = "Lag (years)",
      y = NULL
    ) +
    theme_publication() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.ticks.y = element_blank()
    )
  
  # Save image
  ggsave(file.path(output_dir, "best_lag_visualization_four_factors.pdf"), 
         p, width = 12, height = 6, dpi = 300)
  
  return(p)
}

# 4. Lag distribution heatmap for four factors
create_lag_heatmap_four_factors <- function(summary_table) {
  # Prepare heatmap data
  lag_heatmap_data <- data.frame(
    Bacteria = rep(summary_table$Display_Name, 4),
    Variable = c(
      rep("Temperature", nrow(summary_table)),
      rep("Precipitation", nrow(summary_table)),
      rep("Humidity", nrow(summary_table)),
      rep("Wet Days", nrow(summary_table))
    ),
    Lag = c(
      summary_table$TMP_lag,
      summary_table$PREC_lag,
      summary_table$HUM_lag,
      summary_table$WET_lag
    ),
    Bacteria_Code = rep(summary_table$Bacteria, 4)
  )
  
  # Ensure order
  lag_heatmap_data$Variable <- factor(lag_heatmap_data$Variable, 
                                     levels = c("Temperature", "Precipitation", "Humidity", "Wet Days"))
  lag_heatmap_data$Bacteria <- factor(lag_heatmap_data$Bacteria, 
                                     levels = rev(summary_table$Display_Name))
  
  # Create heatmap
  p <- ggplot(lag_heatmap_data, aes(x = Variable, y = Bacteria, fill = factor(Lag))) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(
      name = "Lag (years)",
      values = lag_colors
    ) +
    geom_text(aes(label = Lag), color = "white", size = 5, fontface = "bold") +
    labs(
      title = "Optimal Climate Lag Combinations for Six Bacterial Species",
      subtitle = "Four Climate Factors Model",
      x = NULL,
      y = NULL
    ) +
    theme_publication() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      panel.grid = element_blank()
    )
  
  # Save image
  ggsave(file.path(output_dir, "lag_distribution_heatmap_four_factors.pdf"), 
         p, width = 8, height = 6, dpi = 300)
  
  return(p)
}

# 5. AIC rank comparison plot for four factors
create_aic_rank_plot_four_factors <- function(all_metrics, summary_table) {
  # For four factors, we'll show top 20 combinations for each bacteria
  plot_data <- all_metrics %>%
    left_join(summary_table %>% select(Bacteria, Display_Name), by = c("bacteria_name" = "Bacteria")) %>%
    group_by(bacteria_name) %>%
    mutate(
      AIC_diff = AIC - min(AIC),
      AIC_diff_capped = pmin(AIC_diff, 20), # Cap max value for better display
      combo_rank = rank(AIC),
      lag_combo = paste0(TMP_lag, ",", PREC_lag, ",", HUM_lag, ",", WET_lag),
      # Create categories based on AIC difference
      AIC_diff_cat = cut(
        AIC_diff,
        breaks = c(-Inf, 2, 5, 10, 15, 20, Inf),
        labels = c("0-2", "2-5", "5-10", "10-15", "15-20", "20+"),
        include.lowest = TRUE
      )
    ) %>%
    # Take only top 20 combinations for each bacteria to avoid overcrowding
    filter(combo_rank <= 20) %>%
    ungroup()
  
  # Create AIC ranking plot
  p <- ggplot(plot_data, aes(x = reorder_within(lag_combo, combo_rank, bacteria_name), 
                            y = AIC_diff_capped, fill = AIC_diff_cat)) +
    # Use category for fill color
    geom_col(width = 0.7) +
    # Use gradient color scheme
    scale_fill_manual(
      values = aic_diff_colors,
      name = "AIC Difference",
      drop = FALSE
    ) +
    # Use proper labels
    scale_x_reordered() +
    facet_wrap(~ Display_Name, scales = "free_x") +
    labs(
      title = "AIC Difference Comparison of Lag Combinations for Six Bacterial Species",
      subtitle = "Top 20 combinations shown per species (Four Climate Factors Model)",
      x = "Lag Combination (Temperature, Precipitation, Humidity, Wet Days)",
      y = "AIC Difference (capped at 20)"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      panel.spacing = unit(1, "lines"),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
  
  # Save image
  ggsave(file.path(output_dir, "aic_rank_comparison_four_factors.pdf"), 
         p, width = 12, height = 8, dpi = 300)
  
  return(p)
}

# 6. Create detailed summary table
create_summary_table_with_metrics_four_factors <- function(summary_table) {
  # Organize table data
  table_data <- summary_table %>%
    select(Display_Name, TMP_lag, PREC_lag, HUM_lag, WET_lag, AIC, Dev_explained, CV_RMSE, Max_VIF) %>%
    rename(
      "Bacterial Species" = Display_Name,
      "Temperature Lag (years)" = TMP_lag,
      "Precipitation Lag (years)" = PREC_lag,
      "Humidity Lag (years)" = HUM_lag,
      "Wet Days Lag (years)" = WET_lag,
      "Explained Deviance (%)" = Dev_explained,
      "Max VIF" = Max_VIF
    )
  
  # Beautify table and export as HTML, if kableExtra available
  if(requireNamespace("kableExtra", quietly = TRUE)) {
    html_table <- kableExtra::kable(table_data, format = "html", 
                                  digits = c(0, 0, 0, 0, 0, 2, 2, 4, 2)) %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                              full_width = TRUE) %>%
      kableExtra::row_spec(0, bold = TRUE, background = "lightgray") %>%
      kableExtra::column_spec(1, bold = TRUE) %>%
      kableExtra::add_header_above(c(" " = 5, "Model Evaluation Metrics" = 4))
    
    # Save HTML table
    writeLines(html_table, file.path(output_dir, "lag_summary_table_four_factors.html"))
  }
  
  # Also save as CSV
  write.csv(table_data, file.path(output_dir, "lag_summary_table_four_factors.csv"), row.names = FALSE)
  
  return(table_data)
}

# 7. Create lag distribution bar plot
create_lag_dist_barplot_four_factors <- function(summary_table) {
  # Calculate lag year distribution
  lag_counts <- data.frame(
    Variable = c(rep("TMP", 3), rep("PREC", 3), rep("HUM", 3), rep("WET", 3)),
    Lag = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
    Count = c(
      sum(summary_table$TMP_lag == 1),
      sum(summary_table$TMP_lag == 2),
      sum(summary_table$TMP_lag == 3),
      sum(summary_table$PREC_lag == 1),
      sum(summary_table$PREC_lag == 2),
      sum(summary_table$PREC_lag == 3),
      sum(summary_table$HUM_lag == 1),
      sum(summary_table$HUM_lag == 2),
      sum(summary_table$HUM_lag == 3),
      sum(summary_table$WET_lag == 1),
      sum(summary_table$WET_lag == 2),
      sum(summary_table$WET_lag == 3)
    )
  ) %>%
    mutate(
      Variable = factor(Variable, levels = c("TMP", "PREC", "HUM", "WET"),
                       labels = c("Temperature", "Precipitation", "Humidity", "Wet Days")),
      Proportion = Count / nrow(summary_table),
      Lag = factor(Lag)
    )
  
  # Create stacked bar chart
  p <- ggplot(lag_counts, aes(x = Variable, y = Count, fill = Lag)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.7), 
              vjust = -0.5, fontface = "bold") +
    scale_fill_manual(values = lag_colors, name = "Lag (years)") +
    labs(
      title = "Distribution of Climate Lag Years Across Six Bacterial Species",
      subtitle = "Four Climate Factors Model",
      x = "Climate Variable",
      y = "Number of Bacterial Species"
    ) +
    theme_publication() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # Save image
  ggsave(file.path(output_dir, "lag_distribution_barplot_four_factors.pdf"), 
         p, width = 9, height = 5, dpi = 300)
  
  return(p)
}

# 8. Create criteria consistency plot for four factors
create_criteria_consistency_plot_four_factors <- function(all_metrics, summary_table) {
  # For each bacteria, find optimal combination based on different criteria
  consistency_data <- all_metrics %>%
    group_by(bacteria_name) %>%
    summarize(
      AIC_best_combo = lag_combo[which.min(AIC)],
      BIC_best_combo = lag_combo[which.min(BIC)],
      DevExp_best_combo = lag_combo[which.max(Dev_explained)],
      CVRMSE_best_combo = lag_combo[which.min(CV_RMSE)]
    ) %>%
    left_join(summary_table %>% select(Bacteria, Display_Name), by = c("bacteria_name" = "Bacteria")) %>%
    pivot_longer(cols = ends_with("best_combo"), names_to = "criterion", values_to = "best_combo") %>%
    mutate(criterion = gsub("_best_combo", "", criterion),
           criterion = factor(criterion, levels = c("AIC", "BIC", "DevExp", "CVRMSE"),
                             labels = c("AIC", "BIC", "Explained Deviance", "CV-RMSE")))
  
  p <- ggplot(consistency_data, aes(x = criterion, y = Display_Name, fill = best_combo)) +
    geom_tile(color = "white") +
    geom_text(aes(label = best_combo), color = "white", fontface = "bold") +
    labs(
      title = "Consistency Analysis: Optimal Lag Combinations Across Evaluation Criteria",
      subtitle = "Values represent lag combinations (Temperature,Precipitation,Humidity,Wet Days) in years",
      x = "Evaluation Criterion",
      y = NULL,
      fill = "Optimal Combination"
    ) +
    scale_fill_viridis_d(option = "turbo") +
    theme_publication() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none"
    )
  
  # Save image
  ggsave(file.path(output_dir, "criteria_consistency_plot_four_factors.pdf"), 
         p, width = 10, height = 6, dpi = 300)
  
  return(p)
}

# 9. Create partial effects plot for the best model
create_partial_effects_plot_four_factors <- function(all_results, summary_table) {
  # Create directory for partial effects plots
  effects_dir <- file.path(output_dir, "partial_effects_four_factors")
  dir.create(effects_dir, showWarnings = FALSE, recursive = TRUE)
  
  for(i in 1:nrow(summary_table)) {
    bacteria_name <- summary_table$Bacteria[i]
    display_name <- summary_table$Display_Name[i]
    
    # Get best model
    if(!bacteria_name %in% names(all_results)) {
      debug_print(sprintf("No results for %s, skipping partial effects plot", bacteria_name), "WARNING")
      next
    }
    
    best_model <- all_results[[bacteria_name]]$best_model
    if(is.null(best_model)) {
      debug_print(sprintf("No best model for %s, skipping partial effects plot", bacteria_name), "WARNING")
      next
    }
    
    # Create partial effects plots
    tryCatch({
      # Create filenames with full path
      climate_effects_file <- file.path(effects_dir, paste0(bacteria_name, "_climate_effects.pdf"))
      
      # First check if model smooths include the expected terms
      model_terms <- names(best_model$smooth)
      climate_terms <- grep("(TMP|PREC|HUM|WET)_scaled_lag", model_terms, value = TRUE)
      
      if(length(climate_terms) > 0) {
        # Create climate effects plot
        pdf(climate_effects_file, width = 10, height = 10, onefile = TRUE)
        
        # Set up layout for multi-panel plot
        par(mfrow = c(2, 2))
        
        # Plot climate variable effects
        term_labels <- c(
          "TMP_scaled_lag" = "Temperature",
          "PREC_scaled_lag" = "Precipitation",
          "HUM_scaled_lag" = "Humidity",
          "WET_scaled_lag" = "Wet Days"
        )
        
        for(term in climate_terms) {
          term_name <- term_labels[grep(paste0("^", term, "$"), names(term_labels))]
          if(length(term_name) == 0) term_name <- term
          
          main_title <- paste0(display_name, " - ", term_name, " Effect")
          x_label <- paste0("Standardized ", term_name)
          
          plot(best_model, select = which(names(best_model$smooth) == term),
               shade = TRUE, main = main_title, xlab = x_label,
               cex.lab = 1.2, cex.main = 1.3)
        }
        
        dev.off()
        
        debug_print(sprintf("Climate effects plot for %s saved to %s", 
                           display_name, climate_effects_file), "INFO")
      } else {
        debug_print(sprintf("No climate terms found for %s", display_name), "WARNING")
      }
      
    }, error = function(e) {
      debug_print(sprintf("Error creating partial effects plot for %s: %s", 
                         display_name, e$message), "ERROR")
    })
  }
}

# 10. Create a correlation analysis plot of climate variables
create_climate_correlation_plot <- function(bacteria_name, pls_components) {
  # Set data file path
  data_file <- file.path(
    base_path,
    paste0(
      bacteria_name,
      "_data_with_PLS_components_",
      pls_components,
      "_with_MERRA2_and_ERA5_climate.csv"
    )
  )
  
  # Load data
  if(!file.exists(data_file)) {
    debug_print(paste("Data file not found:", data_file), "ERROR")
    return(NULL)
  }
  
  data <- read.csv(data_file)
  
  # Create correlation plot if corrplot package is available
  if(requireNamespace("corrplot", quietly = TRUE)) {
    # Extract climate variables
    climate_data <- data[, c("TMP", "PREC", "HUM", "WET")]
    climate_data <- climate_data[complete.cases(climate_data), ]
    
    # Calculate correlation matrix
    cor_matrix <- cor(climate_data, use = "pairwise.complete.obs")
    
    # Create correlation plot
    pdf(file.path(output_dir, "climate_variables_correlation.pdf"), width = 8, height = 8)
    corrplot::corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
                     addCoef.col = "black", tl.col = "black", tl.srt = 45,
                     diag = FALSE)
    dev.off()
    
    # Also output the correlation values
    write.csv(cor_matrix, file.path(output_dir, "climate_variables_correlation.csv"))
    
    return(cor_matrix)
  } else {
    debug_print("corrplot package not available for correlation visualization", "WARNING")
    return(NULL)
  }
}

# 11. Create lag analysis report for four factors
create_lag_report_four_factors <- function(summary_table, output_dir) {
  # Calculate lag distribution
  tmp_counts <- table(summary_table$TMP_lag)
  prec_counts <- table(summary_table$PREC_lag)
  hum_counts <- table(summary_table$HUM_lag)
  wet_counts <- table(summary_table$WET_lag)
  
  # Create report file
  report_file <- file.path(output_dir, "lag_analysis_report_four_factors.md")
  con <- file(report_file, "w", encoding = "UTF-8")
  
  writeLines("# Climate Lag Analysis Report for Six Bacterial Species (Four Climate Factors)\n", con)
  
  writeLines("## Research Background\n", con)
  writeLines("When studying the relationship between climate factors and antimicrobial resistance (AMR), identifying appropriate climate variable lag structures is crucial for revealing true temporal relationships. This analysis extends previous work by including wet days frequency as a fourth climate variable, alongside temperature, precipitation, and relative humidity.\n", con)
  
  writeLines("## Research Methods\n", con)
  writeLines("For each bacterial species, we systematically tested all combinations of 1-3 year lags for temperature, precipitation, humidity, and wet days variables (81 total combinations), using Generalized Additive Mixed Models (GAMM). To handle correlation between variables (particularly between wet days and other humidity variables), we employed select=TRUE parameter in mgcv which applies additional penalization to smooth terms. Models were evaluated based on AIC, explained deviance, cross-validation RMSE, and multicollinearity diagnostics.\n", con)
  
  writeLines("## Optimal Lag Combination Summary\n", con)
  writeLines("| Bacterial Species | Temp. Lag | Precip. Lag | Humidity Lag | Wet Days Lag | AIC | Expl. Dev. (%) | CV-RMSE | Max VIF |\n", con)
  writeLines("|----------|-------|-------|-------|-------|-----|-----------|--------|--------|", con)
  
  for(i in 1:nrow(summary_table)) {
    writeLines(sprintf("| %s | %d | %d | %d | %d | %.2f | %.2f | %.4f | %.2f |", 
                      summary_table$Display_Name[i], 
                      summary_table$TMP_lag[i], 
                      summary_table$PREC_lag[i], 
                      summary_table$HUM_lag[i],
                      summary_table$WET_lag[i],
                      summary_table$AIC[i], 
                      summary_table$Dev_explained[i], 
                      summary_table$CV_RMSE[i],
                      summary_table$Max_VIF[i]), con)
  }
  
  writeLines("\n## Lag Pattern Analysis\n", con)
  
  writeLines("### Temperature Lag Distribution\n", con)
  for(i in 1:3) {
    count <- tmp_counts[as.character(i)]
    if(is.na(count)) count <- 0
    writeLines(sprintf("- %d year lag: %d bacterial species", i, count), con)
  }
  
  writeLines("\n### Precipitation Lag Distribution\n", con)
  for(i in 1:3) {
    count <- prec_counts[as.character(i)]
    if(is.na(count)) count <- 0
    writeLines(sprintf("- %d year lag: %d bacterial species", i, count), con)
  }
  
  writeLines("\n### Humidity Lag Distribution\n", con)
  for(i in 1:3) {
    count <- hum_counts[as.character(i)]
    if(is.na(count)) count <- 0
    writeLines(sprintf("- %d year lag: %d bacterial species", i, count), con)
  }
  
  writeLines("\n### Wet Days Lag Distribution\n", con)
  for(i in 1:3) {
    count <- wet_counts[as.character(i)]
    if(is.na(count)) count <- 0
    writeLines(sprintf("- %d year lag: %d bacterial species", i, count), con)
  }
  
  writeLines("\n## Multicollinearity Assessment\n", con)
  writeLines("The Maximum Variance Inflation Factor (VIF) was calculated for each model to assess potential multicollinearity issues between the four climate variables. Values below 5 generally indicate acceptable levels of correlation between predictors.\n", con)
  
  max_vif_table <- summary_table$Max_VIF
  names(max_vif_table) <- summary_table$Display_Name
  
  writeLines("| Bacterial Species | Maximum VIF |\n", con)
  writeLines("|------------------|------------|", con)
  for(i in 1:nrow(summary_table)) {
    writeLines(sprintf("| %s | %.2f |", summary_table$Display_Name[i], summary_table$Max_VIF[i]), con)
  }
  
  writeLines("\n## Conclusions\n", con)
  writeLines("Our analysis indicates that different bacterial species exhibit different temporal lag patterns in response to all four climate factors. The addition of wet days frequency as a fourth climate variable provides complementary information to the other climate factors, despite modest correlation with precipitation and humidity. The optimal lag combinations identified will be used in subsequent AMR modeling to more accurately capture climate-resistance relationships.", con)
  
  close(con)
  
  # Try rendering to PDF if rmarkdown is available
  if(requireNamespace("rmarkdown", quietly = TRUE)) {
    try(
      rmarkdown::render(
        report_file,
        output_format = "pdf_document",
        output_file = file.path(output_dir, "lag_analysis_report_four_factors.pdf")
      ),
      silent = TRUE
    )
  } else {
    debug_print("rmarkdown package not installed, cannot convert to PDF", "WARNING")
  }
}

#------------------------------------------------------------------------
# Main function - Run analysis and create all visualizations
#------------------------------------------------------------------------
main <- function() {
  # Display startup information
  cat("\n====================================================\n")
  cat(" Climate Lag Combination Analysis (Four Climate Factors) for Six Bacterial Species\n")
  cat("====================================================\n\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check for correlation between climate variables
  debug_print("Analyzing correlation between climate variables...")
  create_climate_correlation_plot("3GCREC", 3)
  
  # Run all analyses
  debug_print("Running lag combination analysis for all bacteria with four climate factors...")
  analysis_results <- analyze_all_bacteria_four_factors(max_lag = max_lag_to_run)
  all_results <- analysis_results$all_results
  all_metrics <- analysis_results$all_metrics
  summary_table <- analysis_results$summary_table
  
  # Create all visualizations
  debug_print("Creating faceted AIC heatmap for four factors...")
  create_facet_aic_heatmap_four_factors(all_metrics, summary_table)
  
  debug_print("Creating multi-criteria comparison plot for four factors...")
  create_multicriteria_plot_four_factors(all_metrics, summary_table)
  
  debug_print("Creating best lag visualization for four factors...")
  create_best_lag_visualization_four_factors(summary_table)
  
  debug_print("Creating lag distribution heatmap for four factors...")
  create_lag_heatmap_four_factors(summary_table)
  
  debug_print("Creating AIC rank comparison plot for four factors...")
  create_aic_rank_plot_four_factors(all_metrics, summary_table)
  
  debug_print("Creating criteria consistency plot for four factors...")
  create_criteria_consistency_plot_four_factors(all_metrics, summary_table)
  
  debug_print("Creating lag distribution bar plot for four factors...")
  create_lag_dist_barplot_four_factors(summary_table)
  
  debug_print("Creating best model partial effects plots for four factors...")
  create_partial_effects_plot_four_factors(all_results, summary_table)
  
  debug_print("Creating summary table for four factors...")
  create_summary_table_with_metrics_four_factors(summary_table)
  
  debug_print("Creating analysis report for four factors...")
  create_lag_report_four_factors(summary_table, output_dir)
  
  debug_print("All analyses and visualizations for four climate factors completed!")
  cat(paste("\nAnalysis complete!\nAll results saved to:", output_dir, "\n"))
  
  # Return main result objects
  return(list(
    all_results = all_results,
    all_metrics = all_metrics,
    summary_table = summary_table
  ))
}

# Run main function
results_four_factors <- main()
