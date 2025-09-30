# ============================================================================
# GAMM Analysis for Climate-Bacteria Resistance Relationships

# Required Libraries --------------------------------------------------------
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

# Configuration --------------------------------------------------------------

# Color scheme for climate variables
climate_colors <- c(
  "Temperature" = "#DD5F60",
  "Precipitation" = "#9AC0CD", 
  "Humidity" = "#3CB371",
  "WetDays" = "#8A2BE2"
)

# Lag time settings for different bacteria-climate combinations
# Based on model selection and cross-validation results
lag_settings <- list(
  "3GCR_Ec" = list(temp_lag = 2, precip_lag = 3, humid_lag = 2, wetdays_lag = 1),
  "3GCR_Kp" = list(temp_lag = 3, precip_lag = 3, humid_lag = 1, wetdays_lag = 1),
  "CR_Ab" = list(temp_lag = 2, precip_lag = 3, humid_lag = 3, wetdays_lag = 1),
  "CR_Ec" = list(temp_lag = 3, precip_lag = 3, humid_lag = 1, wetdays_lag = 2),
  "CR_Kp" = list(temp_lag = 3, precip_lag = 3, humid_lag = 2, wetdays_lag = 3),
  "CR_Pa" = list(temp_lag = 2, precip_lag = 1, humid_lag = 3, wetdays_lag = 1)
)

# Default lag settings for new bacteria
default_lag_settings <- list(
  temp_lag = 2, 
  precip_lag = 3, 
  humid_lag = 2, 
  wetdays_lag = 1
)

# Range settings for visualization [min, max, step]
range_settings <- list(
  temp = c(-10, 40, 10),
  humid = c(30, 100, 10),
  precip = c(0, 3200, 500),
  wetdays = c(0, 300, 50)
)

# Effect size thresholds for threshold detection
effect_thresholds <- list(
  base_threshold = 0.06,
  temperature_modifier = 0.9,
  humidity_modifier = 1.0,
  precipitation_modifier = 1.2,
  wetdays_modifier = 0.95
)

# Helper Functions -----------------------------------------------------------

check_climate_correlations <- function(data) {
  climate_vars <- data %>%
    select(TMP_orig, PREC_orig, HUM_orig, WET_orig) %>%
    na.omit()
  
  if(ncol(climate_vars) == 0 || nrow(climate_vars) == 0) {
    cat("No climate variables found or insufficient data\n")
    return(NULL)
  }
  
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

# Data Preparation Functions ------------------------------------------------

prepare_data <- function(file_path, bacteria_name, lag_config = NULL) {
  
  # Get lag configuration
  if(is.null(lag_config)) {
    if(bacteria_name %in% names(lag_settings)) {
      lag_config <- lag_settings[[bacteria_name]]
    } else {
      lag_config <- default_lag_settings
      cat("Using default lag settings for", bacteria_name, "\n")
    }
  }
  
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  data <- read.csv(file_path)
  
  # Check required columns
  required_cols <- c("year", "lat", "lon", "NAME", "Region", "TMP", "PREC", "HUM", "WET", "logit_R")
  missing_cols <- setdiff(required_cols, names(data))
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
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
  
  # Ensure humidity does not exceed 100%
  data_processed <- data_processed %>%
    mutate(HUM = pmin(HUM, 100))
  
  # Calculate scaling parameters
  scale_params <- data_processed %>%
    summarise(
      across(c(TMP, PREC, HUM, WET), list(
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE)
      ))
    )
  
  # Process data with lag settings
  data_final <- data_processed %>%
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
    group_by(location_id) %>%
    arrange(year) %>%
    mutate(
      TMP_scaled_lag = lag(TMP_scaled, lag_config$temp_lag),
      PREC_scaled_lag = lag(PREC_scaled, lag_config$precip_lag),
      HUM_scaled_lag = lag(HUM_scaled, lag_config$humid_lag),
      WET_scaled_lag = lag(WET_scaled, lag_config$wetdays_lag)
    ) %>%
    filter(!is.na(TMP_scaled_lag) & !is.na(PREC_scaled_lag) & 
             !is.na(HUM_scaled_lag) & !is.na(WET_scaled_lag)) %>%
    ungroup()
  
  cat("\nData preparation completed for", bacteria_name, "\n")
  cat("Sample size after processing:", nrow(data_final), "\n")
  cat("Lag configuration used:\n")
  cat("  Temperature lag:", lag_config$temp_lag, "years\n")
  cat("  Precipitation lag:", lag_config$precip_lag, "years\n")
  cat("  Humidity lag:", lag_config$humid_lag, "years\n")
  cat("  Wet days lag:", lag_config$wetdays_lag, "years\n")
  
  cat("\nChecking climate variable correlations for", bacteria_name, "...\n")
  cor_matrix <- check_climate_correlations(data_final)
  
  return(list(
    data = data_final, 
    scale_params = scale_params, 
    cor_matrix = cor_matrix,
    lag_config = lag_config
  ))
}

# Model Building Functions ---------------------------------------------------

build_gamm_model <- function(data, bacteria_name, output_path) {
  ctrl <- gam.control(nthreads = 4, maxit = 1000, mgcv.tol = 1e-7, mgcv.half = 15)
  
  cat("Building GAMM model for", bacteria_name, "...\n")
  
  # Check for PLS components
  pls_terms <- c()
  for(i in 1:4) {
    comp_name <- paste0("PLS_Comp", i)
    if(comp_name %in% names(data)) {
      pls_terms <- c(pls_terms, paste0("s(", comp_name, ", k = 10, bs = \"cr\")"))
    }
  }
  
  # Build formula
  base_formula <- "logit_R ~ s(TMP_scaled_lag, k = 5, bs = \"cr\") +
    s(PREC_scaled_lag, k = 10, bs = \"cr\") +
    s(HUM_scaled_lag, k = 10, bs = \"cr\") +
    s(WET_scaled_lag, k = 10, bs = \"cr\") +
    s(lat, lon, bs = \"sos\", k = 20) +
    s(year, bs = \"cr\", k = 8) +
    s(Region, bs = \"re\") +
    climate_zone"
  
  if(length(pls_terms) > 0) {
    full_formula <- paste(base_formula, "+", paste(pls_terms, collapse = " + "))
  } else {
    full_formula <- base_formula
  }
  
  tryCatch({
    model <- bam(
      as.formula(full_formula),
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE,
      control = ctrl
    )
  }, error = function(e) {
    warning(paste("bam() failed for", bacteria_name, ", trying gam() instead:", e$message))
    model <- gam(
      as.formula(full_formula),
      data = data,
      family = gaussian(),
      method = "REML",
      select = TRUE
    )
  })
  
  # Save model summary
  summary_file <- file.path(output_path, paste0(bacteria_name, "_model_summary.txt"))
  sink(summary_file)
  cat("====================================\n")
  cat("GAMM Model Summary for", bacteria_name, "\n")
  cat("====================================\n\n")
  cat("Model formula:\n", full_formula, "\n\n")
  cat("Sample size:", nrow(data), "\n\n")
  cat("Model summary:\n")
  print(summary(model))
  cat("\nModel validation statistics:\n")
  print(gam.check(model))
  cat("\n====================================\n")
  sink()
  
  cat("Model summary saved to:", summary_file, "\n")
  
  return(model)
}

# Threshold Detection Functions ----------------------------------------------

detect_thresholds <- function(model, data, scale_params, var_name, bacteria_name) {
  
  orig_var <- str_extract(var_name, "^[A-Z]+")
  var_mean <- scale_params[[paste0(orig_var, "_mean")]]
  var_sd <- scale_params[[paste0(orig_var, "_sd")]]
  
  if(is.na(var_mean) || is.na(var_sd) || var_sd == 0) {
    stop(paste("Invalid scaling parameters for variable", orig_var))
  }
  
  # Determine variable type
  is_precipitation <- grepl("PREC", var_name)
  is_humidity <- grepl("HUM", var_name)
  is_temperature <- grepl("TMP", var_name)
  is_wetdays <- grepl("WET", var_name)
  
  climate_variable <- ifelse(
    is_temperature, "Temperature",
    ifelse(is_humidity, "Humidity",
           ifelse(is_precipitation, "Precipitation", "WetDays"))
  )
  
  # Get range settings
  range_key <- ifelse(is_precipitation, "precip",
                      ifelse(is_humidity, "humid",
                             ifelse(is_wetdays, "wetdays", "temp")))
  var_range <- range_settings[[range_key]]
  var_min <- var_range[1]
  var_max <- var_range[2]
  
  # Set prediction range
  max_pred_value <- if(is_humidity) {
    100
  } else if(is_precipitation) {
    3200
  } else if(is_wetdays) {
    300
  } else {
    var_max
  }
  
  # Calculate scaled value range
  min_scaled_value <- (var_min - var_mean) / var_sd
  max_scaled_value <- (max_pred_value - var_mean) / var_sd
  
  # Generate prediction points
  pred_x <- seq(min_scaled_value, max_scaled_value, length.out = 1200)
  orig_x <- pred_x * var_sd + var_mean
  
  # Handle non-negative constraints
  if (is_precipitation || is_wetdays) {
    orig_x[orig_x <= 0] <- 0.01
    if(is_precipitation) {
      orig_x <- pmin(orig_x, 3200)
    }
    if(is_wetdays) {
      orig_x <- pmin(orig_x, 300)
    }
  }
  
  # Handle humidity constraint
  if (is_humidity) {
    orig_x <- pmin(orig_x, 100)
  }
  
  # Prepare prediction data
  pred_data <- data.frame(
    TMP_scaled_lag = mean(data$TMP_scaled_lag, na.rm = TRUE),
    PREC_scaled_lag = mean(data$PREC_scaled_lag, na.rm = TRUE),
    HUM_scaled_lag = mean(data$HUM_scaled_lag, na.rm = TRUE),
    WET_scaled_lag = mean(data$WET_scaled_lag, na.rm = TRUE),
    lat = mean(data$lat, na.rm = TRUE),
    lon = mean(data$lon, na.rm = TRUE),
    year = mean(data$year, na.rm = TRUE),
    climate_zone = factor(names(which.max(table(data$climate_zone))), 
                          levels = levels(data$climate_zone)),
    Region = factor(names(which.max(table(data$Region))), 
                    levels = levels(data$Region))
  )
  
  # Add PLS components if available
  for(i in 1:4) {
    comp_name <- paste0("PLS_Comp", i)
    if(comp_name %in% names(data)) {
      pred_data[[comp_name]] <- mean(data[[comp_name]], na.rm = TRUE)
    }
  }
  
  pred_data <- pred_data[rep(1, length(pred_x)), ]
  pred_data[[var_name]] <- pred_x
  
  # Make predictions
  suppressWarnings({
    pred <- predict(model, pred_data, type = "terms", se.fit = TRUE)
  })
  
  var_col <- grep(var_name, colnames(pred$fit))
  if(length(var_col) == 0) {
    stop(paste("Variable", var_name, "not found in model terms"))
  }
  
  fit_values <- pred$fit[, var_col]
  
  # Calculate confidence intervals
  lower_ci_99 <- fit_values - 2.576 * pred$se.fit[, var_col]
  upper_ci_99 <- fit_values + 2.576 * pred$se.fit[, var_col]
  lower_ci_95 <- fit_values - 1.96 * pred$se.fit[, var_col]
  upper_ci_95 <- fit_values + 1.96 * pred$se.fit[, var_col]
  lower_ci_90 <- fit_values - 1.645 * pred$se.fit[, var_col]
  upper_ci_90 <- fit_values + 1.645 * pred$se.fit[, var_col]
  
  # Calculate OR values and CIs
  or_values <- exp(fit_values)
  or_lower_99 <- exp(lower_ci_99)
  or_upper_99 <- exp(upper_ci_99)
  or_lower_95 <- exp(lower_ci_95)
  or_upper_95 <- exp(upper_ci_95)
  or_lower_90 <- exp(lower_ci_90)
  or_upper_90 <- exp(upper_ci_90)
  
  # Test for linearity
  is_linear <- FALSE
  lm_slope <- NA
  lm_p <- NA
  lm_r2 <- NA
  
  tryCatch({
    lm_model <- lm(or_values ~ orig_x)
    lm_slope <- coef(lm_model)[2]
    lm_p <- summary(lm_model)$coefficients[2, 4]
    lm_r2 <- summary(lm_model)$r.squared
    
    gam_fit <- gam(or_values ~ s(orig_x, k = 10, bs = "cr"))
    gam_r2 <- summary(gam_fit)$r.sq
    
    # Linearity test
    is_linear <- (gam_r2 - lm_r2) < 0.03 && lm_p < 0.05 && lm_r2 > 0.15
  }, error = function(e) {
    is_linear <- FALSE
  })
  
  # Apply LOESS smoothing
  loess_fit <- function(y, span = 0.25) {
    x <- 1:length(y)
    loess_model <- loess(y ~ x, span = span)
    return(predict(loess_model))
  }
  
  # Calculate adaptive smoothing parameter
  rough_second_deriv <- diff(diff(or_values))
  deriv_variance <- var(rough_second_deriv, na.rm = TRUE)
  span_value <- min(0.5, max(0.1, 0.3 - 0.2 * sqrt(deriv_variance)))
  
  # Apply smoothing
  smooth_or <- loess_fit(or_values, span = span_value)
  smooth_lower_99 <- loess_fit(or_lower_99, span = span_value)
  smooth_upper_99 <- loess_fit(or_upper_99, span = span_value)
  smooth_lower_95 <- loess_fit(or_lower_95, span = span_value)
  smooth_upper_95 <- loess_fit(or_upper_95, span = span_value)
  smooth_lower_90 <- loess_fit(or_lower_90, span = span_value)
  smooth_upper_90 <- loess_fit(or_upper_90, span = span_value)
  
  relationship_type <- ifelse(is_linear, "Linear", "Threshold")
  
  # Create empty threshold tibble
  empty_threshold_tibble <- tibble(
    Bacteria = character(),
    Climate_Variable = character(),
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
  
  # Return linear relationship data
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
      relationship_type = relationship_type,
      linear_info = list(
        slope = lm_slope,
        p_value = lm_p,
        r_squared = lm_r2
      )
    )
    return(result)
  }
  
  # Threshold detection for non-linear relationships
  first_deriv <- diff(smooth_or) / diff(orig_x)
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
  
  # Set effect threshold based on variable type
  base_threshold <- effect_thresholds$base_threshold
  effect_threshold <- base_threshold * if(is_temperature) {
    effect_thresholds$temperature_modifier
  } else if(is_humidity) {
    effect_thresholds$humidity_modifier
  } else if(is_precipitation) {
    effect_thresholds$precipitation_modifier
  } else if(is_wetdays) {
    effect_thresholds$wetdays_modifier
  } else {
    1.0
  }
  
  # Find threshold points
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
    
    # Find global extrema
    max_idx <- which.max(smooth_or)
    min_idx <- which.min(smooth_or)
    
    # Find local extrema and inflection points
    abs_second_deriv <- abs(second_deriv_df$deriv)
    abs_first_deriv <- abs(first_deriv_df$deriv)
    
    potential_turning_points <- which(diff(sign(first_deriv)) != 0)
    sensitivity_threshold <- 0.3
    
    important_turning_points <- potential_turning_points[
      abs_first_deriv[potential_turning_points] > quantile(abs_first_deriv, sensitivity_threshold) |
        abs_second_deriv[potential_turning_points-1] > quantile(abs_second_deriv, sensitivity_threshold)
    ]
    
    # Function to add threshold point
    add_threshold_point <- function(idx, point_type_base) {
      if (idx >= 1 && idx <= length(smooth_or)) {
        effect_size <- abs(smooth_or[idx] - 1)
        
        # Check significance levels
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
          return()
        }
        
        # Calculate priority score
        type_weight <- case_when(
          grepl("Global", point_type_base) ~ 4,
          grepl("Maximum|Minimum", point_type_base) ~ 3,
          grepl("Rapid change", point_type_base) ~ 2,
          grepl("Stable", point_type_base) ~ 1,
          TRUE ~ 0
        )
        
        effect_weight <- min(5, max(1, effect_size * 10))
        priority_score <- sig_weight * 4 + type_weight * 2 + effect_weight
        
        point_type <- if (ci_significant_95) point_type_base else paste("Trend", point_type_base)
        
        # Add to candidates
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
    
    # Add different types of points
    add_threshold_point(max_idx, "Global maximum")
    add_threshold_point(min_idx, "Global minimum")
    
    # Add local extrema
    if (length(important_turning_points) > 0) {
      for (i in important_turning_points) {
        if (i > 1 && i < length(first_deriv)) {
          if (first_deriv[i-1] > 0 && first_deriv[i+1] < 0) {
            add_threshold_point(i, "Local maximum")
          } else if (first_deriv[i-1] < 0 && first_deriv[i+1] > 0) {
            add_threshold_point(i, "Local minimum")
          }
        }
      }
    }
    
    # Add rapid change points
    rapid_changes <- first_deriv_df %>%
      filter(abs(deriv) > quantile(abs_first_deriv, 0.90)) %>%
      arrange(desc(abs(deriv)))
    
    if (nrow(rapid_changes) > 0) {
      for (i in 1:min(3, nrow(rapid_changes))) {
        rc_idx <- rapid_changes$idx[i]
        add_threshold_point(rc_idx, "Rapid change")
      }
    }
    
    # Add stable regions
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
      return(empty_threshold_tibble %>% select(-c(Bacteria, Climate_Variable, Relationship_Type, Value, OR, Lower_CI, Upper_CI, Smooth_Parameter)))
    }
    
    # Sort by priority and apply distance constraints
    sorted_candidates <- all_threshold_candidates %>%
      arrange(desc(priority_score))
    
    max_points <- 5
    min_dist_fraction <- 0.1
    
    final_points <- tibble()
    x_range <- max(orig_x) - min(orig_x)
    
    # Select final points with distance constraints
    for (i in 1:nrow(sorted_candidates)) {
      current_x <- sorted_candidates$x_orig[i]
      
      if (nrow(final_points) == 0 || 
          !any(abs(final_points$x_orig - current_x) < (x_range * min_dist_fraction))) {
        final_points <- bind_rows(final_points, sorted_candidates[i,])
      }
      
      if (nrow(final_points) >= max_points) break
    }
    
    # Ensure minimum points
    if (nrow(final_points) < 2 && nrow(sorted_candidates) >= 2) {
      final_points <- sorted_candidates[1:min(2, nrow(sorted_candidates)),]
    }
    
    return(final_points %>% select(-priority_score))
  }
  
  # Find threshold points
  threshold_points <- find_threshold_points()
  
  # Handle constraints for specific variable types
  if ((is_precipitation || is_wetdays) && nrow(threshold_points) > 0) {
    threshold_points <- threshold_points %>%
      mutate(x_orig = ifelse(x_orig <= 0, 0.01, x_orig))
  }
  
  if (is_humidity && nrow(threshold_points) > 0) {
    threshold_points <- threshold_points %>%
      mutate(x_orig = pmin(x_orig, 100))
  }
  
  # Add required information for output
  if (nrow(threshold_points) > 0) {
    threshold_points <- threshold_points %>%
      mutate(
        Bacteria = bacteria_name,
        Climate_Variable = climate_variable,
        Relationship_Type = relationship_type,
        Value = x_orig,
        OR = y,
        Lower_CI = lower_ci_95,
        Upper_CI = upper_ci_95,
        Smooth_Parameter = span_value
      )
  }
  
  # Create result
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
    relationship_type = relationship_type,
    linear_info = list(
      slope = lm_slope,
      p_value = lm_p,
      r_squared = lm_r2
    )
  )
  
  return(result)
}

# Plotting Functions ---------------------------------------------------------

create_climate_effect_plot <- function(model, data, scale_params, var_name, 
                                       x_lab, title, color, bacteria_name, 
                                       threshold_data, lag_config = NULL) {
  
  if(is.null(lag_config)) {
    if(bacteria_name %in% names(lag_settings)) {
      lag_config <- lag_settings[[bacteria_name]]
    } else {
      lag_config <- default_lag_settings
    }
  }
  
  curve_data <- threshold_data$curve_data
  threshold_points <- threshold_data$threshold_points
  is_linear <- threshold_data$is_linear
  
  orig_var <- str_extract(var_name, "^[A-Z]+")
  
  # Get lag information for axis label
  lag_text <- if(grepl("TMP", var_name)) {
    paste("(lag", lag_config$temp_lag, "yr)")
  } else if(grepl("PREC", var_name)) {
    paste("(lag", lag_config$precip_lag, "yr)")
  } else if(grepl("HUM", var_name)) {
    paste("(lag", lag_config$humid_lag, "yr)")
  } else if(grepl("WET", var_name)) {
    paste("(lag", lag_config$wetdays_lag, "yr)")
  } else {
    ""
  }
  
  x_lab_with_lag <- paste(x_lab, lag_text)
  
  # Determine variable type
  is_precipitation <- grepl("PREC", var_name)
  is_humidity <- grepl("HUM", var_name)
  is_temperature <- grepl("TMP", var_name)
  is_wetdays <- grepl("WET", var_name)
  
  # Get range settings
  range_key <- ifelse(is_precipitation, "precip",
                      ifelse(is_humidity, "humid",
                             ifelse(is_wetdays, "wetdays", "temp")))
  var_range <- range_settings[[range_key]]
  var_min <- var_range[1]
  var_max <- var_range[2]
  step <- var_range[3]
  
  # Add margin for visualization
  x_range <- var_max - var_min
  margin_size <- x_range * 0.04
  var_min <- var_min - margin_size
  var_max <- var_max + margin_size
  
  # Apply variable-specific constraints
  if (is_humidity) {
    curve_data <- curve_data %>%
      filter(x_orig <= 100) %>%
      mutate(x_orig = pmin(x_orig, 100))
    var_max <- min(104, var_max)
  }
  
  if (is_precipitation || is_wetdays) {
    var_min <- max(-margin_size, var_min)
    if (is_precipitation) {
      var_max <- 3200 + margin_size
      curve_data <- curve_data %>% filter(x_orig <= 3200)
    } else if (is_wetdays) {
      var_max <- 300 + margin_size
      curve_data <- curve_data %>% filter(x_orig <= 300)
    }
  }
  
  # Create x-axis breaks
  x_breaks <- if (is_precipitation) {
    seq(0, 3200, by = 800)
  } else if (is_humidity) {
    seq(30, 100, by = 10)
  } else if (is_wetdays) {
    seq(0, 300, by = 50)
  } else {
    seq(round(var_min/10)*10, round(var_max/10)*10, by = step)
  }
  
  # Prepare density data for visualization
  if (is_precipitation || is_wetdays) {
    pos_data <- data[[paste0(orig_var, "_orig")]]
    pos_data <- pos_data[pos_data > 0 & !is.na(pos_data)]
    
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
    orig_data <- data[[paste0(orig_var, "_orig")]]
    orig_data <- orig_data[!is.na(orig_data)]
    
    if (length(orig_data) > 10) {
      dens <- density(orig_data, na.rm = TRUE, adjust = 1.1)
      density_data <- tibble(x = dens$x, density = dens$y)
      
      density_data <- density_data %>%
        filter(x >= var_min & x <= var_max)
      
      if (is_humidity) {
        density_data <- density_data %>% filter(x <= 100)
      }
    } else {
      breaks <- seq(var_min, var_max, length.out = 30)
      hist_data <- hist(orig_data, breaks = breaks, plot = FALSE)
      density_data <- tibble(x = hist_data$mids, density = hist_data$density)
    }
  }
  
  # Set y-axis range
  y_min <- floor(min(c(curve_data$lower_ci, 0.5)) * 4) / 4
  y_max <- ceiling(max(c(curve_data$upper_ci, 1.5)) * 4) / 4
  y_breaks <- seq(y_min, y_max, length.out = 5)
  
  # Filter curve data to plotting range
  curve_data_filtered <- curve_data %>%
    filter(x_orig >= var_min, x_orig <= var_max)
  
  # Scale density for visualization
  if(nrow(density_data) > 0) {
    max_density <- max(density_data$density, na.rm = TRUE)
    if(max_density > 0) {
      density_data <- density_data %>%
        mutate(scaled_density = density / max_density * 0.46)
    } else {
      density_data$scaled_density <- 0
    }
  }
  
  # Create main effect plot
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
  
  # Add y-axis label for temperature plots only
  if (is_temperature) {
    main_plot <- main_plot +
      labs(y = paste0(bacteria_name, " OR (95% CI)"))
  } else {
    main_plot <- main_plot +
      labs(y = "")
  }
  
  # Add threshold points for non-linear relationships
  if (nrow(threshold_points) > 0 && !is_linear) {
    # Add threshold vertical lines
    for (i in 1:nrow(threshold_points)) {
      main_plot <- main_plot +
        geom_vline(xintercept = threshold_points$x_orig[i], 
                   color = "#FF0000", linetype = "dashed", linewidth = 0.4)
    }
    
    # Prepare threshold labels
    threshold_points <- threshold_points %>%
      mutate(
        label_short = case_when(
          grepl("Global maximum", type) ~ "GMax",
          grepl("Global minimum", type) ~ "GMin",
          grepl("Local maximum|Maximum", type) ~ "Max",
          grepl("Local minimum|Minimum", type) ~ "Min",
          grepl("Rapid change", type) ~ "RC",
          grepl("Stable", type) ~ "Stable",
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
    
    # Calculate label positions to avoid overlap
    if (nrow(threshold_points) > 0) {
      threshold_points <- threshold_points %>%
        arrange(x_orig)
      
      vjust_values <- numeric(nrow(threshold_points))
      hjust_values <- numeric(nrow(threshold_points))
      
      for (i in 1:nrow(threshold_points)) {
        # Vertical positioning based on OR value
        if (threshold_points$y[i] > 1) {
          vjust_values[i] <- -0.3  # Above curve
        } else {
          vjust_values[i] <- 1.3   # Below curve
        }
        
        # Horizontal positioning based on x position
        if (threshold_points$x_orig[i] < (var_min + x_range * 0.15)) {
          hjust_values[i] <- 0     # Left align
        } else if (threshold_points$x_orig[i] > (var_max - x_range * 0.15)) {
          hjust_values[i] <- 1     # Right align
        } else {
          hjust_values[i] <- 0.5   # Center align
        }
      }
      
      # Prevent label overlap by alternating positions
      if (nrow(threshold_points) > 1) {
        for (i in 2:nrow(threshold_points)) {
          if (abs(threshold_points$x_orig[i] - threshold_points$x_orig[i-1]) < (x_range * 0.2)) {
            if (vjust_values[i] * vjust_values[i-1] > 0) {
              vjust_values[i] <- -vjust_values[i-1]
            }
          }
        }
      }
      
      threshold_points$vjust <- vjust_values
      threshold_points$hjust <- hjust_values
      
      # Add threshold labels to plot
      main_plot <- main_plot +
        geom_text(
          data = threshold_points,
          aes(x = x_orig, y = y, label = label, hjust = hjust, vjust = vjust),
          size = 2.3,
          color = "black"
        )
    }
  }
  
  # Create density plot (bottom panel)
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
  
  # Combine main plot and density plot
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

# Main Analysis Functions ----------------------------------------------------

combine_all_bacterial_plots <- function(data_configs, output_path) {
  all_plots <- list()
  threshold_summary <- list()
  linear_results <- list()
  
  # Define climate factors for analysis
  climate_factors <- list(
    list(var = "TMP_scaled_lag", x_lab = "Temperature (°C)", color = climate_colors["Temperature"]),
    list(var = "HUM_scaled_lag", x_lab = "Relative humidity (%)", color = climate_colors["Humidity"]),
    list(var = "PREC_scaled_lag", x_lab = "Precipitation (mm)", color = climate_colors["Precipitation"]),
    list(var = "WET_scaled_lag", x_lab = "Wet Days (d)", color = climate_colors["WetDays"])
  )
  
  # Process each bacteria
  for (i in seq_along(data_configs)) {
    cat("Processing", data_configs[[i]]$bacteria_name, "...\n")
    
    # Prepare data
    data_result <- prepare_data(
      data_configs[[i]]$file_path, 
      data_configs[[i]]$bacteria_name,
      data_configs[[i]]$lag_config
    )
    processed_data <- data_result$data
    scale_params <- data_result$scale_params
    lag_config <- data_result$lag_config
    
    # Build model
    model <- build_gamm_model(processed_data, data_configs[[i]]$bacteria_name, output_path)
    
    # Process each climate factor
    for (factor in climate_factors) {
      cat(" Processing", data_configs[[i]]$bacteria_name, factor$var, "...\n")
      
      # Detect thresholds
      threshold_data <- detect_thresholds(
        model, processed_data, scale_params, factor$var, data_configs[[i]]$bacteria_name
      )
      
      # Store threshold summary
      threshold_summary[[paste(data_configs[[i]]$bacteria_name, factor$var, sep = "_")]] <- 
        threshold_data$threshold_points
      
      # Record linear relationship results
      climate_variable <- threshold_data$climate_variable
      linear_info <- threshold_data$linear_info
      if (is.null(linear_info)) {
        linear_info <- list(slope = NA, p_value = NA, r_squared = NA)
      }
      
      linear_results[[length(linear_results) + 1]] <- data.frame(
        Bacteria = data_configs[[i]]$bacteria_name,
        Climate_Variable = climate_variable,
        Relationship_Type = threshold_data$relationship_type,
        Slope = linear_info$slope,
        P_Value = linear_info$p_value,
        R_Squared = linear_info$r_squared,
        stringsAsFactors = FALSE
      )
      
      # Create plot
      plot <- create_climate_effect_plot(
        model, processed_data, scale_params, factor$var, factor$x_lab,
        data_configs[[i]]$bacteria_name, factor$color, data_configs[[i]]$bacteria_name, 
        threshold_data, lag_config
      )
      
      # Determine clean variable name and lag for file naming
      var_name_clean <- ifelse(
        factor$var == "TMP_scaled_lag", "Temperature",
        ifelse(factor$var == "HUM_scaled_lag", "Humidity",
               ifelse(factor$var == "PREC_scaled_lag", "Precipitation", "WetDays"))
      )
      
      lag_var <- paste0(tolower(str_extract(factor$var, "^[A-Z]+")), "_lag")
      lag <- lag_config[[lag_var]]
      
      # Save individual plots
      plot_path <- file.path(output_path, paste0(data_configs[[i]]$bacteria_name, "_", 
                                                 var_name_clean, "_lag", lag, ".pdf"))
      png_path <- file.path(output_path, paste0(data_configs[[i]]$bacteria_name, "_", 
                                                var_name_clean, "_lag", lag, ".png"))
      
      ggsave(plot_path, plot, width = 5, height = 4, dpi = 300)
      ggsave(png_path, plot, width = 5, height = 4, dpi = 300)
      
      cat(" Saved plot to:", plot_path, "\n")
      
      # Store plot for combined visualization
      all_plots[[paste(data_configs[[i]]$bacteria_name, var_name_clean, sep = "_")]] <- plot
    }
    
    # Create combined plot for each bacteria (2x2 grid)
    combined_bacteria_plot <- plot_grid(
      all_plots[[paste(data_configs[[i]]$bacteria_name, "Temperature", sep = "_")]],
      all_plots[[paste(data_configs[[i]]$bacteria_name, "Humidity", sep = "_")]],
      all_plots[[paste(data_configs[[i]]$bacteria_name, "Precipitation", sep = "_")]],
      all_plots[[paste(data_configs[[i]]$bacteria_name, "WetDays", sep = "_")]],
      ncol = 2,
      nrow = 2,
      align = "hv",
      labels = c("A", "B", "C", "D"),
      label_size = 12
    )
    
    # Save combined bacteria plot
    combined_path <- file.path(output_path, paste0(data_configs[[i]]$bacteria_name, "_combined.pdf"))
    ggsave(combined_path, combined_bacteria_plot, width = 10, height = 8, dpi = 300)
    cat("Saved combined plot for", data_configs[[i]]$bacteria_name, "to:", combined_path, "\n")
  }
  
  # Create final combined plot with all bacteria and climate factors
  title_row <- plot_grid(
    ggplot() + theme_void() + ggtitle("Temperature") + 
      theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ggplot() + theme_void() + ggtitle("Humidity") + 
      theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ggplot() + theme_void() + ggtitle("Precipitation") + 
      theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ggplot() + theme_void() + ggtitle("Wet Days") + 
      theme(plot.title = element_text(hjust = 0.5, size = 14)),
    ncol = 4,
    rel_heights = c(0.2)
  )
  
  # Create plot matrix
  plot_list <- list(title_row)
  bacteria_names <- sapply(data_configs, function(x) x$bacteria_name)
  
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
    plot_list[[i+1]] <- row_plots
  }
  
  # Combine all plots
  final_plot <- plot_grid(
    plotlist = plot_list,
    ncol = 1,
    rel_heights = c(0.08, rep(1, length(bacteria_names)))
  )
  
  # Save final combined plots
  final_path <- file.path(output_path, "all_bacteria_climate_factors.pdf")
  ggsave(final_path, final_plot, width = 16, height = 16, dpi = 300)
  
  png_path <- file.path(output_path, "all_bacteria_climate_factors.png")
  ggsave(png_path, final_plot, width = 16, height = 16, dpi = 300)
  
  cat("Saved final combined plot to:", final_path, "\n")
  
  # Save analysis summaries
  threshold_summary_df <- bind_rows(threshold_summary, .id = "bacteria_variable")
  write.csv(threshold_summary_df, file.path(output_path, "threshold_summary.csv"), row.names = FALSE)
  cat("Saved threshold summary to:", file.path(output_path, "threshold_summary.csv"), "\n")
  
  linear_summary_df <- bind_rows(linear_results)
  write.csv(linear_summary_df, file.path(output_path, "linear_relationships_summary.csv"), row.names = FALSE)
  cat("Saved linear relationships summary to:", file.path(output_path, "linear_relationships_summary.csv"), "\n")
  
  return(list(
    plots = all_plots,
    thresholds = threshold_summary,
    linear_relationships = linear_results
  ))
}

# Main Execution Function ---------------------------------------------------

run_gamm_analysis <- function(base_path, data_configs, output_dir = "gamm_results") {
  
  # Set output path
  output_path <- file.path(base_path, output_dir)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    cat("Created output directory:", output_path, "\n")
  }
  
  cat("Starting GAMM analysis...\n")
  cat("Base path:", base_path, "\n")
  cat("Output path:", output_path, "\n")
  cat("Number of bacteria to process:", length(data_configs), "\n")
  
  # Execute main analysis
  results <- combine_all_bacterial_plots(data_configs, output_path)
  
  cat("\n====================================\n")
  cat("Analysis completed successfully!\n")
  cat("Results saved to:", output_path, "\n")
  cat("====================================\n")
  
  return(results)
}

# Script Status Message ------------------------------------------------------

cat("GAMM Analysis Script Loaded Successfully!\n")
cat("Available bacteria lag configurations:\n")
for(bacteria in names(lag_settings)) {
  lag_config <- lag_settings[[bacteria]]
  cat(sprintf("  %s: Temp=%d, Precip=%d, Humid=%d, WetDays=%d\n", 
              bacteria, lag_config$temp_lag, lag_config$precip_lag, 
              lag_config$humid_lag, lag_config$wetdays_lag))
}
cat("\nUse run_gamm_analysis(base_path, data_configs) to start analysis.\n")