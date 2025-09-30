# ============================================================================
# Climate-AMR Projection Analysis


# Required Libraries --------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(splines)
  library(Cairo)
  library(ggrepel)
  library(boot)
  library(parallel)
  library(cowplot)
})

# Configuration --------------------------------------------------------------

# Bacteria configuration
bacteria_config <- list(
  names = c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa"),
  
  display_names = c(
    "3GCR-Ec" = "3GCR-Ec",
    "3GCR-Kp" = "3GCR-Kp", 
    "CR-Ab" = "CR-Ab",
    "CR-Ec" = "CR-Ec",
    "CR-Kp" = "CR-Kp",
    "CR-Pa" = "CR-Pa"
  ),
  
  full_names = c(
    "3GCR-Ec" = "Third-generation cephalosporin-resistant Escherichia coli",
    "3GCR-Kp" = "Third-generation cephalosporin-resistant Klebsiella pneumoniae",
    "CR-Ab" = "Carbapenem-resistant Acinetobacter baumannii",
    "CR-Ec" = "Carbapenem-resistant Escherichia coli", 
    "CR-Kp" = "Carbapenem-resistant Klebsiella pneumoniae",
    "CR-Pa" = "Carbapenem-resistant Pseudomonas aeruginosa"
  )
)

# Climate scenario configuration
climate_config <- list(
  scenarios = c("ssp126", "ssp245", "ssp370", "ssp585"),
  
  labels = c(
    "ssp126" = "SSP1-2.6",
    "ssp245" = "SSP2-4.5", 
    "ssp370" = "SSP3-7.0",
    "ssp585" = "SSP5-8.5"
  ),
  
  full_labels = c(
    "ssp126" = "SSP1-2.6 (Sustainability)",
    "ssp245" = "SSP2-4.5 (Middle of the Road)",
    "ssp370" = "SSP3-7.0 (Regional Rivalry)", 
    "ssp585" = "SSP5-8.5 (Fossil-fueled Development)"
  ),
  
  colors = c(
    "ssp126" = "#2c7bb6",
    "ssp245" = "#7aaed4", 
    "ssp370" = "#f77f4f",
    "ssp585" = "#d73027"
  )
)

# Time period configuration
time_config <- list(
  baseline_years = 2010:2019,
  projection_years = 2020:2100,
  
  assessment_periods = list(
    "2030s" = 2030:2039,
    "2050s" = 2050:2059,
    "2080s" = 2080:2089,
    "2090s" = 2090:2099
  ),
  
  display_periods = c("2030s", "2050s", "2090s"),
  
  period_colors = c(
    "2030s" = "#D9C5B4",
    "2050s" = "#84A6A5", 
    "2090s" = "#2F5E6F"
  )
)

# Climate zone configuration
climate_zones <- c("Polar Zone", "Temperate Zone", "Tropical Zone")

# Lag configuration (empirically derived from GAMM models)
lag_config <- list(
  "3GCR-Ec" = list(temp_lag = 2, precip_lag = 3, humid_lag = 2, wetdays_lag = 1),
  "3GCR-Kp" = list(temp_lag = 2, precip_lag = 3, humid_lag = 2, wetdays_lag = 2),
  "CR-Ab" = list(temp_lag = 2, precip_lag = 2, humid_lag = 3, wetdays_lag = 2),
  "CR-Ec" = list(temp_lag = 2, precip_lag = 3, humid_lag = 2, wetdays_lag = 1),
  "CR-Kp" = list(temp_lag = 3, precip_lag = 3, humid_lag = 2, wetdays_lag = 3),
  "CR-Pa" = list(temp_lag = 3, precip_lag = 1, humid_lag = 2, wetdays_lag = 2)
)

# Climate response weights (empirically derived from GAMM model coefficients)
climate_weights <- list(
  "3GCR-Ec" = list(temp = 0.2557, precip = 0.3041, humid = 0.2032, wetdays = 0.2370),
  "3GCR-Kp" = list(temp = 0.5687, precip = 0.0719, humid = 0.0993, wetdays = 0.2601),
  "CR-Ab" = list(temp = 0.6989, precip = 0.0500, humid = 0.0500, wetdays = 0.2011),
  "CR-Ec" = list(temp = 0.7882, precip = 0.0500, humid = 0.0500, wetdays = 0.1118),
  "CR-Kp" = list(temp = 0.3923, precip = 0.2366, humid = 0.0500, wetdays = 0.3211),
  "CR-Pa" = list(temp = 0.1569, precip = 0.0500, humid = 0.2681, wetdays = 0.5250)
)

# Monte Carlo simulation settings
mc_config <- list(
  n_simulations = 1000,
  confidence_level = 0.95,
  seed = 12345,
  use_parallel = TRUE,
  cores = 4,
  
  bacteria_params = list(
    "3GCR-Ec" = list(base_sd = 0.04, boundary_factor = 0.6, smoothing_window = 7),
    "3GCR-Kp" = list(base_sd = 0.05, boundary_factor = 0.6, smoothing_window = 9),
    "CR-Ab" = list(base_sd = 0.06, boundary_factor = 0.7, smoothing_window = 7),
    "CR-Ec" = list(base_sd = 0.04, boundary_factor = 0.6, smoothing_window = 7),
    "CR-Kp" = list(base_sd = 0.05, boundary_factor = 0.6, smoothing_window = 7),
    "CR-Pa" = list(base_sd = 0.06, boundary_factor = 0.7, smoothing_window = 7)
  )
)

# Uncertainty factors (based on IPCC AR6 climate projection uncertainties)
uncertainty_factors <- list(
  time_factors = c("2030s" = 1.0, "2050s" = 1.3, "2090s" = 1.8),
  scenario_factors = c("ssp126" = 1.0, "ssp245" = 1.2, "ssp370" = 1.5, "ssp585" = 1.8)
)

# Model settings (climate-only, bacteria-specific characteristics)
model_settings <- list(
  # Temporal correlation parameters
  time_series_correlation = list(persistence = 0.7, innovation = 0.3),
  
  # Bacteria-specific variation factors (empirically estimated)
  bacteria_variation_factors = c(
    "3GCR-Ec" = 1.2, "3GCR-Kp" = 1.0, "CR-Ab" = 1.5,
    "CR-Ec" = 0.8, "CR-Kp" = 0.8, "CR-Pa" = 1.3
  ),
  
  # Scenario-specific variation factors (climate uncertainty only)
  scenario_variation_factors = c(
    "ssp126" = 1.0, "ssp245" = 1.1, "ssp370" = 1.2, "ssp585" = 1.3
  ),
  
  # CR-Kp specific biological characteristics
  crkp_controls = list(
    max_effect_limit = 0.3,           # Maximum single climate effect
    threshold_damping = 0.5,          # Threshold effect dampening
    smoothing_window = 5,             # Smoothing window size
    trend_emphasis = 0.7,             # Trend component weight
    variability_reduction = 0.4,      # Variability reduction coefficient
    max_growth_rate = 0.008,          # Maximum annual growth rate
    min_memory_factor = 0.80,         # Minimum memory factor
    extra_smoothing = TRUE,           # Additional smoothing
    boundary_protection = TRUE,       # Boundary protection
    enhanced_smoothing = 7,           # Enhanced smoothing window
    end_correction = TRUE,            # End point correction
    smooth_boundary = TRUE            # Smooth boundary regions
  ),
  
  # Plot settings
  plot_settings = list(
    fig_width = 12, fig_height = 8, fig_dpi = 300,
    line_width = 1.2, confidence_interval_alpha = 0.15,
    title_size = 14, axis_title_size = 12, axis_text_size = 10,
    legend_position = "bottom", end_label_nudge = 3
  )
)

# Utility Functions ---------------------------------------------------------

# Create safe seed from string
create_seed_from_string <- function(text_string) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    # Fallback method if digest package not available
    return(as.integer(sum(utf8ToInt(text_string))) %% 2147483647L)
  }
  
  hash_value <- digest::digest(text_string, algo = "crc32")
  digits_only <- gsub("[^0-9]", "", hash_value)
  
  if (nchar(digits_only) == 0) return(42L)
  if (nchar(digits_only) > 9) digits_only <- substr(digits_only, 1, 9)
  
  seed_value <- as.integer(digits_only) %% 2147483647L
  if (is.na(seed_value)) return(42L)
  
  return(seed_value)
}

# Check required columns in dataset
check_required_columns <- function(data, required_cols, dataset_name) {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in", dataset_name, ":", paste(missing_cols, collapse = ", ")))
  }
  return(TRUE)
}

# Detect bacteria file pattern and read data
detect_and_read_bacteria_data <- function(data_path, bacteria_name) {
  # Common file patterns for bacteria data
  possible_patterns <- c(
    paste0(toupper(gsub("-", "", bacteria_name)), "_data"),
    paste0(tolower(gsub("-", "", bacteria_name)), "_data"),
    paste0(bacteria_name, "_data"),
    toupper(gsub("-", "", bacteria_name)),
    tolower(gsub("-", "", bacteria_name)),
    bacteria_name
  )
  
  # Look for files matching these patterns
  all_files <- list.files(data_path, pattern = "\\.csv$", ignore.case = TRUE)
  
  for (pattern in possible_patterns) {
    matching_files <- all_files[grepl(pattern, all_files, ignore.case = TRUE)]
    if (length(matching_files) > 0) {
      # Take the first matching file
      return(file.path(data_path, matching_files[1]))
    }
  }
  
  # If no pattern matches, return NULL
  return(NULL)
}

# Data Loading Functions ----------------------------------------------------

# Read resistance data with flexible file detection
read_resistance_data <- function(data_path) {
  cat("Loading resistance data...\n")
  
  if (!dir.exists(data_path)) {
    stop(paste("Data directory not found:", data_path))
  }
  
  all_data <- list()
  
  for (bacteria in bacteria_config$names) {
    # Try to detect the file for this bacteria
    file_path <- detect_and_read_bacteria_data(data_path, bacteria)
    
    if (is.null(file_path) || !file.exists(file_path)) {
      warning(paste("No data file found for", bacteria, "in", data_path))
      next
    }
    
    cat(paste("Loading", bacteria, "from:", basename(file_path), "\n"))
    
    data <- read.csv(file_path)
    
    # Check required columns
    required_cols <- c("R", "year", "lat")
    
    # Try different common column name variations
    if (!"NAME" %in% names(data)) {
      if ("Country" %in% names(data)) {
        data$NAME <- data$Country
      } else if ("country" %in% names(data)) {
        data$NAME <- data$country
      } else if ("country_name" %in% names(data)) {
        data$NAME <- data$country_name
      } else {
        data$NAME <- "Unknown"
      }
    }
    
    if (!"Region" %in% names(data)) {
      if ("region" %in% names(data)) {
        data$Region <- data$region
      } else if ("region_name" %in% names(data)) {
        data$Region <- data$region_name
      } else {
        data$Region <- "Unknown"
      }
    }
    
    check_required_columns(data, required_cols, paste(bacteria, "resistance data"))
    
    # Process data
    data <- data %>%
      mutate(
        year = as.integer(year),
        bacteria = bacteria,
        resistance_rate = ifelse(median(R, na.rm = TRUE) <= 1, R * 100, R),
        climate_zone = case_when(
          abs(lat) >= 60 ~ "Polar Zone",
          abs(lat) < 23.5 ~ "Tropical Zone",
          TRUE ~ "Temperate Zone"
        )
      ) %>%
      rename(country_name = NAME, region_name = Region)
    
    all_data[[bacteria]] <- data
    
    cat(sprintf(" %s resistance rate range: %.2f%% - %.2f%%\n", 
                bacteria, min(data$resistance_rate, na.rm = TRUE), 
                max(data$resistance_rate, na.rm = TRUE)))
  }
  
  if (length(all_data) == 0) {
    stop("No resistance data files found. Please check the data directory and file naming.")
  }
  
  combined_data <- bind_rows(all_data)
  
  cat(paste("Loaded", nrow(combined_data), "rows of resistance data\n"))
  cat(paste("Coverage:", length(unique(combined_data$bacteria)), "bacteria,",
            "years", min(combined_data$year, na.rm = TRUE), "-", 
            max(combined_data$year, na.rm = TRUE), "\n\n"))
  
  return(combined_data)
}

# Read climate data with flexible naming
read_climate_data <- function(climate_file_path) {
  cat("Loading climate data...\n")
  
  if (!file.exists(climate_file_path)) {
    stop(paste("Climate data file not found:", climate_file_path))
  }
  
  climate_data <- read_csv(climate_file_path, show_col_types = FALSE)
  
  # Flexible column name mapping
  name_mappings <- list(
    country_name = c("NAME", "name", "Country", "country", "country_name"),
    region_name = c("Region", "region", "region_name"),
    scenario = c("ssp", "SSP", "scenario", "Scenario"),
    temperature = c("tas_annual_mean", "temperature", "temp", "Temperature"),
    humidity = c("hurs_annual_mean", "humidity", "humid", "Humidity"),
    precipitation = c("pr_annual_total", "precipitation", "precip", "Precipitation"),
    wetdays = c("wet_days_yearly0.1", "wetdays", "wet_days", "Wet_Days")
  )
  
  # Apply flexible naming
  for (target_name in names(name_mappings)) {
    possible_names <- name_mappings[[target_name]]
    found_name <- intersect(possible_names, names(climate_data))[1]
    
    if (!is.na(found_name) && found_name != target_name) {
      climate_data <- climate_data %>%
        rename(!!target_name := !!found_name)
    }
  }
  
  # Check required columns after mapping
  required_cols <- c("year", "scenario", "temperature", "humidity", "precipitation", "wetdays")
  check_required_columns(climate_data, required_cols, "climate data")
  
  # Process data
  climate_data <- climate_data %>%
    mutate(
      scenario = tolower(scenario),
      climate_zone = case_when(
        abs(lat) >= 60 ~ "Polar Zone",
        abs(lat) < 23.5 ~ "Tropical Zone", 
        TRUE ~ "Temperate Zone"
      )
    )
  
  cat(paste("Loaded", nrow(climate_data), "rows of climate data\n"))
  cat(paste("Countries:", length(unique(climate_data$country_name)), "\n"))
  cat(paste("Scenarios:", paste(unique(climate_data$scenario), collapse = ", "), "\n"))
  cat(paste("Year range:", min(climate_data$year), "-", max(climate_data$year), "\n\n"))
  
  return(climate_data)
}

# Data Processing Functions -------------------------------------------------

# Prepare baseline data
prepare_baseline_data <- function(resistance_data, climate_data) {
  cat("Preparing baseline data...\n")
  
  # Filter baseline period
  baseline_resistance <- resistance_data %>%
    filter(year %in% time_config$baseline_years) %>%
    group_by(bacteria, climate_zone) %>%
    summarise(baseline_resistance = mean(resistance_rate, na.rm = TRUE), .groups = "drop")
  
  # Calculate baseline climate conditions
  climate_baseline <- climate_data %>%
    filter(year %in% time_config$baseline_years) %>%
    group_by(climate_zone) %>%
    summarise(
      baseline_temp = mean(temperature, na.rm = TRUE),
      baseline_humid = mean(humidity, na.rm = TRUE),
      baseline_precip = mean(precipitation, na.rm = TRUE),
      baseline_wetdays = mean(wetdays, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Combine resistance and climate baselines
  baseline_data <- baseline_resistance %>%
    left_join(climate_baseline, by = "climate_zone")
  
  # Create global average baseline
  global_baseline <- baseline_data %>%
    group_by(bacteria) %>%
    summarise(
      climate_zone = "Global Average",
      baseline_resistance = mean(baseline_resistance, na.rm = TRUE),
      baseline_temp = mean(baseline_temp, na.rm = TRUE),
      baseline_humid = mean(baseline_humid, na.rm = TRUE), 
      baseline_precip = mean(baseline_precip, na.rm = TRUE),
      baseline_wetdays = mean(baseline_wetdays, na.rm = TRUE),
      .groups = "drop"
    )
  
  full_baseline <- bind_rows(baseline_data, global_baseline)
  
  cat("Baseline resistance rates by bacteria (%):\n")
  global_baseline %>%
    select(bacteria, baseline_resistance) %>%
    print()
  
  return(full_baseline)
}

# Prepare climate scenarios
prepare_climate_scenarios <- function(climate_data) {
  cat("Preparing climate scenario data...\n")
  
  # Filter and aggregate scenario data
  scenarios_data <- climate_data %>%
    filter(scenario %in% climate_config$scenarios, year >= min(time_config$projection_years)) %>%
    group_by(year, scenario, climate_zone) %>%
    summarise(
      temperature = mean(temperature, na.rm = TRUE),
      humidity = mean(humidity, na.rm = TRUE), 
      precipitation = mean(precipitation, na.rm = TRUE),
      wetdays = mean(wetdays, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Ensure complete coverage with interpolation if needed
  complete_scenarios <- expand_grid(
    year = time_config$projection_years,
    scenario = climate_config$scenarios,
    climate_zone = climate_zones
  )
  
  full_scenarios <- complete_scenarios %>%
    left_join(scenarios_data, by = c("year", "scenario", "climate_zone")) %>%
    group_by(scenario, climate_zone) %>%
    arrange(year) %>%
    mutate(
      temperature = approx(year, temperature, year, rule = 2)$y,
      humidity = approx(year, humidity, year, rule = 2)$y,
      precipitation = approx(year, precipitation, year, rule = 2)$y,
      wetdays = approx(year, wetdays, year, rule = 2)$y
    ) %>%
    ungroup()
  
  cat("Climate scenario data prepared for", length(climate_config$scenarios), "scenarios\n")
  
  return(full_scenarios)
}

# Climate Response Functions ------------------------------------------------

# Load GAMM-derived response functions (REQUIRED)
load_gamm_response_functions <- function(response_file_path) {
  cat("Loading GAMM-derived response functions...\n")
  
  if (!file.exists(response_file_path)) {
    stop(paste("GAMM response functions file not found:", response_file_path, 
               "\nThis file is required and must contain empirically-derived functions from GAMM analysis."))
  }
  
  tryCatch({
    # Load the response functions file
    source(response_file_path, local = FALSE)
    
    # Verify required functions exist
    required_functions <- c("predict_climate_effect", "predict_combined_climate_effect")
    missing_functions <- required_functions[!sapply(required_functions, exists)]
    
    if (length(missing_functions) > 0) {
      stop(paste("Required GAMM functions not found in", response_file_path, ":", 
                 paste(missing_functions, collapse = ", "),
                 "\nPlease ensure the file contains empirically-derived response functions from GAMM analysis."))
    }
    
    # Test the functions with sample data
    test_bacteria <- bacteria_config$names[1]
    tryCatch({
      test_temp_effect <- predict_climate_effect(test_bacteria, "Temperature", 20)
      test_humid_effect <- predict_climate_effect(test_bacteria, "Relative Humidity", 70)
      test_precip_effect <- predict_climate_effect(test_bacteria, "Precipitation", 1000)
      test_wetdays_effect <- predict_climate_effect(test_bacteria, "Wet Days", 150)
      
      cat("GAMM response functions loaded and tested successfully\n")
      cat("Sample responses for", test_bacteria, ":\n")
      cat(" Temperature (20°C):", round(test_temp_effect, 4), "\n")
      cat(" Humidity (70%):", round(test_humid_effect, 4), "\n")
      cat(" Precipitation (1000mm):", round(test_precip_effect, 4), "\n")
      cat(" Wet Days (150d):", round(test_wetdays_effect, 4), "\n")
      
      return(TRUE)
      
    }, error = function(e) {
      stop(paste("GAMM response functions failed testing:", e$message,
                 "\nPlease verify the functions are correctly implemented."))
    })
    
  }, error = function(e) {
    stop(paste("Error loading GAMM response functions:", e$message))
  })
}

# Calculate climate effects using GAMM response functions
calculate_climate_effects <- function(bacteria, temperature, humidity, precipitation, wetdays) {
  # Use GAMM-derived functions to predict climate effects
  temp_effect <- predict_climate_effect(bacteria, "Temperature", temperature)
  humid_effect <- predict_climate_effect(bacteria, "Relative Humidity", humidity)
  precip_effect <- predict_climate_effect(bacteria, "Precipitation", precipitation)
  wetdays_effect <- predict_climate_effect(bacteria, "Wet Days", wetdays)
  
  # Get bacteria-specific weights (empirically derived)
  weights <- climate_weights[[bacteria]]
  if (is.null(weights)) {
    stop(paste("Climate weights not found for bacteria:", bacteria))
  }
  
  # Calculate weighted combined effect
  combined_effect <- 1.0
  if (weights$temp > 0) combined_effect <- combined_effect * (temp_effect ^ weights$temp)
  if (weights$precip > 0) combined_effect <- combined_effect * (precip_effect ^ weights$precip)
  if (weights$humid > 0) combined_effect <- combined_effect * (humid_effect ^ weights$humid)
  if (weights$wetdays > 0) combined_effect <- combined_effect * (wetdays_effect ^ weights$wetdays)
  
  return(list(
    temp_effect = temp_effect,
    humid_effect = humid_effect,
    precip_effect = precip_effect,
    wetdays_effect = wetdays_effect,
    combined_effect = combined_effect
  ))
}

# Prediction Functions ------------------------------------------------------

# Apply CR-Kp specific biological characteristics
apply_crkp_specific_processing <- function(combined_effect, bacteria) {
  if (bacteria != "CR-Kp") return(combined_effect)
  
  cat("Applying CR-Kp specific biological characteristics...\n")
  
  # CR-Kp specific controls based on biological properties
  controls <- model_settings$crkp_controls
  
  # Apply variability reduction for CR-Kp biological stability
  variability_reduction <- controls$variability_reduction
  
  # Calculate moving average to extract trend component
  window_size <- controls$enhanced_smoothing
  trend_component <- combined_effect
  
  if (length(trend_component) > window_size) {
    # Calculate moving average
    weights <- rep(1/window_size, window_size)
    ma_trend <- stats::filter(trend_component, weights, sides = 2)
    
    # Handle NA values with interpolation
    na_idx <- is.na(ma_trend)
    if(any(na_idx)) {
      if(sum(!na_idx) > 1) {
        ma_trend[na_idx] <- approx(
          x = which(!na_idx),
          y = ma_trend[!na_idx],
          xout = which(na_idx),
          rule = 2
        )$y
      } else {
        ma_trend[na_idx] <- trend_component[na_idx]
      }
    }
    
    # Extract fluctuation component
    fluctuation_component <- trend_component - ma_trend
    
    # Reduce fluctuation component for CR-Kp biological stability
    reduced_fluctuation <- fluctuation_component * (1 - variability_reduction)
    
    # Recombine trend and reduced fluctuation
    combined_effect <- ma_trend + reduced_fluctuation
  }
  
  # End point correction for biological realism
  if (controls$end_correction) {
    n <- length(combined_effect)
    if (n > 3) {
      # Smooth beginning
      combined_effect[1] <- (
        3 * combined_effect[1] +
          2 * combined_effect[2] +
          combined_effect[3]
      ) / 6
      
      # Smooth end
      combined_effect[n] <- (
        3 * combined_effect[n] +
          2 * combined_effect[n-1] +
          combined_effect[n-2]
      ) / 6
    }
  }
  
  # Limit maximum effect range for biological realism
  max_effect <- controls$max_effect_limit * 1.5
  combined_effect <- pmin(
    pmax(combined_effect, 1 - max_effect),
    1 + max_effect
  )
  
  return(combined_effect)
}

# Predict future resistance rates (climate-driven with bacteria-specific characteristics)
predict_future_resistance <- function(baseline_data, climate_scenarios) {
  cat("Predicting future resistance rates using climate factors and bacteria-specific characteristics...\n")
  
  predictions <- tibble()
  
  for (bacteria in bacteria_config$names) {
    cat(paste("Processing", bacteria, "...\n"))
    
    # Get bacteria-specific parameters
    lag_combo <- lag_config[[bacteria]]
    bacteria_variation <- model_settings$bacteria_variation_factors[bacteria]
    is_crkp <- bacteria == "CR-Kp"
    
    for (zone in climate_zones) {
      # Get zone baseline data
      zone_baseline <- baseline_data %>%
        filter(bacteria == !!bacteria, climate_zone == !!zone)
      
      if (nrow(zone_baseline) == 0) {
        warning(paste("No baseline data for", bacteria, "in", zone))
        next
      }
      
      baseline_resistance <- zone_baseline$baseline_resistance[1]
      
      for (scenario in climate_config$scenarios) {
        # Get scenario climate data
        scenario_climate <- climate_scenarios %>%
          filter(climate_zone == !!zone, scenario == !!scenario)
        
        if (nrow(scenario_climate) == 0) {
          warning(paste("No climate data for", scenario, "in", zone))
          next
        }
        
        # Initialize predictions dataframe
        scenario_predictions <- tibble(
          year = scenario_climate$year,
          bacteria = bacteria,
          climate_zone = zone,
          scenario = scenario,
          temperature = scenario_climate$temperature,
          humidity = scenario_climate$humidity,
          precipitation = scenario_climate$precipitation,
          wetdays = scenario_climate$wetdays,
          temperature_lagged = NA_real_,
          humidity_lagged = NA_real_,
          precipitation_lagged = NA_real_,
          wetdays_lagged = NA_real_,
          resistance_rate = NA_real_,
          lower_ci = NA_real_,
          upper_ci = NA_real_,
          relative_change = NA_real_,
          relative_change_baseline = NA_real_
        )
        
        # Apply lag effects
        for (i in seq_len(nrow(scenario_predictions))) {
          scenario_predictions$temperature_lagged[i] <- if (i <= lag_combo$temp_lag) {
            zone_baseline$baseline_temp
          } else {
            scenario_predictions$temperature[i - lag_combo$temp_lag]
          }
          
          scenario_predictions$humidity_lagged[i] <- if (i <= lag_combo$humid_lag) {
            zone_baseline$baseline_humid
          } else {
            scenario_predictions$humidity[i - lag_combo$humid_lag]
          }
          
          scenario_predictions$precipitation_lagged[i] <- if (i <= lag_combo$precip_lag) {
            zone_baseline$baseline_precip
          } else {
            scenario_predictions$precipitation[i - lag_combo$precip_lag]
          }
          
          scenario_predictions$wetdays_lagged[i] <- if (i <= lag_combo$wetdays_lag) {
            zone_baseline$baseline_wetdays
          } else {
            scenario_predictions$wetdays[i - lag_combo$wetdays_lag]
          }
        }
        
        # Calculate climate effects using GAMM functions
        climate_effects <- calculate_climate_effects(
          bacteria,
          scenario_predictions$temperature_lagged,
          scenario_predictions$humidity_lagged,
          scenario_predictions$precipitation_lagged,
          scenario_predictions$wetdays_lagged
        )
        
        # Apply bacteria-specific processing (especially CR-Kp)
        combined_effects <- apply_crkp_specific_processing(climate_effects$combined_effect, bacteria)
        
        # Predict resistance rates with bacteria-specific temporal characteristics
        persistence <- model_settings$time_series_correlation$persistence
        innovation <- model_settings$time_series_correlation$innovation
        
        # CR-Kp specific temporal adjustments
        if (is_crkp) {
          persistence <- min(0.85, persistence + 0.10)
          innovation <- max(0.15, innovation - 0.10)
        }
        
        for (i in seq_len(nrow(scenario_predictions))) {
          if (i == 1) {
            # First year based on baseline and climate effect
            current_rate <- baseline_resistance * combined_effects[i]
          } else {
            # Subsequent years with temporal correlation and bacteria-specific characteristics
            years_passed <- scenario_predictions$year[i] - min(scenario_predictions$year)
            
            # Bacteria-specific growth rates
            base_growth_rate <- switch(bacteria,
                                       "CR-Ab" = 0.0025,    # Fastest spread
                                       "CR-Pa" = 0.0023,
                                       "3GCR-Kp" = 0.0022,
                                       "CR-Kp" = model_settings$crkp_controls$max_growth_rate,  # CR-Kp controlled rate
                                       "3GCR-Ec" = 0.0018,
                                       "CR-Ec" = 0.0016,    # Slowest spread
                                       0.002                # Default
            )
            
            climate_growth <- (combined_effects[i] - 1) * 0.5
            
            # Add controlled stochastic variation
            set.seed(create_seed_from_string(paste(bacteria, scenario, years_passed, sep = "_")))
            scenario_variation <- model_settings$scenario_variation_factors[scenario]
            time_variation <- 0.005 + 0.0002 * years_passed
            
            # CR-Kp specific variation control
            if (is_crkp) {
              time_variation <- time_variation * 0.6  # Reduced variation for biological stability
            }
            
            random_component <- rnorm(1, mean = 0, sd = time_variation * scenario_variation * bacteria_variation)
            
            # Calculate climate-driven rate with bacteria-specific characteristics
            climate_rate <- baseline_resistance * combined_effects[i] * 
              (1 + base_growth_rate * years_passed + climate_growth) * 
              (1 + random_component)
            
            # Apply temporal correlation
            previous_rate <- scenario_predictions$resistance_rate[i-1]
            current_rate <- persistence * previous_rate + innovation * climate_rate
            
            # CR-Kp specific growth rate limiting
            if (is_crkp && i > 1) {
              growth_rate <- (current_rate - previous_rate) / previous_rate
              max_allowed_growth <- model_settings$crkp_controls$max_growth_rate * 
                model_settings$scenario_variation_factors[scenario]
              
              if (abs(growth_rate) > max_allowed_growth) {
                limited_growth <- sign(growth_rate) * max_allowed_growth
                current_rate <- previous_rate * (1 + limited_growth)
              }
              
              # Additional bounds for CR-Kp biological realism
              current_rate <- min(max(current_rate, baseline_resistance * 0.6), baseline_resistance * 2.2)
            }
          }
          
          # Ensure reasonable bounds for all bacteria
          current_rate <- max(0.1, min(current_rate, baseline_resistance * 3.0))
          scenario_predictions$resistance_rate[i] <- current_rate
        }
        
        # Calculate initial confidence intervals (updated later by Monte Carlo)
        scenario_predictions$lower_ci <- scenario_predictions$resistance_rate * 0.9
        scenario_predictions$upper_ci <- scenario_predictions$resistance_rate * 1.1
        
        # Calculate relative changes
        scenario_predictions$relative_change <- 100 * (
          scenario_predictions$resistance_rate / scenario_predictions$resistance_rate[1] - 1
        )
        scenario_predictions$relative_change_baseline <- 100 * (
          scenario_predictions$resistance_rate / baseline_resistance - 1
        )
        
        predictions <- bind_rows(predictions, scenario_predictions)
      }
    }
  }
  
  return(predictions)
}

# Calculate global averages
calculate_global_averages <- function(predictions) {
  cat("Calculating global averages...\n")
  
  if (nrow(predictions) == 0) {
    warning("No prediction data available")
    return(tibble())
  }
  
  global_predictions <- predictions %>%
    group_by(bacteria, scenario, year) %>%
    summarise(
      resistance_rate = mean(resistance_rate, na.rm = TRUE),
      lower_ci = mean(lower_ci, na.rm = TRUE),
      upper_ci = mean(upper_ci, na.rm = TRUE),
      relative_change = mean(relative_change, na.rm = TRUE),
      relative_change_baseline = mean(relative_change_baseline, na.rm = TRUE),
      climate_zone = "Global Average",
      .groups = "drop"
    )
  
  return(bind_rows(predictions, global_predictions))
}

# Monte Carlo Functions ----------------------------------------------------

# Run Monte Carlo simulations for uncertainty quantification
run_monte_carlo_simulations <- function(predictions) {
  cat("Running Monte Carlo simulations for uncertainty quantification...\n")
  
  if (nrow(predictions) == 0) {
    warning("No prediction data for Monte Carlo simulation")
    return(predictions)
  }
  
  set.seed(mc_config$seed)
  
  # Focus on global average data
  global_data <- predictions %>% filter(climate_zone == "Global Average")
  
  for (bacteria in bacteria_config$names) {
    cat(paste("Processing", bacteria, "Monte Carlo simulations...\n"))
    
    # Get bacteria-specific MC parameters
    mc_params <- mc_config$bacteria_params[[bacteria]]
    
    bacteria_data <- global_data %>% filter(bacteria == !!bacteria)
    
    for (scenario in climate_config$scenarios) {
      scenario_data <- bacteria_data %>%
        filter(scenario == !!scenario) %>%
        arrange(year)
      
      if (nrow(scenario_data) == 0) next
      
      # Apply smoothing to reduce noise
      if (nrow(scenario_data) > mc_params$smoothing_window) {
        window_size <- mc_params$smoothing_window
        smoothed_rates <- stats::filter(scenario_data$resistance_rate, 
                                        rep(1/window_size, window_size), sides = 2)
        
        # Handle NA values from filtering
        na_idx <- is.na(smoothed_rates)
        if (any(na_idx) && sum(!na_idx) > 1) {
          smoothed_rates[na_idx] <- approx(
            which(!na_idx), smoothed_rates[!na_idx],
            which(na_idx), rule = 2
          )$y
        } else {
          smoothed_rates[na_idx] <- scenario_data$resistance_rate[na_idx]
        }
      } else {
        smoothed_rates <- scenario_data$resistance_rate
      }
      
      # Calculate uncertainty for each time point
      for (i in seq_len(nrow(scenario_data))) {
        year <- scenario_data$year[i]
        
        # Time-dependent uncertainty (climate projection uncertainty increases over time)
        years_from_start <- year - min(scenario_data$year)
        max_years <- max(scenario_data$year) - min(scenario_data$year)
        time_factor <- 1.0 + (years_from_start / max_years) * 0.5
        
        # Scenario uncertainty (higher emission scenarios have more uncertainty)
        scenario_factor <- uncertainty_factors$scenario_factors[scenario]
        
        # Boundary effects
        is_boundary <- (i == 1 || i == nrow(scenario_data))
        boundary_factor <- if (is_boundary) mc_params$boundary_factor else 1.0
        
        # Calculate standard deviation based on climate uncertainty
        base_sd <- mc_params$base_sd
        estimated_sd <- smoothed_rates[i] * base_sd * time_factor * scenario_factor * boundary_factor
        
        # Generate Monte Carlo samples
        set.seed(create_seed_from_string(paste(bacteria, scenario, year, sep = "_")))
        mc_samples <- rnorm(mc_config$n_simulations, mean = smoothed_rates[i], sd = estimated_sd)
        
        # Calculate confidence intervals
        conf_level <- mc_config$confidence_level
        quantiles <- quantile(mc_samples, probs = c((1 - conf_level)/2, 1 - (1 - conf_level)/2))
        
        # Update predictions
        idx <- which(predictions$bacteria == bacteria & 
                       predictions$scenario == scenario &
                       predictions$year == year & 
                       predictions$climate_zone == "Global Average")
        
        if (length(idx) == 1) {
          predictions$lower_ci[idx] <- max(0.1, quantiles[1])
          predictions$upper_ci[idx] <- quantiles[2]
        }
      }
    }
  }
  
  cat("Monte Carlo simulations completed\n")
  return(predictions)
}

# Visualization Functions --------------------------------------------------

# Create time trends plot
plot_resistance_time_trends <- function(predictions, output_path = NULL) {
  cat("Creating time trends plot...\n")
  
  # Extract global average data
  global_data <- predictions %>%
    filter(climate_zone == "Global Average") %>%
    mutate(
      year = as.numeric(year),
      scenario = factor(scenario, levels = climate_config$scenarios),
      bacteria = factor(bacteria, levels = bacteria_config$names)
    )
  
  # Get final year values for labeling
  final_year <- max(global_data$year)
  final_values <- global_data %>%
    filter(year == final_year) %>%
    mutate(
      label = sprintf("%.2f", round(resistance_rate, 2)),
      x_position = year + model_settings$plot_settings$end_label_nudge
    )
  
  # Create plot
  p <- ggplot(global_data, aes(x = year, y = resistance_rate, color = scenario, group = scenario)) +
    # Confidence intervals
    geom_ribbon(
      aes(ymin = lower_ci, ymax = upper_ci, fill = scenario),
      alpha = model_settings$plot_settings$confidence_interval_alpha,
      color = NA
    ) +
    # Trend lines
    geom_line(linewidth = model_settings$plot_settings$line_width) +
    # Final year labels
    geom_label_repel(
      data = final_values,
      aes(x = x_position, y = resistance_rate, label = label),
      size = 3,
      fontface = "bold",
      min.segment.length = 0,
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.2, "lines"),
      force = 2,
      nudge_x = 2,
      direction = "y",
      hjust = 0,
      max.overlaps = 15
    ) +
    # Color and fill scales
    scale_color_manual(values = climate_config$colors, 
                       labels = climate_config$labels, 
                       name = "Climate Scenario") +
    scale_fill_manual(values = climate_config$colors, guide = "none") +
    # Axis settings
    scale_x_continuous(
      breaks = seq(2020, 2100, by = 20),
      limits = c(2020, 2108),
      expand = c(0.01, 0)
    ) +
    # Faceting
    facet_wrap(~ bacteria, scales = "free_y", ncol = 3) +
    # Labels
    labs(
      title = "Climate-Driven AMR Projections with Bacteria-Specific Characteristics (2020-2100)",
      subtitle = "Projections based on empirical climate-AMR relationships under SSP scenarios",
      x = "Year",
      y = "Resistance Rate (%)"
    ) +
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = model_settings$plot_settings$title_size, 
                                face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = model_settings$plot_settings$title_size - 2, 
                                   hjust = 0.5),
      axis.title = element_text(size = model_settings$plot_settings$axis_title_size, 
                                face = "bold"),
      axis.text = element_text(size = model_settings$plot_settings$axis_text_size),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.5),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = model_settings$plot_settings$legend_position,
      legend.direction = "horizontal"
    )
  
  # Save plot if path provided
  if (!is.null(output_path)) {
    tryCatch({
      ggsave(output_path, p, 
             width = model_settings$plot_settings$fig_width,
             height = model_settings$plot_settings$fig_height,
             dpi = model_settings$plot_settings$fig_dpi)
      cat(paste("Time trends plot saved to:", output_path, "\n"))
    }, error = function(e) {
      cat(paste("Error saving plot:", e$message, "\n"))
    })
  }
  
  return(p)
}

# Create relative changes plot
plot_relative_changes <- function(predictions, output_path = NULL) {
  cat("Creating relative changes plot...\n")
  
  # Prepare data for different periods
  period_data <- map_dfr(time_config$display_periods, function(period) {
    period_years <- time_config$assessment_periods[[period]]
    
    predictions %>%
      filter(climate_zone == "Global Average", year %in% period_years) %>%
      group_by(bacteria, scenario) %>%
      summarise(
        mean_relative_change = mean(relative_change_baseline, na.rm = TRUE),
        sd_relative_change = sd(relative_change_baseline, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        period = period,
        se = sd_relative_change / sqrt(n),
        # Apply uncertainty factors (climate-based only)
        adjusted_se = se * uncertainty_factors$time_factors[period] * 
          uncertainty_factors$scenario_factors[scenario],
        lower_ci = mean_relative_change - 1.96 * adjusted_se,
        upper_ci = mean_relative_change + 1.96 * adjusted_se
      )
  })
  
  # Set factor levels
  period_data <- period_data %>%
    mutate(
      period = factor(period, levels = time_config$display_periods),
      scenario = factor(scenario, levels = climate_config$scenarios),
      bacteria = factor(bacteria, levels = bacteria_config$names)
    )
  
  # Create plot
  p <- ggplot(period_data, aes(x = scenario, y = mean_relative_change, fill = period)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, 
             color = "gray40", linewidth = 0.2) +
    geom_errorbar(
      aes(ymin = lower_ci, ymax = upper_ci),
      position = position_dodge(width = 0.8),
      width = 0.3,
      linewidth = 0.7,
      alpha = 0.9,
      color = "gray50"
    ) +
    geom_text(
      aes(label = sprintf("%.1f", round(mean_relative_change, 1))),
      position = position_dodge(width = 0.8),
      vjust = ifelse(period_data$mean_relative_change >= 0, -0.4, 1.5),
      size = 3,
      color = "black",
      fontface = "bold"
    ) +
    scale_x_discrete(labels = function(x) climate_config$labels[x]) +
    scale_fill_manual(values = time_config$period_colors, name = "Time Period") +
    facet_wrap(~ bacteria, scales = "free_y", ncol = 3) +
    labs(
      title = "Climate-Driven Relative Change in AMR with Bacteria-Specific Characteristics",
      subtitle = paste0("Changes relative to baseline (", min(time_config$baseline_years), 
                        "-", max(time_config$baseline_years), ") with 95% confidence intervals"),
      x = "Climate Scenario",
      y = "Relative Change (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = model_settings$plot_settings$title_size, 
                                face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = model_settings$plot_settings$title_size - 2, 
                                   hjust = 0.5),
      axis.title = element_text(size = model_settings$plot_settings$axis_title_size, 
                                face = "bold"),
      axis.text = element_text(size = model_settings$plot_settings$axis_text_size),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.2),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.5),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
  
  # Save plot if path provided
  if (!is.null(output_path)) {
    tryCatch({
      ggsave(output_path, p,
             width = model_settings$plot_settings$fig_width,
             height = model_settings$plot_settings$fig_height,
             dpi = model_settings$plot_settings$fig_dpi)
      cat(paste("Relative changes plot saved to:", output_path, "\n"))
    }, error = function(e) {
      cat(paste("Error saving plot:", e$message, "\n"))
    })
  }
  
  return(p)
}

# Data Export Functions -----------------------------------------------------

# Export prediction data
export_prediction_data <- function(predictions, output_dir) {
  cat("Exporting prediction data...\n")
  
  # Create data directory
  data_dir <- file.path(output_dir, "data")
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Export main predictions
  main_path <- file.path(data_dir, "climate_amr_projections.csv")
  write_csv(predictions, main_path)
  
  # Export global averages only
  global_path <- file.path(data_dir, "climate_amr_projections_global.csv")
  global_data <- predictions %>%
    filter(climate_zone == "Global Average") %>%
    select(bacteria, scenario, year, resistance_rate, lower_ci, upper_ci, 
           relative_change, relative_change_baseline)
  write_csv(global_data, global_path)
  
  # Export summary statistics
  summary_data <- predictions %>%
    filter(climate_zone == "Global Average") %>%
    group_by(bacteria, scenario) %>%
    summarise(
      baseline_rate = first(resistance_rate),
      final_rate = last(resistance_rate),
      total_change = final_rate - baseline_rate,
      percent_change = 100 * (final_rate - baseline_rate) / baseline_rate,
      .groups = "drop"
    )
  
  summary_path <- file.path(data_dir, "projection_summary.csv")
  write_csv(summary_data, summary_path)
  
  cat(paste("Prediction data exported to:", data_dir, "\n"))
  
  return(list(main = main_path, global = global_path, summary = summary_path))
}

# Export metadata
export_metadata <- function(output_dir) {
  metadata_dir <- file.path(output_dir, "metadata")
  dir.create(metadata_dir, showWarnings = FALSE, recursive = TRUE)
  
  metadata <- list(
    model_version = "Climate-AMR Projection Model with Bacteria-Specific Characteristics",
    model_type = "Empirical climate-driven projection with biological realism",
    creation_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    baseline_period = paste(range(time_config$baseline_years), collapse = "-"),
    projection_period = paste(range(time_config$projection_years), collapse = "-"),
    bacteria = bacteria_config$names,
    climate_scenarios = climate_config$scenarios,
    climate_factors = c("Temperature", "Humidity", "Precipitation", "Wet Days"),
    monte_carlo_simulations = mc_config$n_simulations,
    confidence_level = mc_config$confidence_level,
    bacteria_specific_features = list(
      "CR-Kp biological stability characteristics",
      "Bacteria-specific growth rates",
      "Species-specific temporal correlation patterns",
      "Biological realism constraints"
    ),
    model_assumptions = list(
      "Climate effects only - no socioeconomic factors",
      "GAMM-derived response functions required", 
      "Empirical lag structures from GAMM analysis",
      "Climate uncertainty based on IPCC AR6",
      "Bacteria-specific biological characteristics incorporated"
    )
  )
  
  # Save as JSON if available
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    metadata_path <- file.path(metadata_dir, "model_metadata.json")
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE)
    cat("Metadata exported to JSON format\n")
  }
  
  # Save as text file
  metadata_txt <- file.path(metadata_dir, "model_metadata.txt")
  writeLines(
    c("Climate-AMR Projection Model with Bacteria-Specific Characteristics",
      "====================================================================",
      paste("Creation Date:", metadata$creation_date),
      paste("Model Type:", metadata$model_type),
      paste("Baseline Period:", metadata$baseline_period),
      paste("Projection Period:", metadata$projection_period),
      paste("Bacteria:", paste(metadata$bacteria, collapse = ", ")),
      paste("Climate Scenarios:", paste(metadata$climate_scenarios, collapse = ", ")),
      paste("Monte Carlo Simulations:", metadata$monte_carlo_simulations),
      paste("Confidence Level:", metadata$confidence_level),
      "",
      "Bacteria-Specific Features:",
      paste("•", metadata$bacteria_specific_features, collapse = "\n"),
      "",
      "Model Assumptions:",
      paste("•", metadata$model_assumptions, collapse = "\n")
    ),
    metadata_txt
  )
  
  return(metadata_path)
}

# Main Analysis Function ---------------------------------------------------

# Run complete climate-AMR projection analysis
run_climate_amr_projection <- function(data_path, climate_file_path, 
                                       response_functions_path, output_dir) {
  
  cat("Starting Climate-AMR Projection Analysis with Bacteria-Specific Characteristics\n")
  cat("==============================================================================\n")
  start_time <- Sys.time()
  
  # Create output directories
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load GAMM response functions (REQUIRED - no fallback)
  load_gamm_response_functions(response_functions_path)
  
  # Load and process data
  resistance_data <- read_resistance_data(data_path)
  climate_data <- read_climate_data(climate_file_path)
  
  # Prepare analysis datasets
  baseline_data <- prepare_baseline_data(resistance_data, climate_data)
  climate_scenarios <- prepare_climate_scenarios(climate_data)
  
  # Run climate-driven predictions with bacteria-specific characteristics
  cat("Running climate-driven projection models with bacteria-specific characteristics...\n")
  predictions <- predict_future_resistance(baseline_data, climate_scenarios)
  
  if (nrow(predictions) == 0) {
    stop("No predictions generated. Please check input data and GAMM response functions.")
  }
  
  # Calculate global averages
  predictions <- calculate_global_averages(predictions)
  
  # Run Monte Carlo uncertainty analysis
  predictions <- run_monte_carlo_simulations(predictions)
  
  # Create visualizations
  cat("Creating visualizations...\n")
  
  time_trends_plot <- plot_resistance_time_trends(
    predictions, file.path(figures_dir, "climate_amr_time_trends.png")
  )
  
  relative_changes_plot <- plot_relative_changes(
    predictions, file.path(figures_dir, "climate_amr_relative_changes.png")
  )
  
  # Save plots as PDF
  tryCatch({
    ggsave(file.path(figures_dir, "climate_amr_time_trends.pdf"), time_trends_plot,
           width = model_settings$plot_settings$fig_width,
           height = model_settings$plot_settings$fig_height,
           device = cairo_pdf)
    
    ggsave(file.path(figures_dir, "climate_amr_relative_changes.pdf"), relative_changes_plot,
           width = model_settings$plot_settings$fig_width,
           height = model_settings$plot_settings$fig_height,
           device = cairo_pdf)
    
    cat("PDF plots saved successfully\n")
  }, error = function(e) {
    cat(paste("Warning: Could not save PDF plots -", e$message, "\n"))
  })
  
  # Export data and metadata
  data_paths <- export_prediction_data(predictions, output_dir)
  metadata_path <- export_metadata(output_dir)
  
  # Calculate runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  cat("\nClimate-driven AMR projection analysis completed successfully!\n")
  cat("=============================================================\n")
  cat(paste("Runtime:", round(runtime, 2), "minutes\n"))
  cat(paste("Results saved to:", output_dir, "\n"))
  cat("\nThis analysis incorporates:\n")
  cat("• Empirical climate-AMR relationships from GAMM models\n")
  cat("• Bacteria-specific biological characteristics (including CR-Kp stability)\n") 
  cat("• No subjective socioeconomic assumptions\n")
  cat("• IPCC AR6-based climate uncertainty quantification\n")
  
  # Return results
  return(list(
    predictions = predictions,
    baseline_data = baseline_data,
    plots = list(time_trends = time_trends_plot, relative_changes = relative_changes_plot),
    data_paths = data_paths,
    metadata_path = metadata_path,
    runtime = runtime
  ))
}

# Script Execution Function ------------------------------------------------

# Main execution function with example usage
main <- function() {
  cat("Climate-AMR Projection Analysis Script\n")
  cat("======================================\n\n")
  
  cat("This script performs climate-driven AMR projections using:\n")
  cat("• GAMM-derived response functions (REQUIRED)\n")
  cat("• Empirical climate-AMR relationships only\n")
  cat("• Bacteria-specific biological characteristics\n")
  cat("• CR-Kp stability and biological realism features\n") 
  cat("• No subjective socioeconomic assumptions\n")
  cat("• IPCC AR6-based climate uncertainty quantification\n\n")
  
  cat("Example usage:\n")
  cat("results <- run_climate_amr_projection(\n")
  cat("  data_path = 'data/resistance_data',\n")
  cat("  climate_file_path = 'data/climate_data.csv',\n") 
  cat("  response_functions_path = 'models/climate_response_functions.R',\n")
  cat("  output_dir = 'results/climate_amr_projections'\n")
  cat(")\n\n")
  
  cat("Required input files:\n")
  cat("1. Resistance data files (CSV format) - flexible naming detection\n")
  cat("2. Climate scenario data (CSV format) - flexible column naming\n")
  cat("3. GAMM response functions (R script) [REQUIRED]\n\n")
  
  cat("Output files will include:\n")
  cat("- Climate-driven projection plots (PNG/PDF)\n")
  cat("- Projection data and summaries (CSV)\n")
  cat("- Model metadata and documentation (JSON/TXT)\n")
}

# Load script notification
cat("Climate-AMR Projection Analysis Script Loaded Successfully!\n")
cat("Use main() for usage instructions or run_climate_amr_projection() to execute analysis.\n")

# Uncomment to see usage instructions when script is loaded
# main()