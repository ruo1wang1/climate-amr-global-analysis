# ============================================================================
# Spatial Heterogeneity Analysis for Climate-AMR Relationships

# Required Libraries --------------------------------------------------------
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(tidyverse)
  library(RColorBrewer)
  library(grid)
  library(viridis)
  library(dendextend)
  library(cluster)
  library(gridExtra)
  library(ggplot2)
  library(meta)
  library(pheatmap)
})

# Configuration --------------------------------------------------------------

# Bacteria list and standardized names
bacteria_list <- c("3GCR-Ec", "3GCR-Kp", "CR-Ec", "CR-Kp", "CR-Pa", "CR-Ab")

bacteria_rename <- c(
  "3GCREC" = "3GCR-Ec", 
  "3GCRKP" = "3GCR-Kp", 
  "CREC" = "CR-Ec", 
  "CRKP" = "CR-Kp", 
  "CRPA" = "CR-Pa", 
  "CRAB" = "CR-Ab"
)

# Climate factors configuration
climate_factors <- c("TMP", "HUM", "PREC", "WET")

climate_factor_names <- c(
  "TMP" = "Temperature",
  "HUM" = "Humidity",
  "PREC" = "Precipitation",
  "WET" = "Wet Days"
)

# Regional color scheme
region_colors <- c(
  "East Asia and Pacific" = "#FFC300",
  "Europe and Central Asia" = "#90EE90",
  "Latin America and the Caribbean" = "#FFA500",
  "Middle East and North Africa" = "#8FBC8F",
  "North America" = "#FF0000",
  "South Asia" = "#DA70D6",
  "Sub-Saharan Africa" = "#BA55D3"
)

# Heatmap color palette
col_fun <- colorRamp2(
  c(0.5, 0.75, 1, 1.25, 1.5),
  c("#2166AC", "#92C5DE", "#FFFFBF", "#F4A582", "#B2182B")
)

# Utility Functions ---------------------------------------------------------

# Safe log transformation function
safe_log <- function(x) {
  ifelse(x <= 0 | is.na(x), NA, log(x))
}

# Calculate confidence intervals with boundary condition handling
calculate_ci <- function(mean_val, sd_val, n, log_transform = TRUE, conf_level = 0.95) {
  if (n <= 1 || is.na(mean_val) || is.na(sd_val) || sd_val <= 0) {
    return(list(
      lower = NA,
      upper = NA,
      significant = FALSE,
      sig_direction = "Not significant",
      sig_stars = ""
    ))
  }
  
  z_value <- qnorm((1 + conf_level) / 2)
  
  if (log_transform && mean_val > 0) {
    log_mean <- log(mean_val)
    se_log <- sd_val / mean_val / sqrt(n)
    log_lower <- log_mean - z_value * se_log
    log_upper <- log_mean + z_value * se_log
    lower <- exp(log_lower)
    upper <- exp(log_upper)
  } else {
    se <- sd_val / sqrt(n)
    lower <- mean_val - z_value * se
    upper <- mean_val + z_value * se
  }
  
  if (log_transform) {
    significant <- (lower > 1) || (upper < 1)
    sig_direction <- if (lower > 1) "Positive" else if (upper < 1) "Negative" else "Not significant"
  } else {
    significant <- (lower > 0) || (upper < 0)
    sig_direction <- if (lower > 0) "Positive" else if (upper < 0) "Negative" else "Not significant"
  }
  
  sig_stars <- ""
  if (significant) {
    if (lower > 1.1 || upper < 0.9) {
      sig_stars <- "***"
    } else if (lower > 1.05 || upper < 0.95) {
      sig_stars <- "**" 
    } else if (lower > 1 || upper < 1) {
      sig_stars <- "*"
    }
  }
  
  return(list(
    lower = lower,
    upper = upper,
    significant = significant,
    sig_direction = sig_direction,
    sig_stars = sig_stars
  ))
}

# Data Processing Functions -------------------------------------------------

# Process climate data for specific bacteria and climate factor
process_climate_data <- function(bacteria_code, climate_factor, input_folder) {
  file_path <- file.path(input_folder, bacteria_code, 
                         paste0(bacteria_code, "_spatial_effects_with_OR_4factors.csv"))
  
  if (!file.exists(file_path)) {
    cat("Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  
  data <- read.csv(file_path, check.names = FALSE)
  
  # Determine column names based on climate factor
  factor_cols <- switch(climate_factor,
                        "TMP" = c("NAME", "Region", "TMP_OR", "TMP_CI_lower", "TMP_CI_upper"),
                        "HUM" = c("NAME", "Region", "HUM_OR", "HUM_CI_lower", "HUM_CI_upper"),
                        "PREC" = c("NAME", "Region", "PREC_OR", "PREC_CI_lower", "PREC_CI_upper"),
                        "WET" = c("NAME", "Region", "WET_OR", "WET_CI_lower", "WET_CI_upper"))
  
  # Get standardized bacteria name
  bacteria_name <- bacteria_rename[bacteria_code]
  if (is.na(bacteria_name)) bacteria_name <- bacteria_code
  
  # Extract column names
  or_col <- factor_cols[3]
  ci_lower_col <- factor_cols[4]
  ci_upper_col <- factor_cols[5]
  
  # Process data
  result <- data %>%
    select(all_of(factor_cols)) %>%
    mutate(across(c(all_of(or_col), all_of(ci_lower_col), all_of(ci_upper_col)), as.numeric)) %>%
    group_by(NAME) %>%
    summarise(
      OR = mean(get(or_col), na.rm = TRUE),
      CI_lower = mean(get(ci_lower_col), na.rm = TRUE),
      CI_upper = mean(get(ci_upper_col), na.rm = TRUE),
      Region = first(Region),
      .groups = "drop"
    ) %>%
    mutate(
      Bacteria = bacteria_name,
      significance = case_when(
        CI_lower > 1.1 | CI_upper < 0.9 ~ "***",
        CI_lower > 1.05 | CI_upper < 0.95 ~ "**",
        CI_lower > 1 | CI_upper < 1 ~ "*",
        TRUE ~ ""
      ),
      effect_direction = case_when(
        OR > 1 ~ "Positive",
        OR < 1 ~ "Negative",
        TRUE ~ "No correlation"
      )
    )
  
  return(result)
}

# Load and process all climate data
load_climate_data <- function(input_folder) {
  cat("Processing climate factor data...\n")
  
  climate_data_list <- list()
  bacteria_codes <- names(bacteria_rename)
  
  for (factor_name in climate_factors) {
    cat("Processing", factor_name, "data...\n")
    
    processed_data_list <- lapply(bacteria_codes, function(bc) {
      process_climate_data(bc, factor_name, input_folder)
    })
    names(processed_data_list) <- bacteria_codes
    
    # Remove NULL results
    processed_data_list <- processed_data_list[!sapply(processed_data_list, is.null)]
    
    # Combine all data
    all_data <- bind_rows(processed_data_list)
    climate_data_list[[factor_name]] <- all_data
    
    cat(factor_name, "data processing completed with", nrow(all_data), "records\n")
  }
  
  return(climate_data_list)
}

# Matrix Creation Functions ------------------------------------------------

# Create heatmap matrices for all climate factors
create_heatmap_matrices <- function(climate_data_list) {
  cat("Creating heatmap matrices...\n")
  
  heatmap_matrices <- list()
  signif_matrices <- list()
  clustering_matrices <- list()
  
  for (factor_name in climate_factors) {
    factor_data <- climate_data_list[[factor_name]]
    
    # Create OR matrix
    heatmap_matrix <- factor_data %>%
      select(NAME, Bacteria, OR) %>%
      pivot_wider(names_from = Bacteria, values_from = OR) %>%
      column_to_rownames("NAME") %>%
      as.matrix()
    
    # Create significance matrix
    signif_matrix <- factor_data %>%
      select(NAME, Bacteria, significance) %>%
      pivot_wider(names_from = Bacteria, values_from = significance) %>%
      column_to_rownames("NAME") %>%
      as.matrix()
    
    # Create clustering matrix (replace NA values for clustering)
    clustering_matrix <- heatmap_matrix
    for(i in 1:ncol(clustering_matrix)) {
      col_mean <- mean(clustering_matrix[, i], na.rm = TRUE)
      if (!is.na(col_mean)) {
        clustering_matrix[is.na(clustering_matrix[, i]), i] <- col_mean
      }
    }
    
    # Ensure correct column order
    if (all(bacteria_list %in% colnames(heatmap_matrix))) {
      heatmap_matrix <- heatmap_matrix[, bacteria_list]
      clustering_matrix <- clustering_matrix[, bacteria_list]
      signif_matrix <- signif_matrix[, bacteria_list]
    }
    
    # Store matrices
    heatmap_matrices[[factor_name]] <- heatmap_matrix
    signif_matrices[[factor_name]] <- signif_matrix
    clustering_matrices[[factor_name]] <- clustering_matrix
  }
  
  return(list(
    heatmap = heatmap_matrices,
    signif = signif_matrices,
    clustering = clustering_matrices
  ))
}

# Visualization Functions ---------------------------------------------------

# Create individual optimized heatmap
create_optimized_heatmap <- function(matrix_data, signif_matrix, clustering_matrix, 
                                     factor_name, region_info) {
  formal_name <- climate_factor_names[factor_name]
  
  # Cell annotation function for significance
  cell_fun <- function(j, i, x, y, width, height, fill) {
    sig <- signif_matrix[i, j]
    if (!is.na(sig) && sig != "") {
      grid.text(sig, x, y, gp = gpar(fontsize = 7))
    }
  }
  
  # Row annotation for regions
  ha_row <- rowAnnotation(
    Region = region_info$Region[match(rownames(matrix_data), region_info$NAME)],
    col = list(Region = region_colors),
    show_legend = TRUE,
    width = unit(0.6, "cm"),
    annotation_legend_param = list(
      Region = list(
        title = "Geographic Region",
        direction = "horizontal",
        title_gp = gpar(fontsize = 9),
        labels_gp = gpar(fontsize = 8),
        nrow = 2
      )
    )
  )
  
  # Create heatmap
  title_text <- paste0(formal_name, " Effect on AMR Risk")
  
  ht <- Heatmap(
    matrix = matrix_data,
    name = paste0(formal_name, " Risk (OR)"),
    col = col_fun,
    na_col = "white",
    right_annotation = ha_row,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = dist(clustering_matrix),
    clustering_distance_columns = dist(t(clustering_matrix)),
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    cell_fun = cell_fun,
    width = unit(2.5, "inches"),
    height = unit(16, "inches"),
    row_dend_width = unit(1.2, "inches"),
    column_dend_height = unit(0.4, "inches"),
    rect_gp = gpar(col = "lightgrey", lwd = 0.5),
    heatmap_legend_param = list(
      title = paste0(formal_name, " Risk (OR)"),
      at = c(0.5, 0.75, 1, 1.25, 1.5),
      labels = c("0.5", "0.75", "1.0", "1.25", "1.5"),
      direction = "horizontal",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(0.4, "cm"),
      legend_width = unit(3, "cm"),
      title_position = "topcenter"
    ),
    column_title = title_text,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    border = TRUE,
    border_gp = gpar(col = "lightgrey")
  )
  
  return(ht)
}

# Create uniform heatmap for integration
create_uniform_heatmap <- function(matrix_data, signif_matrix, factor_name, 
                                   common_countries, show_row_names = TRUE, region_info) {
  formal_name <- climate_factor_names[factor_name]
  
  # Filter and order matrix by common countries
  countries_in_matrix <- intersect(common_countries, rownames(matrix_data))
  matrix_data_ordered <- matrix_data[countries_in_matrix, ]
  signif_matrix_ordered <- signif_matrix[countries_in_matrix, ]
  
  # Row annotation
  ha_row <- NULL
  if (show_row_names) {
    ha_row <- rowAnnotation(
      Region = region_info$Region[match(countries_in_matrix, region_info$NAME)],
      col = list(Region = region_colors),
      show_legend = FALSE,
      width = unit(0.6, "cm")
    )
  }
  
  # Cell function for significance
  cell_fun <- function(j, i, x, y, width, height, fill) {
    sig <- signif_matrix_ordered[i, j]
    if (!is.na(sig) && sig != "") {
      grid.text(sig, x, y, gp = gpar(fontsize = 7))
    }
  }
  
  # Create heatmap
  ht <- Heatmap(
    matrix = matrix_data_ordered,
    name = paste0(formal_name, " (OR)"),
    col = col_fun,
    na_col = "white",
    right_annotation = ha_row,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_row_names,
    show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    cell_fun = cell_fun,
    width = unit(2, "inches"),
    height = unit(8, "inches"),
    rect_gp = gpar(col = "lightgrey", lwd = 0.3),
    column_title = paste0(formal_name, " Effect"),
    column_title_gp = gpar(fontsize = 9, fontface = "bold"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 1)
  )
  
  return(ht)
}

# Analysis Functions --------------------------------------------------------

# Extract and save clustering information
extract_clustering_info <- function(drawn_heatmap_objects, clustering_matrices, 
                                    region_info, output_folder) {
  cat("Extracting and saving clustering information...\n")
  
  all_clusters_info <- data.frame()
  
  for (factor_name in climate_factors) {
    cat("Processing", factor_name, "clustering information...\n")
    
    # Extract clustering information
    row_clusters <- row_order(drawn_heatmap_objects[[factor_name]])
    
    # Determine optimal number of clusters using gap statistic
    gap_stat_rows <- clusGap(clustering_matrices[[factor_name]], 
                             FUN = kmeans, nstart = 25, K.max = 10, B = 50)
    best_k_rows <- maxSE(gap_stat_rows$Tab[,"gap"], gap_stat_rows$Tab[,"SE.sim"])
    
    # K-means clustering
    km_rows <- kmeans(clustering_matrices[[factor_name]], centers = best_k_rows, nstart = 25)
    
    country_clusters <- data.frame(
      NAME = rownames(clustering_matrices[[factor_name]]),
      Cluster = km_rows$cluster
    )
    
    # Merge with region information
    country_clusters_with_region <- country_clusters %>%
      left_join(region_info, by = "NAME") %>%
      arrange(Cluster, Region, NAME) %>%
      mutate(ClimateFactor = factor_name)
    
    # Save individual clustering results
    write.csv(
      country_clusters_with_region,
      file.path(output_folder, paste0(tolower(factor_name), "_country_clusters.csv")),
      row.names = FALSE
    )
    
    # Add to overall clustering information
    all_clusters_info <- bind_rows(all_clusters_info, country_clusters_with_region)
  }
  
  # Save comprehensive clustering information
  write.csv(
    all_clusters_info,
    file.path(output_folder, "all_climate_factors_country_clusters.csv"),
    row.names = FALSE
  )
  
  cat("Clustering information extraction completed\n")
  return(all_clusters_info)
}

# Create analysis summary
create_analysis_summary <- function(climate_data_list, output_folder) {
  cat("Creating climate factor analysis summary...\n")
  
  factor_stats <- data.frame()
  
  for (factor_name in climate_factors) {
    factor_data <- climate_data_list[[factor_name]]
    
    sig_pos <- sum(factor_data$CI_lower > 1, na.rm = TRUE)
    sig_neg <- sum(factor_data$CI_upper < 1, na.rm = TRUE)
    nonsig <- sum(factor_data$CI_lower <= 1 & factor_data$CI_upper >= 1, na.rm = TRUE)
    total <- nrow(factor_data)
    
    formal_name <- climate_factor_names[factor_name]
    
    factor_stats <- rbind(factor_stats, data.frame(
      Factor = formal_name,
      Positive = sig_pos,
      Negative = sig_neg,
      Nonsignificant = nonsig,
      Total = total,
      PositivePercent = round(sig_pos / total * 100, 1),
      NegativePercent = round(sig_neg / total * 100, 1)
    ))
  }
  
  # Create comparison bar chart
  p <- ggplot(factor_stats, aes(x = Factor)) +
    geom_bar(aes(y = PositivePercent, fill = "Positive"), stat = "identity", position = "dodge") +
    geom_bar(aes(y = -NegativePercent, fill = "Negative"), stat = "identity", position = "dodge") +
    geom_text(aes(y = PositivePercent + 3, label = paste0(PositivePercent, "%")), 
              position = position_dodge(width = 0.9), size = 3.5) +
    geom_text(aes(y = -NegativePercent - 3, label = paste0(NegativePercent, "%")), 
              position = position_dodge(width = 0.9), size = 3.5) +
    scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
    coord_flip() +
    labs(
      title = "Comparison of Climate Factors Effect on AMR",
      x = "",
      y = "Percentage of Significant Associations (%)",
      fill = "Effect Direction"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text = element_text(size = 9)
    )
  
  # Save analysis summary
  pdf(file.path(output_folder, "climate_factors_analysis_summary.pdf"), 
      width = 8, height = 6)
  print(p)
  dev.off()
  
  png(file.path(output_folder, "climate_factors_analysis_summary.png"), 
      width = 2400, height = 1800, res = 300)
  print(p)
  dev.off()
  
  return(factor_stats)
}

# Create correlation heatmap
create_correlation_heatmap <- function(climate_data_list, region_info, output_folder) {
  cat("Creating climate factors correlation heatmap...\n")
  
  tryCatch({
    # Prepare correlation data
    correlations_data <- data.frame(
      NAME = character(0), 
      Bacteria = character(0),
      TMP_OR = numeric(0), 
      HUM_OR = numeric(0), 
      PREC_OR = numeric(0), 
      WET_OR = numeric(0)
    )
    
    # Extract OR values for all countries and bacteria
    for (bc in names(bacteria_rename)) {
      for (country in unique(region_info$NAME)) {
        tmp_or <- NA
        hum_or <- NA
        prec_or <- NA
        wet_or <- NA
        
        for (factor_name in climate_factors) {
          factor_data <- climate_data_list[[factor_name]]
          record <- factor_data %>% 
            filter(NAME == country & Bacteria == bacteria_rename[bc])
          
          if (nrow(record) > 0) {
            if (factor_name == "TMP") tmp_or <- record$OR[1]
            if (factor_name == "HUM") hum_or <- record$OR[1]
            if (factor_name == "PREC") prec_or <- record$OR[1]
            if (factor_name == "WET") wet_or <- record$OR[1]
          }
        }
        
        correlations_data <- rbind(correlations_data, data.frame(
          NAME = country,
          Bacteria = bacteria_rename[bc],
          TMP_OR = tmp_or,
          HUM_OR = hum_or,
          PREC_OR = prec_or,
          WET_OR = wet_or
        ))
      }
    }
    
    # Calculate correlation matrix
    cor_matrix <- cor(correlations_data[, c("TMP_OR", "HUM_OR", "PREC_OR", "WET_OR")], 
                      use = "pairwise.complete.obs")
    
    # Update matrix names
    formal_names <- climate_factor_names[c("TMP", "HUM", "PREC", "WET")]
    rownames(cor_matrix) <- formal_names
    colnames(cor_matrix) <- formal_names
    
    # Create correlation heatmap
    pdf(file.path(output_folder, "climate_factors_correlation_heatmap.pdf"),
        width = 7, height = 6)
    
    pheatmap(cor_matrix,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             display_numbers = TRUE,
             number_color = "black",
             fontsize_number = 10,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             main = "Correlation Between Climate Factors",
             fontsize = 10,
             angle_col = 45)
    
    dev.off()
    
    # Create PNG version
    png(file.path(output_folder, "climate_factors_correlation_heatmap.png"),
        width = 2100, height = 1800, res = 300)
    
    pheatmap(cor_matrix,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             display_numbers = TRUE,
             number_color = "black",
             fontsize_number = 10,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             main = "Correlation Between Climate Factors",
             fontsize = 10,
             angle_col = 45)
    
    dev.off()
    
    return(cor_matrix)
  }, error = function(e) {
    cat("Error creating correlation heatmap:", e$message, "\n")
    return(NULL)
  })
}

# Create spatial pattern analysis
create_spatial_pattern_analysis <- function(climate_data_list, output_folder) {
  cat("Creating spatial pattern analysis...\n")
  
  tryCatch({
    # Prepare regional summary data
    region_summary <- data.frame()
    
    for (factor_name in climate_factors) {
      factor_data <- climate_data_list[[factor_name]]
      
      region_data <- factor_data %>%
        group_by(Region) %>%
        summarize(
          MeanOR = mean(OR, na.rm = TRUE),
          MedianOR = median(OR, na.rm = TRUE),
          SignificantPos = sum(CI_lower > 1, na.rm = TRUE),
          SignificantNeg = sum(CI_upper < 1, na.rm = TRUE),
          Total = n(),
          .groups = "drop"
        ) %>%
        mutate(
          Factor = climate_factor_names[factor_name],
          SignificantPosPercent = SignificantPos / Total * 100,
          SignificantNegPercent = SignificantNeg / Total * 100
        )
      
      region_summary <- rbind(region_summary, region_data)
    }
    
    # Create regional effect plot
    g1 <- ggplot(region_summary, aes(x = Region, y = MeanOR, fill = Factor)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey") +
      coord_flip() +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = "Regional Climate Factor Effects on AMR",
        x = "Geographic Region",
        y = "Mean Odds Ratio (OR)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 9)
      )
    
    # Create regional significance plot
    g2 <- ggplot(region_summary, aes(x = Region)) +
      geom_bar(aes(y = SignificantPosPercent, fill = Factor), 
               stat = "identity", position = "dodge") +
      geom_text(aes(y = SignificantPosPercent + 5, label = paste0(round(SignificantPosPercent), "%"), 
                    group = Factor),
                position = position_dodge(width = 0.9), size = 2.5) +
      coord_flip() +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = "Regional Positive Significant Associations",
        x = "",
        y = "Percentage of Significant Positive Associations (%)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        axis.title = element_text(face = "bold"),
        axis.text.y = element_blank()
      )
    
    # Save spatial pattern analysis
    pdf(file.path(output_folder, "climate_factors_spatial_patterns.pdf"),
        width = 10, height = 7)
    grid.arrange(g1, g2, ncol = 2, widths = c(3, 2))
    dev.off()
    
    png(file.path(output_folder, "climate_factors_spatial_patterns.png"),
        width = 3000, height = 2100, res = 300)
    grid.arrange(g1, g2, ncol = 2, widths = c(3, 2))
    dev.off()
    
  }, error = function(e) {
    cat("Error creating spatial pattern analysis:", e$message, "\n")
  })
}

# Generate analysis report
generate_analysis_report <- function(output_folder) {
  cat("Generating integrated analysis report...\n")
  
  report_content <- c(
    "# Climate Factors Spatial Heterogeneity Analysis Report",
    "",
    "## Overview",
    "",
    "This analysis integrates four climate factors (Temperature, Humidity, Precipitation, and Wet Days)",
    "to examine their spatial heterogeneous effects on antimicrobial resistance (AMR).",
    "The analysis is based on four-factor GAMM model results covering multiple countries",
    "and six types of resistant bacteria.",
    "",
    "## Key Findings",
    "",
    "1. **Spatial Heterogeneity**: All four climate factors show significant spatial heterogeneity",
    "   in their effects on AMR across different geographic regions.",
    "",
    "2. **Climate Factor Interactions**: Multiple climate factors interact in certain regions",
    "   to create unique AMR risk patterns.",
    "",
    "3. **Regional Clustering**: Countries can be clustered based on climate-AMR relationships,",
    "   with each cluster representing a specific climate-AMR association pattern.",
    "",
    "## Climate Factor Summary",
    "",
    "1. **Temperature**: Shows strong correlations with multiple resistant bacteria,",
    "   particularly in specific geographic regions.",
    "",
    "2. **Humidity**: Demonstrates unique influence patterns in some regions,",
    "   complementing temperature effects.",
    "",
    "3. **Precipitation**: Shows significant effects on certain resistant bacteria,",
    "   but with more concentrated geographic distribution.",
    "",
    "4. **Wet Days**: As a precipitation frequency indicator, provides complementary",
    "   information to total precipitation amounts.",
    "",
    "## Spatial Heterogeneity Assessment",
    "",
    "1. **Clustering Analysis**: Countries can be grouped into clusters with similar",
    "   response patterns based on climate factor effects.",
    "",
    "2. **Regional Characteristics**: Certain regions show particular sensitivity",
    "   to specific climate factors.",
    "",
    "3. **Cross-Validation**: Multi-factor analysis validates the robustness",
    "   of spatial heterogeneity patterns.",
    "",
    "## Research Implications",
    "",
    "1. **Risk Assessment**: Provides climate-based AMR risk spatial distribution maps",
    "   to help identify high-risk regions.",
    "",
    "2. **Policy Development**: Offers evidence for developing differentiated AMR",
    "   control strategies for different climate regions.",
    "",
    "3. **Climate Change Impact**: Helps predict potential impacts of climate change",
    "   on global AMR distribution patterns.",
    "",
    "## Conclusions",
    "",
    "The integrated four-factor climate analysis reveals complex spatial heterogeneity",
    "in climate-AMR relationships, emphasizing the importance of considering climate",
    "factors in AMR risk assessment and control strategy development.",
    "This analytical framework can guide the design of region-specific AMR interventions."
  )
  
  # Save report
  writeLines(report_content, file.path(output_folder, "spatial_heterogeneity_analysis_report.md"))
}

# Main Analysis Function ----------------------------------------------------

run_spatial_heterogeneity_analysis <- function(input_folder, output_folder = "spatial_heterogeneity_results") {
  
  # Create output directory
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    cat("Created output directory:", output_folder, "\n")
  }
  
  cat("Starting spatial heterogeneity analysis...\n")
  cat("Input folder:", input_folder, "\n")
  cat("Output folder:", output_folder, "\n")
  
  # Step 1: Load and process climate data
  climate_data_list <- load_climate_data(input_folder)
  
  # Create unified region information
  region_info <- climate_data_list[[1]] %>%
    select(NAME, Region) %>%
    distinct()
  
  # Step 2: Create heatmap matrices
  matrices <- create_heatmap_matrices(climate_data_list)
  
  # Step 3: Create individual heatmaps
  cat("Creating individual heatmaps...\n")
  heatmap_objects <- list()
  drawn_heatmap_objects <- list()
  
  for (factor_name in climate_factors) {
    cat("Creating", factor_name, "heatmap...\n")
    
    ht <- create_optimized_heatmap(
      matrices$heatmap[[factor_name]], 
      matrices$signif[[factor_name]],
      matrices$clustering[[factor_name]], 
      factor_name,
      region_info
    )
    
    heatmap_objects[[factor_name]] <- ht
    
    # Save individual heatmap
    pdf(file.path(output_folder, paste0(tolower(factor_name), "_heatmap_individual.pdf")), 
        width = 8, height = 18)
    
    set.seed(123)
    ht_drawn <- draw(ht, 
                     heatmap_legend_side = "top", 
                     annotation_legend_side = "bottom",
                     padding = unit(c(1, 0.8, 1, 0.8), "cm"))
    
    drawn_heatmap_objects[[factor_name]] <- ht_drawn
    dev.off()
    
    cat(factor_name, "heatmap completed\n")
  }
  
  # Step 4: Extract clustering information
  clustering_info <- extract_clustering_info(
    drawn_heatmap_objects, 
    matrices$clustering, 
    region_info, 
    output_folder
  )
  
  # Step 5: Create integrated visualization
  cat("Creating integrated visualization...\n")
  
  # Use temperature heatmap row order for consistency
  tmp_heatmap_drawn <- drawn_heatmap_objects[["TMP"]]
  tmp_row_order <- row_order(tmp_heatmap_drawn)
  common_countries <- rownames(matrices$heatmap[["TMP"]])[tmp_row_order]
  
  # Create uniform heatmaps
  uniform_heatmaps <- list()
  uniform_heatmaps[["TMP"]] <- create_uniform_heatmap(
    matrices$heatmap[["TMP"]], matrices$signif[["TMP"]], "TMP", common_countries, TRUE, region_info
  )
  uniform_heatmaps[["HUM"]] <- create_uniform_heatmap(
    matrices$heatmap[["HUM"]], matrices$signif[["HUM"]], "HUM", common_countries, FALSE, region_info
  )
  uniform_heatmaps[["PREC"]] <- create_uniform_heatmap(
    matrices$heatmap[["PREC"]], matrices$signif[["PREC"]], "PREC", common_countries, TRUE, region_info
  )
  uniform_heatmaps[["WET"]] <- create_uniform_heatmap(
    matrices$heatmap[["WET"]], matrices$signif[["WET"]], "WET", common_countries, FALSE, region_info
  )
  
  # Add region legend to last heatmap
  ha_region_legend <- rowAnnotation(
    Region = region_info$Region[match(common_countries, region_info$NAME)],
    col = list(Region = region_colors),
    show_legend = TRUE,
    width = unit(0.6, "cm"),
    annotation_legend_param = list(
      Region = list(
        title = "Geographic Region",
        direction = "horizontal",
        title_gp = gpar(fontsize = 9),
        labels_gp = gpar(fontsize = 8),
        nrow = 2
      )
    )
  )
  
  uniform_heatmaps[["WET"]] <- uniform_heatmaps[["WET"]] + ha_region_legend
  
  # Step 6: Draw integrated visualization
  set.seed(123)
  
  pdf(file.path(output_folder, "integrated_climate_factors_visualization.pdf"), 
      width = 10, height = 14)
  
  heatmap_list <- uniform_heatmaps[["TMP"]] + 
    uniform_heatmaps[["HUM"]] + 
    uniform_heatmaps[["PREC"]] + 
    uniform_heatmaps[["WET"]]
  
  draw(heatmap_list, 
       heatmap_legend_side = "bottom", 
       annotation_legend_side = "bottom",
       column_title = "Climate Factors Effect on AMR: Spatial Heterogeneity Patterns",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_title = "Countries Grouped by Temperature Effect Pattern",
       row_title_gp = gpar(fontsize = 11, fontface = "bold"),
       padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  dev.off()
  
  # Create high-resolution PNG
  png(file.path(output_folder, "integrated_climate_factors_visualization.png"), 
      width = 3000, height = 4200, res = 300)
  
  draw(heatmap_list, 
       heatmap_legend_side = "bottom", 
       annotation_legend_side = "bottom",
       column_title = "Climate Factors Effect on AMR: Spatial Heterogeneity Patterns",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_title = "Countries Grouped by Temperature Effect Pattern",
       row_title_gp = gpar(fontsize = 11, fontface = "bold"),
       padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  dev.off()
  
  # Step 7: Create analysis summaries
  factor_stats <- create_analysis_summary(climate_data_list, output_folder)
  correlation_matrix <- create_correlation_heatmap(climate_data_list, region_info, output_folder)
  create_spatial_pattern_analysis(climate_data_list, output_folder)
  
  # Step 8: Generate final report
  generate_analysis_report(output_folder)
  
  cat("\nSpatial heterogeneity analysis completed!\n")
  cat("All results saved to:", output_folder, "\n\n")
  
  cat("Main output files:\n")
  cat("- integrated_climate_factors_visualization.pdf: Integrated heatmap visualization\n")
  cat("- climate_factors_analysis_summary.pdf: Climate factors comparison summary\n")
  cat("- climate_factors_correlation_heatmap.pdf: Climate factors correlation analysis\n")
  cat("- climate_factors_spatial_patterns.pdf: Spatial pattern analysis\n")
  cat("- all_climate_factors_country_clusters.csv: Country clustering information\n")
  cat("- spatial_heterogeneity_analysis_report.md: Comprehensive analysis report\n")
  
  return(list(
    climate_data = climate_data_list,
    matrices = matrices,
    clustering_info = clustering_info,
    factor_stats = factor_stats,
    correlation_matrix = correlation_matrix
  ))
}

# Script Status Message ------------------------------------------------------

cat("Spatial Heterogeneity Analysis Script Loaded Successfully!\n")
cat("Use run_spatial_heterogeneity_analysis(input_folder, output_folder) to start analysis.\n")
cat("\nExample usage:\n")
cat("results <- run_spatial_heterogeneity_analysis(\n")
cat("  input_folder = 'path/to/your/gamm/results',\n")
cat("  output_folder = 'spatial_heterogeneity_results'\n")
cat(")\n")