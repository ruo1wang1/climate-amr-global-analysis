###############################################################################
# Figure 2 spatial-heterogeneity workflow with optimal-k clustering
###############################################################################

suppressPackageStartupMessages({
  required_packages <- c(
    "tidyverse",
    "mgcv",
    "ComplexHeatmap",
    "cowplot",
    "grid",
    "openxlsx",
    "scales",
    "cluster",
    "RColorBrewer",
    "circlize"
  )

  cran_packages <- setdiff(required_packages, "ComplexHeatmap")
  missing_cran_packages <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran_packages) > 0) {
    install.packages(missing_cran_packages, repos = "https://cloud.r-project.org")
  }

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
  }

  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Unable to install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }

  invisible(lapply(required_packages, library, character.only = TRUE))
})

sanitize_file_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

significance_stars <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
}

relabel_clusters_by_centroid <- function(cluster_matrix, clusters) {
  clusters <- as.integer(clusters)
  centroid_tbl <- as_tibble(cluster_matrix, rownames = "NAME") %>%
    mutate(Cluster = clusters) %>%
    group_by(Cluster) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  ordered_clusters <- centroid_tbl$Cluster[order(rowMeans(centroid_tbl[, -1, drop = FALSE]), decreasing = TRUE)]
  remap <- setNames(seq_along(ordered_clusters), ordered_clusters)
  unname(remap[as.character(clusters)])
}

compute_gap_best_k <- function(cluster_matrix, min_k, max_k, nstart = 25L, b = 50L) {
  max_allowed_k <- min(as.integer(max_k), nrow(cluster_matrix) - 1L)
  min_allowed_k <- max(2L, as.integer(min_k))

  if (max_allowed_k < min_allowed_k) {
    stop("Insufficient rows to evaluate bounded k range for clustering.", call. = FALSE)
  }

  set.seed(123)
  gap_stat <- cluster::clusGap(
    cluster_matrix,
    FUN = kmeans,
    nstart = nstart,
    K.max = max_allowed_k,
    B = b
  )

  gap_tab <- as_tibble(gap_stat$Tab, rownames = "k_raw") %>%
    transmute(
      k = as.integer(k_raw),
      log_w = logW,
      e_log_w = E.logW,
      gap = gap,
      se = SE.sim
    )

  restricted_tab <- gap_tab %>%
    filter(k >= min_allowed_k, k <= max_allowed_k)

  selected_idx <- cluster::maxSE(restricted_tab$gap, restricted_tab$se)
  selected_k <- restricted_tab$k[selected_idx]

  km_fit <- kmeans(cluster_matrix, centers = selected_k, nstart = nstart)
  relabeled_clusters <- relabel_clusters_by_centroid(cluster_matrix, km_fit$cluster)

  list(
    selected_k = as.integer(selected_k),
    gap_table = gap_tab,
    restricted_gap_table = restricted_tab,
    clusters = setNames(relabeled_clusters, rownames(cluster_matrix)),
    kmeans = km_fit
  )
}

extract_country_level_effects <- function(data, model, bacteria_name) {
  pred <- predict(model, type = "terms", se.fit = TRUE)

  term_map <- c(
    "Temperature" = "s(TMP_scaled_lag)",
    "Relative Humidity" = "s(HUM_scaled_lag)",
    "Precipitation" = "s(PREC_scaled_lag)",
    "Wet Days" = "s(WET_scaled_lag)"
  )

  out <- vector("list", length(term_map))
  names(out) <- names(term_map)

  for (factor_name in names(term_map)) {
    term_name <- term_map[[factor_name]]
    term_idx <- match(term_name, colnames(pred$fit))
    if (is.na(term_idx)) {
      next
    }

    out[[factor_name]] <- tibble(
      NAME = data$NAME,
      Region = as.character(data$Region),
      climate_zone = as.character(data$climate_zone),
      effect = pred$fit[, term_idx],
      se = pred$se.fit[, term_idx]
    ) %>%
      group_by(NAME, Region) %>%
      summarise(
        n_observations = dplyr::n(),
        Climate_Zone_Majority = names(which.max(table(climate_zone))),
        effect = mean(effect, na.rm = TRUE),
        se = sqrt(mean(se^2, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        Bacteria = bacteria_name,
        Climate_Factor = factor_name,
        OR = exp(effect),
        CI_Lower = exp(effect - 1.96 * se),
        CI_Upper = exp(effect + 1.96 * se),
        Z_value = ifelse(se > 0, effect / se, NA_real_),
        P_Value = ifelse(se > 0, 2 * pnorm(-abs(Z_value)), NA_real_),
        Significance = significance_stars(P_Value),
        Effect_Direction = case_when(
          OR > 1 ~ "Positive",
          OR < 1 ~ "Negative",
          TRUE ~ "Null"
        )
      ) %>%
      select(
        Bacteria, Climate_Factor, NAME, Region, Climate_Zone_Majority,
        n_observations, effect, se, OR, CI_Lower, CI_Upper,
        Z_value, P_Value, Significance, Effect_Direction
      )
  }

  bind_rows(out)
}

extract_climate_zone_effects <- function(model, bacteria_name) {
  p_table <- as.data.frame(summary(model)$p.table)
  p_table$Term <- rownames(p_table)
  rownames(p_table) <- NULL

  p_table %>%
    filter(grepl("^climate_zone", Term)) %>%
    transmute(
      Bacteria = bacteria_name,
      Term = Term,
      Climate_Zone = gsub(" Zone$", "", gsub("^climate_zone", "", Term)),
      Estimate = Estimate,
      Std_Error = `Std. Error`,
      Test_Statistic = `t value`,
      P_Value = `Pr(>|t|)`,
      OR = exp(Estimate),
      CI_Lower = exp(Estimate - 1.96 * Std_Error),
      CI_Upper = exp(Estimate + 1.96 * Std_Error),
      Reference = "Polar",
      Significance = significance_stars(P_Value)
    ) %>%
    mutate(Climate_Zone = factor(Climate_Zone, levels = c("Temperate", "Tropical")))
}

create_factor_heatmap_inputs <- function(country_effects, factor_name, bacteria_order) {
  df <- country_effects %>%
    filter(Climate_Factor == factor_name) %>%
    mutate(Bacteria = factor(Bacteria, levels = bacteria_order)) %>%
    arrange(Bacteria)

  value_mat <- df %>%
    select(NAME, Bacteria, OR) %>%
    pivot_wider(names_from = Bacteria, values_from = OR) %>%
    tibble::column_to_rownames("NAME") %>%
    as.matrix()

  sig_mat <- df %>%
    select(NAME, Bacteria, Significance) %>%
    pivot_wider(names_from = Bacteria, values_from = Significance) %>%
    tibble::column_to_rownames("NAME") %>%
    as.matrix()

  cluster_mat <- value_mat
  for (j in seq_len(ncol(cluster_mat))) {
    col_mean <- mean(cluster_mat[, j], na.rm = TRUE)
    if (is.finite(col_mean)) {
      cluster_mat[is.na(cluster_mat[, j]), j] <- col_mean
    }
  }

  region_df <- df %>%
    select(NAME, Region) %>%
    distinct() %>%
    arrange(match(NAME, rownames(value_mat)))

  list(
    value_mat = value_mat,
    sig_mat = sig_mat,
    cluster_mat = cluster_mat,
    region_df = region_df
  )
}

build_factor_heatmap <- function(heatmap_bundle, factor_name, cluster_info, region_colors, cluster_distance_method, cluster_linkage_method) {
  cluster_levels <- paste0("Cluster ", seq_len(cluster_info$selected_k))
  cluster_palette <- setNames(
    grDevices::colorRampPalette(c("#5E3C99", "#4C78A8", "#36A7A6", "#7BC87C", "#FDE725"))(cluster_info$selected_k),
    cluster_levels
  )

  region_vec <- factor(
    heatmap_bundle$region_df$Region[match(rownames(heatmap_bundle$value_mat), heatmap_bundle$region_df$NAME)],
    levels = names(region_colors)
  )
  cluster_vec <- factor(
    paste0("Cluster ", cluster_info$clusters[rownames(heatmap_bundle$value_mat)]),
    levels = cluster_levels
  )

  display_numbers <- heatmap_bundle$sig_mat
  display_numbers[is.na(display_numbers)] <- ""

  heat_colors <- c("#2166AC", "#92C5DE", "#FFFFBF", "#F4A582", "#B2182B")
  col_fun <- circlize::colorRamp2(c(0.5, 0.75, 1.0, 1.25, 1.5), heat_colors)

  row_dist <- dist(heatmap_bundle$cluster_mat, method = cluster_distance_method)
  col_dist <- dist(t(heatmap_bundle$cluster_mat), method = cluster_distance_method)
  row_dend <- as.dendrogram(hclust(row_dist, method = cluster_linkage_method))
  col_dend <- as.dendrogram(hclust(col_dist, method = cluster_linkage_method))

  left_annotation <- ComplexHeatmap::rowAnnotation(
    Optimal_Cluster = cluster_vec,
    col = list(Optimal_Cluster = cluster_palette),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    border = TRUE,
    gp = grid::gpar(col = "white", lwd = 0.2),
    width = unit(3.6, "mm")
  )

  right_annotation <- ComplexHeatmap::rowAnnotation(
    Region = region_vec,
    col = list(Region = region_colors),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    border = TRUE,
    gp = grid::gpar(col = "white", lwd = 0.2),
    width = unit(3.6, "mm")
  )

  ht <- ComplexHeatmap::Heatmap(
    heatmap_bundle$value_mat,
    name = "Odds Ratio (OR)",
    col = col_fun,
    na_col = "white",
    cluster_rows = row_dend,
    cluster_columns = col_dend,
    left_annotation = left_annotation,
    right_annotation = right_annotation,
    show_heatmap_legend = FALSE,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 7.2),
    column_names_gp = grid::gpar(fontsize = 8.2),
    column_names_rot = 45,
    rect_gp = grid::gpar(col = "grey86", lwd = 0.55),
    row_dend_width = unit(22, "mm"),
    column_dend_height = unit(16, "mm"),
    border = TRUE,
    column_title = paste0(factor_name, " (k = ", cluster_info$selected_k, ")"),
    column_title_gp = grid::gpar(fontsize = 11.5, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      label <- display_numbers[i, j]
      if (nzchar(label)) {
        grid::grid.text(label, x = x, y = y, gp = grid::gpar(fontsize = 6.0, col = "black"))
      }
    }
  )

  heatmap_grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ht,
      newpage = FALSE,
      padding = unit(c(3, 3, 3, 3), "mm")
    ),
    wrap = TRUE
  )

  cluster_df <- tibble(
    NAME = names(cluster_info$clusters),
    Cluster = unname(cluster_info$clusters),
    Optimal_K = cluster_info$selected_k,
    Selection_Method = "Gap statistic (restricted maxSE)"
  ) %>%
    left_join(heatmap_bundle$region_df, by = "NAME") %>%
    mutate(Climate_Factor = factor_name)

  list(
    grob = heatmap_grob,
    clusters = cluster_df,
    cluster_palette = cluster_palette,
    col_fun = col_fun
  )
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
revision_root <- Sys.getenv("CLIMATE_AMR_WORKSPACE_ROOT", unset = repo_root)
code_dir <- file.path(revision_root, "code")
source_main_script <- file.path(code_dir, "01_historical_associations", "analysis_historical_associations_main_model.R")

existing_results_root <- file.path(
  revision_root,
  "outputs/historical_associations",
  "figure2_spatial_heterogeneity"
)

results_root <- file.path(
  revision_root,
  "outputs/historical_associations",
  "figure2_spatial_heterogeneity_optimal_k"
)

results_dirs <- list(
  model_summaries = file.path(results_root, "00_model_summaries"),
  country_effects = file.path(results_root, "01_country_effect_tables"),
  fig2a = file.path(results_root, "02_Figure2A_spatial_heatmaps"),
  fig2b = file.path(results_root, "03_Figure2B_climate_zone_associations"),
  fig2c = file.path(results_root, "04_Figure2C_gap_statistic_selection"),
  combined = file.path(results_root, "05_combined_figure2"),
  metadata = file.path(results_root, "06_metadata")
)

source_data_root <- file.path(
  revision_root,
  "data/source_data",
  "figure2_spatial_heterogeneity_profiles"
)

source_data_dirs <- list(
  csv = file.path(source_data_root, "01_csv"),
  workbook = file.path(source_data_root, "02_workbook")
)

invisible(lapply(c(results_dirs, source_data_dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

force_refit <- identical(tolower(Sys.getenv("MODELC_FORCE_REFIT", unset = "0")), "1")
cluster_distance_method <- trimws(Sys.getenv("MODELC_FIG2_CLUSTER_DISTANCE", unset = "euclidean"))
cluster_linkage_method <- trimws(Sys.getenv("MODELC_FIG2_CLUSTER_METHOD", unset = "complete"))
cluster_min_k <- suppressWarnings(as.integer(Sys.getenv("MODELC_FIG2_CLUSTER_MIN_K", unset = "2")))
cluster_max_k <- suppressWarnings(as.integer(Sys.getenv("MODELC_FIG2_CLUSTER_MAX_K", unset = "10")))
gap_b <- suppressWarnings(as.integer(Sys.getenv("MODELC_FIG2_GAP_B", unset = "50")))
kmeans_nstart <- suppressWarnings(as.integer(Sys.getenv("MODELC_FIG2_KMEANS_NSTART", unset = "25")))
run_only <- trimws(Sys.getenv("MODELC_RUN_ONLY", unset = ""))

if (!cluster_distance_method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
  cluster_distance_method <- "euclidean"
}
if (!cluster_linkage_method %in% c("complete", "ward.D2", "single", "average", "mcquitty", "median", "centroid")) {
  cluster_linkage_method <- "complete"
}
if (!is.finite(cluster_min_k) || is.na(cluster_min_k) || cluster_min_k < 2) {
  cluster_min_k <- 2L
}
if (!is.finite(cluster_max_k) || is.na(cluster_max_k) || cluster_max_k < cluster_min_k) {
  cluster_max_k <- 10L
}
if (!is.finite(gap_b) || is.na(gap_b) || gap_b < 10) {
  gap_b <- 50L
}
if (!is.finite(kmeans_nstart) || is.na(kmeans_nstart) || kmeans_nstart < 5) {
  kmeans_nstart <- 25L
}

bacteria_order <- c("3GCR-Ec", "3GCR-Kp", "CR-Ab", "CR-Ec", "CR-Kp", "CR-Pa")
climate_factor_order <- c("Temperature", "Relative Humidity", "Precipitation", "Wet Days")

climate_factor_colors <- c(
  "Temperature" = "#DD5F60",
  "Relative Humidity" = "#3CB371",
  "Precipitation" = "#9AC0CD",
  "Wet Days" = "#8A2BE2"
)

zone_colors <- c("Temperate" = "#F8766D", "Tropical" = "#00BFC4")

region_colors <- c(
  "East Asia and Pacific" = "#FFC300",
  "Europe and Central Asia" = "#90EE90",
  "Latin America and the Caribbean" = "#FFA500",
  "Middle East and North Africa" = "#8FBC8F",
  "North America" = "#FF0000",
  "South Asia" = "#DA70D6",
  "Sub-Saharan Africa" = "#BA55D3"
)

existing_input_candidates <- list(
  legacy_results_aliases = c(
    file.path(existing_results_root, "01_country_effect_tables", "ModelC_all_bacteria_country_level_climate_effects_with_or_ci.csv"),
    file.path(existing_results_root, "03_Figure2B_climate_zone_associations", "ModelC_Figure2B_climate_zone_effects.csv"),
    file.path(existing_results_root, "05_metadata", "ModelC_Figure2_model_metadata.csv")
  ),
  current_optimal_k = c(
    file.path(results_dirs$country_effects, "ModelC_all_bacteria_country_level_climate_effects_with_or_ci.csv"),
    file.path(results_dirs$fig2b, "ModelC_Figure2B_climate_zone_effects.csv"),
    file.path(results_dirs$metadata, "ModelC_Figure2_model_metadata.csv")
  )
)

selected_existing_inputs <- NULL
if (!force_refit) {
  for (candidate_name in names(existing_input_candidates)) {
    candidate_paths <- existing_input_candidates[[candidate_name]]
    if (all(file.exists(candidate_paths))) {
      selected_existing_inputs <- candidate_paths
      message("Reusing existing country-level inputs for the optimal-k Figure 2 suite from: ", candidate_name)
      break
    }
  }
}

if (!is.null(selected_existing_inputs)) {
  all_country_effects <- read.csv(selected_existing_inputs[[1]], stringsAsFactors = FALSE)
  all_climate_zone_effects <- read.csv(selected_existing_inputs[[2]], stringsAsFactors = FALSE)
  all_model_metadata <- read.csv(selected_existing_inputs[[3]], stringsAsFactors = FALSE)
} else {
  if (!file.exists(source_main_script)) {
    stop("Source model script not found: ", source_main_script, call. = FALSE)
  }

  Sys.setenv(MODELC_SKIP_MAIN = "1")
  source(source_main_script, local = FALSE)

  if (!exists("bacteria_specs") || !exists("prepare_data") || !exists("build_gamm_model")) {
    stop("Failed to load required model functions from source_main_script.", call. = FALSE)
  }

  if (nzchar(run_only)) {
    run_only_values <- trimws(strsplit(run_only, ",")[[1]])
    bacteria_specs <- Filter(function(x) x$code %in% run_only_values || x$title %in% run_only_values, bacteria_specs)
  }
  if (length(bacteria_specs) == 0) {
    stop("No bacteria matched MODELC_RUN_ONLY for optimal-k Figure 2 analysis.", call. = FALSE)
  }

  country_effect_list <- list()
  climate_zone_list <- list()
  model_metadata <- list()

  for (spec in bacteria_specs) {
    file_path <- file.path(input_data_dir, spec$file_name)
    if (!file.exists(file_path)) {
      stop("Input data file not found: ", file_path, call. = FALSE)
    }

    message("Preparing spatial inputs for ", spec$title, " ...")
    data_result <- prepare_data(file_path, spec$title)
    model <- build_gamm_model(data_result$data, spec$title, results_dirs)

    country_effects <- extract_country_level_effects(data_result$data, model, spec$title)
    climate_zone_effects <- extract_climate_zone_effects(model, spec$title)

    country_effect_list[[spec$title]] <- country_effects
    climate_zone_list[[spec$title]] <- climate_zone_effects
    model_metadata[[spec$title]] <- tibble(
      Bacteria = spec$title,
      Code = spec$code,
      Input_File = spec$file_name,
      Sample_Size = nrow(data_result$data),
      Available_PLS = paste(get_available_pls_components(data_result$data), collapse = ", ")
    )
  }

  all_country_effects <- bind_rows(country_effect_list)
  all_climate_zone_effects <- bind_rows(climate_zone_list)
  all_model_metadata <- bind_rows(model_metadata)
}

all_country_effects <- all_country_effects %>%
  mutate(
    Bacteria = factor(Bacteria, levels = bacteria_order),
    Climate_Factor = factor(Climate_Factor, levels = climate_factor_order)
  ) %>%
  arrange(Climate_Factor, NAME, Bacteria)

all_climate_zone_effects <- all_climate_zone_effects %>%
  mutate(
    Bacteria = factor(Bacteria, levels = bacteria_order),
    Climate_Zone = factor(as.character(Climate_Zone), levels = c("Temperate", "Tropical"))
  ) %>%
  arrange(Bacteria, Climate_Zone)

write.csv(
  all_country_effects,
  file.path(results_dirs$country_effects, "ModelC_all_bacteria_country_level_climate_effects_with_or_ci.csv"),
  row.names = FALSE
)

write.csv(
  all_climate_zone_effects,
  file.path(results_dirs$fig2b, "ModelC_Figure2B_climate_zone_effects.csv"),
  row.names = FALSE
)

write.csv(
  all_model_metadata,
  file.path(results_dirs$metadata, "ModelC_Figure2_model_metadata.csv"),
  row.names = FALSE
)

heatmap_gtables <- list()
cluster_assignments <- list()
optimal_k_summary <- list()
gap_plot_tables <- list()
cluster_palettes <- list()
heatmap_col_fun <- NULL

build_cluster_legend_grob <- function(cluster_palette) {
  cluster_legend <- ComplexHeatmap::Legend(
    title = "Clusters",
    labels = names(cluster_palette),
    legend_gp = grid::gpar(fill = unname(cluster_palette), col = "grey35"),
    title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 7.5),
    ncol = 1,
    grid_width = unit(3.6, "mm")
  )

  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      cluster_legend,
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      just = c("center", "center")
    ),
    wrap = TRUE
  )
}

for (factor_name in climate_factor_order) {
  bundle <- create_factor_heatmap_inputs(all_country_effects, factor_name, bacteria_order)
  cluster_info <- compute_gap_best_k(bundle$cluster_mat, cluster_min_k, cluster_max_k, kmeans_nstart, gap_b)

  gap_tbl <- cluster_info$restricted_gap_table %>%
    mutate(
      Climate_Factor = factor_name,
      Selected = k == cluster_info$selected_k,
      K_Min = cluster_min_k,
      K_Max = cluster_max_k,
      Selection_Method = "Gap statistic (restricted maxSE)"
    )

  gap_plot_tables[[factor_name]] <- gap_tbl
  optimal_k_summary[[factor_name]] <- tibble(
    Climate_Factor = factor_name,
    Optimal_K = cluster_info$selected_k,
    K_Min = cluster_min_k,
    K_Max = cluster_max_k,
    Gap_B = gap_b,
    KMeans_NStart = kmeans_nstart,
    Distance_Method = "kmeans_gap_on_raw_OR",
    Heatmap_Distance_Method = cluster_distance_method,
    Heatmap_Linkage_Method = cluster_linkage_method,
    Selection_Method = "Gap statistic (restricted maxSE)"
  )

  heatmap_result <- build_factor_heatmap(
    bundle,
    factor_name,
    cluster_info,
    region_colors,
    cluster_distance_method,
    cluster_linkage_method
  )

  heatmap_gtables[[factor_name]] <- heatmap_result$grob
  cluster_assignments[[factor_name]] <- heatmap_result$clusters
  cluster_palettes[[factor_name]] <- heatmap_result$cluster_palette
  if (is.null(heatmap_col_fun)) {
    heatmap_col_fun <- heatmap_result$col_fun
  }
}

all_cluster_assignments <- bind_rows(cluster_assignments) %>%
  arrange(Climate_Factor, Cluster, Region, NAME)

all_optimal_k_summary <- bind_rows(optimal_k_summary) %>%
  mutate(Climate_Factor = factor(Climate_Factor, levels = climate_factor_order))

all_gap_plot_data <- bind_rows(gap_plot_tables) %>%
  mutate(Climate_Factor = factor(Climate_Factor, levels = climate_factor_order))

write.csv(
  all_cluster_assignments,
  file.path(results_dirs$fig2a, "ModelC_Figure2A_country_cluster_assignments_optimal_k.csv"),
  row.names = FALSE
)

write.csv(
  all_optimal_k_summary,
  file.path(results_dirs$fig2a, "ModelC_Figure2A_optimal_k_summary.csv"),
  row.names = FALSE
)

write.csv(
  all_gap_plot_data,
  file.path(results_dirs$fig2c, "ModelC_Figure2C_gap_statistic_plot_data.csv"),
  row.names = FALSE
)

or_legend <- ComplexHeatmap::Legend(
  col_fun = heatmap_col_fun,
  title = "Odds Ratio (OR)",
  at = c(0.50, 0.75, 1.00, 1.25, 1.50),
  labels = sprintf("%.2f", c(0.50, 0.75, 1.00, 1.25, 1.50)),
  title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 8),
  legend_width = unit(36, "mm")
)

region_legend <- ComplexHeatmap::Legend(
  title = "Region",
  labels = names(region_colors),
  legend_gp = grid::gpar(fill = unname(region_colors), col = "grey40"),
  title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 8),
  ncol = 1
)

legend_pack <- ComplexHeatmap::packLegend(or_legend, region_legend, direction = "horizontal", gap = unit(8, "mm"))
legend_grob <- grid::grid.grabExpr(
  ComplexHeatmap::draw(legend_pack, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center", "center")),
  wrap = TRUE
)

legend_row <- cowplot::ggdraw() +
  cowplot::draw_grob(legend_grob, x = 0.08, y = 0.12, width = 0.84, height = 0.76)

heatmap_panels <- lapply(climate_factor_order, function(factor_name) {
  cowplot::ggdraw() + cowplot::draw_grob(heatmap_gtables[[factor_name]])
})
names(heatmap_panels) <- climate_factor_order

for (factor_name in climate_factor_order) {
  factor_plot <- cowplot::plot_grid(
    heatmap_panels[[factor_name]],
    cowplot::ggdraw() + cowplot::draw_grob(build_cluster_legend_grob(cluster_palettes[[factor_name]]), x = 0.10, y = 0.08, width = 0.80, height = 0.84),
    ncol = 2,
    rel_widths = c(1, 0.16)
  )
  factor_prefix <- file.path(
    results_dirs$fig2a,
    paste0("ModelC_Figure2A_", sanitize_file_stem(factor_name), "_spatial_heatmap_optimalk_nature")
  )

  ggsave(
    paste0(factor_prefix, ".pdf"),
    plot = factor_plot,
    width = 8.6,
    height = 10.8,
    dpi = 320
  )

  ggsave(
    paste0(factor_prefix, ".png"),
    plot = factor_plot,
    width = 8.6,
    height = 10.8,
    dpi = 320
  )
}

figure2a_panel_grid <- cowplot::plot_grid(
  plotlist = heatmap_panels,
  nrow = 1,
  align = "h",
  rel_widths = c(1.02, 0.95, 0.98, 1.02)
)

figure2a_combined <- cowplot::plot_grid(
  legend_row,
  figure2a_panel_grid,
  ncol = 1,
  rel_heights = c(0.10, 1)
)

ggsave(
  file.path(results_dirs$fig2a, "ModelC_Figure2A_spatial_heterogeneity_heatmaps_optimal_k.pdf"),
  plot = figure2a_combined,
  width = 10,
  height = 14,
  dpi = 320
)

ggsave(
  file.path(results_dirs$fig2a, "ModelC_Figure2A_spatial_heterogeneity_heatmaps_optimal_k.png"),
  plot = figure2a_combined,
  width = 10,
  height = 14,
  dpi = 320
)

fig2b_plot_data <- all_climate_zone_effects %>%
  mutate(
    OR = as.numeric(OR),
    OR_lower = as.numeric(CI_Lower),
    OR_upper = as.numeric(CI_Upper),
    significance = significance_stars(P_Value),
    label = sprintf("%.2f (%.2f-%.2f)%s", OR, OR_lower, OR_upper, significance)
  )

figure2b_plot <- ggplot(fig2b_plot_data, aes(x = Climate_Zone, y = OR, fill = Climate_Zone)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.48) +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.12, linewidth = 0.45) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray55", linewidth = 0.5) +
  geom_text(
    aes(
      label = label,
      y = ifelse(OR >= 1, OR_upper * 1.05, pmax(OR_lower * 0.95, 0.12))
    ),
    size = 2.6,
    hjust = 0.5,
    vjust = ifelse(fig2b_plot_data$OR >= 1, 0, 1)
  ) +
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = scales::number_format(accuracy = 0.01),
    limits = c(0.20, 8)
  ) +
  scale_fill_manual(values = zone_colors) +
  facet_wrap(~ Bacteria, ncol = 6, strip.position = "top") +
  labs(
    x = "Climate Zone",
    y = "Odds Ratio (OR)",
    fill = "Climate Zone",
    title = "Figure 2B. Climate zone associations with AMR strains"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size = unit(0.8, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(0.15, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(4, 6, 4, 6, unit = "pt")
  )

ggsave(
  file.path(results_dirs$fig2b, "ModelC_Figure2B_climate_zone_associations_optimal_k.pdf"),
  plot = figure2b_plot,
  width = 12,
  height = 4.8,
  dpi = 320
)

ggsave(
  file.path(results_dirs$fig2b, "ModelC_Figure2B_climate_zone_associations_optimal_k.png"),
  plot = figure2b_plot,
  width = 12,
  height = 4.8,
  dpi = 320
)

fig2c_label_data <- all_optimal_k_summary %>%
  left_join(
    all_gap_plot_data %>%
      group_by(Climate_Factor) %>%
      summarise(
        label_y = max(gap + se, na.rm = TRUE) + 0.04,
        .groups = "drop"
      ),
    by = "Climate_Factor"
  )

fig2c_selected_points <- all_gap_plot_data %>%
  filter(Selected)

figure2c_plot <- ggplot(
  all_gap_plot_data,
  aes(x = k, y = gap, group = 1)
) +
  geom_ribbon(
    aes(ymin = gap - se, ymax = gap + se, fill = Climate_Factor),
    alpha = 0.14,
    color = NA
  ) +
  geom_errorbar(
    aes(ymin = gap - se, ymax = gap + se),
    width = 0.12,
    linewidth = 0.35,
    color = "grey30"
  ) +
  geom_line(aes(color = Climate_Factor), linewidth = 0.8) +
  geom_point(aes(color = Climate_Factor), size = 2.2) +
  geom_point(
    data = fig2c_selected_points,
    aes(x = k, y = gap),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.4,
    stroke = 0.8,
    fill = "white",
    color = "black"
  ) +
  geom_vline(
    data = fig2c_label_data,
    aes(xintercept = Optimal_K),
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.45
  ) +
  geom_label(
    data = fig2c_label_data,
    aes(x = Optimal_K, y = label_y, label = paste0("k=", Optimal_K)),
    inherit.aes = FALSE,
    size = 3.1,
    label.size = 0.2,
    fill = "white"
  ) +
  facet_wrap(~ Climate_Factor, ncol = 2, scales = "free_y") +
  scale_color_manual(values = climate_factor_colors, guide = "none") +
  scale_fill_manual(values = climate_factor_colors, guide = "none") +
  scale_x_continuous(
    breaks = seq(cluster_min_k, cluster_max_k, by = 1),
    limits = c(cluster_min_k - 0.2, cluster_max_k + 0.2)
  ) +
  labs(
    x = "Candidate cluster number (k)",
    y = "Gap statistic",
    title = "Figure 2C. Bounded gap-statistic selection of optimal cluster number",
    subtitle = sprintf(
      "Points and ribbons show clusGap estimate ± 1 SE; selection rule = restricted maxSE over k = %d to %d",
      cluster_min_k,
      cluster_max_k
    )
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#F3F5F7", color = "grey70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 9),
    plot.margin = margin(6, 8, 6, 8, unit = "pt")
  )

ggsave(
  file.path(results_dirs$fig2c, "ModelC_Figure2C_gap_statistic_selection.pdf"),
  plot = figure2c_plot,
  width = 10.5,
  height = 6.8,
  dpi = 320
)

ggsave(
  file.path(results_dirs$fig2c, "ModelC_Figure2C_gap_statistic_selection.png"),
  plot = figure2c_plot,
  width = 10.5,
  height = 6.8,
  dpi = 320
)

ggsave(
  file.path(results_dirs$fig2b, "ModelC_Figure2B_climate_zone_associations_optimal_k_singlepanel.pdf"),
  plot = figure2b_plot,
  width = 11.2,
  height = 4.4,
  dpi = 320
)

ggsave(
  file.path(results_dirs$fig2b, "ModelC_Figure2B_climate_zone_associations_optimal_k_singlepanel.png"),
  plot = figure2b_plot,
  width = 11.2,
  height = 4.4,
  dpi = 320
)

ggsave(
  file.path(results_dirs$fig2c, "ModelC_Figure2C_gap_statistic_selection_singlepanel.pdf"),
  plot = figure2c_plot,
  width = 8.6,
  height = 6.4,
  dpi = 320
)

ggsave(
  file.path(results_dirs$fig2c, "ModelC_Figure2C_gap_statistic_selection_singlepanel.png"),
  plot = figure2c_plot,
  width = 8.6,
  height = 6.4,
  dpi = 320
)

figure2a_panel <- cowplot::ggdraw() +
  cowplot::draw_plot(figure2a_combined) +
  cowplot::draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 26, fontface = "bold")

figure2b_panel <- cowplot::ggdraw() +
  cowplot::draw_plot(figure2b_plot) +
  cowplot::draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 22, fontface = "bold")

figure2c_panel <- cowplot::ggdraw() +
  cowplot::draw_plot(figure2c_plot) +
  cowplot::draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 22, fontface = "bold")

figure2_combined <- cowplot::plot_grid(
  figure2a_panel,
  figure2b_panel,
  ncol = 1,
  rel_heights = c(1.22, 0.78)
)

ggsave(
  file.path(results_dirs$combined, "ModelC_Figure2_spatial_heterogeneity_panels_optimal_k.pdf"),
  plot = figure2_combined,
  width = 18,
  height = 24,
  dpi = 320
)

ggsave(
  file.path(results_dirs$combined, "ModelC_Figure2_spatial_heterogeneity_panels_optimal_k.png"),
  plot = figure2_combined,
  width = 18,
  height = 24,
  dpi = 320
)

figure2_combined_with_gap <- cowplot::plot_grid(
  figure2a_panel,
  cowplot::plot_grid(
    figure2b_panel,
    figure2c_panel,
    ncol = 2,
    rel_widths = c(1.15, 1)
  ),
  ncol = 1,
  rel_heights = c(1.42, 0.90)
)

ggsave(
  file.path(results_dirs$combined, "ModelC_Figure2_spatial_heterogeneity_panels_with_gap_optimal_k.pdf"),
  plot = figure2_combined_with_gap,
  width = 18,
  height = 26,
  dpi = 320
)

ggsave(
  file.path(results_dirs$combined, "ModelC_Figure2_spatial_heterogeneity_panels_with_gap_optimal_k.png"),
  plot = figure2_combined_with_gap,
  width = 18,
  height = 26,
  dpi = 320
)

write.csv(
  all_country_effects,
  file.path(source_data_dirs$csv, "Figure2_ModelC_Fig2A_country_effects.csv"),
  row.names = FALSE
)

write.csv(
  all_cluster_assignments,
  file.path(source_data_dirs$csv, "Figure2_ModelC_Fig2A_optimal_k_cluster_assignments.csv"),
  row.names = FALSE
)

write.csv(
  all_optimal_k_summary,
  file.path(source_data_dirs$csv, "Figure2_ModelC_Fig2A_optimal_k_summary.csv"),
  row.names = FALSE
)

write.csv(
  fig2b_plot_data,
  file.path(source_data_dirs$csv, "Figure2_ModelC_Fig2B_climate_zone_effects.csv"),
  row.names = FALSE
)

write.csv(
  all_gap_plot_data,
  file.path(source_data_dirs$csv, "Figure2_ModelC_Fig2C_gap_statistics.csv"),
  row.names = FALSE
)

write.csv(
  all_model_metadata,
  file.path(source_data_dirs$csv, "Figure2_ModelC_metadata.csv"),
  row.names = FALSE
)

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "Fig2A_country_effects")
openxlsx::writeData(workbook, "Fig2A_country_effects", all_country_effects)
openxlsx::addWorksheet(workbook, "Fig2A_cluster_assignments")
openxlsx::writeData(workbook, "Fig2A_cluster_assignments", all_cluster_assignments)
openxlsx::addWorksheet(workbook, "Fig2A_optimal_k_summary")
openxlsx::writeData(workbook, "Fig2A_optimal_k_summary", all_optimal_k_summary)
openxlsx::addWorksheet(workbook, "Fig2B_climate_zone_effects")
openxlsx::writeData(workbook, "Fig2B_climate_zone_effects", fig2b_plot_data)
openxlsx::addWorksheet(workbook, "Fig2C_gap_statistics")
openxlsx::writeData(workbook, "Fig2C_gap_statistics", all_gap_plot_data)
openxlsx::addWorksheet(workbook, "Metadata")
openxlsx::writeData(workbook, "Metadata", all_model_metadata)

openxlsx::saveWorkbook(
  workbook,
  file.path(source_data_dirs$workbook, "SourceData_figure2_spatial_heterogeneity_profiles.xlsx"),
  overwrite = TRUE
)

message("Optimal-k Figure 2 suite completed successfully.")
message("Results root: ", results_root)
message("Source data root: ", source_data_root)
