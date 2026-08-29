# Spatiotemporal dynamics between climatic factors and antimicrobial resistance

This repository contains the code, derived source-data files, and figure assets associated with the manuscript:

**Spatiotemporal dynamics between climatic factors and antimicrobial resistance**

It includes materials for the descriptive global landscape, historical multi-country analysis, spatial-heterogeneity analyses, sensitivity analyses and future climate-driven projections.

This repository is intended as a reproducibility-oriented analysis companion rather than as a packaged end-user software product. It provides the R scripts used for the main analytical modules, together with manuscript-oriented derived source data and figure assets.

## Repository structure

- `code/`
  - `00_descriptive_landscape/`: Figure 1 global geographical and temporal summaries.
  - `01_historical_associations/`: primary historical Model C scripts, lag selection, Figure 2 turning-point analysis, observation-level influence diagnostics and source-data checks.
  - `02_spatial_heterogeneity/`: scripts for Figure 3, cluster-based spatial heterogeneity analyses, and clustering-sensitivity analyses.
  - `03_sensitivity_analyses/`: scripts for BH–FDR control, leave-one-country-out refitting and summaries, basis-dimension, climate-variable-specification, detrending, reanalysis-product and covariate-adjustment analyses.
  - `04_future_projections/`: scripts for simplified historical-model preparation, projection framework, variance decomposition, and Figure 4 generation.
- `data/source_data/`
  - Publication-oriented derived source data used for main figures and selected additional analyses.
- `figures/`
  - Current main-figure image files. Supplementary figure and table assets are not included.

## Data-sharing scope

This repository provides **derived, publication-oriented source data** and analysis code for reproducibility of the manuscript. It includes selected processed figure- and table-level source data suitable for public release, but it does not distribute the Figure 1 country-year input, the full internal modelling panels or all upstream datasets. The repository should therefore be interpreted as a reproducibility-oriented code and derived-data repository rather than a redistribution point for upstream inputs.

The repository includes code modules for the descriptive, historical, spatial, sensitivity and projection analyses, together with curated derived source data and four main-figure image assets. Supplementary analyses are represented by code and derived source data rather than duplicated figure or table files.

## System requirements

### Hardware

- No non-standard hardware is required.
- The code was designed for a standard desktop or laptop computer.
- GPU acceleration is not required.

### Operating system

- The repository scripts were checked on `macOS 26.1`.
- The scripts should also run on other operating systems supported by R, provided that the listed package dependencies are installed and file paths are adjusted appropriately.

### Software

- `R 4.3.3` was used for the tested local run documented in this repository.
- Base or recommended R components used by the scripts include `grid` and `splines`.
- Additional R package dependencies used across the repository include:
  - `tidyverse`
  - `mgcv`
  - `ggplot2`
  - `ggh4x`
  - `patchwork`
  - `scales`
  - `zoo`
  - `gridExtra`
  - `cowplot`
  - `viridis`
  - `kableExtra`
  - `tidytext`
  - `cluster`
  - `digest`
  - `openxlsx`
  - `readxl`
  - `ggrepel`
  - `htmltools`
  - `pagedown`
  - `png`
  - `svglite`
  - `ragg`
  - `sf`
  - `rnaturalearth`
  - `systemfonts`
  - `pdftools`
  - `writexl`
  - `forcats`
  - `RColorBrewer`
  - `circlize`
  - `extrafont`
  - `flextable`
  - `officer`
  - `tinytex`
  - `ComplexHeatmap`
- Selected export scripts use `pagedown` for HTML-to-PDF rendering and may require a local Chrome-, Chromium-, or Microsoft Edge-based browser.
- Selected table-export scripts use `tinytex`; TeX-based PDF rendering may require a local TinyTeX or compatible TeX installation.

## Installation guide

### 1. Download the repository

Clone or download the repository to a local directory.

```bash
git clone https://github.com/ruo1wang1/climate-amr-global-analysis.git
cd climate-amr-global-analysis
```

### 2. Install R package dependencies

The example below installs the main non-base dependencies used across the repository.

```r
cran_packages <- c(
  "tidyverse", "mgcv", "ggh4x", "patchwork", "scales", "zoo", "gridExtra",
  "cowplot", "viridis", "kableExtra", "tidytext", "cluster", "digest",
  "openxlsx", "readxl", "ggrepel", "htmltools", "pagedown", "png",
  "svglite", "ragg", "pdftools", "sf", "rnaturalearth", "systemfonts",
  "writexl", "forcats", "RColorBrewer", "circlize", "extrafont",
  "flextable", "officer", "tinytex"
)
install.packages(cran_packages)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")
```

### 3. Optional environment variables

Most scripts detect the repository root automatically. If needed, the following environment variables can be set before running scripts:

- `CLIMATE_AMR_WORKSPACE_ROOT`
  - Overrides the repository root used for reading inputs and writing outputs.
- `MODELC_RUN_ONLY`
  - Restricts the primary historical Model C script to selected phenotypes.
- `MODELC_OUTPUT_SUFFIX`
  - Appends a suffix to selected Model C output folders.
- `CLIMATE_AMR_CMIP6_ROOT`
  - Overrides the CMIP6 archive path for projection scripts.
- `CLIMATE_AMR_WET_DAYS_INPUT`
  - Overrides the yearly wet-days panel path for projection scripts.
- `FIGURE1_INPUT_FILE`
  - Points the Figure 1 workflow to the local country-year input CSV; this input is not distributed in the repository.
- `FIGURE1_OUTPUT_ROOT`
  - Overrides the repository-local Figure 1 output directory.
- `FIGURE2_INPUT_ROOT`
  - Points the Figure 2 workflow to the fitted models, model-ready panels and required analysis inputs.
- `FIGURE2_OUTPUT_ROOT`
  - Overrides the repository-local Figure 2 output directory.
- `MODEL_C_LOCO_INPUT_ROOT`
  - Points the full LOCO refitting workflow to the Figure 2 analysis output containing the frozen Model C objects, model-ready panels and turning-point inputs.
- `MODEL_C_LOCO_REFIT_OUTPUT_ROOT`
  - Overrides the output directory for the full country-deletion refits.
- `MODEL_C_LOCO_CORES`
  - Sets the number of parallel workers used by the LOCO refitting workflow.
- `MODEL_C_LOCO_SOURCE_ROOT`
  - Points the LOCO summary workflow to the validated country-deletion refit outputs.
- `MODEL_C_LOCO_OUTPUT_ROOT`
  - Overrides the repository-local LOCO output directory.

### Typical installation time

- On a normal desktop computer with a working R installation and internet access for package installation, setup typically takes about `10-20 minutes`.
- If all R packages are already installed, repository setup is effectively immediate after download.

## Instructions for use

### Reproducing the provided manuscript outputs

The repository is organised by analytical module:

1. `code/00_descriptive_landscape`
   - Global geographical distributions and temporal summaries for Figure 1.
2. `code/01_historical_associations`
   - Historical lag selection, primary GAMM fitting, Figure 2 turning-point analysis and observation-level influence diagnostics.
3. `code/02_spatial_heterogeneity`
   - Country-level OR profiling, clustering, and climate-zone analyses for Figure 3.
4. `code/03_sensitivity_analyses`
   - BH–FDR control, country-level influence refitting and summaries, basis-dimension, model-specification, detrending, reanalysis-product and covariate-adjustment analyses.
5. `code/04_future_projections`
   - Projection preparation, lag uncertainty, variance decomposition, and Figure 4 future trajectories.

Scripts write their outputs under repository-local `outputs/` folders. Source-data folders under `data/source_data/` provide the manuscript-exported tables used for reporting and figure source data.

### Running the software on your own data

This repository is not a general-purpose packaged application. To adapt the workflow to your own data:

1. Prepare country-year analytical inputs that follow the manuscript pipeline structure.
2. For historical Model C analyses, place phenotype-specific model-ready CSV files in:
   - `outputs/historical_associations/model_ready_inputs/`
3. Follow the filename conventions expected by the primary historical script:
   - `3GCREC_model_ready_data.csv`
   - `3GCRKP_model_ready_data.csv`
   - `CRAB_model_ready_data.csv`
   - `CREC_model_ready_data.csv`
   - `CRKP_model_ready_data.csv`
   - `CRPA_model_ready_data.csv`
4. Ensure that accompanying lag-summary inputs and any required derived source-data tables are available in the repository paths used by the scripts.
5. Set `CLIMATE_AMR_WORKSPACE_ROOT` if you want inputs or outputs to be resolved relative to a directory other than the repository root.
6. Run the relevant module script with `Rscript`.

Because the manuscript workflow uses model-ready analytical panels generated from multi-source surveillance, climate, and covariate preprocessing, adapting the workflow to external data may require additional preprocessing to match the expected input structure.

## Reproduction instructions

### Broader manuscript reproduction

This repository provides the code for the main analytical modules and the manuscript source-data exports. Full end-to-end reproduction of every intermediate result also requires access to the surveillance and climate inputs described in the manuscript and Supplementary Tables S1 and S2, together with selected large intermediate analytical inputs that are not distributed in this repository.

Within those constraints, the recommended execution order is:

Before running Figure 1, set `FIGURE1_INPUT_FILE` to the local country-year input CSV described in the script. The workflow validates the expected rows, keys, phenotypes, countries and year range before analysis.

1. `code/00_descriptive_landscape/figure1_descriptive_landscape.R`
2. `code/01_historical_associations/analysis_lag_selection.R`
3. `code/01_historical_associations/analysis_historical_associations_main_model.R`
4. `code/01_historical_associations/figure2_turning_point_analysis.R`
5. `code/02_spatial_heterogeneity/analysis_spatial_heterogeneity_profiles.R`
6. `code/02_spatial_heterogeneity/figure3_spatial_heterogeneity_optimal_k.R`
7. `code/02_spatial_heterogeneity/figure3_climate_zone_associations.R`
8. `code/03_sensitivity_analyses/analysis_basis_dimension_sensitivity.R`
9. `code/03_sensitivity_analyses/analysis_climate_variable_specification.R`
10. `code/03_sensitivity_analyses/analysis_climate_specification_checks.R`
11. `code/03_sensitivity_analyses/analysis_covariate_adjustment_robustness.R`
12. `code/03_sensitivity_analyses/analysis_detrended_time_series.R`
13. `code/03_sensitivity_analyses/analysis_reanalysis_product_configurations.R`
14. `code/03_sensitivity_analyses/analysis_multiple_testing_bh_fdr.R`
15. `code/03_sensitivity_analyses/analysis_country_level_influence_loco.R`
16. `code/03_sensitivity_analyses/summarise_country_level_influence_loco.R`
17. `code/04_future_projections/analysis_projection_preparation.R`
18. `code/04_future_projections/analysis_projection_lag_selection.R`
19. `code/04_future_projections/analysis_projection_variance_decomposition.R`
20. `code/04_future_projections/figure4_future_projections_bootstrap.R`

The LOCO refitting workflow uses the phenotype-specific frozen Model C objects and model-ready country–year panels generated by the Figure 2 workflow. These inputs are supplied through `MODEL_C_LOCO_INPUT_ROOT`; the downstream summary script reads the resulting refit outputs through `MODEL_C_LOCO_SOURCE_ROOT`.

## Versions tested

The updated scripts were checked with:

- `macOS 26.1`
- `R 4.3.3`

## Code availability

The repository corresponding to the manuscript is:

- https://github.com/ruo1wang1/climate-amr-global-analysis

## License

This repository is released under the MIT License. See the `LICENSE` file for the full license text.
