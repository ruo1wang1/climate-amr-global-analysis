# Spatiotemporal dynamics between climatic factors and antimicrobial resistance

This repository contains the code, derived source-data files, and figure assets associated with the manuscript:

**Spatiotemporal dynamics between climatic factors and antimicrobial resistance**

It includes materials for the historical multi-country analysis, spatial-heterogeneity analyses, robustness and sensitivity analyses, and future climate-driven projection analyses.

This repository is intended as a reproducibility-oriented analysis companion rather than as a packaged end-user software product. It provides the R scripts used for the main analytical modules, together with manuscript-oriented derived source data and figure assets. Some large upstream climate archives, internal model-ready inputs, and restricted-access surveillance inputs are not redistributed here.

## Repository structure

- `code/`
  - `01_historical_associations/`: primary historical Model C scripts, lag-selection scripts, and source-data/table-generation scripts.
  - `02_spatial_heterogeneity/`: scripts for Figure 2, cluster-based spatial heterogeneity analyses, and clustering-sensitivity analyses.
  - `03_sensitivity_analyses/`: scripts for detrending, reanalysis-product comparison, covariate-robustness, model-specification, and country-effect robustness analyses.
  - `04_future_projections/`: scripts for simplified historical-model preparation, projection framework, variance decomposition, and Figure 3 generation.
  - `05_helpers/`: auxiliary plotting and helper scripts used for source-data or additional analysis outputs.
- `data/source_data/`
  - Publication-oriented derived source data used for main figures and selected additional analyses.
- `figures/`
  - Main-figure image files, additional figure exports, and editable figure source files.
- `docs/`
  - Release inventory and documentation related to repository preparation.

## Data-sharing scope

This repository provides **derived, publication-oriented source data** and analysis code for reproducibility of the manuscript.

It does **not** redistribute the full merged country-year analytical panel datasets used internally to fit the models, because these files require an additional review of data-sharing constraints, upstream source permissions, and release format. Some upstream AMR surveillance and covariate sources were obtained from public or research-access platforms with their own access conditions.

Accordingly:

- the repository includes processed figure- and table-level source data suitable for public release;
- it excludes the internal merged phenotype-level modelling panels at this stage;
- it should be interpreted as a reproducibility-oriented code and derived-data repository rather than a redistribution point for all upstream input datasets.

The repository includes code modules for the historical, spatial, sensitivity, and projection analyses, together with the curated source-data package for the main figures and additional robustness analyses, and the corresponding figure exports and editable figure assets.

## System requirements

### Hardware

- No non-standard hardware is required.
- The code was designed for a standard desktop or laptop computer.
- GPU acceleration is not required.

### Operating system

- The quick-start demo described below was tested on `macOS 26.1`.
- The scripts should also run on other operating systems supported by R, provided that the listed package dependencies are installed and file paths are adjusted appropriately.

### Software

- `R 4.3.3` was used for the tested local run documented in this repository.
- Base or recommended R components used by the scripts include `grid` and `splines`.
- Additional R package dependencies used across the repository include:
  - `tidyverse`
  - `mgcv`
  - `ggplot2`
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
  - `png`
  - `writexl`
  - `forcats`
  - `RColorBrewer`
  - `circlize`
  - `extrafont`
  - `flextable`
  - `officer`
  - `ComplexHeatmap`

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
  "tidyverse", "mgcv", "patchwork", "scales", "zoo", "gridExtra",
  "cowplot", "viridis", "kableExtra", "tidytext", "cluster", "digest",
  "openxlsx", "readxl", "ggrepel", "png", "writexl", "forcats",
  "RColorBrewer", "circlize", "extrafont", "flextable", "officer"
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

### Typical installation time

- On a normal desktop computer with a working R installation and internet access for package installation, setup typically takes about `10-20 minutes`.
- If all R packages are already installed, repository setup is effectively immediate after download.

## Demo dataset and quick start

This repository contains a small real derived dataset suitable for a quick demonstration run:

- Demo input file:
  - `data/source_data/figure1_historical_associations/01_csv/Figure1_ModelC_Fig1B_summary.csv`
- Demo script:
  - `code/01_historical_associations/figure1_statistical_summary.R`

### Demo command

Run the following command from the repository root:

```bash
Rscript code/01_historical_associations/figure1_statistical_summary.R
```

### Expected demo output

The script reads the provided Figure 1B source-data summary table and writes:

- `outputs/historical_associations/analysis_historical_associations_main_model/05_figure1B_statistical_summary/ModelC_GAMM_Figure1B_statistical_summary.pdf`
- `outputs/historical_associations/analysis_historical_associations_main_model/05_figure1B_statistical_summary/ModelC_GAMM_Figure1B_statistical_summary.png`
- `outputs/historical_associations/analysis_historical_associations_main_model/05_figure1B_statistical_summary/ModelC_GAMM_Figure1B_plot_ready_data.csv`

### Expected demo run time

- The demo run was verified locally on a normal desktop computer and completed in about `2-3 seconds` after package installation.

## Instructions for use

### Reproducing the provided manuscript outputs

The repository is organised by analytical module:

1. `code/01_historical_associations`
   - Historical lag selection, primary GAMM fitting, and Figure 1 source-data workflows.
2. `code/02_spatial_heterogeneity`
   - Country-level OR profiling, clustering, and climate-zone analyses for Figure 2.
3. `code/03_sensitivity_analyses`
   - Basis-dimension, model-specification, detrending, reanalysis-product, and covariate-adjustment sensitivity analyses.
4. `code/04_future_projections`
   - Projection preparation, lag uncertainty, variance decomposition, and Figure 3 future trajectories.

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

Because the manuscript workflow uses model-ready analytical panels generated from multi-source surveillance, climate, and covariate preprocessing, adapting the workflow to external data may require additional preprocessing beyond what is redistributed here.

## Reproduction instructions

### Minimal reproduction

The easiest verified reproduction step is the Figure 1B demo described above.

### Broader manuscript reproduction

This repository provides the code for the main analytical modules and the manuscript source-data exports. Full end-to-end reproduction of every intermediate result requires access to:

- the upstream surveillance inputs described in the manuscript and Supplementary Tables S1 and S2,
- the climate archives used to construct the model-ready country-year panels,
- and selected large intermediate analytical inputs that are not redistributed in this repository.

Within those constraints, the recommended execution order is:

1. `code/01_historical_associations/analysis_lag_selection.R`
2. `code/01_historical_associations/analysis_historical_associations_main_model.R`
3. `code/01_historical_associations/figure1_source_data_export.R`
4. `code/01_historical_associations/figure1_statistical_summary.R`
5. `code/02_spatial_heterogeneity/analysis_spatial_heterogeneity_profiles.R`
6. `code/02_spatial_heterogeneity/figure2_spatial_heterogeneity_optimal_k.R`
7. `code/02_spatial_heterogeneity/figure2_climate_zone_associations.R`
8. `code/03_sensitivity_analyses/analysis_basis_dimension_sensitivity.R`
9. `code/03_sensitivity_analyses/analysis_climate_specification_checks.R`
10. `code/03_sensitivity_analyses/analysis_covariate_adjustment_robustness.R`
11. `code/03_sensitivity_analyses/analysis_detrended_time_series.R`
12. `code/03_sensitivity_analyses/analysis_reanalysis_product_configurations.R`
13. `code/03_sensitivity_analyses/analysis_model_structure_comparison.R`
14. `code/03_sensitivity_analyses/analysis_country_effect_robustness.R`
15. `code/04_future_projections/analysis_projection_preparation.R`
16. `code/04_future_projections/analysis_projection_lag_selection.R`
17. `code/04_future_projections/analysis_projection_variance_decomposition.R`
18. `code/04_future_projections/figure3_future_projections_bootstrap.R`

## Versions tested

The quick-start demo in this repository update was tested with:

- `macOS 26.1`
- `R 4.3.3`

## Code availability

The repository corresponding to the manuscript is:

- https://github.com/ruo1wang1/climate-amr-global-analysis

## License

This repository is released under the MIT License. See the `LICENSE` file for the full license text.
