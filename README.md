# Spatiotemporal dynamics between climatic factors and antimicrobial resistance

This folder is a curated GitHub upload package prepared from the larger internal release-staging workspace. It retains the code and derived source-data components that are most directly needed for public presentation of the revised manuscript while omitting duplicated working files, older alternative outputs, and internal merged modeling panels.

## Included components

- `code/01_historical_associations/`
  - Primary historical Model C analysis
  - Lag-selection analysis
  - Figure 1 source-data export and summary scripts
- `code/02_spatial_heterogeneity/`
  - Figure 2 spatial-heterogeneity profile analysis
  - Climate-zone association analysis
  - Optimal-k clustering workflow
- `code/03_sensitivity_analyses/`
  - Alternative climatic-variable specification analysis
  - Climate-specification diagnostic checks
  - Detrended time-series analysis
  - Reanalysis-product comparison
  - Covariate-adjustment robustness analysis
- `code/04_future_projections/`
  - Projection preparation
  - Projection lag-selection analysis
  - Projection variance decomposition
  - Figure 3 future-projection workflow
- `data/source_data/`
  - Lag-selection source data for historical and projection workflows
  - Main figure source data for Figures 1-3
  - Climate-variable specification source data
  - Basis-dimension sensitivity source data
  - Detrended time-series source data
  - Reanalysis-product sensitivity source data
  - Covariate-adjustment robustness source data
- `figures/`
  - Main figure image files

## Intentionally omitted in this upload package

- Internal merged country-year modeling panels
- Redundant working versions and intermediate staging folders
- Additional appendix-only plotting suites that are not needed for the initial public release
- Alternative projection source-data branches not used for the current public upload set
- Editable design files that are not required for the initial public code-and-source-data release

## Data-sharing note

This package is designed as a curated, reproducibility-oriented public release. It provides code and publication-style derived source data aligned with the revised manuscript, but it is not intended to function as a one-command full rerun of the entire internal analysis workspace. The package does not redistribute all upstream raw inputs or all internally merged analytical panels.
