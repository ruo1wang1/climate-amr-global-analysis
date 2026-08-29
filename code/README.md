# Analysis scripts

This folder groups the analysis scripts by study component:

- `00_descriptive_landscape`
- `01_historical_associations`
- `02_spatial_heterogeneity`
- `03_sensitivity_analyses`
- `04_future_projections`

Current revision-specific scripts include:

- `00_descriptive_landscape/figure1_descriptive_landscape.R`: Figure 1 global geographical and temporal summaries.
- `01_historical_associations/figure2_turning_point_analysis.R`: Figure 2 turning-point and observation-level influence workflow.
- `03_sensitivity_analyses/analysis_multiple_testing_bh_fdr.R`: BH–FDR adjustment for the 24 primary climate smooth-term tests.
- `03_sensitivity_analyses/analysis_country_level_influence_loco.R`: full country-deletion Model C refitting and diagnostic export.
- `03_sensitivity_analyses/summarise_country_level_influence_loco.R`: downstream aggregation and visualization of validated LOCO refit outputs.

For system requirements, package installation, demo instructions, expected outputs, and broader reproduction guidance, see the repository-level [README.md](../README.md).
