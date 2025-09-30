# Spatiotemporal Dynamics Between Climatic Factors and Antimicrobial Resistance: A Global Analysis

## Overview
This repository contains the core analysis code for our global study examining climate-antimicrobial resistance associations across 101 countries from 1999-2022.

## Repository Contents

### Code Files
- **`GAMM_analysis.R`** - Generalized Additive Mixed Models for nonlinear climate-AMR associations
- **`Spatial_heterogeneity.R`** - Spatial clustering and heterogeneity analysis  
- **`Projection_analysis.R`** - Future scenario projections under SSP pathways

## Key Findings
- Temperature exhibits latitude-dependent effects: protective in high-latitude regions, harmful in low-latitude areas
- Wet days frequency more predictive than total precipitation volume
- Future AMR burden strongly depends on emission scenarios: 51% increase under high-emission vs. 1% decrease under sustainable pathways by 2090s

## Methods Summary
- **Statistical Framework**: Generalized Additive Mixed Models (GAMMs) with lag optimization
- **Spatial Analysis**: Hierarchical clustering with heterogeneity assessment (I² statistics)
- **Future Projections**: Recursive prediction framework with Monte Carlo uncertainty quantification
- **Data Sources**: ResistanceMap, ERA5 climate reanalysis, CMIP6 projections

## System Requirements
- R (≥ 4.0.0)
- Required packages: mgcv, dplyr, ggplot2, cluster, metafor, parallel, forecast

## Acknowledgments
This research contributes to understanding global patterns of antimicrobial resistance in the context of climate change.
