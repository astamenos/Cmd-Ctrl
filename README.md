# DomesticAbuseLGCP

Bayesian spatiotemporal analysis of domestic assault point pattern data from Minneapolis (2020–2022) using a marked Log Gaussian Cox Process (LGCP) implemented in NIMBLE.

## Overview

The model jointly estimates:
- **Intensity** `λ(s, t)` — where and when aggravated assaults occur, as a function of spatial covariates (population density), a low-rank Gaussian Process spatial random effect, and year random effects
- **Mark** `P(domestic | s, t)` — the probability that an assault is domestic in nature, modeled as a logistic regression on call type features and year

Crime events are filtered to aggravated assault calls from the Minneapolis open crime dataset and classified as domestic or non-domestic based on offense type and initial call type.

## Data

| Source | Content |
|--------|---------|
| Minneapolis ArcGIS Crime API | Assault incidents 2020–2022 |
| OpenStreetMap (`osmdata`) | Police station locations |
| US Census ACS 2019 (`tidycensus`) | Population and poverty by tract |
| TIGER/Line (`tigris`) | Water bodies, primary roads |
| `data/Minneapolis_Neighborhoods/` | Neighborhood boundary shapefile (checked in) |

Processed data is written to `data/MPLS_crime.csv` by `data_processing.R`.

## Usage

All scripts share globals set in `master.R` and must be run through it:

```r
source('master.R')
```

For a short diagnostic MCMC run (after data is loaded):

```r
result <- test_run(model_code, data_list, constants_list, cov_beta, cov_gamma, n_iter = 5000)
```

Production MCMC: 150,000 iterations, 50% burn-in, thinning by 10, 3 chains.

## Requirements

- R with packages: `nimble`, `spatstat`, `sf`, `tidyverse`, `tidycensus`, `tigris`, `osmdata`, `spdep`, `gstat`, `FNN`, `lme4`, `patchwork`, `cowplot`, `ggnewscale`, `scales`
- Census API key in `.Renviron` as `CENSUS_API_KEY` (required for ACS data)
- Internet access for live data pulls (crime API, OSM, TIGER)

## Structure

```
master.R      # Entry point; sets globals, sources all scripts in order
00_utils.R    # Helper functions, custom NIMBLE distributions
01_data.R     # Data ingestion, cleaning, spatial joins
02_model.R    # GP setup, NIMBLE model, MCMC
data/         # Raw inputs and processed CSV
plots/        # Output figures
```
