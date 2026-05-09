# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bayesian spatiotemporal analysis of domestic assault point pattern data from Minneapolis (2020–2022) using a marked Log Gaussian Cox Process (LGCP) implemented in NIMBLE.

## Running the Analysis

All scripts must be run through `master.R`, which sets globals and sources the others in order:

```r
# In R or RStudio
source('master.R')
```

Never source `00_utils.R`, `01_data.R`, or `02_model.R` directly — they depend on globals defined in `master.R` (`crs_use`, `years`, `years_fac`, `call_types`, `domestic_terms`).

**Test a short MCMC run** (uses `test_run()` from `utils.R`):
```r
# After sourcing master.R through data_processing.R:
result <- test_run(model_code, data_list, constants_list, cov_beta, cov_gamma, n_iter = 5000)
```

## Architecture

### Script Execution Order

```
master.R
├── sets globals: crs_use (UTM 15N), years (2020-2022), call_types, domestic_terms
├── source('00_utils.R')   — helper functions and NIMBLE custom distributions
├── source('01_data.R')    — data ingestion, cleaning, spatial joins
└── source('02_model.R')   — model definition, GP setup, MCMC
```

### Coordinate System

All spatial objects are projected to `"+proj=utm +zone=15 +datum=WGS84"` then scaled to **kilometers** via `mutate(geometry = geometry/1000)`. This is a deliberate design choice — NIMBLE operates on these km-scale coordinates directly.

### Data Flow

1. **`data_processing.R`** pulls live data (Minneapolis Crime ArcGIS API, OSM police stations, ACS Census tracts, tigris water/roads), filters to aggravated assaults 2020–2022, deduplicates by datetime+geometry, derives `domestic_flag` from call type text, jitters duplicate locations, and writes `data/MPLS_crime.csv`.

2. **`spatiotemp_mark_LGCP.R`** reads `data/MPLS_crime.csv`, constructs the spatstat point process objects, fits a variogram to estimate the GP range (`phi_E0`), places adaptive knots (denser in high-crime areas), builds the Cholesky-factored GP basis matrices (`L1`, `C_cross_C_inv_*`), and runs the NIMBLE MCMC.

### Model Structure

The NIMBLE model (`model_code`) is a **marked LGCP**:
- **Intensity** (`lambda`): log-linear in fixed effects (`beta`), population density offset, zero-mean GP spatial random effect (`W1star`), and centered year random effects (`delta_1`)
- **Mark** (`Y`, binary domestic/non-domestic): logistic regression with `gamma` coefficients and year RE (`delta_2`)
- **GP**: low-rank approximation using `n_K` knots; non-centered parameterization (`z1 ~ N(0,1)`, then `W1star = L1 @ z1`)
- **MCMC**: block RW samplers — `beta` and `gamma` use empirical proposal covariances updated between runs; `(log_sigma, z1)` sampled jointly; `(u1, u2)` use AF_slice

### Key Functions (`utils.R`)

| Function | Purpose |
|----------|---------|
| `test_run()` | Short diagnostic MCMC run; returns updated proposal covariances |
| `posterior_intensity()` | Plots posterior mean intensity surface by offense type × year |
| `check_gp_covariate_collinearity()` | Flags knots with posterior correlation > threshold vs. fixed effects |
| `visualize_problematic_knots()` | Map of collinear knots |
| `inits_fn()` | Random initial values from MLE estimates |
| `d_lgcp` / `r_lgcp` | Custom NIMBLE distribution for Poisson LGCP likelihood |
| `jsd_im()` | Jensen-Shannon divergence between two spatstat intensity images |

## External Data Dependencies

`data_processing.R` fetches live from:
- Minneapolis Crime API (ArcGIS FeatureServer) — requires internet
- OpenStreetMap via `osmdata` — requires internet
- US Census via `tidycensus` — requires a Census API key in `.Renviron` as `CENSUS_API_KEY`
- `tigris` for TIGER water/roads shapefiles — cached locally after first run

The Minneapolis Neighborhoods shapefile is checked in at `data/Minneapolis_Neighborhoods/`.

## MCMC Configuration

Production run settings (in `spatiotemp_mark_LGCP.R`):
- `niter = 150000`, `nburnin = 75000`, `thin = 10`, `nchains = 3`
- Yields `(niter - nburnin) / thin = 7500` posterior samples per chain

Convergence is assessed via `coda::gelman.diag()` on `beta` and `gamma` parameters.
