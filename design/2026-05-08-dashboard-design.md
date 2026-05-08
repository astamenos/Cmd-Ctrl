# Dashboard Design Spec
**Project:** DomesticAbuseLGCP  
**Date:** 2026-05-08  
**Status:** Approved

---

## Decisions Log

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Hosting | GitHub Pages | Already in use for firm website |
| Framework | Quarto multi-page site (`format: html`) | Stays in R, fits existing Quarto workflow, reads as credibility artifact |
| Interactivity | Fixed layout (htmlwidgets: leaflet, plotly) | Simpler first version; user filtering deferred as future feature |
| Live data | GitHub Actions nightly cron → re-render → deploy | Scheduled re-render, not browser-side fetching |
| LGCP artifacts | Pre-computed RDS committed by separate process | NIMBLE too intensive for Actions; dashboard reads existing RDS only |
| Audience | Public; positions author as expert and informed source | |

---

## Site Structure

Five pages in `_quarto.yml` navbar:

1. **Overview** — value boxes, monthly trend chart
2. **Hotspots** — KDE leaflet map, three year tabs
3. **LGCP Results** — posterior intensity PNGs, covariate table
4. **Domestic Violence** — education content, local resources
5. **Methods** — model spec, data sources, GitHub link

The Education page is a top-level navbar item, not nested under Methods.

---

## Data Flow

```
GitHub Actions (nightly cron)
│
├── R: source('master.R') → source('01_data.R')
│   └── ArcGIS FeatureServer → data/MPLS_crime.csv
│
├── load output/lgcp_posterior.rds  [read-only; written by separate LGCP job]
│
├── quarto render dashboard.qmd --output-dir docs/
│
└── git commit docs/ → gh-pages branch → Pages deploys
```

LGCP artifacts (`output/lgcp_posterior.rds`) are committed by a separate process on its own schedule. The Actions render job never re-runs NIMBLE.

---

## Script Architecture

```
master.R                  — globals, sources helpers
00_utils.R                — unchanged; helper functions
01_data.R                 — unchanged; live data fetch → data/MPLS_crime.csv
02_model.R                — unchanged; NIMBLE MCMC (offline, separate schedule)
03_dashboard.R            — NEW; all data prep for dashboard render:
                            loads MPLS_crime.csv, computes KDE,
                            formats tables, builds ggplot/leaflet objects
dashboard.qmd             — NEW; Quarto site entry point, sources 03_dashboard.R
_quarto.yml               — NEW; site config, navbar, theme
```

`dashboard.qmd` contains no heavy computation — all prep is in `03_dashboard.R`.

---

## Page Components

### Overview
- Value box: total incidents (labeled "as of [render date]")
- Value box: % flagged domestic
- Value box: most recent incident date
- `plotly` line chart: monthly counts 2020–2022, colored domestic vs. non-domestic
- Source: `MPLS_crime.csv`, computed at render time

### Hotspots
- `leaflet` map with KDE raster overlay (`density.ppp` → raster → `addRasterImage`)
- Neighborhood boundaries and police station markers
- Three tabs: 2020 | 2021 | 2022 (fixed, no user filtering)
- Palette: `viridis` (colorblind-safe)

### LGCP Results
- Pre-saved PNGs from `posterior_intensity()` via `knitr::include_graphics`
- `gt` table: beta/gamma posterior means + 95% CIs per covariate
- Plain-language gloss per parameter (one sentence each)

### Domestic Violence
- Static Markdown prose
- Sections: prevalence, underreporting, institutional betrayal, Minneapolis context, local resources
- No R computation

### Methods
- Model equation (LaTeX)
- Data provenance table
- GitHub repo link
- Render timestamp

---

## Deferred Features

- User-controlled filtering (year, call type, neighborhood) via `crosstalk`
- Animated temporal playback
- Downloadable data export

---

## Files to Create

| File | Action |
|------|--------|
| `_quarto.yml` | Create |
| `dashboard.qmd` | Create |
| `03_dashboard.R` | Create |
| `output/lgcp_posterior.rds` | Written by LGCP job (not this scope) |
| `.github/workflows/render-dashboard.yml` | Create |

## Files to Modify

| File | Change |
|------|--------|
| `master.R` | None — sourcing unchanged |
| `01_data.R` | None — already writes `data/MPLS_crime.csv` |
