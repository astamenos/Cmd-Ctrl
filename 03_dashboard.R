library(tidyverse)
library(sf)
library(leaflet)
#library(leaflet.extras)
library(plotly)
library(gt)
library(spatstat.geom)
library(spatstat.explore)
library(terra)
library(here)

# ── Crime data ────────────────────────────────────────────────────────────────

load_crime_data <- function(path = "data/MPLS_crime.csv") {
  required_cols <- c("year", "Offense", "Y", "domestic_flag",
                     "x", "y", "pop_dens", "pct_pov",
                     "Precinct", "Ward", "occurred_date")
  df <- read_csv(here::here(path), show_col_types = FALSE) |>
    mutate(year = as.character(year))
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))
  df
}

# ── Annual trend ──────────────────────────────────────────────────────────────

compute_annual_trend <- function(crime_df) {
  crime_df |>
    group_by(year, Offense) |>
    summarise(
      n            = n(),
      pct_domestic = 100 * mean(Y == 1, na.rm = TRUE),
      .groups      = "drop"
    )
}

# ── LGCP posterior ────────────────────────────────────────────────────────────

load_lgcp_posterior <- function(path = "output/lgcp_posterior.rds") {
  full_path <- here::here(path)
  if (!file.exists(full_path)) stop("LGCP artifact not found at: ", full_path)
  post <- readRDS(full_path)
  required <- c("samples_comb", "ips_sf", "model_summary",
                "p_int", "p_mark", "n_years",
                "covariates_int", "covariates_mark")
  missing <- setdiff(required, names(post))
  if (length(missing) > 0) stop("RDS missing elements: ", paste(missing, collapse = ", "))
  post
}
