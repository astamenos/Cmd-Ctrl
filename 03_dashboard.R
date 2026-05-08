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

# ── KDE heatmap raster ────────────────────────────────────────────────────────
# Returns a WGS84 SpatRaster for use with leaflet::addRasterImage().
# crime_df:  data frame with columns x, y (km, UTM Zone 15N)
# window:    spatstat owin in matching km coordinates
# sigma_km:  KDE bandwidth in km
# offense:   "all" | "domestic" | "nondomestic"

kde_crime_raster <- function(crime_df, window, sigma_km = 0.5,
                              offense = c("all", "domestic", "nondomestic")) {
  offense <- match.arg(offense)
  pts <- switch(offense,
    "domestic"    = crime_df[crime_df$Y == 1, ],
    "nondomestic" = crime_df[crime_df$Y == 0, ],
    crime_df
  )

  ppp_obj <- spatstat.geom::ppp(x = pts$x, y = pts$y, window = window,
                                 check = FALSE)
  dens    <- spatstat.explore::density.ppp(ppp_obj, sigma = sigma_km,
                                           positive = TRUE)

  # spatstat im: v[row, col] where row 1 = bottom y (ascending).
  # terra rast fills top-to-bottom, left-to-right.
  mat_flipped <- dens$v[nrow(dens$v):1, ]  # flip y so row 1 = top

  r <- terra::rast(
    nrows = nrow(mat_flipped),
    ncols = ncol(mat_flipped),
    xmin  = dens$xrange[1], xmax = dens$xrange[2],
    ymin  = dens$yrange[1], ymax = dens$yrange[2],
    crs   = "+proj=utm +zone=15 +datum=WGS84 +units=km"
  )
  terra::values(r) <- as.vector(t(mat_flipped))
  terra::project(r, "EPSG:4326")
}
