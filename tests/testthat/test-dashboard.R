library(testthat)
library(dplyr)
library(here)

# Source the file under test (functions exist once sourced)
source(here::here("03_dashboard.R"))

# ── load_crime_data ───────────────────────────────────────────────────────────

test_that("load_crime_data returns tibble with required columns", {
  df <- load_crime_data("data/MPLS_crime.csv")
  expect_s3_class(df, "tbl_df")
  required <- c("year", "Offense", "Y", "domestic_flag", "x", "y",
                "pop_dens", "pct_pov", "Precinct", "Ward", "occurred_date")
  expect_true(all(required %in% names(df)))
})

test_that("load_crime_data filters to expected years", {
  df <- load_crime_data("data/MPLS_crime.csv")
  expect_setequal(unique(df$year), c("2020", "2021", "2022"))
})

test_that("load_crime_data has no zero-N rows", {
  df <- load_crime_data("data/MPLS_crime.csv")
  expect_true(all(!is.na(df$x) & !is.na(df$y)))
})

# ── compute_annual_trend ──────────────────────────────────────────────────────

test_that("compute_annual_trend returns one row per year × offense", {
  df <- load_crime_data("data/MPLS_crime.csv")
  trend <- compute_annual_trend(df)
  expect_equal(nrow(trend), 6L)  # 3 years × 2 offense types
  expect_true(all(c("year", "Offense", "n", "pct_domestic") %in% names(trend)))
})

test_that("compute_annual_trend pct_domestic is in [0, 100]", {
  df <- load_crime_data("data/MPLS_crime.csv")
  trend <- compute_annual_trend(df)
  expect_true(all(trend$pct_domestic >= 0 & trend$pct_domestic <= 100))
})

# ── load_lgcp_posterior ───────────────────────────────────────────────────────

test_that("load_lgcp_posterior returns list with required elements", {
  skip_if_not(file.exists("output/lgcp_posterior.rds"),
              "LGCP artifact not present — run 02_model.R first")
  post <- load_lgcp_posterior("output/lgcp_posterior.rds")
  required <- c("samples_comb", "ips_sf", "model_summary",
                "p_int", "p_mark", "n_years",
                "covariates_int", "covariates_mark")
  expect_true(all(required %in% names(post)))
})

test_that("load_lgcp_posterior samples_comb has expected beta/gamma columns", {
  skip_if_not(file.exists("output/lgcp_posterior.rds"))
  post <- load_lgcp_posterior("output/lgcp_posterior.rds")
  cols <- colnames(post$samples_comb)
  expect_true("beta[1]"  %in% cols)
  expect_true("beta[2]"  %in% cols)
  expect_true("gamma[1]" %in% cols)
  expect_true("gamma[3]" %in% cols)
})

# ── kde_crime_raster ──────────────────────────────────────────────────────────

test_that("kde_crime_raster returns a SpatRaster in WGS84", {
  df    <- load_crime_data("data/MPLS_crime.csv")
  mpls  <- sf::st_read(here::here("data/Minneapolis_Neighborhoods/Minneapolis_Neighborhoods.shp"),
                       quiet = TRUE) |>
    sf::st_transform("+proj=utm +zone=15 +datum=WGS84") |>
    dplyr::mutate(geometry = geometry / 1000) |>
    sf::st_set_crs("+proj=utm +zone=15 +datum=WGS84") |>
    sf::st_union()
  win   <- spatstat.geom::as.owin(mpls)
  r     <- kde_crime_raster(df, win, sigma_km = 0.5, offense = "all")
  expect_s4_class(r, "SpatRaster")
  expect_equal(as.character(terra::crs(r, describe = TRUE)$authority),
               "EPSG")
  expect_equal(terra::crs(r, describe = TRUE)$code, "4326")
})

test_that("kde_crime_raster domestic-only has fewer non-NA cells than all", {
  df   <- load_crime_data("data/MPLS_crime.csv") |> dplyr::filter(year == "2021")
  mpls <- sf::st_read(here::here("data/Minneapolis_Neighborhoods/Minneapolis_Neighborhoods.shp"),
                      quiet = TRUE) |>
    sf::st_transform("+proj=utm +zone=15 +datum=WGS84") |>
    dplyr::mutate(geometry = geometry / 1000) |>
    sf::st_set_crs("+proj=utm +zone=15 +datum=WGS84") |>
    sf::st_union()
  win  <- spatstat.geom::as.owin(mpls)
  r_all <- kde_crime_raster(df, win, sigma_km = 0.5, offense = "all")
  r_dom <- kde_crime_raster(df, win, sigma_km = 0.5, offense = "domestic")
  expect_lte(sum(!is.na(terra::values(r_dom))),
             sum(!is.na(terra::values(r_all))))
})

# ── format_posterior_table ────────────────────────────────────────────────────

test_that("format_posterior_table returns gt table with correct rows", {
  skip_if_not(file.exists("output/lgcp_posterior.rds"))
  post <- load_lgcp_posterior("output/lgcp_posterior.rds")
  tbl  <- format_posterior_table(post$samples_comb, post$p_int, post$p_mark)
  expect_s3_class(tbl, "gt_tbl")
  # p_int + p_mark + 3 variance components (sigma, tau[1], tau[2])
  n_expected <- post$p_int + post$p_mark + 3L
  expect_equal(nrow(tbl$`_data`), n_expected)
})
