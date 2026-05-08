# refresh_data.R
# Called by GitHub Actions nightly to update data/MPLS_crime.csv.
# Mirrors the globals from master.R, then sources 00_utils.R and 01_data.R only.
# Never sources 02_model.R — LGCP is run separately and artifacts are checked in.

library(scales)
library(patchwork)
library(tigris)
library(tidycensus)
require(tigris)
require(spdep)
require(osmdata)
library(spdep)
library(ggnewscale)
library(parallel)
suppressPackageStartupMessages(library(spatstat))
library(sf)
library(tidyverse)
library(nimble)
library(lme4)
require(FNN)
require(gstat)

set.seed(51225)

crs_use <- "+proj=utm +zone=15 +datum=WGS84"
d_jitt  <- 0.5
years   <- as.character(seq(2020, 2022, 1))
years_fac <- factor(years)

call_types <- c("no_info", "weapon", "unknown_trouble",
                "unknown_phone", "welfare", "emo_dist")

domestic_terms <- c(
  "domestic",
  "Person with a Weapon",
  "Unknown Trouble",
  "Unknown Wireless/Cell Phone",
  "Check the Welfare",
  "Emotionally Disturb Person"
)

source("00_utils.R")
source("01_data.R")

cat("data/MPLS_crime.csv updated:", format(Sys.time(), "%Y-%m-%d %H:%M UTC"), "\n")
