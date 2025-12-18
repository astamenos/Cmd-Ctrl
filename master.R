# Libraries
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

# For reproducibility
set.seed(51225)

# Common coordinate reference system to use for points and shapefile
crs_use <- "+proj=utm +zone=15 +datum=WGS84"
d_jitt <- 0.5
years <- as.character(seq(2020, 2022, 1))
years_fac <- factor(years)

# Call types
call_types <- c("no_info", "weapon", "unknown_trouble",
                "unknown_phone", "welfare", "emo_dist")

# Call types that may indicate domestic
domestic_terms <- c(
  "domestic",
  "Person with a Weapon",
  "Unknown Trouble",
  "Unknown Wireless/Cell Phone",
  "Check the Welfare",
  "Emotionally Disturb Person"
)

# Sourcing helper functions
source('utils.R')

# Sourcing data processing script
source('data_processing.R')

# Sourcing analysis script
source('spatiotemp_mark_LGCP.R')

