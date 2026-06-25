#!/usr/bin/env Rscript

# Entry point for PMR harmonisation: clean raw PMR, attach geography, restrict to the cohort, write the modelling-ready file.

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(splines)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_path <- if (length(script_path)) normalizePath(script_path) else ""
script_dir <- if (nzchar(script_path)) dirname(script_path) else getwd()
setwd(script_dir)

source(file.path(script_dir, "PMR_harmonise.R"))
source(file.path(script_dir, "merge_geography.R"))
source(file.path(script_dir, "row_harmonisation_PMR.R"))
source(file.path(dirname(dirname(dirname(script_dir))), "framework_config.R"))

cfg <- read_framework_config(script_path)
project_root <- file.path(cfg$framework_root, "mortality_risk_prediction")
results_dir <- file.path(project_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

pmr_input_path <- file.path(project_root, "data", "PMR", "PMR.csv")
output_path <- file.path(results_dir, "PMR_harmonised_with_ukb.csv")

pmr <- read.csv(pmr_input_path, stringsAsFactors = FALSE) %>%
  pmr_prepare(verbose = TRUE)

pmr <- merge_geography(pmr, verbose = TRUE)

pmr <- pmr_filter_age_ukb_eligible(
  pmr,
  age_col = "age_at_baseline",
  verbose = TRUE
)
pmr <- pmr_filter_england_only(
  pmr,
  col = "RGN11CD",
  verbose = TRUE
)

pmr <- pmr %>% select(-any_of(c("LSOA11CD_join", "uresindpuk11_census")))
validate_harmonised_ethnicity5(pmr$ethnicity5, context = "PMR harmonised file")

required_vars <- c(
  "age_at_baseline",
  "sex",
  "ethnicity5",
  "tenure",
  "household_size",
  "econstatus",
  "education",
  "ruralurban",
  "health",
  "disability",
  "imd_decile",
  "RGN11CD"
)
required_vars <- intersect(required_vars, names(pmr))
if (length(required_vars) == 0L) {
  stop("No phase-2 modelling variables found in the harmonised PMR dataset.", call. = FALSE)
}

pmr_out <- pmr %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(required_vars), ~ !is.na(.)))

write.csv(pmr_out, output_path, row.names = FALSE)

message("Saved harmonised PMR dataset to ", output_path)
