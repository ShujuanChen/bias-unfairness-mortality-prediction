#!/usr/bin/env Rscript

# Build the harmonised UK Biobank analysis file aligned to the PMR covariates.

suppressPackageStartupMessages({
  library(dplyr)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_path <- if (length(script_path)) normalizePath(script_path) else ""
script_dir <- if (nzchar(script_path)) dirname(script_path) else getwd()
setwd(script_dir)

source(file.path(script_dir, "UKB_harmonise.R"))
source(file.path(script_dir, "merge_geography_UKB.R"))
source(file.path(dirname(dirname(dirname(script_dir))), "framework_config.R"))

cfg <- read_framework_config(script_path)
project_root <- file.path(cfg$framework_root, "mortality_risk_prediction")
results_dir <- file.path(project_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

ukb_input_path <- file.path(project_root, "data", "UKB", "UKB_for_harmonisation_with_PMR.csv")
output_path <- file.path(results_dir, "UKB_harmonised_with_pmr.csv")

ukb <- read.csv(ukb_input_path, stringsAsFactors = FALSE)
ukb <- ukb_prepare(ukb, verbose = TRUE)
ukb <- ukb %>% filter(!is.na(age_at_baseline), age_at_baseline >= 40, age_at_baseline <= 69)

# Drop Scotland (LSOA codes starting "S"); England/Wales only at this stage.
ukb <- ukb %>% filter(substr(LSOA11CD, 1, 1) != "S")

# Census baseline; exclude anyone who died before it so all enter follow-up alive.
t0 <- as.Date("2011-03-27")
ukb <- ukb %>% filter(is.na(date_of_death) | as.Date(date_of_death) >= t0)

ukb <- merge_geography_ukb(ukb, verbose = TRUE)
ukb <- ukb %>% filter(!is.na(RGN11CD), startsWith(RGN11CD, "E"))
validate_harmonised_ethnicity5(ukb$ethnicity5, context = "UKB harmonised file")

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
required_vars <- intersect(required_vars, names(ukb))
if (length(required_vars) == 0L) {
  stop("No phase-2 modelling variables found in the harmonised UKB dataset.", call. = FALSE)
}

ukb_out <- ukb %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(required_vars), ~ !is.na(.)))

write.csv(ukb_out, output_path, row.names = FALSE)

message("Saved harmonised UKB dataset to ", output_path)
