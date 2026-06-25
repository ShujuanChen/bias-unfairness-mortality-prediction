#!/usr/bin/env Rscript

# Generate UK Biobank participation weights against the 2011 Census 5% microdata.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(SuperLearner)
  library(glmnet)
  library(ranger)
  library(xgboost)
})

script_match <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("--file=", "", script_match[1]), mustWork = TRUE)
script_dir <- dirname(script_path)
project_root <- dirname(script_dir)

source(file.path(script_dir, "census_weights.R"))

predictors <- census_predictors()
ukb_input <- file.path(project_root, "data", "UKB", "UKB_for_harmonisation_with_census.csv")
census_input <- file.path(project_root, "data", "census", "recodev12.csv")

message("Harmonising UKB + Census for Census-based weighting")
ukb_prepared <- prepare_ukb_for_census(ukb_input, predictors)

census_prepared <- prepare_census_microdata(census_input)
# Census ids sit above the UKB range so the stacked frame has unique keys.
census_prepared$eid <- max(ukb_prepared$eid, na.rm = TRUE) + seq_len(nrow(census_prepared))
census_prepared$source <- "Census"
census_prepared$sample_weight <- 20  # 5% sample => each record ~ 20 people
census_prepared <- census_prepared[, c("eid", "source", "sample_weight", predictors), drop = FALSE]

stacked <- rbind(ukb_prepared, census_prepared)
model_input <- stacked[complete.cases(stacked[, predictors, drop = FALSE]), , drop = FALSE]

model_input_path <- census_model_input_path(project_root)
write.csv(model_input, model_input_path, row.names = FALSE)
write_note_file(
  census_model_input_note_path(project_root),
  c(
    "Census + UKB participation-model input",
    paste("Ethnicity protocol:", census_ethnicity_protocol()),
    "Sample weights: UKB = 1, Census = 20.",
    paste("UKB rows:", sum(model_input$source == "UKB")),
    paste("Census rows:", sum(model_input$source == "Census")),
    paste("Predictors:", paste(predictors, collapse = ", "))
  )
)

# Read the model input back before fitting.
model_input <- read.csv(model_input_path, stringsAsFactors = FALSE)
fit <- fit_census_superlearner(model_input, predictors)

weight_path <- census_weight_path(project_root)
write.csv(fit$ukb_weights, weight_path, row.names = FALSE)
write_note_file(
  census_weight_note_path(project_root),
  c(
    "Census + UKB SuperLearner weights",
    paste("Rows used:", nrow(model_input)),
    paste("Predictors:", paste(fit$active_predictors, collapse = ", ")),
    paste("Learner library:", paste(fit$sl_library, collapse = ", ")),
    "UKB rows carry the normalised inverse-odds weight; Census rows keep their design weight.",
    paste("Output:", weight_path)
  )
)

message("Saved Census+UKB SuperLearner weights to ", weight_path)
