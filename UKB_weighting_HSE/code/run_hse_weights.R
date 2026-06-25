#!/usr/bin/env Rscript

# Generate UK Biobank participation weights against the Health Survey for England.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(haven)
  library(plyr)
  library(fastDummies)
  library(Matrix)
  library(glmnet)
  library(SuperLearner)
  library(ranger)
  library(xgboost)
})

script_match <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("--file=", "", script_match[1]), mustWork = TRUE)
script_dir <- dirname(script_path)
project_root <- dirname(script_dir)

source(file.path(script_dir, "hse_weights.R"))

model <- commandArgs(trailingOnly = TRUE)
model <- if (length(model)) model[1] else "superlearner"
if (!model %in% c("superlearner", "lassologit")) {
  stop("Unknown model '", model, "'. Use 'superlearner' or 'lassologit'.")
}

cfg <- load_hse_config(project_root)

fit <- if (identical(model, "superlearner")) {
  fit_superlearner_weights(project_root, cfg)
} else {
  fit_lasso_weights(project_root, cfg)
}

model_input_path <- hse_model_input_path(project_root)
weight_path <- hse_weight_path(project_root, model)

data.table::fwrite(fit$model_input, model_input_path, na = "NA")
data.table::fwrite(fit$ukb_weights, weight_path, na = "NA")

write_note_file(
  hse_model_input_note_path(project_root),
  c(
    "HSE + UKB participation-model input",
    paste("Config:", file.path(project_root, "config", "hse_weighting.json")),
    paste("Rows:", nrow(fit$model_input)),
    paste("HSE rows:", sum(fit$model_input$source == "HSE")),
    paste("UKB rows:", sum(fit$model_input$source == "UKB")),
    paste("Predictors:", paste(fit$prediction_labels, collapse = ", "))
  )
)

write_note_file(
  hse_weight_note_path(project_root, model),
  c(
    paste0("HSE participation weights (", model, ")"),
    paste("Rows used:", nrow(fit$model_input)),
    paste("Predictors:", paste(fit$prediction_labels, collapse = ", ")),
    "UKB rows carry the normalised inverse-odds weight; HSE rows keep their sample weight.",
    paste("Output:", weight_path)
  )
)

message("Saved HSE ", model, " weights to ", weight_path)
