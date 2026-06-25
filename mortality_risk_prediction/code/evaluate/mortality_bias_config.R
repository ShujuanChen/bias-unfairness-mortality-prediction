# Config accessors for the linear_cox sensitivity evaluation.

script_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path_for_config <- if (length(script_arg)) sub("--file=", "", script_arg[1]) else ""

find_framework_config_r <- function(start_path = getwd()) {
  path <- normalizePath(start_path, mustWork = TRUE)
  if (file.exists(path) && !dir.exists(path)) {
    path <- dirname(path)
  }

  repeat {
    candidate <- file.path(path, "framework_config.R")
    if (file.exists(candidate)) {
      return(candidate)
    }

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate framework_config.R from: ", start_path, call. = FALSE)
    }
    path <- parent
  }
}

source(find_framework_config_r(if (nzchar(script_path_for_config) && file.exists(script_path_for_config)) script_path_for_config else getwd()))

read_mortality_bias_config <- function(start_path = current_script_path()) {
  read_framework_config(start_path)
}

get_mortality_risk_prediction_repo_dir <- function(cfg) {
  path_from_framework_root(cfg, "mortality_risk_prediction")
}

get_mortality_results_dir <- function(cfg) {
  results_dir <- file.path(get_mortality_risk_prediction_repo_dir(cfg), "results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  results_dir
}

get_mortality_harmonised_dataset_paths <- function(cfg) {
  results_dir <- get_mortality_results_dir(cfg)

  list(
    ukb = file.path(results_dir, "UKB_harmonised_with_pmr.csv"),
    pmr = file.path(results_dir, "PMR_harmonised_with_ukb.csv")
  )
}

get_mortality_write_model_fitting_outputs <- function(cfg) {
  isTRUE(cfg_get(cfg, c("phase2", "write_model_fitting_outputs"), default = FALSE, required = FALSE))
}

normalise_mortality_transfer_model_family <- function(family) {
  family <- tolower(trimws(as.character(family)))

  if (!nzchar(family)) {
    stop("`phase2.transfer_model_family` must be a non-empty string.", call. = FALSE)
  }

  if (family %in% c("linear_cox", "linearcox", "cox", "cox_spline")) {
    return("linear_cox")
  }

  if (family %in% c("deep_surv", "deepsurv", "deep_survival", "deepcox")) {
    stop(
      "The deep survival model is not dispatched through the R wrapper. ",
      "Run the standalone pipeline at mortality_risk_prediction/code/evaluate/deep_surv.py ",
      "(run it per cause/source/fold, then ensemble per cause).",
      call. = FALSE
    )
  }

  stop(
    "Unsupported `phase2.transfer_model_family`: ", family,
    ". Use `linear_cox` in framework_config.json. ",
    "The deep survival model runs as a standalone Python pipeline.",
    call. = FALSE
  )
}

get_mortality_transfer_model_family <- function(cfg) {
  normalise_mortality_transfer_model_family(
    cfg_get(cfg, c("phase2", "transfer_model_family"), default = "linear_cox", required = FALSE)
  )
}

normalise_mortality_outcome_key <- function(key) {
  key <- tolower(trimws(as.character(key)))
  if (!nzchar(key)) {
    stop("`phase2.outcome_key` must be a non-empty string.", call. = FALSE)
  }
  key
}

get_mortality_outcome_key <- function(cfg) {
  normalise_mortality_outcome_key(cfg_get(cfg, c("phase2", "outcome_key")))
}

get_mortality_outcome_config <- function(cfg, key = NULL) {
  if (is.null(key) || !nzchar(key)) {
    key <- get_mortality_outcome_key(cfg)
  } else {
    key <- normalise_mortality_outcome_key(key)
  }

  outcome_cfg <- cfg_get(cfg, c("phase2", "outcomes", key))
  type <- cfg_get(outcome_cfg, c("type"))
  label <- cfg_get(outcome_cfg, c("label"))
  match_scope <- cfg_get(outcome_cfg, c("match_scope"), default = "underlying", required = FALSE)

  icd10 <- cfg_get(outcome_cfg, c("icd10"), default = list(), required = FALSE)
  icd10 <- unname(as.character(unlist(icd10)))

  if (!type %in% c("all_cause", "cause_specific")) {
    stop(
      "Unsupported mortality outcome type for `phase2.outcomes.", key, "`: ", type,
      call. = FALSE
    )
  }

  if (identical(type, "cause_specific") && !length(icd10)) {
    stop(
      "Cause-specific outcome `", key,
      "` must define at least one ICD-10 range in framework_config.json.",
      call. = FALSE
    )
  }

  list(
    key = key,
    label = label,
    type = type,
    match_scope = "underlying",
    icd10 = icd10
  )
}

get_mortality_weight_source_keys <- function(cfg) {
  keys <- cfg_get(cfg, c("phase2", "weight_source_keys"), default = NULL, required = FALSE)

  if (!is.null(keys)) {
    keys <- unname(unlist(keys))
    keys <- keys[nzchar(keys)]
    if (length(keys) > 0L) {
      return(keys)
    }
  }

  key <- cfg_get(cfg, c("phase2", "weight_source_key"))
  key <- unname(unlist(key))
  key[nzchar(key)]
}

get_mortality_weight_source <- function(cfg, key = NULL) {
  if (is.null(key) || !nzchar(key)) {
    key <- cfg_get(cfg, c("phase2", "weight_source_key"))
  }

  source_cfg <- cfg_get(cfg, c("phase2", "weight_sources", key))
  list(
    key = key,
    label = cfg_get(source_cfg, c("label")),
    path = path_from_framework_root(cfg, cfg_get(source_cfg, c("path"))),
    column = cfg_get(source_cfg, c("column"))
  )
}

get_mortality_evaluation_dir <- function(cfg) {
  eval_dir <- file.path(get_mortality_results_dir(cfg), "evaluation")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  eval_dir
}

get_mortality_evaluation_model_dir <- function(cfg) {
  model_dir <- file.path(
    get_mortality_evaluation_dir(cfg),
    get_mortality_transfer_model_family(cfg)
  )
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  model_dir
}

get_mortality_outcome_dir <- function(cfg, outcome_key = NULL) {
  outcome_cfg <- get_mortality_outcome_config(cfg, key = outcome_key)
  outcome_dir <- file.path(get_mortality_evaluation_model_dir(cfg), outcome_cfg$key)
  if (identical(outcome_cfg$type, "cause_specific")) {
    outcome_dir <- file.path(outcome_dir, outcome_cfg$match_scope)
  }
  dir.create(outcome_dir, recursive = TRUE, showWarnings = FALSE)
  outcome_dir
}

normalise_mortality_run_slug_component <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

get_mortality_fit_dir <- function(cfg, fit_key, outcome_key = NULL, weight_source_key = NULL) {
  fit_key <- normalise_mortality_run_slug_component(fit_key)
  outcome_dir <- get_mortality_outcome_dir(cfg, outcome_key = outcome_key)

  if (identical(fit_key, "ukbw_fit")) {
    if (is.null(weight_source_key) || !nzchar(weight_source_key)) {
      stop("`weight_source_key` is required for `ukbw_fit`.", call. = FALSE)
    }
    fit_dir <- file.path(outcome_dir, fit_key, normalise_mortality_run_slug_component(weight_source_key))
  } else {
    fit_dir <- file.path(outcome_dir, fit_key)
  }

  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  fit_dir
}

get_mortality_summary_dir <- function(cfg, outcome_key = NULL) {
  summary_dir <- file.path(get_mortality_outcome_dir(cfg, outcome_key = outcome_key), "summaries")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  summary_dir
}
