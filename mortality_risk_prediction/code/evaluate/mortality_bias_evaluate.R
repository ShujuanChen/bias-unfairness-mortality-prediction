#!/usr/bin/env Rscript
# linear_cox sensitivity evaluation: Cox vs deep survival comparison (Fig S13).

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_path <- if (length(script_path)) normalizePath(script_path) else ""
script_dir <- if (nzchar(script_path)) dirname(script_path) else getwd()
setwd(script_dir)

source(file.path("..", "PMR", "clean_outcomes.R"))
source(file.path(script_dir, "models.R"))
source(file.path(script_dir, "mortality_bias_config.R"))

t0 <- as.Date("2011-03-27")
t_admin_end <- as.Date("2023-02-15")
horizons_years <- c(5, 10)
horizons_days <- as.integer(round(365.25 * horizons_years))

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

safe_weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & (w > 0)
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

extract_rhs_covariates <- function() {
  c(
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
    "imd_decile"
  )
}

extract_categorical_rhs_covariates <- function() {
  setdiff(extract_rhs_covariates(), "age_at_baseline")
}

filter_to_common_support <- function(ukb_df, pmr_df, vars, max_iter = 10L) {
  support_levels <- list()

  for (iter in seq_len(max_iter)) {
    changed <- FALSE

    for (v in vars) {
      ukb_vals <- sort(unique(as.character(stats::na.omit(ukb_df[[v]]))))
      pmr_vals <- sort(unique(as.character(stats::na.omit(pmr_df[[v]]))))
      common_vals <- intersect(ukb_vals, pmr_vals)

      if (length(common_vals) == 0L) {
        stop("No common support for variable: ", v)
      }

      ukb_keep <- !is.na(ukb_df[[v]]) & as.character(ukb_df[[v]]) %in% common_vals
      pmr_keep <- !is.na(pmr_df[[v]]) & as.character(pmr_df[[v]]) %in% common_vals

      if (sum(!ukb_keep) > 0L || sum(!pmr_keep) > 0L) {
        changed <- TRUE
        ukb_df <- ukb_df[ukb_keep, , drop = FALSE]
        pmr_df <- pmr_df[pmr_keep, , drop = FALSE]
      }

      support_levels[[v]] <- common_vals
    }

    if (!changed) {
      return(list(ukb = ukb_df, pmr = pmr_df, support_levels = support_levels))
    }
  }

  stop("Common-support filtering did not converge within ", max_iter, " iterations.")
}

make_age_band_10y <- function(age) {
  age_num <- suppressWarnings(as.numeric(age))
  out <- rep(NA_character_, length(age_num))

  ok <- !is.na(age_num)
  lower <- floor(age_num[ok] / 10) * 10
  upper <- lower + 9
  out[ok] <- paste0(lower, "-", upper)

  factor(out, levels = sort(unique(out[ok])))
}

prepare_horizon_rows <- function(df, horizon_days) {
  tt <- as.numeric(horizon_days)
  keep <- !(df$status == 0 & df$time_days < tt)
  d <- df[keep, , drop = FALSE]
  d$event_h <- as.numeric(d$status == 1 & d$time_days <= tt)
  d
}

summarise_weighted_overall <- function(df, weight_col, value_cols) {
  w <- as.numeric(df[[weight_col]])
  out <- lapply(value_cols, function(col) safe_weighted_mean(as.numeric(df[[col]]), w))
  as.data.frame(out, stringsAsFactors = FALSE)
}

summarise_weighted_by_group <- function(df, strata_var, weight_col, value_cols) {
  strata_level <- as.character(df[[strata_var]])
  keep <- !is.na(strata_level) & strata_level != ""

  if (!any(keep)) {
    return(data.frame(strata_level = character(0), stringsAsFactors = FALSE))
  }

  strata_level <- strata_level[keep]
  w <- as.numeric(df[[weight_col]][keep])
  strata_index <- split(seq_along(strata_level), strata_level)

  out <- data.frame(
    strata_level = names(strata_index),
    stringsAsFactors = FALSE
  )

  for (out_name in names(value_cols)) {
    x <- as.numeric(df[[value_cols[[out_name]]]][keep])
    grouped_mean <- vapply(
      strata_index,
      function(idx) safe_weighted_mean(x[idx], w[idx]),
      numeric(1)
    )
    out[[out_name]] <- as.numeric(grouped_mean[out$strata_level])
  }

  out
}

extract_model_predictions <- function(pred_df, id_col, model_name, prefix) {
  out <- pred_df[pred_df$model == model_name, c(id_col, "risk_5y", "risk_10y"), drop = FALSE]
  names(out) <- c(id_col, paste0(prefix, "_5y"), paste0(prefix, "_10y"))
  out
}

attach_prediction_pattern_id <- function(df, covars) {
  df_with_id <- df %>%
    group_by(across(all_of(covars))) %>%
    mutate(pred_pattern_id = dplyr::cur_group_id()) %>%
    ungroup()

  pred_input <- df_with_id %>%
    select(pred_pattern_id, all_of(covars)) %>%
    distinct() %>%
    arrange(pred_pattern_id)

  list(
    df = df_with_id,
    pred_input = pred_input
  )
}

build_risk_summary <- function(pmr_rows, ukb_rows) {
  strata_vars <- c(
    "age_band_10y",
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

  res_list <- list()

  for (i in seq_along(horizons_years)) {
    hz_years <- horizons_years[i]
    hz_days <- horizons_days[i]
    suffix <- if (hz_years == 5) "5y" else "10y"

    pmr_h <- prepare_horizon_rows(pmr_rows, hz_days)
    ukb_h <- prepare_horizon_rows(ukb_rows, hz_days)

    pmr_value_cols <- c(
      pmr_obs = "event_h",
      pmr_pred = paste0("pmr_pred_", suffix),
      ukb_pred = paste0("ukb_pred_", suffix),
      ukbw_pred = paste0("ukbw_pred_", suffix)
    )

    overall_joined <- data.frame(
      strata_variable = "Overall",
      strata_level = "All",
      horizon_years = hz_years,
      horizon_days = hz_days,
      stringsAsFactors = FALSE
    ) %>%
      bind_cols(summarise_weighted_overall(pmr_h, "w", pmr_value_cols)) %>%
      bind_cols(summarise_weighted_overall(ukb_h, "w_unweighted", c(ukb_obs = "event_h"))) %>%
      bind_cols(summarise_weighted_overall(ukb_h, "ukb_eval_weight", c(ukbw_obs = "event_h"))) %>%
      mutate(
        ref_min_pmr_obs = pmr_pred - pmr_obs,
        ref_min_ukb_pred = pmr_pred - ukb_pred,
        ref_min_ukbw_pred = pmr_pred - ukbw_pred
      )

    res_list[[paste0("Overall__", hz_years)]] <- overall_joined

    for (v in strata_vars) {
      pmr_sum <- summarise_weighted_by_group(pmr_h, v, "w", pmr_value_cols)
      ukb_obs <- summarise_weighted_by_group(ukb_h, v, "w_unweighted", c(ukb_obs = "event_h"))
      ukbw_obs <- summarise_weighted_by_group(ukb_h, v, "ukb_eval_weight", c(ukbw_obs = "event_h"))

      joined <- Reduce(
        function(x, y) merge(x, y, by = "strata_level", all = TRUE),
        list(pmr_sum, ukb_obs, ukbw_obs)
      )

      joined$strata_variable <- v
      joined$horizon_years <- hz_years
      joined$horizon_days <- hz_days
      joined$ref_min_pmr_obs <- joined$pmr_pred - joined$pmr_obs
      joined$ref_min_ukb_pred <- joined$pmr_pred - joined$ukb_pred
      joined$ref_min_ukbw_pred <- joined$pmr_pred - joined$ukbw_pred

      res_list[[paste(v, hz_years, sep = "__")]] <- joined
    }
  }

  bind_rows(res_list) %>%
    arrange(
      dplyr::if_else(strata_variable == "Overall", 0L, 1L),
      strata_variable,
      horizon_years,
      strata_level
    ) %>%
    distinct() %>%
    select(
      strata_variable,
      strata_level,
      horizon_years,
      horizon_days,
      pmr_obs,
      ukb_obs,
      ukbw_obs,
      pmr_pred,
      ukb_pred,
      ukbw_pred,
      ref_min_pmr_obs,
      ref_min_ukb_pred,
      ref_min_ukbw_pred
    )
}

get_model_suffix <- function(transfer_model_family) {
  switch(
    transfer_model_family,
    linear_cox = "linear_cox"
  )
}

make_model_key <- function(transfer_model_family, dataset_name) {
  paste(dataset_name, get_model_suffix(transfer_model_family), sep = "__")
}

make_model_fns <- function(transfer_model_family) {
  switch(
    transfer_model_family,
    linear_cox = list(
      linear_cox = function(df, w_col) {
        train_cox_spline(df, time_col = "time_days", status_col = "status", w_col = w_col)
      }
    )
  )
}

write_model_object_rds <- function(component_dir, model_payload, enabled = TRUE) {
  if (!isTRUE(enabled)) return(invisible(NULL))
  saveRDS(model_payload, file.path(component_dir, "model_object.rds"))
}

make_model_payload <- function(transfer_model_family, fit_component, fit_obj, metadata = list()) {
  list(
    model_family = transfer_model_family,
    fit_component = fit_component,
    model_name = fit_obj$name %||% transfer_model_family,
    fit = fit_obj$fit,
    metadata = metadata
  )
}

# This R driver handles only the linear_cox family; the deep survival model is
# the standalone Python pipeline at mortality_risk_prediction/code/evaluate/deep_surv.py.

prepare_base_evaluation_data <- function(cfg, dataset_paths, outcome_cfg, transfer_model_family) {
  stop_if_missing(dataset_paths$ukb)
  stop_if_missing(dataset_paths$pmr)

  ukb <- read.csv(dataset_paths$ukb, stringsAsFactors = FALSE) %>%
    mutate(participant_id = as.numeric(participant_id))

  pmr <- read.csv(dataset_paths$pmr, stringsAsFactors = FALSE) %>%
    mutate(pmr_row_id = dplyr::row_number())

  ukb_rows_raw <- nrow(ukb)
  pmr_rows_raw <- nrow(pmr)

  ukb <- ukb %>%
    mutate(
      hse_age_ok = !is.na(age_at_baseline) & age_at_baseline >= 40 & age_at_baseline <= 69,
      hse_england_ok = !is.na(RGN11CD) & grepl("^E", RGN11CD)
    ) %>%
    filter(hse_age_ok, hse_england_ok)

  pmr <- pmr %>%
    mutate(
      hse_age_ok = !is.na(age_at_baseline) & age_at_baseline >= 40 & age_at_baseline <= 69,
      hse_england_ok = !is.na(RGN11CD) & grepl("^E", RGN11CD)
    ) %>%
    filter(hse_age_ok, hse_england_ok)

  ukb_rows_after_cohort_filter <- nrow(ukb)
  pmr_rows_after_cohort_filter <- nrow(pmr)

  ukb <- ukb %>%
    mutate(date_of_death = as.Date(date_of_death)) %>%
    filter(is.na(date_of_death) | date_of_death >= t0) %>%
    bind_cols(build_survival_outcome_from_config(
      .,
      death_date_col = "date_of_death",
      outcome_cfg = outcome_cfg,
      t0 = t0,
      t_admin_end = t_admin_end
    )) %>%
    filter(!is.na(time_days) & time_days >= 0)

  pmr <- pmr %>%
    mutate(dod_deaths = as.Date(dod_deaths)) %>%
    filter(is.na(dod_deaths) | dod_deaths >= t0) %>%
    bind_cols(build_survival_outcome_from_config(
      .,
      death_date_col = "dod_deaths",
      outcome_cfg = outcome_cfg,
      t0 = t0,
      t_admin_end = t_admin_end
    )) %>%
    mutate(w = as.numeric(w)) %>%
    filter(!is.na(time_days) & time_days >= 0)

  ukb_rows_after_survival <- nrow(ukb)
  pmr_rows_after_survival <- nrow(pmr)

  ukbX <- coerce_harmonised_covariates(ukb)
  pmrX <- coerce_harmonised_covariates(pmr)

  rhs_covs <- extract_rhs_covariates()

  ukb_complete <- complete.cases(ukbX[, rhs_covs, drop = FALSE], ukbX$time_days, ukbX$status)
  pmr_complete <- complete.cases(pmrX[, rhs_covs, drop = FALSE], pmrX$time_days, pmrX$status, pmrX$w)

  ukbX <- ukbX[ukb_complete, , drop = FALSE]
  pmrX <- pmrX[pmr_complete, , drop = FALSE]

  ukb_rows_after_modelling_filter <- nrow(ukbX)
  pmr_rows_after_modelling_filter <- nrow(pmrX)

  support_filtered <- filter_to_common_support(
    ukb_df = ukbX,
    pmr_df = pmrX,
    vars = extract_categorical_rhs_covariates()
  )

  ukbX <- droplevels(support_filtered$ukb)
  pmrX <- droplevels(support_filtered$pmr)

  ukb_rows_after_common_support <- nrow(ukbX)
  pmr_rows_after_common_support <- nrow(pmrX)

  prediction_covars <- extract_rhs_covariates()
  pmr_prediction_layout <- attach_prediction_pattern_id(pmrX, prediction_covars)
  pmrX <- pmr_prediction_layout$df

  common_ethnicity_levels <- support_filtered$support_levels$ethnicity5
  outcome_icd10 <- if (length(outcome_cfg$icd10)) paste(outcome_cfg$icd10, collapse = ", ") else "none"

  summary_lines <- c(
    paste0("transfer_model_family=", transfer_model_family),
    paste0("outcome_key=", outcome_cfg$key),
    paste0("outcome_label=", outcome_cfg$label),
    paste0("outcome_type=", outcome_cfg$type),
    paste0("outcome_match_scope=", outcome_cfg$match_scope),
    paste0("outcome_icd10=", outcome_icd10),
    "Cohort restriction: England and age 40-69, then non-missing modelling fields and common categorical support.",
    "Cause-specific mortality handling: non-target deaths are treated as right-censoring at the death date.",
    paste0("Common-support ethnicity levels retained: ", paste(common_ethnicity_levels, collapse = ", ")),
    "",
    paste0("UKB rows in harmonised input: ", ukb_rows_raw),
    paste0("UKB rows after England/age cohort filter: ", ukb_rows_after_cohort_filter),
    paste0("UKB rows after survival filter: ", ukb_rows_after_survival),
    paste0("UKB target-cause deaths after survival filter: ", sum(ukb$target_death_in_followup == 1L, na.rm = TRUE)),
    paste0("UKB rows after non-missing modelling-fields filter: ", ukb_rows_after_modelling_filter),
    paste0("UKB target-cause deaths after non-missing modelling-fields filter: ", sum(ukbX$target_death_in_followup == 1L, na.rm = TRUE)),
    paste0("UKB rows after common categorical support filter: ", ukb_rows_after_common_support),
    paste0("UKB target-cause deaths after common categorical support filter: ", sum(ukbX$target_death_in_followup == 1L, na.rm = TRUE)),
    "",
    paste0("PMR rows in harmonised input: ", pmr_rows_raw),
    paste0("PMR rows after England/age cohort filter: ", pmr_rows_after_cohort_filter),
    paste0("PMR rows after survival filter: ", pmr_rows_after_survival),
    paste0("PMR target-cause deaths after survival filter: ", sum(pmr$target_death_in_followup == 1L, na.rm = TRUE)),
    paste0("PMR rows after non-missing modelling-fields filter: ", pmr_rows_after_modelling_filter),
    paste0("PMR target-cause deaths after non-missing modelling-fields filter: ", sum(pmrX$target_death_in_followup == 1L, na.rm = TRUE)),
    paste0("PMR rows after common categorical support filter: ", pmr_rows_after_common_support),
    paste0("PMR target-cause deaths after common categorical support filter: ", sum(pmrX$target_death_in_followup == 1L, na.rm = TRUE))
  )

  list(
    ukbX = ukbX,
    pmrX = pmrX,
    pmr_pred_input = pmr_prediction_layout$pred_input,
    prediction_covars = prediction_covars,
    shared_summary_lines = summary_lines
  )
}

write_component_outputs <- function(component_dir, cohort_lines) {
  writeLines(cohort_lines, file.path(component_dir, "cohort_filter_summary.txt"))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

run_evaluation_workflow <- function(cfg, dataset_paths, weight_sources) {
  outcome_cfg <- get_mortality_outcome_config(cfg)
  transfer_model_family <- get_mortality_transfer_model_family(cfg)
  write_model_fitting_outputs <- get_mortality_write_model_fitting_outputs(cfg)
  model_suffix <- get_model_suffix(transfer_model_family)

  log_step(
    "Evaluation config: transfer_model_family=", transfer_model_family,
    "; outcome_key=", outcome_cfg$key,
    "; outcome_label=", outcome_cfg$label,
    "; outcome_type=", outcome_cfg$type,
    "; match_scope=", outcome_cfg$match_scope,
    "; weight_source_keys=", paste(vapply(weight_sources, `[[`, character(1), "key"), collapse = ","),
    "; write_model_fitting_outputs=", write_model_fitting_outputs
  )

  if (identical(transfer_model_family, "deep_surv")) {
    stop(
      "DeepSurv is run via the standalone Python pipeline at mortality_risk_prediction/code/evaluate/deep_surv.py ",
      "(per cause/source/fold, then ensemble per cause); this R wrapper handles only the linear_cox family.",
      call. = FALSE
    )
  }

  base <- prepare_base_evaluation_data(
    cfg = cfg,
    dataset_paths = dataset_paths,
    outcome_cfg = outcome_cfg,
    transfer_model_family = transfer_model_family
  )

  model_fns <- make_model_fns(transfer_model_family)
  model_key <- function(dataset_name) make_model_key(transfer_model_family, dataset_name)
  prediction_covars <- base$prediction_covars

  shared_train_sets <- list(
    UKB = list(df = base$ukbX, w_col = NULL),
    PMR = list(df = base$pmrX, w_col = "w")
  )

  shared_fits <- list()
  shared_pred_pattern_tables <- list()

  for (ds_name in names(shared_train_sets)) {
    ds <- shared_train_sets[[ds_name]]
    for (m_name in names(model_fns)) {
      key <- paste(ds_name, m_name, sep = "__")
      log_step("Training start: ", key)
      shared_fits[[key]] <- model_fns[[m_name]](ds$df, ds$w_col)
      log_step("Training done: ", key)

      log_step("Prediction start on PMR patterns: ", key)
      risks <- shared_fits[[key]]$predict_risk(
        base$pmr_pred_input[, prediction_covars, drop = FALSE],
        horizons_days = horizons_days
      )
      log_step("Prediction done on PMR patterns: ", key)

      shared_pred_pattern_tables[[key]] <- base$pmr_pred_input %>%
        transmute(
          pred_pattern_id,
          risk_5y = risks[[paste0("risk_", horizons_days[1], "d")]],
          risk_10y = risks[[paste0("risk_", horizons_days[2], "d")]]
        )
    }
  }

  shared_pred_needed <- Reduce(
    function(x, y) merge(x, y, by = "pred_pattern_id", all = TRUE),
    list(
      extract_model_predictions(transform(shared_pred_pattern_tables[[model_key("PMR")]], model = model_key("PMR")), "pred_pattern_id", model_key("PMR"), "pmr_pred"),
      extract_model_predictions(transform(shared_pred_pattern_tables[[model_key("UKB")]], model = model_key("UKB")), "pred_pattern_id", model_key("UKB"), "ukb_pred")
    )
  )

  pmr_rows_shared <- base$pmrX %>%
    select(
      pmr_row_id, time_days, status, w,
      age_at_baseline, sex, ethnicity5, tenure, household_size,
      econstatus, education, ruralurban, health, disability,
      imd_decile, RGN11CD, pred_pattern_id
    ) %>%
    left_join(shared_pred_needed, by = "pred_pattern_id")

  pmr_rows_shared$age_band_10y <- make_age_band_10y(pmr_rows_shared$age_at_baseline)

  ukb_rows_base <- base$ukbX %>%
    transmute(
      participant_id,
      time_days,
      status,
      w_unweighted = 1,
      age_at_baseline,
      sex,
      ethnicity5,
      tenure,
      household_size,
      econstatus,
      education,
      ruralurban,
      health,
      disability,
      imd_decile,
      RGN11CD
    )
  ukb_rows_base$age_band_10y <- make_age_band_10y(ukb_rows_base$age_at_baseline)

  ukb_dir <- get_mortality_fit_dir(cfg, "ukb_fit", outcome_key = outcome_cfg$key)
  pmr_dir <- get_mortality_fit_dir(cfg, "pmr_fit", outcome_key = outcome_cfg$key)

  write_component_outputs(
    component_dir = ukb_dir,
    cohort_lines = c(base$shared_summary_lines, "", "fit_component=ukb_fit")
  )
  write_model_object_rds(
    component_dir = ukb_dir,
    model_payload = make_model_payload(
      transfer_model_family = transfer_model_family,
      fit_component = "ukb_fit",
      fit_obj = shared_fits[[model_key("UKB")]],
      metadata = list(outcome_key = outcome_cfg$key)
    ),
    enabled = write_model_fitting_outputs
  )

  write_component_outputs(
    component_dir = pmr_dir,
    cohort_lines = c(base$shared_summary_lines, "", "fit_component=pmr_fit")
  )
  write_model_object_rds(
    component_dir = pmr_dir,
    model_payload = make_model_payload(
      transfer_model_family = transfer_model_family,
      fit_component = "pmr_fit",
      fit_obj = shared_fits[[model_key("PMR")]],
      metadata = list(outcome_key = outcome_cfg$key)
    ),
    enabled = write_model_fitting_outputs
  )

  summary_dir <- get_mortality_summary_dir(cfg, outcome_key = outcome_cfg$key)

  risk_summaries <- list()

  for (weight_source in weight_sources) {
    stop_if_missing(weight_source$path)

    wts <- read.csv(weight_source$path, stringsAsFactors = FALSE) %>%
      mutate(eid = as.numeric(eid)) %>%
      select(eid, all_of(weight_source$column))

    ukb_with_weights <- base$ukbX %>%
      left_join(wts, by = c("participant_id" = "eid")) %>%
      rename(ukb_eval_weight = all_of(weight_source$column))

    ukb_with_weights$ukb_eval_weight <- as.numeric(ukb_with_weights$ukb_eval_weight)

    ukb_rows_missing_weight <- sum(is.na(ukb_with_weights$ukb_eval_weight))
    ukb_rows_nonmissing_weight <- sum(!is.na(ukb_with_weights$ukb_eval_weight))

    ukbwX <- ukb_with_weights %>% filter(!is.na(ukb_eval_weight))

    if (ukb_rows_missing_weight > 0L) {
      warning(
        "Dropping ", ukb_rows_missing_weight,
        " UKB rows with missing attached evaluation weights for weighted fit: ",
        weight_source$key
      )
    }

    ukbw_key <- model_key("UKB_weighted")
    log_step("Training start: ", ukbw_key, " [", weight_source$key, "]")
    ukbw_fit <- model_fns[[model_suffix]](ukbwX, "ukb_eval_weight")
    log_step("Training done: ", ukbw_key, " [", weight_source$key, "]")

    log_step("Prediction start on PMR patterns: ", ukbw_key, " [", weight_source$key, "]")
    ukbw_risks <- ukbw_fit$predict_risk(
      base$pmr_pred_input[, prediction_covars, drop = FALSE],
      horizons_days = horizons_days
    )
    log_step("Prediction done on PMR patterns: ", ukbw_key, " [", weight_source$key, "]")

    ukbw_pred_tbl <- base$pmr_pred_input %>%
      transmute(
        pred_pattern_id,
        ukbw_pred_5y = ukbw_risks[[paste0("risk_", horizons_days[1], "d")]],
        ukbw_pred_10y = ukbw_risks[[paste0("risk_", horizons_days[2], "d")]]
      )

    pmr_rows_weighted <- pmr_rows_shared %>%
      left_join(ukbw_pred_tbl, by = "pred_pattern_id")

    required_cols <- c("pmr_pred_5y", "pmr_pred_10y", "ukb_pred_5y", "ukb_pred_10y", "ukbw_pred_5y", "ukbw_pred_10y")
    if (anyNA(pmr_rows_weighted[, required_cols, drop = FALSE])) {
      stop("PMR evaluation rows contain missing prediction values after merge for weight source: ", weight_source$key)
    }

    ukb_rows_weighted <- ukb_with_weights %>%
      transmute(
        participant_id,
        time_days,
        status,
        ukb_eval_weight,
        w_unweighted = 1,
        age_at_baseline,
        sex,
        ethnicity5,
        tenure,
        household_size,
        econstatus,
        education,
        ruralurban,
        health,
        disability,
        imd_decile,
        RGN11CD
      )
    ukb_rows_weighted$age_band_10y <- make_age_band_10y(ukb_rows_weighted$age_at_baseline)

    risk_summary <- build_risk_summary(
      pmr_rows = pmr_rows_weighted,
      ukb_rows = ukb_rows_weighted
    )

    risk_summary_path <- file.path(summary_dir, paste0("pmr_reference_vs_transferred_risk__", weight_source$key, ".csv"))
    write.csv(risk_summary, risk_summary_path, row.names = FALSE)
    log_step("Saved PMR reference vs transferred risk summary to ", risk_summary_path)

    ukbw_dir <- get_mortality_fit_dir(cfg, "ukbw_fit", outcome_key = outcome_cfg$key, weight_source_key = weight_source$key)
    write_component_outputs(
      component_dir = ukbw_dir,
      cohort_lines = c(
        base$shared_summary_lines,
        "",
        paste0("fit_component=ukbw_fit/", weight_source$key),
        paste0("weight_source_key=", weight_source$key),
        paste0("weight_label=", weight_source$label),
        paste0("weight_path=", weight_source$path),
        paste0("weight_column=", weight_source$column),
        paste0("UKB rows after weight merge onto common-support UKB cohort: ", nrow(ukb_with_weights)),
        paste0("UKB rows with non-missing matched evaluation weights: ", ukb_rows_nonmissing_weight),
        paste0("UKB rows with missing/unmatched evaluation weights: ", ukb_rows_missing_weight),
        paste0("UKBW rows used for weighted model fit: ", nrow(ukbwX)),
        paste0("UKBW target-cause deaths used for weighted model fit: ", sum(ukbwX$target_death_in_followup == 1L, na.rm = TRUE))
      )
    )
    write_model_object_rds(
      component_dir = ukbw_dir,
      model_payload = make_model_payload(
        transfer_model_family = transfer_model_family,
        fit_component = paste0("ukbw_fit/", weight_source$key),
        fit_obj = ukbw_fit,
        metadata = list(
          outcome_key = outcome_cfg$key,
          weight_source_key = weight_source$key,
          weight_path = weight_source$path,
          weight_column = weight_source$column
        )
      ),
      enabled = write_model_fitting_outputs
    )

    risk_summaries[[weight_source$key]] <- risk_summary
  }

  invisible(list(risk_summaries = risk_summaries))
}

cfg <- read_mortality_bias_config(script_path)
dataset_paths <- get_mortality_harmonised_dataset_paths(cfg)
weight_source_keys <- get_mortality_weight_source_keys(cfg)
weight_sources <- lapply(weight_source_keys, function(key) get_mortality_weight_source(cfg, key = key))

results <- run_evaluation_workflow(
  cfg = cfg,
  dataset_paths = dataset_paths,
  weight_sources = weight_sources
)

log_step("Completed mortality-risk evaluation for: ", paste(vapply(weight_sources, function(x) x$key, character(1)), collapse = ", "))
