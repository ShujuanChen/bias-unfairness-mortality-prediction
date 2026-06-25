options(stringsAsFactors = FALSE)

# Core of the Census participation-weighting pipeline.

census_predictors <- function() {
  c("birthyear", "sex", "carsnoc", "Education", "HealthSelfReport",
    "SingleHousehold", "Ethnicity", "Empstat", "Tenure")
}

census_ethnicity_protocol <- function() "HSE5_CHINESE_OTHER"

normalise_weights <- function(x) x / mean(x, na.rm = TRUE)

winsorise_weights <- function(w, probs = c(0.01, 0.99)) {
  stopifnot(length(probs) == 2L, probs[1] >= 0, probs[2] <= 1, probs[1] < probs[2])
  qs <- stats::quantile(w, probs = probs, na.rm = TRUE, names = FALSE)
  pmin(pmax(w, qs[1]), qs[2])
}

write_note_file <- function(path, lines) writeLines(as.character(lines), path, useBytes = TRUE)

# ---- result-tree paths -------------------------------------------------------

census_results_dir <- function(project_root) file.path(project_root, "results")

census_intermediate_dir <- function(project_root) {
  d <- file.path(census_results_dir(project_root), "intermediate")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

census_weights_subdir <- function(project_root) {
  d <- file.path(census_results_dir(project_root), "superlearner", "weights")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

census_model_input_path <- function(project_root) file.path(census_intermediate_dir(project_root), "model_input.csv")
census_model_input_note_path <- function(project_root) file.path(census_intermediate_dir(project_root), "model_input.txt")
census_weight_path <- function(project_root) file.path(census_weights_subdir(project_root), "ukb_weights.csv")
census_weight_note_path <- function(project_root) file.path(census_weights_subdir(project_root), "ukb_weights.txt")

# ---- harmonisation -----------------------------------------------------------

.census_code_dir <- tryCatch({
  arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("--file=", "", arg[1]))) else getwd()
}, error = function(e) getwd())

source(file.path(.census_code_dir, "census_harmonise.R"))

# ---- Super Learner fit -------------------------------------------------------

coerce_census_predictors <- function(df, predictors) {
  for (v in predictors) df[[v]] <- as.factor(df[[v]])
  df
}

census_active_predictors <- function(df, predictors) {
  keep <- vapply(predictors, function(v) {
    length(unique(df[[v]][!is.na(df[[v]])])) > 1
  }, logical(1))
  predictors[keep]
}

# Pin the tree learners to one thread so the fit reproduces from the seed alone.
SL.ranger.1t  <- function(...) SuperLearner::SL.ranger(...,  num.threads = 1)
SL.xgboost.1t <- function(...) SuperLearner::SL.xgboost(..., nthread = 1)

fit_census_superlearner <- function(model_input, predictors = census_predictors()) {
  prepared <- coerce_census_predictors(model_input, predictors)
  active <- census_active_predictors(prepared, predictors)
  if (!length(active)) stop("No varying predictors remain in the Census+UKB model input.")

  x_mm <- stats::model.matrix(
    stats::as.formula(paste("~", paste(active, collapse = "+"))),
    data = prepared[, active, drop = FALSE]
  )[, -1, drop = FALSE]
  colnames(x_mm) <- make.names(colnames(x_mm), unique = TRUE)

  y <- as.integer(prepared$source == "UKB")
  sample_weight <- prepared$sample_weight
  sl_library <- c("SL.glm", "SL.glmnet", "SL.ranger.1t", "SL.xgboost.1t")

  set.seed(8791)
  fit_sl <- SuperLearner::SuperLearner(
    Y = y, X = as.data.frame(x_mm), family = stats::binomial(),
    SL.library = sl_library, method = "method.NNLS",
    cvControl = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE),
    obsWeights = sample_weight
  )

  prob_ukb <- as.numeric(fit_sl$SL.predict)
  ukb_idx <- prepared$source == "UKB"
  if (any(!is.finite(prob_ukb) | prob_ukb <= 0 | prob_ukb >= 1, na.rm = TRUE)) {
    stop("SuperLearner participation probabilities produced invalid inverse-odds weights.")
  }

  inverse_odds <- (1 - prob_ukb) / prob_ukb
  prepared$prob_ukb <- prob_ukb
  prepared$inverse_odds <- inverse_odds
  prepared$w <- NA_real_
  prepared$w[ukb_idx] <- normalise_weights(winsorise_weights(inverse_odds[ukb_idx]))
  prepared$w[!ukb_idx] <- prepared$sample_weight[!ukb_idx]

  list(
    combined = prepared,
    active_predictors = active,
    fit_sl = fit_sl,
    sl_library = sl_library,
    p_ukb_marginal = mean(ukb_idx),
    ukb_weights = prepared[ukb_idx, c("eid", "w", "prob_ukb"), drop = FALSE]
  )
}
