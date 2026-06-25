options(stringsAsFactors = FALSE)

# Core of the HSE participation-weighting pipeline.

# ---- 5-group ethnicity -------------------------------------------------------

HSE_ETHNICITY5_LEVELS <- c("White", "Mixed", "Asian", "Black", "Chinese/Other")

coerce_hse_ethnicity5 <- function(x, context = "ethnicity") {
  chr <- as.character(x)
  chr[is.na(chr) | trimws(chr) == ""] <- NA_character_
  unexpected <- sort(unique(chr[!is.na(chr) & !chr %in% HSE_ETHNICITY5_LEVELS]))
  if (length(unexpected)) {
    stop(context, " has unexpected ethnicity values: ",
         paste(shQuote(unexpected), collapse = ", "),
         ". Expected only: ", paste(HSE_ETHNICITY5_LEVELS, collapse = ", "),
         call. = FALSE)
  }
  factor(chr, levels = HSE_ETHNICITY5_LEVELS)
}

normalise_weights <- function(x) x / mean(x, na.rm = TRUE)

winsorise_weights <- function(w, probs = c(0.01, 0.99)) {
  stopifnot(length(probs) == 2L, probs[1] >= 0, probs[2] <= 1, probs[1] < probs[2])
  qs <- stats::quantile(w, probs = probs, na.rm = TRUE, names = FALSE)
  pmin(pmax(w, qs[1]), qs[2])
}

write_note_file <- function(path, lines) writeLines(as.character(lines), path, useBytes = TRUE)

# ---- configuration -----------------------------------------------------------

load_hse_config <- function(project_root) {
  raw <- jsonlite::read_json(file.path(project_root, "config", "hse_weighting.json"),
                             simplifyVector = FALSE)

  predictors <- data.frame(
    label              = vapply(raw$predictors, function(p) p$label, character(1)),
    ukb_field          = vapply(raw$predictors, function(p) if (is.null(p$ukb_field)) NA_real_ else as.numeric(p$ukb_field), numeric(1)),
    include_prediction = vapply(raw$predictors, function(p) isTRUE(p$include_prediction), logical(1)),
    recode             = vapply(raw$predictors, function(p) p$recode, character(1)),
    type               = vapply(raw$predictors, function(p) p$type, character(1)),
    stringsAsFactors = FALSE
  )

  prediction_labels <- predictors$label[predictors$include_prediction]
  # ethnicity_harmonised is derived (not a UKB field) and appended as the last
  # categorical predictor of the participation model.
  model_info <- rbind(
    data.frame(label = prediction_labels,
               type = predictors$type[predictors$include_prediction],
               stringsAsFactors = FALSE),
    data.frame(label = "ethnicity_harmonised", type = "cat", stringsAsFactors = FALSE)
  )

  list(
    predictors = predictors,
    prediction_labels = prediction_labels,
    model_info = model_info,
    model_labels = model_info$label,
    hse_waves = raw$hse_waves,
    hse_variable_map = raw$hse_variable_map,
    hse_weight_column = raw$hse_weight_column,
    ukb_age_field_override = raw$ukb_age_field_override,
    ukb_england_exclude_centres = as.numeric(unlist(raw$ukb_england_exclude_centres)),
    age_min = raw$age_min,
    age_max = raw$age_max
  )
}

# ---- result-tree paths -------------------------------------------------------

hse_results_dir <- function(project_root) file.path(project_root, "results")

hse_intermediate_dir <- function(project_root) {
  d <- file.path(hse_results_dir(project_root), "intermediate")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

hse_weights_dir <- function(project_root, model) {
  d <- file.path(hse_results_dir(project_root), model, "weights")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

hse_model_input_path <- function(project_root) file.path(hse_intermediate_dir(project_root), "model_input.csv")
hse_model_input_note_path <- function(project_root) file.path(hse_intermediate_dir(project_root), "model_input.txt")
hse_weight_path <- function(project_root, model) file.path(hse_weights_dir(project_root, model), "ukb_weights.csv")
hse_weight_note_path <- function(project_root, model) file.path(hse_weights_dir(project_root, model), "ukb_weights.txt")

# ---- recodes + harmonisation -------------------------------------------------

.hse_code_dir <- tryCatch({
  arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("--file=", "", arg[1]))) else getwd()
}, error = function(e) getwd())

source(file.path(.hse_code_dir, "hse_recodes.R"))
source(file.path(.hse_code_dir, "hse_harmonise.R"))

# ---- stack HSE + UKB ---------------------------------------------------------

prepare_training_data <- function(project_root, cfg) {
  hse_obj <- prepare_hse_cohort(project_root, cfg)
  ukb_obj <- prepare_ukb_cohort(project_root, cfg)

  info <- cfg$model_info
  labels <- info$label

  hse_model <- hse_obj$final_cohort[, c("eid", "weight_individual", labels), drop = FALSE]
  hse_model$source <- "HSE"
  names(hse_model)[names(hse_model) == "weight_individual"] <- "sample_weight"
  hse_model$sample_weight <- normalise_weights(hse_model$sample_weight)

  ukb_model <- ukb_obj$final_cohort[, c("eid", "weight_individual", labels), drop = FALSE]
  ukb_model$source <- "UKB"
  names(ukb_model)[names(ukb_model) == "weight_individual"] <- "sample_weight"

  stacked <- data.table::rbindlist(list(hse_model, ukb_model), fill = TRUE)
  complete <- complete.cases(stacked[, ..labels])
  combined <- as.data.frame(stacked[complete, , drop = FALSE])

  for (label in info$label[info$type != "con"]) {
    if (identical(label, "ethnicity_harmonised")) {
      combined[[label]] <- coerce_hse_ethnicity5(combined[[label]], context = "HSE training ethnicity_harmonised")
    } else {
      combined[[label]] <- as.factor(combined[[label]])
    }
  }

  active <- vapply(info$label, function(label) {
    vals <- combined[[label]]
    length(unique(as.character(vals[!is.na(vals)]))) > 1
  }, logical(1))
  info_active <- info[active, , drop = FALSE]
  if (!nrow(info_active)) stop("No varying predictors remain after the HSE/UKB restrictions.")

  combined$sample <- as.integer(combined$source == "UKB")

  list(
    hse_obj = hse_obj,
    ukb_obj = ukb_obj,
    prediction_info = info_active,
    prediction_labels = info_active$label,
    combined = combined,
    model_input = combined[, c("eid", "source", "sample_weight", info_active$label), drop = FALSE]
  )
}

# ---- design matrix for the lasso (one-hot with pairwise interactions) --------

build_dummy_design <- function(df, info) {
  cat_labels <- info$label[info$type == "cat"]
  bin_labels <- info$label[info$type == "bin"]
  con_labels <- info$label[info$type == "con"]

  bin_part <- if (length(bin_labels)) {
    fastDummies::dummy_cols(subset(df, select = bin_labels), select_columns = bin_labels,
                            remove_first_dummy = TRUE, ignore_na = TRUE)
  } else data.frame()
  cat_part <- if (length(cat_labels)) {
    fastDummies::dummy_cols(subset(df, select = cat_labels), select_columns = cat_labels,
                            remove_first_dummy = FALSE, ignore_na = TRUE)
  } else data.frame()

  expanded <- c(colnames(cat_part), colnames(bin_part))
  con_part <- subset(df, select = c("eid", con_labels))
  design <- cbind(con_part, cat_part, bin_part)
  keep <- c(con_labels, unique(expanded[!expanded %in% c(cat_labels, bin_labels)]))
  list(data = design, vars = keep)
}

# ---- fits --------------------------------------------------------------------

assemble_inverse_odds_weights <- function(combined, prob_ukb) {
  if (any(!is.finite(prob_ukb) | prob_ukb <= 0 | prob_ukb >= 1, na.rm = TRUE)) {
    stop("Participation probabilities produced invalid inverse-odds weights.")
  }
  is_ukb <- combined$source == "UKB"
  combined$prob_ukb <- prob_ukb
  combined$inverse_odds <- (1 - prob_ukb) / prob_ukb
  combined$w <- NA_real_
  combined$w[is_ukb] <- normalise_weights(winsorise_weights(combined$inverse_odds[is_ukb]))
  combined$w[!is_ukb] <- combined$sample_weight[!is_ukb]
  combined
}

fit_lasso_weights <- function(project_root, cfg) {
  prepared <- prepare_training_data(project_root, cfg)
  combined <- prepared$combined
  design <- build_dummy_design(combined, prepared$prediction_info)
  dummy_data <- data.frame(design$data[, design$vars, drop = FALSE], check.names = TRUE)

  x <- Matrix::sparse.model.matrix(stats::as.formula("~ .*."), data = dummy_data)
  y <- as.integer(combined$source == "UKB")

  set.seed(1234)
  cvfit <- glmnet::cv.glmnet(
    x = x, y = y, family = "binomial", type.measure = "class",
    nfolds = 5, weights = combined$sample_weight, nlambda = 100, parallel = FALSE
  )
  prob_ukb <- as.numeric(stats::predict(cvfit, newx = x, s = "lambda.min", type = "response")[, 1])
  combined <- assemble_inverse_odds_weights(combined, prob_ukb)

  ukb_rows <- combined$source == "UKB"
  prepared$ukb_weights <- combined[ukb_rows, c("eid", "w", "prob_ukb"), drop = FALSE]
  prepared$combined <- combined
  prepared
}

SL.ranger.1t  <- function(...) SuperLearner::SL.ranger(...,  num.threads = 1)
SL.xgboost.1t <- function(...) SuperLearner::SL.xgboost(..., nthread = 1)

fit_superlearner_weights <- function(project_root, cfg,
                                     sl_library = c("SL.glm", "SL.glmnet", "SL.ranger.1t", "SL.xgboost.1t")) {
  prepared <- prepare_training_data(project_root, cfg)
  combined <- prepared$combined
  vars <- prepared$prediction_labels

  for (v in vars) {
    if (!is.numeric(combined[[v]]) && !is.factor(combined[[v]])) combined[[v]] <- as.factor(combined[[v]])
  }

  design_df <- combined[, c("sample", vars), drop = FALSE]
  formula <- stats::as.formula(paste("sample ~ (", paste(vars, collapse = "+"), ")"))
  x_mm <- stats::model.matrix(formula, design_df)[, -1, drop = FALSE]
  if (!ncol(x_mm)) stop("The HSE/UKB design matrix has no non-intercept columns.")
  colnames(x_mm) <- make.names(colnames(x_mm), unique = TRUE)

  set.seed(8791)
  fit_sl <- SuperLearner::SuperLearner(
    Y = combined$sample, X = as.data.frame(x_mm), family = stats::binomial(),
    SL.library = sl_library, method = "method.NNLS",
    cvControl = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE),
    obsWeights = combined$sample_weight
  )

  prob_ukb <- as.numeric(fit_sl$SL.predict)
  if (length(prob_ukb) != nrow(combined)) stop("SuperLearner predictions are length-mismatched.")
  combined <- assemble_inverse_odds_weights(combined, prob_ukb)

  ukb_rows <- combined$source == "UKB"
  prepared$ukb_weights <- combined[ukb_rows, c("eid", "w", "prob_ukb"), drop = FALSE]
  prepared$combined <- combined
  prepared$sl_library <- sl_library
  prepared$p_ukb_marginal <- mean(ukb_rows)
  prepared
}
