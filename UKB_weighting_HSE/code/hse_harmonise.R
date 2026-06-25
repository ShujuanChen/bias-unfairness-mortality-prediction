# Harmonise the HSE survey waves and the UK Biobank extract to a shared schema.

# ---- ethnicity collapse ------------------------------------------------------

labelled_to_text <- function(x) as.character(haven::as_factor(x, levels = "labels"))

WHITE_HSE_LABELS <- c("White", "White - British", "White - Irish",
                      "Any other white background")

collapse_hse_ethnicity <- function(detail) {
  group <- rep(NA_character_, length(detail))
  group[detail %in% WHITE_HSE_LABELS] <- "White"
  group[grepl("^Mixed", detail) | detail == "Mixed" | detail == "Any other mixed background"] <- "Mixed"
  group[grepl("^Asian", detail) | detail == "Asian or Asian British"] <- "Asian"
  group[grepl("^Black", detail) | detail == "Black or Black British"] <- "Black"
  group[detail %in% c("Chinese", "Chinese or other ethnic group", "Any other (please describe)")] <- "Chinese/Other"
  coerce_hse_ethnicity5(group, context = "HSE ethnicity broad mapping")
}

collapse_ukb_ethnicity <- function(codes) {
  numeric_codes <- suppressWarnings(as.numeric(as.character(codes)))
  group <- rep(NA_character_, length(numeric_codes))
  group[numeric_codes %in% c(1, 1001, 1002, 1003)] <- "White"
  group[numeric_codes %in% c(2, 2001, 2002, 2003, 2004)] <- "Mixed"
  group[numeric_codes %in% c(3, 3001, 3002, 3003, 3004)] <- "Asian"
  group[numeric_codes %in% c(4, 4001, 4002, 4003)] <- "Black"
  group[numeric_codes %in% c(5, 6)] <- "Chinese/Other"
  coerce_hse_ethnicity5(group, context = "UKBB ethnicity broad mapping")
}

# ---- HSE side ----------------------------------------------------------------

prepare_hse_wave <- function(dat, variable_map, year, ethnicity_column, weight_column) {
  use <- vapply(variable_map, function(v) {
    col <- v$columns[[year]]
    !is.null(col) && nzchar(col)
  }, logical(1))
  used <- variable_map[use]

  source_cols <- vapply(used, function(v) v$columns[[year]], character(1))
  labels <- vapply(used, function(v) v$label, character(1))
  types <- vapply(used, function(v) v$type, character(1))

  wave <- as.data.frame(dat[, source_cols])
  names(wave) <- labels
  for (i in seq_along(labels)) {
    if (identical(types[i], "con")) {
      wave[[labels[i]]] <- suppressWarnings(as.numeric(wave[[labels[i]]]))
    } else {
      wave[[labels[i]]] <- labelled_to_text(wave[[labels[i]]])
    }
  }

  wave$weight_individual <- suppressWarnings(as.numeric(dat[[weight_column]]))
  wave$ethnicity_detail <- labelled_to_text(dat[[ethnicity_column]])
  wave$ethnicity_broad <- collapse_hse_ethnicity(wave$ethnicity_detail)
  wave$ethnicity_harmonised <- wave$ethnicity_broad
  wave$year <- year
  wave
}

recode_hse_columns <- function(df) {
  df$sex <- hse_recode_sex(df$sex)
  df$overallhealth <- hse_recode_health(df$overallhealth)
  df$smoking_status <- hse_recode_smoking(df$smoking_status)
  df$employment_status <- hse_recode_employment(df$employment_status)
  df$income <- hse_recode_income(df$income)
  df$alcfrequency <- hse_recode_alcohol(df$alcfrequency)
  df$urbanisation <- hse_recode_urbanisation(df$urbanisation)
  df$education_age <- hse_recode_education_age(df$education_age)
  df$household_size <- hse_recode_household_size(df$household_size)
  df$age <- suppressWarnings(as.numeric(df$age))
  df$height <- suppressWarnings(as.numeric(df$height))
  df$weight <- suppressWarnings(as.numeric(df$weight))
  df$bmi <- df$weight / (df$height / 100)^2
  df$bmi_cat <- hse_recode_bmi_category(df$bmi)
  df
}

prepare_hse_cohort <- function(project_root, cfg) {
  prediction_labels <- cfg$prediction_labels

  waves <- lapply(cfg$hse_waves, function(w) {
    sav <- haven::read_sav(file.path(project_root, "data", "HSE", w$path))
    prepare_hse_wave(
      dat = sav,
      variable_map = cfg$hse_variable_map,
      year = w$year,
      ethnicity_column = w$ethnicity_column,
      weight_column = cfg$hse_weight_column
    )
  })

  stacked <- as.data.frame(data.table::rbindlist(waves, fill = TRUE))
  stacked <- recode_hse_columns(stacked)
  stacked$eid <- seq_len(nrow(stacked))

  eligible <- !is.na(stacked$age) & stacked$age >= cfg$age_min & stacked$age <= cfg$age_max
  cohort <- stacked[eligible, , drop = FALSE]
  cohort$complete_predictors <- complete.cases(cohort[, prediction_labels, drop = FALSE])

  list(all = stacked, final_cohort = cohort, prediction_labels = prediction_labels)
}

# ---- UKB side ----------------------------------------------------------------

resolve_ukb_column <- function(field, label, available, age_override) {
  if (identical(label, "age") && paste0("p", age_override) %in% available) {
    return(paste0("p", age_override))
  }
  candidates <- c(paste0("p", field, "_i0"), paste0("p", field))
  hit <- candidates[candidates %in% available]
  if (!length(hit)) stop(sprintf("Cannot map '%s' (field %s) to a UKB column.", label, field))
  hit[1]
}

split_codes <- function(x) {
  lapply(strsplit(as.character(x), "|", fixed = TRUE), function(tok) {
    tok <- trimws(tok)
    tok[nzchar(tok)]
  })
}

collapse_employment_codes <- function(x) {
  vapply(split_codes(x), function(tok) {
    tok <- tok[!tok %in% c("-3", "-7")]
    if (!length(tok)) return(NA_character_)
    if ("1" %in% tok) return("1")
    if ("5" %in% tok) return("5")
    if ("2" %in% tok) return("2")
    if (any(tok %in% c("3", "4", "6", "7"))) return("3")
    NA_character_
  }, character(1))
}

collapse_degree_codes <- function(x) {
  vapply(split_codes(x), function(tok) {
    tok <- tok[!tok %in% c("-3")]
    if (!length(tok)) return(NA_character_)
    if ("1" %in% tok) "1" else "0"
  }, character(1))
}

prepare_ukb_cohort <- function(project_root, cfg) {
  predictors <- cfg$predictors
  with_field <- predictors[!is.na(predictors$ukb_field), , drop = FALSE]
  prediction_labels <- cfg$prediction_labels

  raw <- data.table::fread(
    file.path(project_root, "data", "UKB", "UKB_for_harmonisation_with_HSE_participant.csv"),
    na.strings = c("", "NA")
  )

  selected <- data.frame(eid = raw$eid, check.names = FALSE)
  for (i in seq_len(nrow(with_field))) {
    label <- with_field$label[i]
    col <- resolve_ukb_column(with_field$ukb_field[i], label, names(raw), cfg$ukb_age_field_override)
    selected[[label]] <- raw[[col]]
  }

  selected$employment_status <- collapse_employment_codes(selected$employment_status)
  selected$education_degree <- collapse_degree_codes(selected$education_degree)

  recoded <- data.frame(eid = selected$eid, check.names = FALSE)
  for (i in seq_len(nrow(predictors))) {
    label <- predictors$label[i]
    recoded[[label]] <- apply_ukb_recode(predictors$recode[i], selected, label)
  }

  recoded$weight_individual <- 1
  recoded$ethnicity_harmonised <- collapse_ukb_ethnicity(recoded$ethnic_background)

  eligible_age <- !is.na(recoded$age) & recoded$age >= cfg$age_min & recoded$age <= cfg$age_max
  eligible_england <- !is.na(recoded$assessment_center) &
    !recoded$assessment_center %in% cfg$ukb_england_exclude_centres

  cohort <- recoded[eligible_age & eligible_england, , drop = FALSE]
  cohort$complete_predictors <- complete.cases(cohort[, prediction_labels, drop = FALSE])

  list(all = recoded, final_cohort = cohort, prediction_labels = prediction_labels)
}
