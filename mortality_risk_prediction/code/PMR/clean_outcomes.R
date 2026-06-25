# ICD-10 outcome helpers: build all-cause or cause-specific survival outcomes from harmonised death columns.

normalise_icd_code <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[is.na(x)] <- NA_character_
  y[y %in% c("", "NA")] <- NA_character_
  y <- gsub("[^A-Z0-9]", "", y)
  y[!nzchar(y)] <- NA_character_
  y
}

normalise_icd_spec <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[is.na(x)] <- NA_character_
  y[y %in% c("", "NA")] <- NA_character_
  y <- gsub("[^A-Z0-9-]", "", y)
  y[!nzchar(y)] <- NA_character_
  y
}

expand_icd_specs_for_matching <- function(specs) {
  specs <- unique(stats::na.omit(normalise_icd_spec(specs)))

  out <- list(
    icd10_root3 = character(),
    icd10_exact4 = character()
  )

  for (spec in specs) {
    if (grepl("^[A-Z][0-9]{2}$", spec)) {
      out$icd10_root3 <- c(out$icd10_root3, spec)
      next
    }

    if (grepl("^[A-Z][0-9]{2}-[A-Z][0-9]{2}$", spec)) {
      parts <- strsplit(spec, "-", fixed = TRUE)[[1]]
      left <- parts[1]
      right <- parts[2]
      if (substr(left, 1, 1) != substr(right, 1, 1)) {
        stop("ICD-10 root ranges must share the same letter prefix: ", spec, call. = FALSE)
      }
      lo <- as.integer(substr(left, 2, 3))
      hi <- as.integer(substr(right, 2, 3))
      if (is.na(lo) || is.na(hi) || lo > hi) {
        stop("Invalid ICD-10 range: ", spec, call. = FALSE)
      }
      out$icd10_root3 <- c(out$icd10_root3, sprintf("%s%02d", substr(left, 1, 1), lo:hi))
      next
    }

    if (grepl("^[A-Z][0-9]{3}$", spec)) {
      out$icd10_exact4 <- c(out$icd10_exact4, spec)
      next
    }

    if (grepl("^[A-Z][0-9]{3}-[A-Z][0-9]{3}$", spec)) {
      parts <- strsplit(spec, "-", fixed = TRUE)[[1]]
      left <- parts[1]
      right <- parts[2]
      if (substr(left, 1, 1) != substr(right, 1, 1)) {
        stop("ICD-10 exact ranges must share the same letter prefix: ", spec, call. = FALSE)
      }
      lo <- as.integer(substr(left, 2, 4))
      hi <- as.integer(substr(right, 2, 4))
      if (is.na(lo) || is.na(hi) || lo > hi) {
        stop("Invalid ICD-10 exact range: ", spec, call. = FALSE)
      }
      out$icd10_exact4 <- c(out$icd10_exact4, sprintf("%s%03d", substr(left, 1, 1), lo:hi))
      next
    }

    stop("Unsupported ICD spec: ", spec, call. = FALSE)
  }

  out$icd10_root3 <- unique(out$icd10_root3)
  out$icd10_exact4 <- unique(out$icd10_exact4)
  out
}

icd_vector_matches_plan <- function(codes, plan) {
  norm <- normalise_icd_code(codes)
  out <- rep(FALSE, length(norm))

  if (length(plan$icd10_root3) > 0L) {
    idx <- grepl("^[A-Z][0-9]{2}", norm)
    out[idx] <- out[idx] | substr(norm[idx], 1, 3) %in% plan$icd10_root3
  }

  if (length(plan$icd10_exact4) > 0L) {
    idx <- grepl("^[A-Z][0-9]{3}", norm)
    out[idx] <- out[idx] | substr(norm[idx], 1, 4) %in% plan$icd10_exact4
  }

  out[is.na(norm)] <- FALSE
  out
}

get_harmonised_death_icd_columns <- function(df, match_scope = "underlying") {
  underlying_col <- "underlying_cause_of_death_icd"
  secondary_cols <- grep("^secondary_cause_of_death_[0-9]+_icd$", names(df), value = TRUE)

  if (!underlying_col %in% names(df)) {
    stop("Missing underlying cause-of-death column: ", underlying_col, call. = FALSE)
  }

  if (identical(match_scope, "underlying")) {
    return(underlying_col)
  }

  stop(
    "Unsupported outcome match_scope: ", match_scope,
    ". Use `underlying`.",
    call. = FALSE
  )
}

build_outcome_death_match <- function(df, outcome_cfg) {
  if (identical(outcome_cfg$type, "all_cause")) {
    return(rep(TRUE, nrow(df)))
  }

  specs <- outcome_cfg$icd10
  plan <- expand_icd_specs_for_matching(specs)
  icd_cols <- get_harmonised_death_icd_columns(df, match_scope = outcome_cfg$match_scope)

  matched <- rep(FALSE, nrow(df))
  for (col in icd_cols) {
    matched <- matched | icd_vector_matches_plan(df[[col]], plan)
  }

  matched
}

build_survival_outcome_from_config <- function(df, death_date_col, outcome_cfg, t0, t_admin_end) {
  stopifnot(is.data.frame(df))
  if (!death_date_col %in% names(df)) {
    stop("Missing death date column: ", death_date_col, call. = FALSE)
  }

  death_date <- as.Date(df[[death_date_col]])
  matched_target <- build_outcome_death_match(df, outcome_cfg)
  death_any_in_followup <- !is.na(death_date) & death_date <= t_admin_end
  target_death_in_followup <- death_any_in_followup & matched_target

  censor_date <- rep(as.Date(t_admin_end), nrow(df))
  has_death_date <- !is.na(death_date)
  censor_date[has_death_date] <- pmin(death_date[has_death_date], as.Date(t_admin_end))

  data.frame(
    censor_date = as.Date(censor_date, origin = "1970-01-01"),
    time_days = as.integer(difftime(censor_date, as.Date(t0), units = "days")),
    status = as.integer(target_death_in_followup),
    death_any_in_followup = as.integer(death_any_in_followup),
    target_death_in_followup = as.integer(target_death_in_followup),
    stringsAsFactors = FALSE
  )
}
