# Harmonise the 2011 Census 5% microdata and the UK Biobank extract to the Census-weighting scheme.

# ---- helpers -----------------------------------------------------------------

rename_column <- function(df, old, new) {
  if (old %in% names(df)) names(df)[names(df) == old] <- new
  df
}

# Years of education implied by each UKB qualification code; codes 5/6 fall back
# to the age finished education (capped), and anything below 7 floors at 7.
years_of_education <- function(degree_code, age_finished) {
  years <- dplyr::case_when(
    degree_code == -7 ~ 7, degree_code == 3 ~ 10, degree_code == 4 ~ 10,
    degree_code == 5 ~ -10, degree_code == 2 ~ 13, degree_code == 1 ~ 20,
    degree_code == 6 ~ -11
  )
  is5 <- degree_code == 5 & !is.na(degree_code)
  years[is5] <- age_finished[is5] - 5
  years[is5 & years >= 19] <- 19
  is6 <- degree_code == 6 & !is.na(degree_code)
  years[is6] <- age_finished[is6] - 5
  years[is6 & years >= 15] <- 15
  years[years <= 7] <- 7
  years
}

# Collapse the pipe-delimited UKB employment field to a single status, optionally
# overridden by the corrected field; precedence is 1 (employed) > 5 > 2, with the
# residual group folded into 2.
collapse_empstat <- function(x, override = NULL) {
  raw <- as.character(x)
  if (!is.null(override)) {
    override <- as.character(override)
    override[override %in% c("-7", "-3", "")] <- NA_character_
    raw[!is.na(override)] <- override[!is.na(override)]
  }
  raw[is.na(raw)] <- ""

  tokens <- lapply(strsplit(raw, "|", fixed = TRUE), function(v) {
    v <- trimws(v)
    suppressWarnings(as.numeric(v[nzchar(v)]))
  })
  multi <- vapply(tokens, function(v) sum(!is.na(v)) >= 2, logical(1))

  out <- suppressWarnings(as.numeric(raw))
  out[raw == ""] <- NA_real_
  for (i in which(multi)) {
    v <- tokens[[i]][!is.na(tokens[[i]])]
    out[i] <- if (!length(v) || all(v %in% c(-7, -3))) NA_real_
      else if (any(v == 1)) 1 else if (any(v == 5)) 5 else if (any(v == 2)) 2 else 6
  }
  out[out == 6] <- 2
  out[out %in% c(-7, -3)] <- NA_real_
  out
}

# ---- Census microdata --------------------------------------------------------

prepare_census_microdata <- function(path) {
  census <- data.table::fread(
    path,
    select = c("country", "ageh", "carsnoc", "ecopuk11", "ethnicityew",
               "health", "hlqupuk11", "meighuk11", "sex", "tenure"),
    data.table = FALSE
  )

  keep <- !is.na(census$country) & census$country == 1 &
    !is.na(census$ageh) & census$ageh >= 9 & census$ageh <= 14
  census <- census[keep, , drop = FALSE]

  census <- dplyr::mutate(
    census,
    Ethnicity = dplyr::case_when(
      ethnicityew <= 3 ~ 1,
      ethnicityew >= 4 & ethnicityew <= 5 ~ 2,
      ethnicityew >= 6 & ethnicityew <= 10 ~ 3,
      ethnicityew >= 11 & ethnicityew <= 12 ~ 4,
      ethnicityew >= 13 ~ 5
    ),
    YearsEducation = dplyr::case_when(
      hlqupuk11 == 10 ~ 7, hlqupuk11 == 11 ~ 10, hlqupuk11 == 12 ~ 10,
      hlqupuk11 == 13 ~ 12, hlqupuk11 == 14 ~ 13, hlqupuk11 == 15 ~ 20,
      hlqupuk11 == 16 ~ 15
    ),
    Education = cut(YearsEducation, breaks = c(0, 8.5, 11, 17.5, 20),
                    labels = c("1", "2", "3", "4")),
    birthyear = dplyr::case_when(
      ageh == 9 ~ "40-44", ageh == 10 ~ "45-49", ageh == 11 ~ "50-54",
      ageh == 12 ~ "55-59", ageh == 13 ~ "60-64", ageh == 14 ~ "65-69"
    ),
    SingleHousehold = dplyr::case_when(meighuk11 == 1 ~ 1, meighuk11 > 1 ~ 0),
    HealthSelfReport = dplyr::case_when(
      health %in% c(5, 4) ~ 1, health == 3 ~ 2, health %in% c(1, 2) ~ 3
    ),
    Empstat = dplyr::case_when(
      ecopuk11 <= 6 ~ 1, ecopuk11 == 7 ~ 5, ecopuk11 == 8 ~ 1, ecopuk11 == 9 ~ 7,
      ecopuk11 == 11 ~ 7, ecopuk11 == 10 ~ 2, ecopuk11 == 14 ~ 2,
      ecopuk11 == 12 ~ 3, ecopuk11 == 13 ~ 4
    )
  )

  census <- census[!is.na(census$Ethnicity), , drop = FALSE]
  census$carsnoc[census$carsnoc < 0] <- NA_real_
  census$Tenure <- suppressWarnings(as.numeric(census$tenure))
  census$Tenure[census$Tenure < 0] <- NA_real_
  census$ethnicity_protocol <- census_ethnicity_protocol()

  census[, c("birthyear", "sex", "carsnoc", "Education", "HealthSelfReport",
             "SingleHousehold", "Ethnicity", "Empstat", "Tenure", "ethnicity_protocol")]
}

# ---- UKB extract -------------------------------------------------------------

prepare_ukb_for_census <- function(path, predictors, exclude_centres = c(11004, 11005, 11003, 11022, 11023)) {
  ukb <- data.table::fread(path, header = TRUE, sep = ",", data.table = FALSE)

  renames <- c("31-0.0" = "sex", "34-0.0" = "birthyear_detail", "54-0.0" = "assessment_center",
               "21022-0.0" = "age", "6138-0.0" = "Degree", "6142-0.0" = "Empstat",
               "20119-0.0" = "Empstat_corrected", "680-0.0" = "Tenure", "728-0.0" = "carsnoc",
               "2178-0.0" = "HealthSelfReport", "709-0.0" = "NoInHH", "21000-0.0" = "Ethnicity",
               "845-0.0" = "AgeCompletedEduc")
  for (old in names(renames)) ukb <- rename_column(ukb, old, renames[[old]])

  ukb <- ukb[
    !is.na(ukb$age) & ukb$age >= 40 & ukb$age <= 69 &
      !is.na(ukb$assessment_center) & !ukb$assessment_center %in% exclude_centres,
    , drop = FALSE
  ]

  is_subgroup <- ukb$Ethnicity >= 1000 & !is.na(ukb$Ethnicity)
  ukb$Ethnicity[is_subgroup] <- round(ukb$Ethnicity[is_subgroup] / 1000)
  ukb$Ethnicity[ukb$Ethnicity %in% c(5, 6)] <- 5
  ukb$Ethnicity[ukb$Ethnicity <= 0] <- NA
  ukb$AgeCompletedEduc[ukb$AgeCompletedEduc < 0] <- NA

  degree_tokens <- lapply(strsplit(as.character(ukb$Degree), "|", fixed = TRUE), function(v) {
    v <- trimws(v)
    suppressWarnings(as.numeric(v[nzchar(v)]))
  })
  degree_mat <- t(vapply(degree_tokens, function(v) {
    length(v) <- 6L  # pad/truncate to 6 codes
    v
  }, numeric(6)))
  years_cols <- vapply(seq_len(6L), function(j) {
    years_of_education(degree_mat[, j], ukb$AgeCompletedEduc)
  }, numeric(nrow(ukb)))
  ukb$YearsEducation <- apply(years_cols, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  ukb$Education <- cut(ukb$YearsEducation, breaks = c(0, 8.5, 11, 17.5, 20),
                       labels = c("1", "2", "3", "4"))

  ukb$age_2011 <- 2011L - ukb$birthyear_detail
  ukb$birthyear <- dplyr::case_when(
    ukb$age_2011 >= 40 & ukb$age_2011 <= 44 ~ "40-44",
    ukb$age_2011 >= 45 & ukb$age_2011 <= 49 ~ "45-49",
    ukb$age_2011 >= 50 & ukb$age_2011 <= 54 ~ "50-54",
    ukb$age_2011 >= 55 & ukb$age_2011 <= 59 ~ "55-59",
    ukb$age_2011 >= 60 & ukb$age_2011 <= 64 ~ "60-64",
    ukb$age_2011 >= 65 & ukb$age_2011 <= 69 ~ "65-69"
  )

  ukb$SingleHousehold <- NA_real_
  ukb$SingleHousehold[ukb$NoInHH == 1] <- 1
  ukb$SingleHousehold[ukb$NoInHH > 1] <- 0

  ukb$carsnoc[ukb$carsnoc %in% c(-1, -3)] <- NA
  ukb$HealthSelfReport[ukb$HealthSelfReport %in% c(-1, -3)] <- NA
  ukb$HealthSelfReport <- dplyr::case_when(
    ukb$HealthSelfReport == 4 ~ 1, ukb$HealthSelfReport == 3 ~ 2,
    ukb$HealthSelfReport %in% c(1, 2) ~ 3
  )
  ukb$carsnoc <- ukb$carsnoc - 1
  ukb$sex[ukb$sex == 0] <- 2
  ukb$Tenure <- dplyr::case_when(
    ukb$Tenure == 1 ~ 1, ukb$Tenure == 2 ~ 2, ukb$Tenure == 3 ~ 4,
    ukb$Tenure == 4 ~ 4, ukb$Tenure == 5 ~ 3, ukb$Tenure == 6 ~ 5
  )
  ukb$Empstat <- collapse_empstat(ukb$Empstat, ukb$Empstat_corrected)

  prepared <- ukb[, c("eid", predictors), drop = FALSE]
  prepared$source <- "UKB"
  prepared$sample_weight <- 1
  prepared[, c("eid", "source", "sample_weight", predictors), drop = FALSE]
}
