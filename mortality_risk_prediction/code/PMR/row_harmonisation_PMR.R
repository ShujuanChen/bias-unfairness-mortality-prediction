# Row filters restricting PMR to the analysis cohort: England residents aged 40-69 at the 2011 Census.

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})

pmr_filter_age_ukb_eligible <- function(df,
                                        age_col = "age_at_baseline",
                                        min_age_2011 = 40L,
                                        max_age_2011 = 69L,
                                        keep = TRUE,
                                        keep_debug_cols = FALSE,
                                        verbose = TRUE) {
  stopifnot(is.data.frame(df))

  if (!age_col %in% names(df)) stop("Age filter missing column: ", age_col)
  age_2011 <- suppressWarnings(as.integer(df[[age_col]]))
  eligible <- !is.na(age_2011) & age_2011 >= as.integer(min_age_2011) & age_2011 <= as.integer(max_age_2011)

  if (isTRUE(keep_debug_cols)) {
    df$age_2011_03_27 <- age_2011
    df$ukb_eligible_age <- as.integer(eligible)
  }

  if (isTRUE(verbose)) {
    message(
      "Age filter: keeping ", sum(eligible, na.rm = TRUE), " / ", nrow(df),
      " (dropped ", sum(!eligible, na.rm = TRUE), ")."
    )
  }

  if (isTRUE(keep)) df <- df[eligible, , drop = FALSE]
  df
}

pmr_filter_england_only <- function(df,
                                     col = "RGN11CD",
                                     keep = TRUE,
                                     keep_debug_cols = FALSE,
                                     verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop("England filter missing column: ", col)

  x <- as.character(df[[col]])
  ok <- !is.na(x) & startsWith(x, "E")

  if (isTRUE(keep_debug_cols)) {
    df$ukb_eligible_england <- as.integer(ok)
  }

  if (isTRUE(verbose)) {
    message(
      "England filter: keeping ", sum(ok, na.rm = TRUE), " / ", nrow(df),
      " (dropped ", sum(!ok, na.rm = TRUE), ")."
    )
  }

  if (isTRUE(keep)) df <- df[ok, , drop = FALSE]
  df
}
