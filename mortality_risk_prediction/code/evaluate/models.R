# Cox proportional-hazards model library for the linear_cox sensitivity analysis.

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it first.", call. = FALSE)
  }
}

coerce_harmonised_covariates <- function(df, rgn_levels = NULL) {
  df %>%
    dplyr::mutate(
      age_at_baseline = as.numeric(age_at_baseline),
      sex = factor(sex, levels = c("Male", "Female")),
      ethnicity5 = factor(ethnicity5, levels = c("White", "Mixed", "Asian", "Black", "Chinese/Other")),
      tenure = factor(
        tenure,
        levels = c(
          "Owned outright", "Owned with a mortgage or loan", "Shared ownership",
          "Social rented", "Private rented", "Living rent free", "Other"
        )
      ),
      household_size = factor(household_size, levels = c("1", "2+")),
      econstatus = factor(econstatus, levels = c("Employed", "Unemployed", "Retired", "Other")),
      education = factor(education, levels = c("Level 1", "Level 2", "Level 3", "Level 4")),
      ruralurban = factor(
        ruralurban,
        levels = c("Urban", "Town and Fringe", "Village", "Hamlet and Isolated Dwelling")
      ),
      health = factor(health, levels = c("Good", "Fair", "Bad")),
      disability = factor(disability, levels = c("Yes", "No")),
      RGN11CD = if (!is.null(rgn_levels)) {
        factor(RGN11CD, levels = rgn_levels)
      } else {
        factor(RGN11CD)
      },
      imd_decile = factor(as.integer(as.character(imd_decile)), levels = 1:10)
    )
}

get_rhs_term_map <- function(use_spline_age = FALSE) {
  c(
    age_at_baseline = if (use_spline_age) "splines::ns(age_at_baseline, df = 3)" else "age_at_baseline",
    sex = "sex",
    ethnicity5 = "ethnicity5",
    tenure = "tenure",
    household_size = "household_size",
    econstatus = "econstatus",
    education = "education",
    health = "health",
    disability = "disability",
    ruralurban = "ruralurban",
    imd_decile = "imd_decile"
  )
}

get_rhs_terms_for_df <- function(df, use_spline_age = FALSE) {
  term_map <- get_rhs_term_map(use_spline_age = use_spline_age)
  keep <- vapply(
    names(term_map),
    function(v) {
      x <- df[[v]]
      length(unique(x[!is.na(x)])) > 1L
    },
    logical(1)
  )

  terms <- unname(term_map[keep])

  if (length(terms) == 0L) {
    stop("No varying RHS covariates available for model fitting.")
  }

  terms
}

train_cox_spline <- function(df, time_col, status_col, w_col = NULL) {
  require_pkg("survival")
  require_pkg("splines")

  rhs <- paste(get_rhs_terms_for_df(df, use_spline_age = TRUE), collapse = " + ")

  fml <- stats::as.formula(
    paste0("survival::Surv(", time_col, ", ", status_col, ") ~ ", rhs)
  )

  fit_args <- list(
    formula = fml,
    data = df,
    ties = "efron",
    robust = !is.null(w_col),
    x = TRUE,
    model = TRUE
  )

  if (!is.null(w_col)) {
    fit_args$weights <- df[[w_col]]
  }

  fit <- do.call(survival::coxph, fit_args)

  bh <- survival::basehaz(fit, centered = FALSE)

  H0_at <- function(days) {
    # Weighted Breslow cumulative hazard: right-continuous step at event times.
    stats::approx(
      x = bh$time,
      y = bh$hazard,
      xout = days,
      yleft = 0,
      yright = max(bh$hazard),
      method = "constant",
      f = 0
    )$y
  }

  S0_at <- function(days) exp(-H0_at(days))

  list(
    name = "cox_spline",
    fit = fit,
    predict_lp = function(newdf) {
      as.numeric(stats::predict(fit, newdata = newdf, type = "lp", reference = "zero"))
    },
    predict_risk = function(newdf, horizons_days) {
      lp <- as.numeric(stats::predict(fit, newdata = newdf, type = "lp", reference = "zero"))

      out <- sapply(as.numeric(horizons_days), function(tt) {
        1 - (S0_at(tt)^exp(lp))
      })

      out <- as.matrix(out)
      if (length(horizons_days) == 1L) out <- matrix(out, ncol = 1L)
      colnames(out) <- paste0("risk_", horizons_days, "d")
      as.data.frame(out)
    }
  )
}
