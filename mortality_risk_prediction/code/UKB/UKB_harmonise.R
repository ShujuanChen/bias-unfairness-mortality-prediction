# Clean and harmonise raw UK Biobank records into PMR-aligned covariates.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(lubridate)
})

# Harmonised 5-group ethnicity (White, Mixed, Asian, Black, Chinese/Other)
HARMONISED_ETHNICITY5_LEVELS <- c("White", "Mixed", "Asian", "Black", "Chinese/Other")

validate_harmonised_ethnicity5 <- function(x, context = "ethnicity5") {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- NA_character_
  x_chr[trimws(x_chr) == ""] <- NA_character_
  unexpected <- sort(unique(x_chr[!is.na(x_chr) & !x_chr %in% HARMONISED_ETHNICITY5_LEVELS]))
  if (length(unexpected) > 0L) {
    stop(
      context, " contains unexpected ethnicity5 values: ",
      paste(shQuote(unexpected), collapse = ", "),
      ". Expected only: ", paste(HARMONISED_ETHNICITY5_LEVELS, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

coerce_harmonised_ethnicity5 <- function(x, context = "ethnicity5") {
  validate_harmonised_ethnicity5(x, context = context)
  factor(as.character(x), levels = HARMONISED_ETHNICITY5_LEVELS)
}

# Canonical factor levels (kept identical to PMR_harmonise.R)
.tenure_levels <- c(
  "Owned outright",
  "Owned with a mortgage or loan",
  "Shared ownership",
  "Social rented",
  "Private rented",
  "Living rent free",
  "Other"
)

.as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

.report <- function(name, x, verbose) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  n_na <- sum(is.na(x))
  lev  <- if (is.factor(x)) paste(levels(x), collapse = " | ") else ""
  message(sprintf("[%-28s] N=%d, NA=%d%s",
                  name, length(x), n_na,
                  if (lev != "") paste0(", levels: ", lev) else ""))
  invisible(NULL)
}

# IMD (EIMD2010 + WIMD2011) harmonisation to LSOA11 deciles.
# Maps LSOA01 IMD scores to LSOA11 via lookup (unweighted mean for many-to-one
# changes), then country-specific deciles (1 = most deprived, 10 = least).
process_imd_decile <- function(df,
                               lsoa11_col = "LSOA11CD",
                               out = "imd_decile",
                               eimd_path = "../../data/LSOA_external_features/EIMD2010.xls",
                               wimd_path = "../../data/LSOA_external_features/WIMD2011.xls",
                               lookup_path = "../../data/lookup/LSOA01_LSOA11_LAD11_Lookup_EW.csv",
                               verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!lsoa11_col %in% names(df)) stop(sprintf("Column '%s' not found.", lsoa11_col))
  required_files <- c(eimd_path, wimd_path, lookup_path)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    if (isTRUE(verbose)) {
      warning(
        "process_imd_decile: missing IMD/lookup files; setting imd_decile to NA. Missing: ",
        paste(missing_files, collapse = ", ")
      )
    }
    df[[out]] <- NA_integer_
    .report(out, df[[out]], verbose)
    return(df)
  }

  eimd_raw <- readxl::read_xls(eimd_path, sheet = 2)
  wimd_raw <- readxl::read_xls(wimd_path, sheet = 2)
  lu <- read.csv(lookup_path, stringsAsFactors = FALSE)

  need_lu <- c("LSOA01CD", "LSOA11CD")
  if (!all(need_lu %in% names(lu))) {
    stop("Lookup file must contain columns: ", paste(need_lu, collapse = ", "))
  }

  eimd_cols <- c("LSOA CODE", "IMD SCORE")
  wimd_cols <- c("LSOA Code", "WIMD 2011 score")
  if (!all(eimd_cols %in% names(eimd_raw))) {
    stop("EIMD sheet must contain columns: ", paste(eimd_cols, collapse = ", "))
  }
  if (!all(wimd_cols %in% names(wimd_raw))) {
    stop("WIMD sheet must contain columns: ", paste(wimd_cols, collapse = ", "))
  }

  e01 <- eimd_raw %>%
    transmute(
      lsoa01 = as.character(.data[["LSOA CODE"]]),
      imd_score = suppressWarnings(as.numeric(.data[["IMD SCORE"]])),
      country = "E"
    )
  w01 <- wimd_raw %>%
    transmute(
      lsoa01 = as.character(.data[["LSOA Code"]]),
      imd_score = suppressWarnings(as.numeric(.data[["WIMD 2011 score"]])),
      country = "W"
    )

  imd01 <- bind_rows(e01, w01) %>%
    filter(!is.na(lsoa01), nzchar(lsoa01), !is.na(imd_score))

  pairs <- imd01 %>%
    inner_join(unique(lu[, need_lu]), by = c("lsoa01" = "LSOA01CD")) %>%
    transmute(lsoa11 = as.character(LSOA11CD), imd_score, country)

  agg11 <- pairs %>%
    group_by(lsoa11, country) %>%
    summarise(imd_score = mean(imd_score, na.rm = TRUE), .groups = "drop") %>%
    group_by(country) %>%
    mutate(imd_decile = ntile(desc(imd_score), 10L)) %>%
    ungroup()

  out_map <- agg11 %>% select(lsoa11, imd_decile)
  df[[lsoa11_col]] <- as.character(df[[lsoa11_col]])
  df <- df %>% left_join(out_map, by = setNames("lsoa11", lsoa11_col))
  df[[out]] <- df$imd_decile
  if (out != "imd_decile") df$imd_decile <- NULL

  .report(out, df[[out]], verbose)
  df
}


# Sex -> factor {Male, Female}
process_ukb_sex <- function(df, col = "p31", out = "sex") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  
  vals_chr <- str_trim(as.character(df[[col]]))
  vals_low <- str_to_lower(vals_chr)
  out_vec  <- ifelse(vals_low %in% c("male", "female"),
                     str_to_title(vals_low),
                     NA_character_)
  df[[out]] <- factor(out_vec, levels = c("Male", "Female"))
  df
}


# Age at Census baseline (2011-03-27) from year + month of birth.
# Mid-month (15th) is assumed so the result aligns with PMR age_at_baseline.
process_ukb_age_at_baseline <- function(df,
                                        year_col  = "p34",
                                        month_col = "p52",
                                        out       = "age_at_baseline",
                                        day       = 15L,
                                        census_date = as.Date("2011-03-27")) {
  stopifnot(is.data.frame(df))
  if (!year_col %in% names(df))  stop(sprintf("Column '%s' not found.", year_col))
  if (!month_col %in% names(df)) stop(sprintf("Column '%s' not found.", month_col))
  
  x <- df[[month_col]]
  x_trim <- str_trim(as.character(x))
  x_low  <- str_to_lower(x_trim)
  full   <- str_to_lower(month.name)
  mob <- suppressWarnings(as.integer(match(x_low, full)))
  
  yob <- suppressWarnings(as.integer(df[[year_col]]))
  dob <- ifelse(!is.na(yob) & !is.na(mob),
                as.Date(suppressWarnings(make_date(yob, mob, as.integer(day)))),
                as.Date(NA))
  
  age <- as.numeric(census_date - dob, units = "days") / 365.25
  age[!is.finite(age)] <- NA_real_
  
  df[[out]] <- floor(age)
  df
}


# Ethnicity -> 5 harmonised broad classes
process_ukb_ethnicity5 <- function(df, col = "p21000_i0", out = "ethnicity5") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  
  x0     <- str_trim(replace(df[[col]], is.na(df[[col]]), ""))
  x0_low <- str_to_lower(x0)
  
  unknown <- x0_low %in% c("", "prefer not to answer", "do not know")
  white   <- x0_low %in% c("british", "irish", "white", "any other white background")
  mixed   <- x0_low %in% c("mixed", "any other mixed background",
                           "white and black african", "white and black caribbean", "white and asian")
  asian   <- x0_low %in% c("asian or asian british", "any other asian background",
                           "indian", "pakistani", "bangladeshi")
  black   <- x0_low %in% c("black or black british", "any other black background",
                           "african", "caribbean")
  chinese_other <- x0_low %in% c("chinese", "other ethnic group", "other ethnic group or background")
  
  out_vec <- character(length(x0))
  out_vec[unknown] <- NA_character_
  out_vec[white]   <- "White"
  out_vec[mixed]   <- "Mixed"
  out_vec[asian]   <- "Asian"
  out_vec[black]   <- "Black"
  out_vec[chinese_other] <- "Chinese/Other"

  df[[out]] <- coerce_harmonised_ethnicity5(out_vec, context = paste0(out, " mapped from ", col))
  df
}


# Tenure of household -> canonical 7-class tenure (matches PMR).
# Unexpected raw values error out so coding changes are not silently dropped.
process_ukb_tenure <- function(df, col = "p680_i0", out = "tenure") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  
  raw <- as.character(df[[col]])
  raw[is.na(raw)] <- ""
  
  allowed <- c(
    "Own outright (by you or someone in your household)",
    "Own with a mortgage",
    "Rent - from private landlord or letting agency",
    "Rent - from local authority, local council, housing association",
    "Pay part rent and part mortgage (shared ownership)",
    "Live in accommodation rent free",
    "None of the above",
    "Prefer not to answer",
    ""
  )
  
  bad <- setdiff(unique(raw), allowed)
  if (length(bad) > 0) {
    stop(sprintf(
      "process_ukb_tenure_harmonised(): unexpected values in '%s': %s",
      col, paste(shQuote(sort(bad)), collapse = ", ")
    ))
  }
  
  mapped <- dplyr::case_when(
    raw == "Own outright (by you or someone in your household)" ~ "Owned outright",
    raw == "Own with a mortgage" ~ "Owned with a mortgage or loan",
    raw == "Pay part rent and part mortgage (shared ownership)" ~ "Shared ownership",
    raw == "Rent - from local authority, local council, housing association" ~ "Social rented",
    raw == "Rent - from private landlord or letting agency" ~ "Private rented",
    raw == "Live in accommodation rent free" ~ "Living rent free",
    raw == "None of the above" ~ "Other",
    raw == "" | raw == "Prefer not to answer" ~ NA_character_
  )
  
  df[[out]] <- factor(mapped, levels = .tenure_levels)
  df
}


# Household size -> factor {1, 2+}
process_ukb_household_size <- function(df, col = "p709_i0", out = "household_size") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  
  raw <- df[[col]]
  num <- suppressWarnings(as.numeric(as.character(raw)))
  
  mapped <- dplyr::case_when(
    is.na(num)   ~ NA_character_,
    num == 1     ~ "1",
    num > 1      ~ "2+",
    TRUE         ~ NA_character_
  )
  
  df[[out]] <- factor(mapped, levels = c("1", "2+"))
  df
}


# Economic status -> Employed, Unemployed, Retired, Other
process_ukb_econstatus <- function(df,
                                   col_override = "p20119_i0",
                                   col_multi    = "p6142_i0",
                                   out          = "econstatus") {
  stopifnot(is.data.frame(df))
  if (!col_override %in% names(df)) stop("Missing col: ", col_override)
  if (!col_multi    %in% names(df)) stop("Missing col: ", col_multi)

  # Single-answer override field takes priority over the multi-select field.
  df[[col_override]][df[[col_override]] %in% c("", "Prefer not to answer", "None of the above")] <- NA
  idx_override <- !is.na(df[[col_override]])
  col_multi_vec <- as.character(df[[col_multi]])
  col_multi_vec[idx_override] <- df[[col_override]][idx_override]
  col_multi_vec[is.na(col_multi_vec)] <- ""

  # Pipe-separated multi-select; classify with precedence Employed > Unemployed > Retired > Other.
  tokens_list <- str_split(col_multi_vec, "\\|", simplify = FALSE)
  tokens_list <- lapply(tokens_list, function(v) {
    v <- str_trim(v)
    v <- v[nzchar(v)]
    v
  })
  
  is_missing <- (col_multi_vec == "") |
    str_detect(str_to_lower(col_multi_vec), "prefer not to answer")
  
  has_emp <- vapply(tokens_list, function(tok)
    any(str_to_lower(tok) == "in paid employment or self-employed"), logical(1))
  
  has_un  <- vapply(tokens_list, function(tok)
    any(str_to_lower(tok) == "unemployed"), logical(1))
  
  has_ret <- vapply(tokens_list, function(tok)
    any(str_to_lower(tok) == "retired"), logical(1))
  
  mapped <- character(length(col_multi_vec))
  mapped[] <- NA_character_
  
  mapped[!is_missing & has_emp] <- "Employed"
  mapped[!is_missing & !has_emp & has_un]  <- "Unemployed"
  mapped[!is_missing & !has_emp & !has_un & has_ret] <- "Retired"
  mapped[!is_missing & is.na(mapped)] <- "Other"
  
  df[[out]] <- factor(mapped, levels = c("Employed", "Unemployed", "Retired", "Other"))
  df
}


# Education -> years of schooling, then 4 levels (Level 1-4, matches PMR).
# Qualifications are mapped to years and the maximum is taken across answers.
process_ukb_education <- function(df,
                                  col = "p6138_i0",
                                  age_col = "p845_i0",
                                  out_level = "education") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop("Missing col: ", col)
  if (!age_col %in% names(df)) stop("Missing col: ", age_col)
  
  # 'Never went to school' is tracked separately so it can be assigned 0 years.
  age_raw <- as.character(df[[age_col]])
  never_school <- !is.na(age_raw) & trimws(age_raw) == "Never went to school"

  age_clean <- suppressWarnings(as.integer(age_raw))
  bad_ans <- is.na(age_raw) | age_raw %in% c("", "Do not know", "Prefer not to answer", "Never went to school")
  age_clean[bad_ans] <- NA_integer_

  raw <- as.character(df[[col]])
  raw[is.na(raw)] <- ""
  tokens_list <- str_split(raw, "\\|", simplify = FALSE)
  tokens_list <- lapply(tokens_list, function(v) {
    v <- str_trim(v)
    v <- v[nzchar(v)]
    v
  })
  
  # Canonical tokens (lowercase) expected from UKB
  TOK_DEGREE   <- "college or university degree"
  TOK_ALEVEL   <- "a levels/as levels or equivalent"
  TOK_OLEVEL   <- "o levels/gcses or equivalent"
  TOK_CSE      <- "cses or equivalent"
  TOK_CSE_SING <- "cse or equivalent"
  TOK_NVQ_HND  <- "nvq or hnd or hnc or equivalent"
  TOK_PROF     <- "other professional qualifications eg: nursing, teaching"
  TOK_NONE     <- "none of the above"
  TOK_PNA      <- "prefer not to answer"
  
  # Map each row to years (max across tokens); exact matches only to avoid
  # accidental partial-string hits between qualification names.
  years_list <- lapply(seq_along(tokens_list), function(i) {
    if (never_school[i]) return(0)

    toks <- tolower(tokens_list[[i]])
    if (length(toks) == 0) return(NA_real_)

    yrs <- numeric(0)

    if (any(toks == TOK_NONE))   yrs <- c(yrs, 7)
    if (any(toks == TOK_CSE | toks == TOK_CSE_SING)) yrs <- c(yrs, 10)
    if (any(toks == TOK_OLEVEL)) yrs <- c(yrs, 10)
    if (any(toks == TOK_ALEVEL)) yrs <- c(yrs, 13)
    if (any(toks == TOK_DEGREE)) yrs <- c(yrs, 20)
    
    if (any(toks == TOK_NVQ_HND)) {
      yrs <- c(yrs, ifelse(is.na(age_clean[i]), NA_real_, pmin(age_clean[i] - 5, 19)))
    }
    
    if (any(toks == TOK_PROF)) {
      yrs <- c(yrs, ifelse(is.na(age_clean[i]), NA_real_, pmin(age_clean[i] - 5, 15)))
    }
    
    if (any(toks == "" | toks == TOK_PNA)) yrs <- c(yrs, NA_real_)

    if (!length(yrs) || all(is.na(yrs))) return(NA_real_)
    max(yrs, na.rm = TRUE)
  })

  mapped_years <- unlist(years_list)
  mapped_years[is.infinite(mapped_years)] <- NA_real_

  mapped_level <- dplyr::case_when(
    is.na(mapped_years)     ~ NA_character_,
    mapped_years <= 8.5     ~ "Level 1",
    mapped_years <= 11      ~ "Level 2",
    mapped_years <= 17.5    ~ "Level 3",
    mapped_years > 17.5     ~ "Level 4",
    TRUE                    ~ NA_character_
  )
  
  df[[out_level]] <- factor(mapped_level, levels = c("Level 1","Level 2","Level 3","Level 4"))
  df
}


# Rural/urban -> 4 classes (as in PMR). Scotland codes -> NA.
process_ukb_ruralurban <- function(df, col = "p20118_i0", out = "ruralurban") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  
  raw <- as.character(df[[col]])
  raw[is.na(raw)] <- ""
  norm <- str_squish(str_to_lower(raw))
  is_scotland <- str_detect(norm, "^scotland\\s*-")
  
  allowed <- c(
    "england/wales - urban - sparse",
    "england/wales - urban - less sparse",
    "england/wales - town and fringe - sparse",
    "england/wales - town and fringe - less sparse",
    "england/wales - village - sparse",
    "england/wales - village - less sparse",
    "england/wales - hamlet and isolated dwelling - sparse",
    "england/wales - hamlet and isolated dwelling - less sparse",
    "postcode not linkable",
    ""
  )
  bad <- setdiff(unique(norm[!is_scotland]), allowed)
  if (length(bad) > 0) {
    stop(sprintf("Unexpected values in '%s': %s", col, paste(shQuote(sort(bad)), collapse = ", ")))
  }
  
  mapped <- case_when(
    is_scotland ~ NA_character_,
    norm %in% c("england/wales - urban - sparse", "england/wales - urban - less sparse") ~ "Urban",
    norm %in% c("england/wales - town and fringe - sparse", "england/wales - town and fringe - less sparse") ~ "Town and Fringe",
    norm %in% c("england/wales - village - sparse", "england/wales - village - less sparse") ~ "Village",
    norm %in% c("england/wales - hamlet and isolated dwelling - sparse", "england/wales - hamlet and isolated dwelling - less sparse") ~ "Hamlet and Isolated Dwelling",
    norm %in% c("postcode not linkable", "") ~ NA_character_,
    TRUE ~ NA_character_
  )
  
  df[[out]] <- factor(mapped, levels = c("Urban", "Town and Fringe", "Village", "Hamlet and Isolated Dwelling"))
  df
}


# Overall health rating -> Good / Fair / Bad (matches PMR `health` factor)
process_ukb_health <- function(df, col = "p2178_i0", out = "health") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  
  raw <- as.character(df[[col]])
  raw[is.na(raw)] <- ""
  norm <- str_squish(str_to_lower(raw))
  
  allowed <- c("poor","fair","good","excellent","do not know","prefer not to answer","")
  bad <- setdiff(unique(norm), allowed)
  if (length(bad) > 0) {
    stop(sprintf("Unexpected values in '%s': %s", col, paste(shQuote(sort(bad)), collapse = ", ")))
  }
  
  mapped <- case_when(
    norm == "poor"       ~ "Bad",
    norm == "fair"       ~ "Fair",
    norm %in% c("good","excellent") ~ "Good",
    norm %in% c("do not know","prefer not to answer","") ~ NA_character_,
    TRUE ~ NA_character_
  )
  
  df[[out]] <- factor(mapped, levels = c("Good", "Fair", "Bad"))
  df
}

# Disability limiting daily activities -> binary Yes / No
process_ukb_disability_binary <- function(df, col = "p2188_i0", out = "disability") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- as.character(df[[col]])
  raw[is.na(raw)] <- ""
  norm <- str_squish(str_to_lower(raw))

  allowed <- c("yes", "no", "do not know", "prefer not to answer", "")
  bad <- setdiff(unique(norm), allowed)
  if (length(bad) > 0) {
    stop(sprintf("Unexpected values in '%s': %s", col, paste(shQuote(sort(bad)), collapse = ", ")))
  }

  mapped <- case_when(
    norm == "yes" ~ "Yes",
    norm == "no" ~ "No",
    norm %in% c("do not know", "prefer not to answer", "") ~ NA_character_,
    TRUE ~ NA_character_
  )

  df[[out]] <- factor(mapped, levels = c("Yes", "No"))
  df
}

# Geography: copy UKB LSOA 2011 code into the harmonised LSOA11CD column
process_ukb_lsoa11 <- function(df, col = "p20274_i0", out = "LSOA11CD") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  df[[col]][df[[col]] == ""] <- NA
  df[[out]] <- df[[col]]
  df
}

# Death-cause ICD columns -> ICD code only (first token before space).
# Strips the trailing description text; empty/NA -> NA.
process_ukb_icd_code <- function(df, col, out) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
  raw <- as.character(df[[col]])
  raw[is.na(raw) | trimws(raw) == ""] <- NA_character_
  code <- vapply(str_split(raw, "\\s+", n = 2L), function(x) {
    if (length(x) >= 1L && nzchar(trimws(x[1]))) trimws(x[1]) else NA_character_
  }, character(1))
  code[!nzchar(code)] <- NA_character_
  df[[out]] <- code
  df
}

# Entry point: takes a raw UKB dataframe, keeps only the columns needed for
# harmonisation + death/ICD, applies every harmoniser, and returns a cleaned
# dataframe whose covariates mirror the PMR harmonised set.

# Raw UKB columns consumed by the pipeline (harmonisation + death/ICD only).
# By design this UKB extract carries its own field naming and coding (the p*-style
# field identifiers below); the UKB side is harmonised separately against PMR/HSE/Census,
# so these conventions intentionally differ from the other source pipelines.
.ukb_prepare_cols <- c(
  "eid", "p53_i0", "p31", "p34", "p52", "p20274_i0",
  "p21000_i0", "p680_i0", "p709_i0", "p20119_i0", "p6142_i0",
  "p6138_i0", "p845_i0", "p20118_i0", "p2178_i0", "p2188_i0",
  "p40000_i0", "p40001_i0",
  paste0("p40002_i0_a", 0:14)
)

# Final output columns (no raw p* or eid)
.ukb_prepare_out_cols <- c(
  "participant_id", "assessment_date", "date_of_death",
  "underlying_cause_of_death_icd",
  paste0("secondary_cause_of_death_", 0:14, "_icd"),
  "age_at_baseline", "sex", "ethnicity5", "tenure", "household_size",
  "econstatus", "education", "ruralurban", "health", "disability",
  "LSOA11CD", "imd_decile"
)

ukb_prepare <- function(df, verbose = TRUE) {
  stopifnot(is.data.frame(df))

  needed <- intersect(.ukb_prepare_cols, names(df))
  missing <- setdiff(.ukb_prepare_cols, names(df))
  if (length(missing) > 0L) {
    stop("ukb_prepare: missing required columns: ", paste(missing, collapse = ", "))
  }
  out <- df %>% select(all_of(needed))

  out <- out %>%
    mutate(
      participant_id = eid,
      assessment_date = as.Date(p53_i0),
      date_of_death = as.Date(p40000_i0)
    )

  if (isTRUE(verbose)) message("Processing ICD field: p40001_i0 -> underlying_cause_of_death_icd")
  out <- out %>% process_ukb_icd_code("p40001_i0", "underlying_cause_of_death_icd")
  for (k in 0:14) {
    col_k <- paste0("p40002_i0_a", k)
    out_k <- paste0("secondary_cause_of_death_", k, "_icd")
    if (isTRUE(verbose)) message("Processing ICD field: ", col_k, " -> ", out_k)
    out <- out %>% process_ukb_icd_code(col_k, out_k)
  }
  
  out <- out %>%
    process_ukb_sex(col = "p31", out = "sex") %>%
    process_ukb_age_at_baseline(year_col = "p34",
                                month_col = "p52",
                                out = "age_at_baseline") %>%
    process_ukb_ethnicity5(col = "p21000_i0", out = "ethnicity5") %>%
    process_ukb_tenure(col = "p680_i0", out = "tenure") %>%
    process_ukb_household_size(col = "p709_i0", out = "household_size") %>%
    process_ukb_econstatus(col_override = "p20119_i0",
                           col_multi = "p6142_i0",
                           out = "econstatus") %>%
    process_ukb_education(col = "p6138_i0",
                          age_col = "p845_i0",
                          out_level = "education") %>%
    process_ukb_ruralurban(col = "p20118_i0", out = "ruralurban") %>%
    process_ukb_health(col = "p2178_i0", out = "health") %>%
    process_ukb_disability_binary(col = "p2188_i0", out = "disability") %>%
    process_ukb_lsoa11(col = "p20274_i0", out = "LSOA11CD") %>%
    process_imd_decile(lsoa11_col = "LSOA11CD", out = "imd_decile", verbose = verbose)

  out <- out %>% select(all_of(.ukb_prepare_out_cols))

  if (isTRUE(verbose)) {
    message("ukb_prepare: harmonisation summary")
    for (nm in names(out)) {
      .report(nm, out[[nm]], verbose = TRUE)
    }
  }
  
  tibble::as_tibble(out)
}
