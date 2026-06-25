# Clean and harmonise raw PMR to UKB-aligned covariates; defines pmr_prepare().

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(lubridate)
})

# 5-group ethnicity shared with the UKB side: White, Mixed, Asian, Black, Chinese/Other
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
  message(sprintf("[%-20s] N=%d, NA=%d%s",
                  name, length(x), n_na,
                  if (lev != "") paste0(", levels: ", lev) else ""))
  invisible(NULL)
}

# IMD: read England (EIMD2010) and Wales (WIMD2011) LSOA01 scores, map LSOA01 -> LSOA11
# (unweighted mean when several LSOA01 collapse to one LSOA11), then cut into
# country-specific deciles (1 = most deprived, 10 = least). Missing files -> imd_decile NA.
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

# Treat blank strings and the literal "NA" as missing in ICD columns; keep everything else.
process_pmr_icd_missing <- function(df, cols, verbose = TRUE) {
  stopifnot(is.data.frame(df))
  for (col in cols) {
    if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))
    x <- as.character(df[[col]])
    x[is.na(x)] <- NA_character_
    x_trim <- trimws(x)
    x_trim[x_trim == "" | toupper(x_trim) == "NA"] <- NA_character_
    df[[col]] <- x_trim
    .report(col, df[[col]], verbose)
  }
  df
}

pmr_sampling_weight <- function(df, verbose = TRUE) {
  x <- .as_num(df$sampling_weight)
  .report("sampling_weight", x, verbose)
  x
}

pmr_age_at_baseline <- function(df, verbose = TRUE) {
  x <- .as_num(df$age_census)
  .report("age_at_baseline", x, verbose)
  x
}

pmr_dod_deaths <- function(df, verbose = TRUE) {
  x <- suppressWarnings(as.Date(df$dod_deaths))
  .report("dod_deaths", x, verbose)
  x
}

pmr_sex <- function(df, verbose = TRUE) {
  samp <- as.character(df$sample)
  samp[is.na(samp)] <- NA_character_
  use_death <- !is.na(samp) & samp == "Died"

  c_cen <- .as_num(df$sex_census)
  c_dth <- .as_num(df$sex_deaths)
  code <- ifelse(use_death, c_dth, c_cen)

  lab <- dplyr::case_when(
    code == 1 ~ "Male",
    code == 2 ~ "Female",
    TRUE ~ NA_character_
  )
  f <- factor(lab, levels = c("Male", "Female"))
  .report("sex", f, verbose)
  f
}

# Collapse PMR ethnicity labels to the 5 shared classes; unknown/blank -> NA (no Unknown level).
# Indian and Bangladeshi/Pakistani -> Asian; Chinese and Other -> Chinese/Other.
process_pmr_ethnicity5 <- function(df, col = "ethnicity", out = "ethnicity5") {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  x <- as.character(df[[col]])
  xl <- str_to_lower(str_trim(replace(x, is.na(x), "")))

  out_vec <- case_when(
    xl %in% c("white") ~ "White",
    xl %in% c("black") ~ "Black",
    xl %in% c("mixed") ~ "Mixed",
    xl %in% c("indian", "bangladeshi and pakistani") ~ "Asian",
    xl %in% c("chinese", "other") ~ "Chinese/Other",
    TRUE ~ NA_character_
  )

  df[[out]] <- coerce_harmonised_ethnicity5(out_vec, context = paste0(out, " mapped from ", col))
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

# Tenure (tenhuk11_census) -> 7-class harmonised factor.
# Codes: 0 owned outright; 1 owned with mortgage; 2 shared ownership; 3-4 social rented;
#        5 private rented; 6-8 other; 9 living rent free; X/NA -> NA
process_pmr_tenure <- function(df, col = "tenhuk11_census", out = "tenure", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  code <- .as_num(df[[col]])
  lab <- dplyr::case_when(
    code == 0 ~ "Owned outright",
    code == 1 ~ "Owned with a mortgage or loan",
    code == 2 ~ "Shared ownership",
    code %in% 3:4 ~ "Social rented",
    code == 5 ~ "Private rented",
    code %in% 6:8 ~ "Other",
    code == 9 ~ "Living rent free",
    TRUE ~ NA_character_
  )
  f <- factor(lab, levels = .tenure_levels)
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

# Economic status (ecocatpuk11_census) -> Employed/Unemployed/Retired/Other.
# Codes: 1-4 Employed; 5-6 Unemployed; 7 Retired; 8/X Other; NA -> NA
process_pmr_econstatus <- function(df, col = "ecocatpuk11_census", out = "econstatus", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop("Missing col: ", col)

  chr <- as.character(df[[col]])
  chr[is.na(chr)] <- NA_character_
  chr[trimws(chr) == "NA"] <- NA_character_

  allowed <- c("1", "2", "3", "4", "5", "6", "7", "8", "X", NA_character_)
  bad <- setdiff(unique(chr), allowed)
  if (length(bad) > 0) stop("Unexpected codes in ", col, ": ", paste(sort(bad), collapse = ", "))

  mapped <- dplyr::case_when(
    chr %in% c("1", "2", "3", "4") ~ "Employed",
    chr %in% c("5", "6")           ~ "Unemployed",
    chr == "7"                     ~ "Retired",
    chr %in% c("8", "X")      ~ "Other",
    TRUE                          ~ NA_character_
  )

  f <- factor(mapped, levels = c("Employed", "Unemployed", "Retired", "Other"))
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

# Education (hlqpuk11_census): map each qualification to an equivalent years-of-schooling
# value, then cut into the 4 UKB-aligned levels so both sources share one education scale.
# Codes: 10 no quals; 11-12 L1/L2; 13 apprenticeship; 14 L3; 15 degree; 16 other vocational; XX/NA -> NA
process_pmr_education <- function(df, col = "hlqpuk11_census", out_level = "education", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop("Missing col: ", col)

  chr <- as.character(df[[col]])
  chr[is.na(chr)] <- NA_character_
  chr[trimws(chr) == "NA"] <- NA_character_
  chr_up <- toupper(chr)
  chr_up[chr_up == "NA"] <- NA_character_

  mapped_years <- dplyr::case_when(
    chr == "10"              ~ 7,
    chr %in% c("11", "12")   ~ 10,
    chr == "13"              ~ 12,
    chr == "14"              ~ 13,
    chr == "15"              ~ 20,
    chr == "16"              ~ 15,
    chr_up == "XX" | is.na(chr) ~ NA_real_,
    TRUE                     ~ NA_real_
  )

  mapped_level <- dplyr::case_when(
    is.na(mapped_years)     ~ NA_character_,
    mapped_years <= 8.5     ~ "Level 1",
    mapped_years <= 11       ~ "Level 2",
    mapped_years <= 17.5     ~ "Level 3",
    mapped_years > 17.5      ~ "Level 4",
    TRUE                     ~ NA_character_
  )

  f <- factor(mapped_level, levels = c("Level 1", "Level 2", "Level 3", "Level 4"))
  .report(out_level, f, verbose)
  df[[out_level]] <- f
  if (col != out_level) df <- df %>% select(-any_of(col))
  df
}

# Rural-urban (ruralurban_code_census) -> 4 categories.
# A1/A2/C1/C2 Urban; D1/D2 Town and Fringe; E1/E2 Village; F1/F2 Hamlet and Isolated; NA -> NA
process_pmr_ruralurban <- function(df, col = "ruralurban_code_census", out = "ruralurban", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  chr <- as.character(df[[col]])
  chr[is.na(chr)] <- NA_character_
  chr_trim <- trimws(chr)
  chr_trim[chr_trim == ""] <- NA_character_
  chr_up <- toupper(chr_trim)
  chr_up[chr_up == "NA"] <- NA_character_

  allowed_codes <- c("A1", "A2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2", NA_character_)
  bad <- setdiff(unique(chr_up), allowed_codes)
  if (length(bad) > 0) {
    stop(sprintf("Unexpected codes in '%s': %s", col, paste(shQuote(sort(bad)), collapse = ", ")))
  }

  mapped <- case_when(
    chr_up %in% c("A1", "A2", "C1", "C2") ~ "Urban",
    chr_up %in% c("D1", "D2")             ~ "Town and Fringe",
    chr_up %in% c("E1", "E2")             ~ "Village",
    chr_up %in% c("F1", "F2")             ~ "Hamlet and Isolated Dwelling",
    is.na(chr_up)                         ~ NA_character_,
    TRUE                                  ~ NA_character_
  )

  f <- factor(mapped, levels = c("Urban", "Town and Fringe", "Village", "Hamlet and Isolated Dwelling"))
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

# Household size (hhchuk11_census) -> binary 1 vs 2+.
# Codes: 01-02 -> "1"; 03-26 -> "2+"; XX/NA -> NA
process_pmr_household_size <- function(df, col = "hhchuk11_census", out = "household_size", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  raw <- df[[col]]
  chr <- as.character(raw)
  chr[is.na(raw)] <- NA_character_
  chr_trim <- str_trim(chr)
  chr_trim[chr_trim == ""] <- NA_character_
  chr_trim[toupper(chr_trim) == "NA"] <- NA_character_

  allowed <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14",
               "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "XX", NA_character_)
  norm <- ifelse(nchar(chr_trim) == 1L & chr_trim %in% as.character(1:9), paste0("0", chr_trim), chr_trim)
  norm[is.na(chr_trim)] <- NA_character_
  norm[toupper(chr_trim) == "XX" & !is.na(chr_trim)] <- "XX"

  bad <- setdiff(unique(norm), allowed)
  if (length(bad) > 0) {
    stop(sprintf("Unexpected codes in '%s': %s", col, paste(shQuote(sort(bad)), collapse = ", ")))
  }

  mapped <- dplyr::case_when(
    is.na(norm) | norm == "XX"  ~ NA_character_,
    norm %in% c("01", "02")     ~ "1",
    TRUE                        ~ "2+"
  )

  f <- factor(mapped, levels = c("1", "2+"))
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

# Self-rated health (health_census) -> Good/Fair/Bad.
# Codes: 1-2 Good; 3 Fair; 4-5 Bad; X/NA -> NA
process_pmr_health <- function(df, col = "health_census", out = "health", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  code <- .as_num(df[[col]])
  mapped <- dplyr::case_when(
    code %in% c(1, 2) ~ "Good",
    code == 3 ~ "Fair",
    code %in% c(4, 5) ~ "Bad",
    TRUE ~ NA_character_
  )

  f <- factor(mapped, levels = c("Good", "Fair", "Bad"))
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

# Disability limiting daily activities (disability_census) -> Yes/No.
# Codes: 1-2 Yes; 3 No; X/NA -> NA
process_pmr_disability_binary <- function(df, col = "disability_census", out = "disability_census", verbose = TRUE) {
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(sprintf("Column '%s' not found.", col))

  chr <- as.character(df[[col]])
  chr[is.na(chr)] <- NA_character_
  chr_trim <- trimws(chr)
  chr_trim[chr_trim == ""] <- NA_character_
  chr_up <- toupper(chr_trim)
  chr_up[chr_up == "NA"] <- NA_character_

  allowed <- c("1", "2", "3", "X", NA_character_)
  bad <- setdiff(unique(chr_up), allowed)
  if (length(bad) > 0) {
    stop(sprintf("Unexpected codes in '%s': %s", col, paste(shQuote(sort(bad)), collapse = ", ")))
  }

  mapped <- dplyr::case_when(
    chr_up %in% c("1", "2") ~ "Yes",
    chr_up == "3" ~ "No",
    TRUE ~ NA_character_
  )

  f <- factor(mapped, levels = c("Yes", "No"))
  .report(out, f, verbose)
  df[[out]] <- f
  if (col != out) df <- df %>% select(-any_of(col))
  df
}

# Select the raw PMR fields, run each cleaner exactly once, and return the harmonised dataframe.
pmr_prepare <- function(df, verbose = TRUE) {
  stopifnot(is.data.frame(df))

  # sampling_weight == 1 marks the death subset; pmr_sex() then prefers the death record over the census.
  df$sample <- "Alive"
  df$sample[df$sampling_weight == 1] <- "Died"

  out <- tibble::tibble(
    w = pmr_sampling_weight(df, verbose),
    dod_deaths = pmr_dod_deaths(df, verbose),
    age_at_baseline = pmr_age_at_baseline(df, verbose),
    sex = pmr_sex(df, verbose),
    LSOA11CD = df$LSOA11CD,
    ethnicity = df$ethnicity,
    disability_census = df$disability_census,
    health_census = df$health_census,
    tenhuk11_census = df$tenhuk11_census,
    ecocatpuk11_census = df$ecocatpuk11_census,
    hlqpuk11_census = df$hlqpuk11_census,
    hhchuk11_census = df$hhchuk11_census,
    ruralurban_code_census = df$ruralurban_code_census,
    uresindpuk11_census = df$uresindpuk11_census,
    underlying_cause_of_death_icd = df$fic10und_deaths,
    secondary_cause_of_death_0_icd = df$fic10men1_deaths,
    secondary_cause_of_death_1_icd = df$fic10men2_deaths,
    secondary_cause_of_death_2_icd = df$fic10men3_deaths,
    secondary_cause_of_death_3_icd = df$fic10men4_deaths,
    secondary_cause_of_death_4_icd = df$fic10men5_deaths,
    secondary_cause_of_death_5_icd = df$fic10men6_deaths,
    secondary_cause_of_death_6_icd = df$fic10men7_deaths,
    secondary_cause_of_death_7_icd = df$fic10men8_deaths,
    secondary_cause_of_death_8_icd = df$fic10men9_deaths,
    secondary_cause_of_death_9_icd = df$fic10men10_deaths,
    secondary_cause_of_death_10_icd = df$fic10men11_deaths,
    secondary_cause_of_death_11_icd = df$fic10men12_deaths,
    secondary_cause_of_death_12_icd = df$fic10men13_deaths,
    secondary_cause_of_death_13_icd = df$fic10men14_deaths,
    secondary_cause_of_death_14_icd = df$fic10men15_deaths
  )

  out <- out %>%
    process_pmr_icd_missing(
      cols = c(
        "underlying_cause_of_death_icd",
        paste0("secondary_cause_of_death_", 0:14, "_icd")
      ),
      verbose = verbose
    ) %>%
    process_pmr_ethnicity5(col = "ethnicity", out = "ethnicity5") %>%
    process_pmr_health(col = "health_census", out = "health", verbose = verbose) %>%
    process_pmr_disability_binary(col = "disability_census", out = "disability", verbose = verbose) %>%
    process_pmr_tenure(col = "tenhuk11_census", out = "tenure", verbose = verbose) %>%
    process_pmr_econstatus(col = "ecocatpuk11_census", out = "econstatus", verbose = verbose) %>%
    process_pmr_education(col = "hlqpuk11_census", out_level = "education", verbose = verbose) %>%
    process_pmr_household_size(col = "hhchuk11_census", out = "household_size", verbose = verbose) %>%
    process_pmr_ruralurban(col = "ruralurban_code_census", out = "ruralurban", verbose = verbose) %>%
    process_imd_decile(lsoa11_col = "LSOA11CD", out = "imd_decile", verbose = verbose)

  out
}
