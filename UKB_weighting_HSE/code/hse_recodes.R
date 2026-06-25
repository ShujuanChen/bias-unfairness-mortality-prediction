# Recode registry for HSE participation weighting.

suppressPackageStartupMessages(library(plyr))

# Map raw codes to labels (codes with no entry are kept); `lvls` optionally
# fixes the factor levels.
relabel_factor <- function(x, mapping, lvls = NULL) {
  relabelled <- plyr::revalue(as.factor(x), mapping, warn_missing = FALSE)
  if (is.null(lvls)) relabelled else factor(relabelled, levels = lvls)
}

# ---- UKB-side recoders (registry keyed by the JSON `recode` field) -----------

.recode_drop_negative <- function(df, col) {
  values <- df[[col]]
  ifelse(values < 0, NA, values)
}

.recode_sex <- function(df, col) {
  relabel_factor(df[[col]], c("0" = "female", "1" = "male"), lvls = c("male", "female"))
}

.recode_education_age <- function(df, col) {
  has_degree <- as.numeric(df[["education_degree"]]) == 1
  years <- ifelse(has_degree, 20, as.numeric(df[[col]]))
  years <- ifelse(years < 0, NA, years)
  years <- ifelse(years <= 14, 14, years)
  ifelse(years >= 19, 19, years)
}

.recode_alcohol <- function(df, col) {
  relabel_factor(df[[col]], c(
    "-3" = NA, "6" = "never", "5" = "few_times_year", "4" = "monthly",
    "3" = "once_twice_weekly", "2" = "three_four_times_weekly", "1" = "daily"
  ))
}

.recode_smoking <- function(df, col) {
  relabel_factor(df[[col]], c("-3" = NA, "0" = "never", "1" = "previous", "2" = "current"))
}

.recode_income <- function(df, col) {
  relabel_factor(df[[col]], c(
    "-1" = "not_shared", "-3" = "not_shared", "1" = "<18k", "2" = "18k-31k",
    "3" = "31k-52k", "4" = "52k-100k", "5" = ">100k"
  ))
}

.recode_household_size <- function(df, col) {
  values <- df[[col]]
  values <- ifelse(values == "-3" | values == "-1", NA, values)
  ifelse(values >= 7, "7", values)
}

.recode_employment <- function(df, col) {
  relabel_factor(df[[col]], c(
    "-7" = NA, "-3" = NA, "1" = "employed", "2" = "retired",
    "3" = "economically_inactive", "4" = "economically_inactive",
    "5" = "unemployed", "6" = "economically_inactive", "7" = "economically_inactive"
  ))
}

.bmi_from_measures <- function(df) round(df$weight / (df$height / 100)^2, 0)

.recode_bmi_value <- function(df, col) .bmi_from_measures(df)

.recode_bmi_category <- function(df, col) {
  bmi <- .bmi_from_measures(df)
  band <- ifelse(bmi < 18.5, "underweight", NA_character_)
  band <- ifelse(bmi >= 18.5 & bmi < 25, "healthyweight", band)
  band <- ifelse(bmi >= 25 & bmi < 30, "overweight", band)
  band <- ifelse(bmi >= 30, "obese", band)
  factor(band, levels = c("underweight", "healthyweight", "overweight", "obese"))
}

.recode_health <- function(df, col) {
  as.factor(plyr::revalue(
    as.factor(df[[col]]),
    c("-1" = NA, "-3" = NA, "1" = "good", "2" = "good", "3" = "fair", "4" = "poor")
  ))
}

.recode_urbanisation <- function(df, col) {
  relabel_factor(
    df[[col]],
    c(
      "1" = "urban", "2" = "town_fringe", "3" = "village_hamlet", "4" = "village_hamlet",
      "5" = "urban", "6" = "town_fringe", "7" = "village_hamlet", "8" = "village_hamlet",
      "9" = NA, "11" = NA, "12" = NA, "13" = NA, "14" = NA, "15" = NA, "16" = NA,
      "17" = NA, "18" = NA
    ),
    lvls = c("village_hamlet", "town_fringe", "urban")
  )
}

UKB_RECODERS <- list(
  drop_negative  = .recode_drop_negative,
  sex            = .recode_sex,
  education_age  = .recode_education_age,
  alcohol        = .recode_alcohol,
  smoking        = .recode_smoking,
  income         = .recode_income,
  household_size = .recode_household_size,
  employment     = .recode_employment,
  bmi_value      = .recode_bmi_value,
  bmi_category   = .recode_bmi_category,
  health         = .recode_health,
  urbanisation   = .recode_urbanisation
)

apply_ukb_recode <- function(recode_name, df, col) {
  recoder <- UKB_RECODERS[[recode_name]]
  if (is.null(recoder)) stop(sprintf("No UKB recoder registered for '%s'.", recode_name))
  recoder(df, col)
}

# ---- HSE-side recoders: align HSE survey categories to the harmonised levels --
# (HSE columns arrive as value-label strings; the explicit factor levels here are
# the ones that survive the HSE-first stack and drive the design matrix.)

hse_recode_sex <- function(x) {
  relabel_factor(
    as.character(x),
    c("Men" = "male", "Women" = "female", "Male" = "male", "Female" = "female",
      "Refused/not obtained" = NA, "Schedule not obtained" = NA,
      "Schedule not applicable" = NA),
    lvls = c("male", "female")
  )
}

hse_recode_health <- function(x) {
  relabel_factor(
    as.character(x),
    c("...very good," = "good", "Good," = "good", "good," = "good",
      "Fair," = "fair", "fair," = "fair", "Bad, or" = "poor", "bad, or" = "poor",
      "very bad?" = "poor", "Very bad?" = "poor", "Very good/good" = "good",
      "Fair" = "fair", "Bad/very bad" = "poor"),
    lvls = c("poor", "fair", "good")
  )
}

hse_recode_smoking <- function(x) {
  relabel_factor(
    as.character(x),
    c("Never smoked cigarettes at all" = "never",
      "Used to smoke cigarettes occasionally" = "previous",
      "Used to smoke cigarettes regularly" = "previous",
      "Current cigarette smoker" = "current"),
    lvls = c("never", "previous", "current")
  )
}

hse_recode_employment <- function(x) {
  relabel_factor(
    as.character(x),
    c("ILO unemployed" = "unemployed", "In employment" = "employed",
      "Other economically inactive" = "economically_inactive", "Retired" = "retired"),
    lvls = c("unemployed", "employed", "economically_inactive", "retired")
  )
}

hse_recode_income <- function(x) {
  pound <- "£"
  bands <- c(
    paste0("<", pound, "520"),
    paste0(pound, c("520<", "1,600<", "2,600<", "3,600<", "5,200<", "7,800<",
                    "10,400<", "13,000<", "15,600<", "18,200<", "20,800<", "23,400<",
                    "26,000<", "28,600<", "31,200<", "33,800<", "36,400<", "41,600<",
                    "46,800<", "52,000<", "60,000<", "70,000<", "80,000<", "90,000<",
                    "100,000<", "110,000<", "120,000<", "130,000<", "140,000<"),
           pound, c("1,600", "2,600", "3,600", "5,200", "7,800", "10,400", "13,000",
                    "15,600", "18,200", "20,800", "23,400", "26,000", "28,600", "31,200",
                    "33,800", "36,400", "41,600", "46,800", "52,000", "60,000", "70,000",
                    "80,000", "90,000", "100,000", "110,000", "120,000", "130,000",
                    "140,000", "150,000")),
    paste0(">=", pound, "150,000"), "Do not know", "Refused"
  )
  groups <- c(rep("<18k", 10), rep("18k-31k", 5), rep("31k-52k", 5),
              rep("52k-100k", 5), rep(">100k", 6), "not_shared", "not_shared")
  relabel_factor(
    as.character(x), stats::setNames(groups, bands),
    lvls = c("not_shared", "<18k", "18k-31k", "31k-52k", "52k-100k", ">100k")
  )
}

hse_recode_alcohol <- function(x) {
  relabel_factor(
    as.character(x),
    c("Not at all in the last 12 months/Non-drinker" = "never",
      "Once or twice a year" = "few_times_year",
      "Once every couple of months" = "few_times_year",
      "Once or twice a month" = "monthly",
      "Once or twice a week" = "once_twice_weekly",
      "Three or four days a week" = "three_four_times_weekly",
      "Almost every day" = "daily", "Five or six days a week" = "daily"),
    lvls = c("never", "few_times_year", "monthly", "once_twice_weekly",
             "three_four_times_weekly", "daily")
  )
}

hse_recode_urbanisation <- function(x) {
  relabel_factor(
    as.character(x),
    c("Hamlet & Isolated Dwelling" = "village_hamlet",
      "Hamlet and Isolated Dwelling - sparse" = "village_hamlet",
      "Town & Fringe - less sparse" = "town_fringe",
      "Town & Fringe - sparse" = "town_fringe",
      "Urban >= 10k - less sparse" = "urban", "Urban >= 10k - sparse" = "urban",
      "Village - less sparse" = "village_hamlet", "Village - sparse" = "village_hamlet",
      "Town & fringe" = "town_fringe", "Urban" = "urban",
      "Village, hamlet and isolated dwellings" = "village_hamlet"),
    lvls = c("village_hamlet", "town_fringe", "urban")
  )
}

hse_recode_education_age <- function(x) {
  chr <- as.character(x)
  years <- rep(NA_real_, length(chr))
  years[chr == "14 or under"] <- 14
  years[chr == "15"] <- 15
  years[chr == "16"] <- 16
  years[chr == "17"] <- 17
  years[chr == "18"] <- 18
  years[chr == "19 or over"] <- 19
  years
}

hse_recode_household_size <- function(x) {
  size <- suppressWarnings(as.numeric(as.character(x)))
  factor(ifelse(is.na(size), NA_character_, ifelse(size >= 7, "7", as.character(size))))
}

hse_recode_bmi_category <- function(bmi) {
  band <- ifelse(bmi < 18.5, "underweight", NA_character_)
  band <- ifelse(bmi >= 18.5 & bmi < 25, "healthyweight", band)
  band <- ifelse(bmi >= 25 & bmi < 30, "overweight", band)
  band <- ifelse(bmi >= 30, "obese", band)
  factor(band, levels = c("underweight", "healthyweight", "overweight", "obese"))
}
