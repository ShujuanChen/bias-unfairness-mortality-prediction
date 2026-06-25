#!/usr/bin/env Rscript
# Synthetic data generator: UKB (HSE-weighting extract) + HSE survey waves.

suppressPackageStartupMessages(library(haven))
set.seed(20240617)

# resolve this script's data/ directory
.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
data_dir <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()

# continuous: standard normal squashed into [lo,hi] (no real mean/sd used)
gen_con <- function(n, lo, hi, int = FALSE) {
  x <- lo + (hi - lo) * stats::plogis(stats::rnorm(n))
  if (int) round(x) else round(x, 1)
}
# categorical (plain): uniform sample from the recognised value pool
gen_cat <- function(n, pool) sample(as.character(pool), n, replace = TRUE)
# categorical for .sav: SPSS value-labelled (labels = recognised strings)
gen_lab <- function(n, pool) {
  labs <- stats::setNames(as.double(seq_along(pool)), pool)
  haven::labelled(as.double(sample(seq_along(pool), n, replace = TRUE)), labels = labs)
}
mk_df <- function(lst, n) { attr(lst, "row.names") <- .set_row_names(n); class(lst) <- "data.frame"; lst }

# ---- UKB extract (raw UKB field names; age special-cased to p21022) ----------
n_ukb <- 5000
england_centres <- c(11012, 11021, 11011, 11008, 11024, 11020, 11018, 11010, 11016,
                     11001, 11017, 11009, 11013, 11002, 11007, 11014, 11006, 10003,
                     11025, 11026, 11027, 11028)
eth_codes <- c(1, 1001, 1002, 1003, 2, 2001, 2002, 2003, 2004, 3, 3001, 3002, 3003,
               3004, 4, 4001, 4002, 4003, 5, 6)
ukb <- data.frame(
  eid       = seq_len(n_ukb),
  p31_i0    = gen_cat(n_ukb, c(0, 1)),
  p21022    = gen_con(n_ukb, 40, 69, int = TRUE),
  p845_i0   = gen_con(n_ukb, 14, 19, int = TRUE),
  p1558_i0  = gen_cat(n_ukb, 1:6),
  p20116_i0 = gen_cat(n_ukb, 0:2),
  p738_i0   = gen_cat(n_ukb, 1:5),
  p709_i0   = gen_cat(n_ukb, 1:6),
  p6142_i0  = gen_cat(n_ukb, c(1, 2, 3, 5)),
  p2178_i0  = gen_cat(n_ukb, 1:4),
  p50_i0    = gen_con(n_ukb, 150, 190),
  p20118_i0 = gen_cat(n_ukb, 1:8),
  p21002_i0 = gen_con(n_ukb, 55, 110),
  p54_i0    = gen_cat(n_ukb, england_centres),
  p21000_i0 = gen_cat(n_ukb, eth_codes),
  p6138_i0  = gen_cat(n_ukb, 1:6),
  check.names = FALSE, stringsAsFactors = FALSE
)
dir.create(file.path(data_dir, "UKB"), showWarnings = FALSE, recursive = TRUE)
write.csv(ukb, file.path(data_dir, "UKB", "UKB_for_harmonisation_with_HSE_participant.csv"),
          row.names = FALSE)

# ---- HSE waves (SPSS .sav, value-labelled categoricals) ----------------------
income_pool <- c("<£520", "£5,200<£7,800", "£18,200<£20,800",
                 "£31,200<£33,800", "£52,000<£60,000",
                 ">=£150,000", "Do not know")
urb8 <- c("Urban >= 10k - less sparse", "Urban >= 10k - sparse",
          "Town & Fringe - less sparse", "Town & Fringe - sparse",
          "Village - less sparse", "Village - sparse",
          "Hamlet & Isolated Dwelling", "Hamlet and Isolated Dwelling - sparse")
urb3 <- c("Urban", "Town & fringe", "Village, hamlet and isolated dwellings")
eth_pool <- c("White", "White - British", "White - Irish", "Any other white background",
              "Mixed", "Any other mixed background", "Asian or Asian British",
              "Black or Black British", "Chinese", "Chinese or other ethnic group",
              "Any other (please describe)")
sex_pool <- c("Male", "Female")
educ_pool <- c("14 or under", "15", "16", "17", "18", "19 or over")
alc_pool <- c("Not at all in the last 12 months/Non-drinker", "Once or twice a month",
              "Once or twice a week", "Almost every day")
smoke_pool <- c("Never smoked cigarettes at all", "Used to smoke cigarettes occasionally",
                "Used to smoke cigarettes regularly", "Current cigarette smoker")
emp_pool <- c("ILO unemployed", "In employment", "Other economically inactive", "Retired")
health_pool <- c("Very good/good", "Fair", "Bad/very bad")

waves <- list(
  list(dir = "UKDA-5809-spss/spss/spss12", file = "hse06ai.sav", alc = "dnoft2", hh = "hhsize",  urb = "URINDEW", urb_pool = urb8, eth = "ethinda"),
  list(dir = "UKDA-6112-spss/spss/spss12", file = "hse07ai.sav", alc = "dnoft2", hh = "hhsizeD", urb = "URINDEW", urb_pool = urb8, eth = "ethinda"),
  list(dir = "UKDA-6397-spss/spss/spss24", file = "hse08ai.sav", alc = "dnoft3", hh = "hhsize",  urb = "Urban",   urb_pool = urb3, eth = "origin"),
  list(dir = "UKDA-6732-spss/spss/spss19", file = "hse09ai.sav", alc = "dnoft3", hh = "hhsize",  urb = "Urban",   urb_pool = urb3, eth = "origin"),
  list(dir = "UKDA-6986-spss/spss/spss19", file = "hse10ai.sav", alc = "dnoft3", hh = "hhsize",  urb = "urban",   urb_pool = urb3, eth = "origin")
)

n_hse <- 2000
for (w in waves) {
  lst <- list()
  lst[["sex"]]     <- gen_lab(n_hse, sex_pool)
  lst[["age"]]     <- gen_con(n_hse, 40, 69, int = TRUE)
  lst[["educend"]] <- gen_lab(n_hse, educ_pool)
  lst[[w$alc]]     <- gen_lab(n_hse, alc_pool)
  lst[["cigst1"]]  <- gen_lab(n_hse, smoke_pool)
  lst[["totinc"]]  <- gen_lab(n_hse, income_pool)
  lst[[w$hh]]      <- gen_con(n_hse, 1, 6, int = TRUE)
  lst[["econact"]] <- gen_lab(n_hse, emp_pool)
  lst[["bmiok"]]   <- gen_con(n_hse, 18, 35)
  lst[["genhelf"]] <- gen_lab(n_hse, health_pool)
  lst[["estht"]]   <- gen_con(n_hse, 150, 190)
  lst[[w$urb]]     <- gen_lab(n_hse, w$urb_pool)
  lst[["estwt"]]   <- gen_con(n_hse, 55, 110)
  lst[["wt_int"]]  <- gen_con(n_hse, 0.3, 5)
  lst[[w$eth]]     <- gen_lab(n_hse, eth_pool)
  df <- mk_df(lst, n_hse)
  outdir <- file.path(data_dir, "HSE", w$dir)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  haven::write_sav(df, file.path(outdir, w$file))
}

cat(sprintf("Synthetic HSE-weighting data written: %d UKB rows, %d HSE rows x 5 waves.\n",
            n_ukb, n_hse))
