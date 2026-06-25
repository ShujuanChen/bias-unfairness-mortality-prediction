#!/usr/bin/env Rscript
# Synthetic data generator: UKB (mortality extract) + PMR (Census-linked mortality).

set.seed(20240619)

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
data_dir <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()

# continuous: standard normal squashed into [lo,hi] (no real mean/sd used)
gen_con <- function(n, lo, hi, int = FALSE) {
  x <- lo + (hi - lo) * stats::plogis(stats::rnorm(n))
  if (int) round(x) else round(x, 1)
}
gen_cat <- function(n, pool) sample(as.character(pool), n, replace = TRUE)
gen_date <- function(n, start, end) {
  s <- as.Date(start); e <- as.Date(end)
  format(s + round(as.numeric(e - s) * stats::plogis(stats::rnorm(n))), "%Y-%m-%d")
}
fake_lsoa <- function(n) sprintf("E01%06d", sample.int(330000, n, replace = TRUE))
icd_pool <- c("C349", "C509", "C189", "C61", "I64", "I219", "I259", "J449", "J189",
              "K219", "K703", "K929", "A419", "N179", "E119", "G309", "F03")

# ---- UKB mortality extract (35 raw columns; ICD/date as free text) -----------
n_ukb <- 5000
ukb_dead <- runif(n_ukb) < 0.2
ukb <- data.frame(
  eid       = seq_len(n_ukb),
  p53_i0    = gen_date(n_ukb, "2006-01-01", "2010-12-31"),
  p31       = gen_cat(n_ukb, c("Male", "Female")),
  p34       = gen_con(n_ukb, 1942, 1971, int = TRUE),
  p52       = gen_cat(n_ukb, month.name),
  p20274_i0 = fake_lsoa(n_ukb),
  p21000_i0 = gen_cat(n_ukb, c("British", "Irish", "Any other white background",
                               "White and Black African", "White and Asian", "Indian",
                               "Pakistani", "African", "Caribbean", "Chinese",
                               "Other ethnic group")),
  p680_i0   = gen_cat(n_ukb, c("Own outright (by you or someone in your household)",
                               "Own with a mortgage",
                               "Rent - from private landlord or letting agency",
                               "Rent - from local authority, local council, housing association",
                               "Pay part rent and part mortgage (shared ownership)",
                               "Live in accommodation rent free", "None of the above")),
  p709_i0   = gen_con(n_ukb, 1, 8, int = TRUE),
  p20119_i0 = gen_cat(n_ukb, c("In paid employment or self-employed", "Retired",
                               "Unemployed", "Looking after home and/or family")),
  p6142_i0  = gen_cat(n_ukb, c("In paid employment or self-employed", "Retired",
                               "Unemployed", "None of the above")),
  p6138_i0  = gen_cat(n_ukb, c("College or University degree",
                               "A levels/AS levels or equivalent",
                               "O levels/GCSEs or equivalent", "None of the above")),
  p845_i0   = gen_cat(n_ukb, c("16", "17", "18", "15", "19", "20")),
  p20118_i0 = gen_cat(n_ukb, c("England/Wales - Urban - less sparse",
                               "England/Wales - Town and Fringe - less sparse",
                               "England/Wales - Village - less sparse",
                               "England/Wales - Hamlet and Isolated Dwelling - less sparse")),
  p2178_i0  = gen_cat(n_ukb, c("Excellent", "Good", "Fair", "Poor")),
  p2188_i0  = gen_cat(n_ukb, c("Yes", "No")),
  p40000_i0 = "",
  p40001_i0 = "",
  check.names = FALSE, stringsAsFactors = FALSE
)
ukb$p40000_i0[ukb_dead] <- gen_date(sum(ukb_dead), "2011-03-27", "2023-02-15")
ukb$p40001_i0[ukb_dead] <- gen_cat(sum(ukb_dead), icd_pool)
for (a in 0:14) ukb[[sprintf("p40002_i0_a%d", a)]] <- ""   # secondary causes: empty
dir.create(file.path(data_dir, "UKB"), showWarnings = FALSE, recursive = TRUE)
write.csv(ukb, file.path(data_dir, "UKB", "UKB_for_harmonisation_with_PMR.csv"),
          row.names = FALSE)

# ---- PMR (case-cohort: deaths weight 1, 5% survivor sample weight 20) ---------
n_pmr <- 5000
pmr_case <- runif(n_pmr) < 0.3
pmr <- data.frame(
  sampling_weight        = ifelse(pmr_case, 1, 20),
  age_census             = gen_con(n_pmr, 40, 69, int = TRUE),
  dod_deaths             = "",
  sex_census             = gen_cat(n_pmr, 1:2),
  sex_deaths             = gen_cat(n_pmr, 1:2),
  LSOA11CD               = fake_lsoa(n_pmr),
  ethnicity              = gen_cat(n_pmr, c("White", "Black", "Mixed", "Indian",
                                            "Bangladeshi and Pakistani", "Chinese", "Other")),
  disability_census      = gen_cat(n_pmr, 1:3),
  health_census          = gen_cat(n_pmr, 1:5),
  tenhuk11_census        = gen_cat(n_pmr, 0:9),
  ecocatpuk11_census     = gen_cat(n_pmr, 1:8),
  hlqpuk11_census        = gen_cat(n_pmr, 10:16),
  hhchuk11_census        = gen_cat(n_pmr, sprintf("%02d", 1:10)),
  ruralurban_code_census = gen_cat(n_pmr, c("A1", "A2", "C1", "C2", "D1", "D2",
                                            "E1", "E2", "F1", "F2")),
  uresindpuk11_census    = "1",
  fic10und_deaths        = "",
  check.names = FALSE, stringsAsFactors = FALSE
)
pmr$dod_deaths[pmr_case]      <- gen_date(sum(pmr_case), "2011-03-27", "2023-02-15")
pmr$fic10und_deaths[pmr_case] <- gen_cat(sum(pmr_case), icd_pool)
for (m in 1:15) pmr[[sprintf("fic10men%d_deaths", m)]] <- ""   # secondary causes: empty
dir.create(file.path(data_dir, "PMR"), showWarnings = FALSE, recursive = TRUE)
write.csv(pmr, file.path(data_dir, "PMR", "PMR.csv"), row.names = FALSE)

cat(sprintf("Synthetic mortality data written: %d UKB rows, %d PMR rows.\n", n_ukb, n_pmr))
