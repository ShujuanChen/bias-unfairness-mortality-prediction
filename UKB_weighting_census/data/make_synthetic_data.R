#!/usr/bin/env Rscript
# Synthetic data generator: UKB (Census-weighting extract) + 2011 Census 5% microdata.

set.seed(20240618)

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
data_dir <- if (length(.fa)) dirname(normalizePath(.fa)) else getwd()

# continuous: standard normal squashed into [lo,hi] (no real mean/sd used)
gen_con <- function(n, lo, hi, int = FALSE) {
  x <- lo + (hi - lo) * stats::plogis(stats::rnorm(n))
  if (int) round(x) else round(x, 1)
}
gen_cat <- function(n, pool) sample(as.character(pool), n, replace = TRUE)

# ---- UKB extract (raw UKB field-instance names '<field>-0.0') -----------------
n_ukb <- 5000
england_centres <- c(11012, 11021, 11011, 11008, 11024, 11020, 11018, 11010, 11016,
                     11001, 11017, 11009, 11013, 11002, 11007, 11014, 11006,
                     11025, 11026, 11027, 11028)
eth_codes <- c(1, 2, 3, 4, 5, 6, 1001, 1002, 1003, 2001, 2002, 2003, 2004,
               3001, 3002, 3003, 3004, 4001, 4002, 4003)
ukb <- data.frame(
  eid          = seq_len(n_ukb),
  `31-0.0`     = gen_cat(n_ukb, c(0, 1)),
  `34-0.0`     = gen_con(n_ukb, 1942, 1971, int = TRUE),
  `54-0.0`     = gen_cat(n_ukb, england_centres),
  `21022-0.0`  = gen_con(n_ukb, 40, 69, int = TRUE),
  `6138-0.0`   = gen_cat(n_ukb, 1:4),
  `6142-0.0`   = gen_cat(n_ukb, 1:7),
  `20119-0.0`  = gen_cat(n_ukb, 1:7),
  `680-0.0`    = gen_cat(n_ukb, 1:6),
  `728-0.0`    = gen_con(n_ukb, 1, 5, int = TRUE),
  `2178-0.0`   = gen_cat(n_ukb, 1:4),
  `709-0.0`    = gen_con(n_ukb, 1, 6, int = TRUE),
  `21000-0.0`  = gen_cat(n_ukb, eth_codes),
  `845-0.0`    = gen_con(n_ukb, 10, 20, int = TRUE),
  check.names = FALSE, stringsAsFactors = FALSE
)
dir.create(file.path(data_dir, "UKB"), showWarnings = FALSE, recursive = TRUE)
write.csv(ukb, file.path(data_dir, "UKB", "UKB_for_harmonisation_with_census.csv"),
          row.names = FALSE)

# ---- 2011 Census 5% microdata (recodev12.csv) --------------------------------
n_cen <- 8000
census <- data.frame(
  country     = gen_cat(n_cen, 1),               # 1 = England (only England kept)
  ageh        = gen_cat(n_cen, 9:14),            # 9..14 -> age bands 40-69
  carsnoc     = gen_con(n_cen, 0, 4, int = TRUE),
  ecopuk11    = gen_cat(n_cen, 1:14),
  ethnicityew = gen_cat(n_cen, 1:18),
  health      = gen_cat(n_cen, 1:5),
  hlqupuk11   = gen_cat(n_cen, 10:16),
  meighuk11   = gen_cat(n_cen, 1:8),
  sex         = gen_cat(n_cen, 1:2),
  tenure      = gen_cat(n_cen, 1:5),
  check.names = FALSE, stringsAsFactors = FALSE
)
dir.create(file.path(data_dir, "census"), showWarnings = FALSE, recursive = TRUE)
write.csv(census, file.path(data_dir, "census", "recodev12.csv"), row.names = FALSE)

cat(sprintf("Synthetic Census-weighting data written: %d UKB rows, %d Census rows.\n",
            n_ukb, n_cen))
