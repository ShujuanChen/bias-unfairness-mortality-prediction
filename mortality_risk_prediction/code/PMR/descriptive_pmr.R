#!/usr/bin/env Rscript

# Build the PMR column of the descriptive characteristics table (Table S1).

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
script_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("--file=", "", script_arg[1]), mustWork = TRUE)
} else {
  normalizePath("mortality_risk_prediction/code/PMR/descriptive_pmr.R", mustWork = FALSE)
}
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
pmr_path <- file.path(project_root, "mortality_risk_prediction", "results",
                      "PMR_harmonised_with_ukb.csv")
output_path <- file.path(project_root, "mortality_risk_prediction", "results", "tables",
                         "descriptive_pmr.xlsx")

if (!file.exists(pmr_path)) stop("Missing: ", pmr_path)

# ── Load ──────────────────────────────────────────────────────────────────────
pmr <- read.csv(pmr_path, stringsAsFactors = FALSE)
pmr$w <- as.numeric(pmr$w)

# ── Helpers ───────────────────────────────────────────────────────────────────
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_wcount_wpct <- function(df, var, level, w_col = "w") {
  mask <- !is.na(df[[var]]) & df[[var]] == level
  w_level <- sum(df[[w_col]][mask], na.rm = TRUE)
  w_total <- sum(df[[w_col]][!is.na(df[[var]])], na.rm = TRUE)
  wpct <- 100 * w_level / w_total
  paste0(fmt_n(round(w_level)), " (", sprintf("%.1f", wpct), "%)")
}

bold_rows <- integer(0)

add_total <- function(rows, pmr) {
  rows[[length(rows) + 1]] <- data.frame(
    Variable = "Total N",
    PMR = fmt_n(round(sum(pmr$w, na.rm = TRUE))),
    stringsAsFactors = FALSE
  )
  bold_rows <<- c(bold_rows, length(rows))
  rows
}

add_categorical <- function(rows, pmr, var_name, display_name, levels_order) {
  rows[[length(rows) + 1]] <- data.frame(
    Variable = display_name, PMR = "", stringsAsFactors = FALSE
  )
  bold_rows <<- c(bold_rows, length(rows))
  for (lv in levels_order) {
    rows[[length(rows) + 1]] <- data.frame(
      Variable = paste0("    ", lv),
      PMR = fmt_wcount_wpct(pmr, var_name, lv),
      stringsAsFactors = FALSE
    )
  }
  rows
}

# ── Build table ───────────────────────────────────────────────────────────────
rows <- list()
rows <- add_total(rows, pmr)
rows <- add_categorical(rows, pmr, "disability", "Disability", c("Yes", "No"))
rows <- add_categorical(rows, pmr, "education", "Education",
                        c("Level 1", "Level 2", "Level 3", "Level 4"))
rows <- add_categorical(rows, pmr, "tenure", "Tenure",
                        c("Owned outright",
                          "Owned with a mortgage or loan",
                          "Shared ownership",
                          "Social rented",
                          "Private rented",
                          "Living rent free"))
rows <- add_categorical(rows, pmr, "ruralurban", "Rural-urban classification",
                        c("Urban", "Town and Fringe", "Village",
                          "Hamlet and Isolated Dwelling"))
rows <- add_categorical(rows, pmr, "imd_decile", "IMD decile", as.character(1:10))

# ── Output ────────────────────────────────────────────────────────────────────
out <- do.call(rbind, rows)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()
addWorksheet(wb, "PMR")
writeData(wb, "PMR", out)
bold_style <- createStyle(textDecoration = "bold")
for (r in bold_rows) {
  addStyle(wb, "PMR", bold_style, rows = r + 1, cols = 1:2)
}
setColWidths(wb, "PMR", cols = 1:2, widths = c(35, 25))
saveWorkbook(wb, output_path, overwrite = TRUE)
message("Saved: ", output_path)
