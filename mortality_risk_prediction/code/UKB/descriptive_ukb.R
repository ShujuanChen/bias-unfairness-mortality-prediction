#!/usr/bin/env Rscript

# Descriptive statistics for the harmonised UK Biobank cohort (Table S1).

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

# Paths
script_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("--file=", "", script_arg[1]), mustWork = TRUE)
} else {
  normalizePath("mortality_risk_prediction/code/UKB/descriptive_ukb.R", mustWork = FALSE)
}
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
ukb_path <- file.path(project_root, "mortality_risk_prediction", "results",
                      "UKB_harmonised_with_pmr.csv")
output_path <- file.path(project_root, "mortality_risk_prediction", "results", "tables",
                         "descriptive_ukb.xlsx")

if (!file.exists(ukb_path)) stop("Missing: ", ukb_path)

ukb <- read.csv(ukb_path, stringsAsFactors = FALSE)

# Formatting and row-builder helpers
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_pct <- function(n, total) sprintf("%.1f", 100 * n / total)
fmt_count_pct <- function(n, total) paste0(fmt_n(n), " (", fmt_pct(n, total), "%)")

bold_rows <- integer(0)

add_total <- function(rows, n_ukb) {
  rows[[length(rows) + 1]] <- data.frame(
    Variable = "Total N",
    UKB = fmt_n(n_ukb),
    stringsAsFactors = FALSE
  )
  bold_rows <<- c(bold_rows, length(rows))
  rows
}

add_categorical <- function(rows, ukb, var_name, display_name, levels_order) {
  n_ukb <- nrow(ukb)
  rows[[length(rows) + 1]] <- data.frame(
    Variable = display_name, UKB = "", stringsAsFactors = FALSE
  )
  bold_rows <<- c(bold_rows, length(rows))
  for (lv in levels_order) {
    c_ukb <- sum(ukb[[var_name]] == lv, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      Variable = paste0("    ", lv),
      UKB = fmt_count_pct(c_ukb, n_ukb),
      stringsAsFactors = FALSE
    )
  }
  rows
}

# Assemble variable rows
rows <- list()
rows <- add_total(rows, nrow(ukb))
rows <- add_categorical(rows, ukb, "disability", "Disability", c("Yes", "No"))
rows <- add_categorical(rows, ukb, "education", "Education",
                        c("Level 1", "Level 2", "Level 3", "Level 4"))
rows <- add_categorical(rows, ukb, "tenure", "Tenure",
                        c("Owned outright",
                          "Owned with a mortgage or loan",
                          "Shared ownership",
                          "Social rented",
                          "Private rented",
                          "Living rent free"))
rows <- add_categorical(rows, ukb, "ruralurban", "Rural-urban classification",
                        c("Urban", "Town and Fringe", "Village",
                          "Hamlet and Isolated Dwelling"))
rows <- add_categorical(rows, ukb, "imd_decile", "IMD decile", as.character(1:10))

# Write styled workbook
out <- do.call(rbind, rows)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()
addWorksheet(wb, "UKB")
writeData(wb, "UKB", out)
bold_style <- createStyle(textDecoration = "bold")
for (r in bold_rows) {
  addStyle(wb, "UKB", bold_style, rows = r + 1, cols = 1:2)
}
setColWidths(wb, "UKB", cols = 1:2, widths = c(35, 25))
saveWorkbook(wb, output_path, overwrite = TRUE)
message("Saved: ", output_path)
