# Merge public geography lookups onto the harmonised UK Biobank dataset.

suppressPackageStartupMessages({
  library(dplyr)
})

.ukb_read_lookup_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.ukb_dedupe_lsoa <- function(df, key = "LSOA11CD", cols, name_for_errors = "lookup") {
  stopifnot(key %in% names(df))
  stopifnot(all(cols %in% names(df)))

  out <- df %>%
    transmute(
      LSOA11CD = trimws(.data[[key]]),
      across(all_of(cols), ~ trimws(.x))
    ) %>%
    filter(!is.na(LSOA11CD), LSOA11CD != "") %>%
    distinct()

  dup_keys <- out %>%
    count(LSOA11CD, name = "n") %>%
    filter(n > 1)
  if (nrow(dup_keys) > 0) {
    stop(
      sprintf(
        "UKB geography merge aborted: %s has duplicated LSOA11CD after column selection/distinct (showing first 10): %s",
        name_for_errors,
        paste(utils::head(dup_keys$LSOA11CD, 10), collapse = ", ")
      )
    )
  }

  out
}

merge_geography_ukb <- function(
  ukb,
  oa_to_lsoa_msoa_lad_path = "../../data/lookup/Output_Area_to_Lower_layer_Super_Output_Area_to_Middle_layer_Super_Output_Area_to_Local_Authority_District_(December_2011)_Lookup_in_England_and_Wales.csv",
  lsoa_to_region_path = "../../data/lookup/Lower_Layer_Super_Output_Area_(2011)_to_Built-up_Area_Sub-division_to_Built-up_Area_to_Local_Authority_District_to_Region_(December_2011)_Lookup_in_England_and_Wales.csv",
  verbose = TRUE
) {
  stopifnot(is.data.frame(ukb))
  if (!"LSOA11CD" %in% names(ukb)) stop("merge_geography_ukb(): UKB data is missing `LSOA11CD`.")

  if (isTRUE(verbose)) message("UKB geography: reading OA->LSOA/MSOA/LAD lookup…")
  oa_lookup <- .ukb_read_lookup_csv(oa_to_lsoa_msoa_lad_path)
  needed1 <- c("LSOA11CD", "MSOA11CD", "LAD11CD")
  if (!all(needed1 %in% names(oa_lookup))) {
    stop(
      "UKB OA lookup missing required columns. Need: ",
      paste(needed1, collapse = ", "),
      ". Found: ",
      paste(names(oa_lookup), collapse = ", ")
    )
  }

  oa_lsoa <- .ukb_dedupe_lsoa(
    oa_lookup,
    key = "LSOA11CD",
    cols = c("MSOA11CD", "LAD11CD"),
    name_for_errors = "OA->LSOA/MSOA/LAD lookup"
  )

  if (isTRUE(verbose)) message("UKB geography: reading LSOA->Region lookup…")
  reg_lookup <- .ukb_read_lookup_csv(lsoa_to_region_path)
  needed2 <- c("LSOA11CD", "RGN11CD")
  if (!all(needed2 %in% names(reg_lookup))) {
    stop(
      "UKB Region lookup missing required columns. Need: ",
      paste(needed2, collapse = ", "),
      ". Found: ",
      paste(names(reg_lookup), collapse = ", ")
    )
  }

  lsoa_reg <- .ukb_dedupe_lsoa(
    reg_lookup,
    key = "LSOA11CD",
    cols = c("RGN11CD"),
    name_for_errors = "LSOA->Region lookup"
  )

  out <- ukb %>%
    mutate(LSOA11CD = trimws(as.character(LSOA11CD))) %>%
    left_join(oa_lsoa,  by = "LSOA11CD") %>%
    left_join(lsoa_reg, by = "LSOA11CD")

  if (isTRUE(verbose)) {
    message("UKB geography merge complete: added MSOA11CD, LAD11CD, RGN11CD.")
  }

  out
}

