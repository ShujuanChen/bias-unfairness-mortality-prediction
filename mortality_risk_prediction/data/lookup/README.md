# Public geography lookups

Place these public ONS geography lookups here for the Phase 0 harmonisation
(`UKB_prepare.R`, `PMR_prepare.R`):

- `LSOA01_LSOA11_LAD11_Lookup_EW.csv` — LSOA 2001 → LSOA 2011 (+ LAD 2011); used
  to derive `imd_decile`. Read by `UKB/UKB_harmonise.R`, `PMR/PMR_harmonise.R`.
- `Output_Area_to_Lower_layer_Super_Output_Area_to_Middle_layer_Super_Output_Area_to_Local_Authority_District_(December_2011)_Lookup_in_England_and_Wales.csv`
  — Output Area → LSOA → MSOA → LAD (December 2011). Read by
  `PMR/merge_geography.R`, `UKB/merge_geography_UKB.R`.
- `Lower_Layer_Super_Output_Area_(2011)_to_Built-up_Area_Sub-division_to_Built-up_Area_to_Local_Authority_District_to_Region_(December_2011)_Lookup_in_England_and_Wales.csv`
  — LSOA 2011 → Region (December 2011). Read by `PMR/merge_geography.R`,
  `UKB/merge_geography_UKB.R`.

Source: ONS Open Geography Portal (geoportal.statistics.gov.uk).
