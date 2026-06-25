#!/usr/bin/env python3
"""Per-model-variant summary tables behind the in-text bias and correction numbers."""

import gc
from pathlib import Path
import numpy as np
import pandas as pd

from helpers import (
    CAUSES, CAUSE_LABELS, ALL_CAUSE_KEYS, STRATA_SUBSET, LEVEL_DISPLAY_NAMES,
    get_results_root, get_output_dir, load_cause, compute_rates,
)


STRATA_ALL = [("Overall", ["All"])] + STRATA_SUBSET

WUKB_SOURCES_DEEPSURV = [
    ("HSE SL",    "ukbw_hse_superlearner"),
    ("HSE LL",    "ukbw_hse_lassologit"),
    ("Census SL", "ukbw_census_superlearner"),
]


def _col_to_str(series):
    def _fmt(x):
        try:
            if pd.isna(x):
                return ""
            f = float(x)
            if f == int(f):
                return str(int(f))
            return str(x)
        except (TypeError, ValueError):
            return str(x)
    return series.apply(_fmt)


def _row_metrics(obs_rate, ukb_rate, wukb_rate, n, n_events):
    obs_pm  = obs_rate  * 1_000_000 if np.isfinite(obs_rate)  else np.nan
    ukb_pm  = ukb_rate  * 1_000_000 if np.isfinite(ukb_rate)  else np.nan
    wukb_pm = wukb_rate * 1_000_000 if np.isfinite(wukb_rate) else np.nan
    ukb_abs_bias  = ukb_pm  - obs_pm
    wukb_abs_bias = wukb_pm - obs_pm
    ukb_rel  = (ukb_rate  - obs_rate) / obs_rate * 100 if (np.isfinite(obs_rate) and obs_rate > 0) else np.nan
    wukb_rel = (wukb_rate - obs_rate) / obs_rate * 100 if (np.isfinite(obs_rate) and obs_rate > 0) else np.nan
    ukb_abs_abs  = abs(ukb_abs_bias)  if np.isfinite(ukb_abs_bias)  else np.nan
    wukb_abs_abs = abs(wukb_abs_bias) if np.isfinite(wukb_abs_bias) else np.nan
    bc = ((ukb_abs_abs - wukb_abs_abs) / ukb_abs_abs * 100
          if np.isfinite(ukb_abs_abs) and ukb_abs_abs > 0 and np.isfinite(wukb_abs_abs)
          else np.nan)
    return {
        "n":                      int(n) if np.isfinite(n) else np.nan,
        "n_events":               int(n_events) if np.isfinite(n_events) else np.nan,
        "observed_risk_per_mil":  obs_pm,
        "ukb_pred_risk_per_mil":  ukb_pm,
        "wukb_pred_risk_per_mil": wukb_pm,
        "ukb_abs_bias_per_mil":   ukb_abs_bias,
        "wukb_abs_bias_per_mil":  wukb_abs_bias,
        "ukb_rel_bias_pct":       ukb_rel,
        "wukb_rel_bias_pct":      wukb_rel,
        "bias_correction_pct":    bc,
    }


def _display_level(var, level):
    return LEVEL_DISPLAY_NAMES.get(level, level)


# ── DeepSurv table: rates from individual_predictions ─────────────────────────
def build_deepsurv_table(results_root, wukb_prefix, label):
    rows = []
    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            print(f"[FLAG] {label}: individual_predictions missing for {cause}")
            continue
        for hz in [5, 10]:
            event_col = f"event_{hz}y"
            ukb_col   = f"ukb_pred_{hz}y"
            wukb_col  = f"{wukb_prefix}_pred_{hz}y"
            if event_col not in df.columns or ukb_col not in df.columns:
                print(f"[FLAG] {label}: missing columns in {cause} ({event_col}/{ukb_col})")
                continue
            if wukb_col not in df.columns:
                print(f"[FLAG] {label}: missing wUKB column {wukb_col} for {cause}")

            for var, levels in STRATA_ALL:
                for level in levels:
                    if var == "Overall":
                        sub = df
                    else:
                        sub = df[_col_to_str(df[var]) == level]
                        if len(sub) < 50:
                            rows.append({
                                "cause": cause, "cause_label": CAUSE_LABELS[cause],
                                "horizon_years": hz,
                                "strata_var": var, "strata_level": _display_level(var, level),
                                **_row_metrics(np.nan, np.nan, np.nan, len(sub), np.nan),
                            })
                            continue

                    rates = compute_rates(
                        sub, event_col,
                        {"ukb": ukb_col, "wukb": wukb_col},
                    )
                    metrics = _row_metrics(
                        rates["observed_rate"],
                        rates["ukb_pred_rate"],
                        rates["wukb_pred_rate"],
                        rates["n"], rates["n_events"],
                    )
                    rows.append({
                        "cause": cause, "cause_label": CAUSE_LABELS[cause],
                        "horizon_years": hz,
                        "strata_var": var, "strata_level": _display_level(var, level),
                        **metrics,
                    })
        del df; gc.collect()
    return pd.DataFrame(rows)


# ── Linear-Cox table: rates from pmr_reference_vs_transferred_risk__*.csv ──────
def _linearcox_summary_path(results_root, cause, weight_source="hse_superlearner"):
    if cause in ALL_CAUSE_KEYS:
        d = results_root / "evaluation" / "linear_cox" / cause / "summaries"
    else:
        d = results_root / "evaluation" / "linear_cox" / cause / "underlying" / "summaries"
    return d / f"pmr_reference_vs_transferred_risk__{weight_source}.csv"


def build_linearcox_table(results_root, weight_source="hse_superlearner", label="HSE SL (linear Cox)"):
    rows = []
    wanted = [("Overall", ["All"])] + STRATA_SUBSET
    wanted_set = {(v, lv) for v, levels in wanted for lv in levels}

    for cause in CAUSES:
        path = _linearcox_summary_path(results_root, cause, weight_source)
        if not path.exists():
            print(f"[FLAG] {label}: missing summary file {path}")
            continue
        df = pd.read_csv(path)
        df["strata_variable"] = df["strata_variable"].astype(str)
        df["strata_level"]    = df["strata_level"].astype(str)

        for hz in [5, 10]:
            for var, levels in wanted:
                for level in levels:
                    row = df[(df["strata_variable"] == var) &
                             (df["strata_level"]    == level) &
                             (df["horizon_years"]   == hz)]
                    if row.empty:
                        print(f"[FLAG] {label}: {cause} {hz}y {var}={level} not found")
                        metrics = _row_metrics(np.nan, np.nan, np.nan, np.nan, np.nan)
                        rows.append({
                            "cause": cause, "cause_label": CAUSE_LABELS[cause],
                            "horizon_years": hz,
                            "strata_var": var, "strata_level": _display_level(var, level),
                            **metrics,
                        })
                        continue
                    r = row.iloc[0]
                    obs_rate  = float(r["pmr_obs"])       # PMR observed rate is the reference
                    ukb_rate  = float(r["ukb_pred"])
                    wukb_rate = float(r["ukbw_pred"])
                    metrics = _row_metrics(obs_rate, ukb_rate, wukb_rate, np.nan, np.nan)
                    rows.append({
                        "cause": cause, "cause_label": CAUSE_LABELS[cause],
                        "horizon_years": hz,
                        "strata_var": var, "strata_level": _display_level(var, level),
                        **metrics,
                    })
    return pd.DataFrame(rows)


# ── Write: column order and workbook ──────────────────────────────────────────
COL_ORDER = [
    "cause", "cause_label", "horizon_years", "strata_var", "strata_level",
    "n", "n_events",
    "observed_risk_per_mil", "ukb_pred_risk_per_mil", "wukb_pred_risk_per_mil",
    "ukb_abs_bias_per_mil", "wukb_abs_bias_per_mil",
    "ukb_rel_bias_pct", "wukb_rel_bias_pct",
    "bias_correction_pct",
]


def _order_cols(df):
    cols = [c for c in COL_ORDER if c in df.columns]
    extra = [c for c in df.columns if c not in cols]
    return df[cols + extra]


def main():
    results_root = get_results_root()
    out_dir = get_output_dir("results_table")
    out_path = out_dir / "results_table.xlsx"

    print("Building summary 1: HSE SL (DeepSurv, primary)")
    t1 = _order_cols(build_deepsurv_table(results_root, "ukbw_hse_superlearner", "HSE SL"))
    print("Building summary 2: HSE LL (DeepSurv)")
    t2 = _order_cols(build_deepsurv_table(results_root, "ukbw_hse_lassologit",   "HSE LL"))
    print("Building summary 3: Census SL (DeepSurv)")
    t3 = _order_cols(build_deepsurv_table(results_root, "ukbw_census_superlearner", "Census SL"))
    print("Building summary 4: HSE SL (linear Cox)")
    t4 = _order_cols(build_linearcox_table(results_root, "hse_superlearner", "HSE SL (linear Cox)"))

    with pd.ExcelWriter(out_path, engine="openpyxl") as xl:
        t1.to_excel(xl, sheet_name="1_HSE_SL_deepsurv",    index=False)
        t2.to_excel(xl, sheet_name="2_HSE_LL_deepsurv",    index=False)
        t3.to_excel(xl, sheet_name="3_Census_SL_deepsurv", index=False)
        t4.to_excel(xl, sheet_name="4_HSE_SL_linearcox",   index=False)

    t1.to_csv(out_dir / "summary1_HSE_SL_deepsurv.csv",    index=False)
    t2.to_csv(out_dir / "summary2_HSE_LL_deepsurv.csv",    index=False)
    t3.to_csv(out_dir / "summary3_Census_SL_deepsurv.csv", index=False)
    t4.to_csv(out_dir / "summary4_HSE_SL_linearcox.csv",   index=False)

    print(f"Done — wrote {out_path}")


if __name__ == "__main__":
    main()
