#!/usr/bin/env python3
"""Append cause-level (whole-cohort) bootstrap CI rows to the bootstrap_ci CSVs."""

import gc
import numpy as np
import pandas as pd
from pathlib import Path
from helpers import (
    CAUSES, get_results_root, get_output_dir, find_individual_predictions,
)

N_BOOTSTRAP = 200
SEED = 42

PRED_KEYS = {
    "pmr":  "pmr_pred_{hz}y",
    "ukb":  "ukb_pred_{hz}y",
    "wukb": "ukbw_hse_superlearner_pred_{hz}y",
}


def _needed_cols(all_cols):
    need = {"w"}
    for hz in [5, 10]:
        need.add(f"event_{hz}y")
        for tmpl in PRED_KEYS.values():
            need.add(tmpl.replace("{hz}", str(hz)))
    return [c for c in all_cols if c in need]


def load_slim(results_root, cause):
    path = find_individual_predictions(results_root, cause)
    if not path or not path.exists():
        return None
    header = pd.read_csv(path, nrows=0).columns.tolist()
    usecols = _needed_cols(header)
    dtypes = {c: "float32" for c in usecols}
    return pd.read_csv(path, usecols=usecols, dtype=dtypes)


def bootstrap_overall(df, hz, rng):
    event_col = f"event_{hz}y"
    if event_col not in df.columns:
        return None

    w   = df["w"].to_numpy(dtype="float64")
    obs = df[event_col].to_numpy(dtype="float64")
    n   = len(df)

    pred_arrays = {}
    for name, tmpl in PRED_KEYS.items():
        col = tmpl.replace("{hz}", str(hz))
        if col in df.columns:
            pred_arrays[name] = df[col].to_numpy(dtype="float64")

    if "wukb" not in pred_arrays:
        return None

    names = list(pred_arrays.keys())
    mat   = np.empty((N_BOOTSTRAP, 1 + len(names)), dtype="float64")

    for b in range(N_BOOTSTRAP):
        idx  = rng.integers(0, n, size=n)
        w_b  = w[idx]
        wsum = w_b.sum()
        mat[b, 0] = (obs[idx] * w_b).sum() / wsum
        for i, name in enumerate(names):
            mat[b, 1 + i] = (pred_arrays[name][idx] * w_b).sum() / wsum

    obs_vals = mat[:, 0]
    row = {
        "cause": None, "horizon_years": hz,
        "var": "overall", "level": "all",
        "n_bootstraps": N_BOOTSTRAP,
        "obs_rate_mean":  float(np.mean(obs_vals)),
        "obs_rate_ci_lo": float(np.percentile(obs_vals, 2.5)),
        "obs_rate_ci_hi": float(np.percentile(obs_vals, 97.5)),
    }
    for i, name in enumerate(names):
        pred_vals = mat[:, 1 + i]
        bias      = (pred_vals - obs_vals) * 1_000_000  # rate difference expressed per million
        row[f"{name}_pred_rate_mean"]  = float(np.mean(pred_vals))
        row[f"{name}_pred_rate_ci_lo"] = float(np.percentile(pred_vals, 2.5))
        row[f"{name}_pred_rate_ci_hi"] = float(np.percentile(pred_vals, 97.5))
        row[f"{name}_bias_per_mil_mean"]  = float(np.mean(bias))
        row[f"{name}_bias_ci_lo_per_mil"] = float(np.percentile(bias, 2.5))
        row[f"{name}_bias_ci_hi_per_mil"] = float(np.percentile(bias, 97.5))

    return row


def main():
    results_root = get_results_root()
    ci_dir       = get_output_dir("bootstrap_ci")
    rng          = np.random.default_rng(SEED)

    for cause in CAUSES:
        print(f"\n[{cause}]")
        df = load_slim(results_root, cause)
        if df is None:
            print("  [SKIP] file not found")
            continue

        for hz in [5, 10]:
            csv_path = ci_dir / f"bootstrap_ci__{cause}__{hz}y.csv"
            if not csv_path.exists():
                print(f"  [SKIP] {hz}y — CSV not found")
                continue

            existing = pd.read_csv(csv_path)
            existing = existing[existing["var"] != "overall"]  # avoid duplicating on re-run

            row = bootstrap_overall(df, hz, rng)
            if row is None:
                print(f"  [SKIP] {hz}y — missing columns")
                continue

            row["cause"] = cause
            new_row = pd.DataFrame([row])
            for col in existing.columns:  # match the strata-level CSV schema before concat
                if col not in new_row.columns:
                    new_row[col] = np.nan
            new_row = new_row[existing.columns]

            out = pd.concat([existing, new_row], ignore_index=True)
            out.to_csv(csv_path, index=False)
            print(f"  Appended overall row to {csv_path.name}")

        del df; gc.collect()

    print("\nDone: bootstrap_ci_overall")


if __name__ == "__main__":
    main()
