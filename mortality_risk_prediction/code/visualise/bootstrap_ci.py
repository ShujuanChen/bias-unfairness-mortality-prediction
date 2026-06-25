#!/usr/bin/env python3
"""Strata-level bootstrap confidence intervals on PMR-observed and predicted rates."""

import gc
import numpy as np
import pandas as pd
from helpers import (
    CAUSES, STRATA_SUBSET,
    get_results_root, get_output_dir, find_individual_predictions, save_csv,
)

N_BOOTSTRAP = 200
SEED = 42

PRED_KEYS = {
    "pmr":  "pmr_pred_{hz}y",
    "ukb":  "ukb_pred_{hz}y",
    "wukb": "ukbw_hse_superlearner_pred_{hz}y",
}

ALL_STRATA = STRATA_SUBSET
STRATA_VARS = list(dict.fromkeys(var for var, _ in ALL_STRATA))


# ── Loading ───────────────────────────────────────────────────────────────────

def _build_usecols(df_columns):
    need = set(STRATA_VARS) | {"w"}
    for hz in [5, 10]:
        need.add(f"event_{hz}y")
        for tmpl in PRED_KEYS.values():
            need.add(tmpl.replace("{hz}", str(hz)))
    return [c for c in df_columns if c in need]


def load_cause_slim(results_root, cause):
    path = find_individual_predictions(results_root, cause)
    if not path or not path.exists():
        return None
    header = pd.read_csv(path, nrows=0).columns.tolist()
    usecols = _build_usecols(header)
    dtypes = {c: "float32" for c in usecols
              if c not in STRATA_VARS and c != "w"}
    dtypes["w"] = "float32"
    df = pd.read_csv(path, usecols=usecols, dtype=dtypes)
    for v in STRATA_VARS:
        if v in df.columns:
            df[v] = df[v].astype("category")
    return df


# ── Bootstrap ─────────────────────────────────────────────────────────────────

def _bootstrap_horizon(w, obs, pred_arrays, strata_masks, strata_meta, n, rng):
    """Run N_BOOTSTRAP resamples; return compact numpy matrix per stratum.

    Returns dict: (var, level) -> array of shape (N_BOOTSTRAP, 1 + n_preds)
                  columns: [obs_rate, pred0_rate, pred1_rate, ...]
    """
    n_preds = len(pred_arrays)
    pred_names = list(pred_arrays.keys())
    pred_arrs  = [pred_arrays[k] for k in pred_names]

    results = {
        meta: np.empty((N_BOOTSTRAP, 1 + n_preds), dtype="float32")
        for meta in strata_meta
    }

    for b in range(N_BOOTSTRAP):
        idx   = rng.integers(0, n, size=n)
        w_b   = w[idx]
        obs_b = obs[idx]

        for (var, level), mask in zip(strata_meta, strata_masks):
            m_b   = mask[idx]
            wsum  = float(w_b[m_b].sum())
            if wsum <= 0 or m_b.sum() < 10:
                results[(var, level)][b, :] = np.nan
                continue
            results[(var, level)][b, 0] = float(
                (obs_b[m_b] * w_b[m_b]).sum() / wsum
            )
            for i, arr in enumerate(pred_arrs):
                results[(var, level)][b, 1 + i] = float(
                    (arr[idx][m_b] * w_b[m_b]).sum() / wsum
                )

    return results, pred_names


def _summarise_ci(boot_results, pred_names, cause, hz):
    rows = []
    for (var, level), mat in boot_results.items():
        valid = ~np.isnan(mat[:, 0])
        if valid.sum() < 10:
            continue
        m = mat[valid]
        obs_vals  = m[:, 0]
        row = {
            "cause": cause, "horizon_years": hz,
            "var": var, "level": level,
            "n_bootstraps": int(valid.sum()),
            "obs_rate_mean":  float(np.mean(obs_vals)),
            "obs_rate_ci_lo": float(np.percentile(obs_vals, 2.5)),
            "obs_rate_ci_hi": float(np.percentile(obs_vals, 97.5)),
        }
        for i, name in enumerate(pred_names):
            pred_vals = m[:, 1 + i]
            bias = (pred_vals - obs_vals) * 1_000_000  # rate difference expressed per million
            row[f"{name}_pred_rate_mean"]  = float(np.mean(pred_vals))
            row[f"{name}_pred_rate_ci_lo"] = float(np.percentile(pred_vals, 2.5))
            row[f"{name}_pred_rate_ci_hi"] = float(np.percentile(pred_vals, 97.5))
            row[f"{name}_bias_per_mil_mean"]  = float(np.mean(bias))
            row[f"{name}_bias_ci_lo_per_mil"] = float(np.percentile(bias, 2.5))
            row[f"{name}_bias_ci_hi_per_mil"] = float(np.percentile(bias, 97.5))
        rows.append(row)
    return pd.DataFrame(rows) if rows else None


# ── Per-cause driver ──────────────────────────────────────────────────────────

def process_cause(df, cause, output_dir, rng):
    col_vals = {var: df[var].astype(str).to_numpy()
                for var, _ in ALL_STRATA if var in df.columns}

    strata_meta, strata_masks = [], []
    for var, levels in ALL_STRATA:
        if var not in df.columns:
            continue
        cv = col_vals[var]
        for level in levels:
            mask = cv == str(level)
            # Permissive floor (>= 10) for resampling stability; the strata figures
            # display only the well-powered cells (>= 50).
            if mask.sum() >= 10:
                strata_meta.append((var, level))
                strata_masks.append(mask)

    w = df["w"].to_numpy(dtype="float32")
    n = len(df)

    for hz in [5, 10]:
        event_col = f"event_{hz}y"
        if event_col not in df.columns:
            continue

        obs = df[event_col].to_numpy(dtype="float32")
        pred_arrays = {}
        for name, tmpl in PRED_KEYS.items():
            col = tmpl.replace("{hz}", str(hz))
            if col in df.columns:
                pred_arrays[name] = df[col].to_numpy(dtype="float32")

        if "wukb" not in pred_arrays:
            print(f"  [SKIP] {cause} {hz}y — wukb column missing")
            continue

        print(f"  Bootstrapping {cause} {hz}y "
              f"({n:,} rows, {len(strata_meta)} strata)...")

        boot_results, pred_names = _bootstrap_horizon(
            w, obs, pred_arrays, strata_masks, strata_meta, n, rng
        )
        result = _summarise_ci(boot_results, pred_names, cause, hz)
        del boot_results

        if result is not None:
            fname = f"bootstrap_ci__{cause}__{hz}y"
            save_csv(result, output_dir, fname)
            print(f"    Saved {fname}.csv  ({len(result)} strata rows)")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    results_root = get_results_root()
    output_dir   = get_output_dir("bootstrap_ci")
    rng          = np.random.default_rng(SEED)

    for cause in CAUSES:
        print(f"\n[{cause}]")
        df = load_cause_slim(results_root, cause)
        if df is None:
            print("  [SKIP] individual predictions file not found")
            continue
        process_cause(df, cause, output_dir, rng)
        del df
        gc.collect()

    print("\nDone: bootstrap_ci")


if __name__ == "__main__":
    main()
