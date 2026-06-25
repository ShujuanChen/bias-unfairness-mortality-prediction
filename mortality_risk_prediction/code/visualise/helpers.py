"""Shared constants, IO, and plot helpers for the visualisation scripts."""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Constants ────────────────────────────────────────────────────────────────

CAUSES = [
    "all_cause_mortality",
    "cancer_mortality",
    "cardiovascular_mortality",
    "digestive_mortality",
    "respiratory_mortality",
]
CAUSE_LABELS = {
    "all_cause_mortality": "All-cause",
    "cancer_mortality": "Cancer",
    "cardiovascular_mortality": "Cardiovascular",
    "digestive_mortality": "Digestive",
    "respiratory_mortality": "Respiratory",
}
ALL_CAUSE_KEYS = {"all_cause_mortality"}

COLOURS = {
    "pmr": "#888888",
    "ukb": "#640404",
    "ukbw_hse_superlearner": "#2C5D8A",
    "ukbw_hse_lassologit": "#D2691E",
    "ukbw_census_superlearner": "#D2691E",
}

LEVEL_DISPLAY_NAMES = {
    "Owned with a mortgage or loan": "Owned with mortgage/loan",
    "Hamlet and Isolated Dwelling": "Hamlet/Isolated Dwelling",
}

C_PMR = "#888888"
C_UKB = "#640404"
C_SENS_UKB_ALT = "#bd3c3c"
C_SENS_PMR_ALT = "#2F4F4F"

MARKERS = {
    "pmr": "D",
    "ukb": "s",
    "ukbw_hse_superlearner": "o",
    "ukbw_hse_lassologit": "^",
    "ukbw_census_superlearner": "^",
}

# Strata shown in the figures
STRATA_SUBSET = [
    ("disability", ["Yes", "No"]),
    ("education", ["Level 1", "Level 2", "Level 3", "Level 4"]),
    ("tenure", ["Owned outright", "Owned with a mortgage or loan",
                "Shared ownership", "Private rented", "Living rent free", "Social rented"]),
    ("ruralurban", ["Urban", "Town and Fringe", "Village", "Hamlet and Isolated Dwelling"]),
    ("imd_decile", [str(i) for i in range(1, 11)]),
]


# ── Path helpers ─────────────────────────────────────────────────────────────

def get_results_root():
    return Path(__file__).resolve().parent.parent.parent / "results"


def find_individual_predictions(results_root, cause):
    """Path to the individual-predictions CSV:
    evaluation/deep_surv/{cause}/individual_predictions.csv
    """
    return results_root / "evaluation" / "deep_surv" / cause / "individual_predictions.csv"


def get_output_dir(name):
    """Resolve the output directory for a visualisation.

    Honours the optional `VISUALISE_OUTPUT_PREFIX` env var so a single
    run can divert all visualisation outputs into a sibling subdir
    (e.g., `results/visualise/new/<name>/`) without overwriting the
    canonical results. With the env var unset the behaviour is
    unchanged (`results/visualise/<name>/`).
    """
    prefix = os.environ.get("VISUALISE_OUTPUT_PREFIX", "").strip()
    base = get_results_root() / "visualise"
    d = (base / prefix / name) if prefix else (base / name)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ── Data helpers ─────────────────────────────────────────────────────────────

def load_cause(results_root, cause):
    path = find_individual_predictions(results_root, cause)
    if not path.exists():
        return None
    return pd.read_csv(path)


def weighted_mean(vals, w):
    ok = np.isfinite(vals) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return np.nan
    return float(np.sum(vals[ok] * w[ok]) / np.sum(w[ok]))


def compute_rates(df, event_col, pred_cols, w_col="w"):
    """Compute observed rate and predicted rates (weighted)."""
    obs = df[event_col].to_numpy(dtype=float)
    w = df[w_col].to_numpy(dtype=float) if w_col in df.columns else np.ones(len(df))
    obs_rate = weighted_mean(obs, w)
    result = {"observed_rate": obs_rate, "n": len(df), "n_events": int(np.nansum(obs))}
    for name, col in pred_cols.items():
        if col in df.columns:
            pred = df[col].to_numpy(dtype=float)
            result[f"{name}_pred_rate"] = weighted_mean(pred, w)
            result[f"{name}_relative_bias"] = (result[f"{name}_pred_rate"] - obs_rate) / obs_rate if obs_rate > 0 else np.nan
        else:
            result[f"{name}_pred_rate"] = np.nan
            result[f"{name}_relative_bias"] = np.nan
    return result


def build_strata_labels(groups):
    """Build ordered label list with empty-line placeholders between variable groups."""
    labels = []
    for gi, (var, levels) in enumerate(groups):
        if gi > 0:
            labels.append(("__spacer__", ""))
        for level in levels:
            labels.append((var, level))
    return labels


def weighted_percentile_edges(pred, w, n_bins):
    sort_idx = np.argsort(pred)
    pred_sorted = pred[sort_idx]
    w_sorted = w[sort_idx]
    cum_w = np.cumsum(w_sorted)
    cum_w_frac = cum_w / cum_w[-1]
    target_fracs = np.linspace(0, 1, n_bins + 1)
    edges = np.interp(target_fracs, cum_w_frac, pred_sorted)
    edges[-1] += 1e-10
    return edges


# ── Save helpers ─────────────────────────────────────────────────────────────

def save_fig(fig, output_dir, fname_base):
    """Save as both PNG and SVG."""
    fig.savefig(
        output_dir / f"{fname_base}.png",
        dpi=300,
        bbox_inches="tight",
        transparent=True,
    )
    fig.savefig(
        output_dir / f"{fname_base}.svg",
        bbox_inches="tight",
        transparent=True,
    )
    plt.close(fig)


def save_csv(df, output_dir, fname_base):
    df_out = df.copy()
    for col in df_out.columns:
        if pd.api.types.is_float_dtype(df_out[col]):
            df_out[col] = df_out[col].apply(
                lambda x: "" if pd.isna(x) else f"{x:.15g}"
            )
    df_out.to_csv(output_dir / f"{fname_base}.csv", index=False)
