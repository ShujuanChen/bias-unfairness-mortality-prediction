#!/usr/bin/env python3
"""Primary weighted deep survival mortality-prediction pipeline."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

# Make the sibling dl_components module importable regardless of CWD.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import dl_components as dl

OUTCOMES = dl.OUTCOMES
TRAINING_SOURCES = dl.TRAINING_SOURCES
FOLDS = dl.FOLDS
HORIZONS_DAYS = dl.HORIZONS_DAYS


# ────────────────────────────────────────────────────────────────────────────
# Data loading and preparation
# ────────────────────────────────────────────────────────────────────────────

def _load_harmonised(framework_root: Path):
    """Load harmonised UKB and PMR data."""
    results_dir = framework_root / "mortality_risk_prediction" / "results"
    ukb = pd.read_csv(results_dir / "UKB_harmonised_with_pmr.csv")
    pmr = pd.read_csv(results_dir / "PMR_harmonised_with_ukb.csv")
    return ukb, pmr


def _load_ukb_weight(framework_root: Path, weight_source_key: str) -> pd.DataFrame:
    """Load a UKB participation-weight file. Returns DataFrame with columns
    [eid, w]. The framework_config.json drives the source path."""
    cfg = dl.read_framework_config(framework_root)
    src = dl.get_weight_source_from_framework(cfg, weight_source_key)
    rel_path = src["path"]
    abs_path = (framework_root / rel_path).resolve()
    if not abs_path.exists():
        raise RuntimeError(f"Missing weight file for {weight_source_key}: {abs_path}")
    wts = pd.read_csv(abs_path)
    wts["eid"] = pd.to_numeric(wts["eid"], errors="coerce")
    wts = wts.dropna(subset=["eid"]).copy()
    wts["eid"] = wts["eid"].astype(int)
    col = src["column"]
    if col not in wts.columns:
        raise RuntimeError(f"Weight column {col!r} not in {abs_path}")
    return wts[["eid", col]].rename(columns={col: "w"})


def prepare_data(framework_root: Path, cause: str):
    """Load + harmonise + outcome-build for one cause.

    Returns dict with:
        ukb_df      training-ready UKB rows (covariates, time/status/event_5y/10y)
        pmr_df     training-ready PMR rows (with PMR design weights as 'w')
        levels      shared categorical levels across UKB and PMR
    """
    cfg = dl.read_framework_config(framework_root)
    outcome_cfg = dl.get_outcome_cfg_from_framework(cfg, cause)

    ukb_raw, pmr_raw = _load_harmonised(framework_root)
    ukb_raw = dl.coerce_covariates(ukb_raw)
    pmr_raw = dl.coerce_covariates(pmr_raw)

    # Restrict to common categorical support so model fitting and one-hot
    # encoding share an identical level set.
    ukb_raw, pmr_raw, _support = dl.filter_to_common_support(
        ukb_raw, pmr_raw, dl.RHS_CATEGORICAL,
    )

    # Build cause-specific outcome on each side.
    ukb_outcome = dl.build_survival_outcome(ukb_raw, "date_of_death", outcome_cfg)
    pmr_outcome = dl.build_survival_outcome(pmr_raw, "dod_deaths", outcome_cfg)

    ukb = pd.concat([ukb_raw.reset_index(drop=True), ukb_outcome.reset_index(drop=True)], axis=1)
    pmr = pd.concat([pmr_raw.reset_index(drop=True), pmr_outcome.reset_index(drop=True)], axis=1)

    # Drop pre-baseline deaths (impossible time_days < 0) just in case.
    ukb = ukb.loc[ukb["time_days"] >= 0].reset_index(drop=True)
    pmr = pmr.loc[pmr["time_days"] >= 0].reset_index(drop=True)

    # Per-horizon event indicators for downstream evaluation.
    for hz in HORIZONS_DAYS:
        suffix = f"{int(round(hz / 365.25))}y"
        for d in (ukb, pmr):
            d[f"event_{suffix}"] = ((d["status"].to_numpy() == 1) &
                                    (d["time_days"].to_numpy(dtype=float) <= float(hz))).astype(int)

    # PMR carries its case-base sampling weight in column 'w'.
    if "w" not in pmr.columns:
        raise RuntimeError("PMR data is missing the 'w' (sampling weight) column.")

    levels = dl.build_combined_levels(ukb, pmr)
    return {"ukb": ukb, "pmr": pmr, "levels": levels, "outcome_cfg": outcome_cfg}


def attach_ukb_weight(ukb_df: pd.DataFrame, framework_root: Path, source: str) -> pd.DataFrame:
    """Attach a participation weight to UKB. Drops rows with missing weight.
    For source == 'ukb', returns UKB unchanged with w = 1."""
    out = ukb_df.copy()
    if source == "ukb":
        out["w"] = 1.0
        return out
    if not source.startswith("ukbw_"):
        raise RuntimeError(f"Unexpected UKB-derived source: {source}")
    weight_key = source[len("ukbw_"):]
    wts = _load_ukb_weight(framework_root, weight_key)
    n_before = len(out)
    out = out.merge(wts, left_on="participant_id", right_on="eid", how="left")
    out["w"] = pd.to_numeric(out["w"], errors="coerce")
    n_dropped = int(out["w"].isna().sum())
    out = out.loc[out["w"].notna()].copy()
    if n_dropped:
        dl.log_step(f"attach_ukb_weight({source}): dropped {n_dropped}/{n_before} rows with missing weight")
    return out


# ────────────────────────────────────────────────────────────────────────────
# Fold assignment management (cached on disk per (source-group, cause))
# ────────────────────────────────────────────────────────────────────────────

def _fold_assignment_path(framework_root: Path, cause: str, group: str) -> Path:
    """Fold assignments are cached at:
        results/evaluation/deep_surv/<cause>/folds/_assignments/<group>.csv

    `group` is 'ukb' (used by all UKB-derived sources, since UKB and wUKB
    share participant IDs) or 'pmr'. Storing the full assignment ensures
    every (UKB, wUKB-HSE-SL, wUKB-HSE-LL, wUKB-Census-SL) training run uses
    the same fold splits per cause.
    """
    base = dl.get_cause_dir(framework_root, cause) / "folds" / "_assignments"
    base.mkdir(parents=True, exist_ok=True)
    return base / f"{group}.csv"


def get_or_build_fold_assignments(framework_root: Path, cause: str, group: str,
                                  ids: pd.Series, status: np.ndarray,
                                  base_seed: int = 11) -> pd.DataFrame:
    """Return a DataFrame with columns [id, fold] for every row in `ids`.

    On first call for a (cause, group), generates a deterministic stratified
    5-fold split and writes it to disk; subsequent calls read from the cache.
    """
    out_path = _fold_assignment_path(framework_root, cause, group)
    if out_path.exists():
        cached = pd.read_csv(out_path)
        if set(cached["id"].astype(int)) != set(np.asarray(ids).astype(int)):
            raise RuntimeError(
                f"Cached fold assignment at {out_path} does not match current id set "
                f"for cause={cause}, group={group}. Delete and regenerate to align."
            )
        return cached
    # Deterministic seed depends on cause + group so different causes get
    # different fold compositions, while different sources within a (cause,
    # group) share splits.
    seed = int(base_seed) + 1000 * (OUTCOMES.index(cause) + 1) + (1 if group == "pmr" else 0)
    folds = dl.stratified_kfold_indices(status=status, k=5, seed=seed)
    fold_id = np.empty(len(ids), dtype=int)
    for k, idx in enumerate(folds):
        fold_id[idx] = k
    out = pd.DataFrame({"id": np.asarray(ids).astype(int), "fold": fold_id})
    out.to_csv(out_path, index=False)
    dl.log_step(f"Wrote fold assignments to {out_path} (seed={seed})")
    return out


# ────────────────────────────────────────────────────────────────────────────
# Per-fold training
# ────────────────────────────────────────────────────────────────────────────

def _resolve_seed(cause: str, source: str, fold: int, base: int = 11) -> int:
    return (int(base)
            + 1000 * (OUTCOMES.index(cause) + 1)
            + 100 * (TRAINING_SOURCES.index(source) + 1)
            + int(fold))


def train_one_fold_task(framework_root: Path, cause: str, source: str, fold: int,
                        base_seed: int = 11):
    """Run one (cause, source, fold) training task."""
    fold = int(fold)
    if fold not in FOLDS:
        raise RuntimeError(f"Invalid fold {fold}; expected one of {FOLDS}")
    if cause not in OUTCOMES:
        raise RuntimeError(f"Invalid cause {cause}; expected one of {OUTCOMES}")
    if source not in TRAINING_SOURCES:
        raise RuntimeError(f"Invalid source {source}; expected one of {TRAINING_SOURCES}")

    seed = _resolve_seed(cause, source, fold, base_seed)
    dl.log_step(f"=== fold task: cause={cause} source={source} fold={fold} seed={seed} ===")

    # One locked architecture (selected once on all-cause mortality) is applied
    # to every outcome; the LR schedule and other mechanics come from
    # TRAINING_SCHEDULE_DEFAULTS.
    locked_spec = dl.load_locked_spec(framework_root)
    dl.log_step(f"loaded locked spec: hidden_dims={locked_spec.get('hidden_dims')} "
                f"dropout={locked_spec.get('dropout')} optimizer={locked_spec.get('optimizer')} "
                f"weight_decay={locked_spec.get('weight_decay')}")

    data = prepare_data(framework_root, cause)
    ukb, pmr, levels = data["ukb"], data["pmr"], data["levels"]

    if source == "pmr":
        train_df_full = pmr.copy()
        id_col = "pmr_id" if "pmr_id" in pmr.columns else None
        if id_col is None:
            train_df_full["pmr_id"] = np.arange(len(train_df_full))
            id_col = "pmr_id"
        # PMR fold assignments are built on PMR's own ID set.
        fold_assignments = get_or_build_fold_assignments(
            framework_root, cause, "pmr",
            ids=train_df_full[id_col].astype(int),
            status=train_df_full["status"].to_numpy(dtype=int),
            base_seed=base_seed,
        )
    else:
        # Build fold assignments on the FULL UKB cohort (the largest UKB-derived
        # row set) so that wUKB sources, which drop rows where the IPW weight
        # is NA, can simply look up their subset without rebuilding the cache.
        # Status here is the cause-specific event indicator on the full UKB
        # cohort; same value for the corresponding row in any wUKB subset.
        get_or_build_fold_assignments(
            framework_root, cause, "ukb",
            ids=ukb["participant_id"].astype(int),
            status=ukb["status"].to_numpy(dtype=int),
            base_seed=base_seed,
        )

        if source == "ukb":
            train_df_full = ukb.copy()
            train_df_full["w"] = 1.0
        else:
            train_df_full = attach_ukb_weight(ukb, framework_root, source)
        id_col = "participant_id"

        # Re-read the (possibly just-built) cache and filter to this source's IDs.
        cache_path = _fold_assignment_path(framework_root, cause, "ukb")
        fold_assignments = pd.read_csv(cache_path)

    ids = train_df_full[id_col].astype(int)
    status = train_df_full["status"].to_numpy(dtype=int)
    fold_lookup = dict(zip(fold_assignments["id"].astype(int), fold_assignments["fold"].astype(int)))
    missing_ids = ids[~ids.isin(fold_lookup)].unique()
    if len(missing_ids) > 0:
        raise RuntimeError(
            f"{len(missing_ids)} ids in source={source} are not in the cached fold "
            f"assignments for group={'pmr' if source == 'pmr' else 'ukb'}; "
            f"the cache was built on a smaller cohort. Delete the assignments "
            f"file and rerun."
        )
    fold_per_row = ids.map(fold_lookup).to_numpy(dtype=int)

    test_idx_full = np.where(fold_per_row == fold)[0]
    train_idx, val_idx, test_idx = dl.stratified_train_val_test_split_within_fold(
        status=status, fold_test_idx=test_idx_full, seed=seed,
    )
    dl.log_step(f"split sizes: train={len(train_idx)} val={len(val_idx)} test={len(test_idx)}")

    # Preprocessing fitted on TRAIN ONLY to avoid val/test leakage.
    preprocessor = dl.fit_preprocessor(train_df_full.iloc[train_idx], levels)
    X_full, _ = preprocessor.transform(train_df_full)
    time_days = train_df_full["time_days"].to_numpy(dtype=float)
    sample_weight = train_df_full["w"].to_numpy(dtype=np.float32)

    train_design = dl.build_cox_design(X_full[train_idx], time_days[train_idx],
                                       status[train_idx], sample_weight[train_idx])
    val_design = (
        dl.build_cox_design(X_full[val_idx], time_days[val_idx],
                            status[val_idx], sample_weight[val_idx])
        if len(val_idx) else None
    )

    input_dim = X_full.shape[1]
    model = dl.make_torch_model(input_dim, spec=locked_spec, seed=seed)

    t0 = time.time()
    model, history, best_val = dl.train_with_lr_schedule(
        model, train_design, val_design, spec=locked_spec, seed=seed,
    )
    elapsed = time.time() - t0
    dl.log_step(f"training complete: best_val_loss={best_val:.6f}, "
                f"elapsed={elapsed:.1f}s")

    # The 80%-trained model at its best validation checkpoint is the deployed
    # predictor for the PMR predictions; the model is not refit on the full
    # data. The weighted Breslow baseline is computed on the 80% train
    # partition to match the training distribution.
    baseline = dl.build_weighted_breslow_baseline(model, train_design)

    n_test = int(len(test_idx))
    n_test_events = int(status[test_idx].sum())
    test_metrics = {
        "fold": int(fold),
        "cause": cause,
        "source": source,
        "n_train": int(len(train_idx)),
        "n_val": int(len(val_idx)),
        "n_test": n_test,
        "n_test_events": n_test_events,
        "best_val_loss": float(best_val),
        "training_seconds": float(elapsed),
        "input_dim": int(input_dim),
        "hidden_dims": list(locked_spec.get("hidden_dims", [])),
        "dropout": float(locked_spec.get("dropout", 0.0)),
        "seed": int(seed),
    }

    # Predict on the FULL PMR cohort (not just the test fold of PMR — every
    # fold-model is applied to all PMR rows; ensemble step averages).
    X_pmr, _ = preprocessor.transform(pmr)
    lp_pmr = dl.predict_lp_torch(model, X_pmr)
    risk_pmr = dl.predict_risk_from_baseline(lp_pmr, baseline, HORIZONS_DAYS)
    pmr_id_arr = (pmr["pmr_id"].astype(int).to_numpy()
                   if "pmr_id" in pmr.columns
                   else np.arange(len(pmr), dtype=int))
    risk_5y_arr = risk_pmr[HORIZONS_DAYS[0]].astype(float)
    risk_10y_arr = risk_pmr[HORIZONS_DAYS[1]].astype(float)

    fold_dir = dl.get_fold_dir(framework_root, cause, source, fold)
    torch.save(model.state_dict(), fold_dir / "model_state.pt")
    with open(fold_dir / "preprocessor.json", "w") as f:
        json.dump(dl.serialise_preprocessor(preprocessor), f, indent=2)
    np.savez_compressed(fold_dir / "baseline_hazard.npz",
                        event_times=baseline["event_times"],
                        cum_baseline_hazard=baseline["cum_baseline_hazard"])
    np.savez_compressed(fold_dir / "pmr_predictions.npz",
                        pmr_id=pmr_id_arr,
                        risk_5y=risk_5y_arr,
                        risk_10y=risk_10y_arr)
    with open(fold_dir / "test_metrics.json", "w") as f:
        json.dump(test_metrics, f, indent=2)
    with open(fold_dir / "history.json", "w") as f:
        json.dump(history, f, indent=2)
    summary_lines = [
        f"cause={cause}",
        f"source={source}",
        f"fold={fold}",
        f"seed={seed}",
        f"n_train={len(train_idx)}",
        f"n_val={len(val_idx)}",
        f"n_test={len(test_idx)}",
        f"input_dim={input_dim}",
        f"hidden_dims={list(locked_spec.get('hidden_dims', []))}",
        f"dropout={float(locked_spec.get('dropout', 0.0))}",
        f"best_val_loss={best_val:.6f}",
        f"training_seconds={elapsed:.1f}",
    ]
    (fold_dir / "summary.txt").write_text("\n".join(summary_lines) + "\n")
    dl.log_step(f"wrote outputs to {fold_dir}")


# ────────────────────────────────────────────────────────────────────────────
# Ensemble averaging across folds
# ────────────────────────────────────────────────────────────────────────────

def ensemble_for_cause(framework_root: Path, cause: str):
    """For each training source, average the 5 fold predictions on PMR and
    write the per-cause individual_predictions.csv."""
    cause_dir = dl.get_cause_dir(framework_root, cause)
    folds_dir = cause_dir / "folds"

    # Load PMR data once for the row-level metadata in the output CSV.
    data = prepare_data(framework_root, cause)
    pmr = data["pmr"]
    if "pmr_id" not in pmr.columns:
        pmr = pmr.copy()
        pmr["pmr_id"] = np.arange(len(pmr), dtype=int)
    base_cols = ["pmr_id", "w", "dod_deaths"] + dl.RHS_COVARIATES + ["RGN11CD",
                  "event_5y", "event_10y"]
    base_cols = [c for c in base_cols if c in pmr.columns]
    out = pmr[base_cols].copy()
    out["outcome"] = cause

    pred_columns = {
        "ukb":                       ("ukb_pred_5y", "ukb_pred_10y"),
        "ukbw_hse_superlearner":     ("ukbw_hse_superlearner_pred_5y",
                                      "ukbw_hse_superlearner_pred_10y"),
        "ukbw_hse_lassologit":       ("ukbw_hse_lassologit_pred_5y",
                                      "ukbw_hse_lassologit_pred_10y"),
        "ukbw_census_superlearner":  ("ukbw_census_superlearner_pred_5y",
                                      "ukbw_census_superlearner_pred_10y"),
        "pmr":                      ("pmr_pred_5y", "pmr_pred_10y"),
    }

    for source in TRAINING_SOURCES:
        risk5 = []
        risk10 = []
        for fold in FOLDS:
            fdir = folds_dir / source / f"fold_{fold}"
            pred_path = fdir / "pmr_predictions.npz"
            if not pred_path.exists():
                raise RuntimeError(f"Missing fold outputs at {fdir}")
            pred_arrs = np.load(pred_path)
            pred = pd.DataFrame({
                "pmr_id": pred_arrs["pmr_id"].astype(int),
                "risk_5y": pred_arrs["risk_5y"].astype(float),
                "risk_10y": pred_arrs["risk_10y"].astype(float),
            }).set_index("pmr_id")
            risk5.append(pred["risk_5y"].astype(float))
            risk10.append(pred["risk_10y"].astype(float))
        risk5_mat = pd.concat(risk5, axis=1)
        risk10_mat = pd.concat(risk10, axis=1)
        risk5_mean = risk5_mat.mean(axis=1)
        risk10_mean = risk10_mat.mean(axis=1)
        col5, col10 = pred_columns[source]
        # Align by pmr_id.
        out_indexed = out.set_index("pmr_id")
        out_indexed[col5] = risk5_mean
        out_indexed[col10] = risk10_mean
        out = out_indexed.reset_index()

    # Write per-cause individual_predictions.csv.
    pred_cols_all = [c for pair in pred_columns.values() for c in pair]
    final_cols = base_cols + ["outcome"] + pred_cols_all
    final_cols = [c for c in final_cols if c in out.columns]
    out = out[final_cols]
    ind_path = cause_dir / "individual_predictions.csv"
    out.to_csv(ind_path, index=False)
    dl.log_step(f"wrote {ind_path} (rows={len(out)})")


# ────────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────────

def _resolve_framework_root(arg: str | None) -> Path:
    if arg:
        return Path(arg).resolve()
    return dl.find_framework_root(Path(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["train", "ensemble"], help="train one fold, or ensemble all folds for a cause")
    parser.add_argument("--cause", required=True, choices=OUTCOMES)
    parser.add_argument("--source", choices=TRAINING_SOURCES,
                        help="required for mode=train")
    parser.add_argument("--fold", type=int, choices=FOLDS,
                        help="required for mode=train")
    parser.add_argument("--seed", type=int, default=11, help="base seed")
    parser.add_argument("--framework-root", default=None,
                        help="override framework_config.json root (default: auto-detect upward)")
    args = parser.parse_args()

    framework_root = _resolve_framework_root(args.framework_root)

    if args.mode == "train":
        if args.source is None or args.fold is None:
            parser.error("--source and --fold are required when mode=train")
        train_one_fold_task(framework_root, args.cause, args.source, args.fold, args.seed)
    elif args.mode == "ensemble":
        ensemble_for_cause(framework_root, args.cause)


if __name__ == "__main__":
    main()
