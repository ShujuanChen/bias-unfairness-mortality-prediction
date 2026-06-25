#!/usr/bin/env python3
"""Building blocks for the deep survival mortality-prediction pipeline."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn


# ────────────────────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────────────────────

T0 = np.datetime64("2011-03-27")
T_ADMIN_END = np.datetime64("2023-02-15")

HORIZONS_YEARS = [5, 10]
HORIZONS_DAYS = [int(round(365.25 * x)) for x in HORIZONS_YEARS]

RHS_NUMERIC = ["age_at_baseline"]
RHS_CATEGORICAL = [
    "sex",
    "ethnicity5",
    "tenure",
    "household_size",
    "econstatus",
    "education",
    "ruralurban",
    "health",
    "disability",
    "imd_decile",
]
RHS_COVARIATES = RHS_NUMERIC + RHS_CATEGORICAL

CATEGORY_LEVEL_ORDERS = {
    "sex": ["Male", "Female"],
    "ethnicity5": ["White", "Mixed", "Asian", "Black", "Chinese/Other"],
    "tenure": [
        "Owned outright",
        "Owned with a mortgage or loan",
        "Shared ownership",
        "Social rented",
        "Private rented",
        "Living rent free",
        "Other",
    ],
    "household_size": ["1", "2+"],
    "econstatus": ["Employed", "Unemployed", "Retired", "Other"],
    "education": ["Level 1", "Level 2", "Level 3", "Level 4"],
    "ruralurban": ["Urban", "Town and Fringe", "Village", "Hamlet and Isolated Dwelling"],
    "health": ["Good", "Fair", "Bad"],
    "disability": ["Yes", "No"],
    "imd_decile": [str(x) for x in range(1, 11)],
}

OUTCOMES = [
    "all_cause_mortality",
    "cancer_mortality",
    "cardiovascular_mortality",
    "respiratory_mortality",
    "digestive_mortality",
]

# Training-source configurations: UKB unweighted, three weighted-UKB variants,
# and PMR with case-base design weights.
TRAINING_SOURCES = [
    "ukb",
    "ukbw_hse_superlearner",
    "ukbw_hse_lassologit",
    "ukbw_census_superlearner",
    "pmr",
]

FOLDS = list(range(5))

# Network architecture and dropout come from a single locked-spec JSON
# (deep_surv_specs/final_locked_v1.json), selected once on all-cause mortality
# and applied to every outcome. The values below are the training defaults that
# override or extend the locked spec, capturing decisions the spec does not:
#   - LR schedule decays 1e-2 -> 1e-3 -> 1e-4 (x1/10 per stage), 80 epochs and
#     patience 8 per stage; this schedule is defined here, not in the spec
#     (only hidden_dims/dropout/seed are read from the spec).
#   - The best validation checkpoint per fold scores both the held-out test fold
#     and PMR directly; the model is not refit on the full data.
#   - Full-batch gradient descent (one step per epoch over the whole training
#     partition); prediction is mini-batched here, not driven by the spec.
TRAINING_SCHEDULE_DEFAULTS = {
    "optimizer": "adamw",
    "weight_decay": 1e-5,
    "ties": "efron",
    "baseline_estimator": "weighted_breslow",
    "categorical_encoding": "one_hot",
    "input_dropout": 0.0,
    "base_learning_rate": 1e-2,
    "lr_floor": 1e-4,
    "lr_decay_factor": 10.0,
    "epochs_per_stage": 80,
    "early_stop_patience": 8,
    "grad_clip_max_norm": 5.0,
    "seed": 11,
}


# ────────────────────────────────────────────────────────────────────────────
# Logging
# ────────────────────────────────────────────────────────────────────────────

def log_step(*parts):
    from datetime import datetime
    print(f"[{datetime.now().isoformat(timespec='seconds')}] " + "".join(str(p) for p in parts), flush=True)


# ────────────────────────────────────────────────────────────────────────────
# Coercion
# ────────────────────────────────────────────────────────────────────────────

def coerce_text_series(x):
    s = x.astype("string").str.strip()
    s = s.replace({"": pd.NA, "NA": pd.NA, "nan": pd.NA})
    return s


def coerce_covariates(df):
    d = df.copy()
    d["age_at_baseline"] = pd.to_numeric(d["age_at_baseline"], errors="coerce")
    for col in RHS_CATEGORICAL:
        if col == "imd_decile":
            v = pd.to_numeric(d[col], errors="coerce")
            d[col] = v.round().astype("Int64").astype("string")
            d[col] = d[col].replace({"<NA>": pd.NA})
        else:
            d[col] = coerce_text_series(d[col])
    if "RGN11CD" in d.columns:
        d["RGN11CD"] = coerce_text_series(d["RGN11CD"])
    return d


# ────────────────────────────────────────────────────────────────────────────
# ICD matching for outcome construction
# ────────────────────────────────────────────────────────────────────────────

def normalise_icd_code(x):
    s = pd.Series(x, copy=False).astype("string").str.upper().str.strip()
    s = s.replace({"": pd.NA, "NA": pd.NA})
    s = s.str.replace(r"[^A-Z0-9]", "", regex=True)
    s = s.replace({"": pd.NA})
    return s


def normalise_icd_spec(x):
    s = str(x).strip().upper()
    s = "".join(ch for ch in s if ch.isalnum() or ch == "-")
    if not s:
        return None
    return s


def expand_icd_specs_for_matching(specs):
    specs = [normalise_icd_spec(x) for x in specs]
    specs = [x for x in specs if x]
    out = {"icd10_root3": set(), "icd10_exact4": set()}

    for spec in specs:
        if len(spec) == 3 and spec[0].isalpha() and spec[1:].isdigit():
            out["icd10_root3"].add(spec)
            continue
        if len(spec) == 7 and spec[0].isalpha() and spec[3] == "-" and spec[4].isalpha():
            left, right = spec.split("-")
            if left[0] != right[0]:
                raise RuntimeError(f"ICD-10 root ranges must share the same letter prefix: {spec}")
            lo = int(left[1:]); hi = int(right[1:])
            if lo > hi:
                raise RuntimeError(f"Invalid ICD-10 range: {spec}")
            for value in range(lo, hi + 1):
                out["icd10_root3"].add(f"{left[0]}{value:02d}")
            continue
        if len(spec) == 4 and spec[0].isalpha() and spec[1:].isdigit():
            out["icd10_exact4"].add(spec)
            continue
        if len(spec) == 9 and spec[0].isalpha() and spec[4] == "-" and spec[5].isalpha():
            left, right = spec.split("-")
            if left[0] != right[0]:
                raise RuntimeError(f"ICD-10 exact ranges must share the same letter prefix: {spec}")
            lo = int(left[1:]); hi = int(right[1:])
            if lo > hi:
                raise RuntimeError(f"Invalid ICD-10 exact range: {spec}")
            for value in range(lo, hi + 1):
                out["icd10_exact4"].add(f"{left[0]}{value:03d}")
            continue
        raise RuntimeError(f"Unsupported ICD spec: {spec}")

    return out


def icd_vector_matches_plan(codes, plan):
    norm = normalise_icd_code(codes)
    out = pd.Series(False, index=norm.index)

    if plan["icd10_root3"]:
        idx = norm.str.match(r"^[A-Z][0-9]{2}", na=False)
        out.loc[idx] = out.loc[idx] | norm.loc[idx].str.slice(0, 3).isin(plan["icd10_root3"])
    if plan["icd10_exact4"]:
        idx = norm.str.match(r"^[A-Z][0-9]{3}", na=False)
        out.loc[idx] = out.loc[idx] | norm.loc[idx].str.slice(0, 4).isin(plan["icd10_exact4"])
    return out.fillna(False)


def get_harmonised_death_icd_columns(df, match_scope="underlying"):
    underlying = "underlying_cause_of_death_icd"
    if underlying not in df.columns:
        raise RuntimeError(f"Missing underlying cause-of-death column: {underlying}")
    if match_scope != "underlying":
        raise RuntimeError(f"Only 'underlying' match scope is supported in this pipeline; got {match_scope!r}")
    return [underlying]


def build_outcome_death_match(df, outcome_cfg):
    if outcome_cfg["type"] == "all_cause":
        return pd.Series(True, index=df.index)
    plan = expand_icd_specs_for_matching(list(outcome_cfg.get("icd10", [])))
    matched = pd.Series(False, index=df.index)
    for col in get_harmonised_death_icd_columns(df, "underlying"):
        matched = matched | icd_vector_matches_plan(df[col], plan)
    return matched


def build_survival_outcome(df, death_date_col, outcome_cfg):
    death_date = pd.to_datetime(df[death_date_col], errors="coerce").values.astype("datetime64[D]")
    matched_target = build_outcome_death_match(df, outcome_cfg).to_numpy(dtype=bool)
    death_any_in_followup = (~pd.isna(death_date)) & (death_date <= T_ADMIN_END)
    target_death_in_followup = death_any_in_followup & matched_target

    censor_date = np.full(len(df), T_ADMIN_END, dtype="datetime64[D]")
    has_date = ~pd.isna(death_date)
    censor_date[has_date] = np.minimum(death_date[has_date], T_ADMIN_END)

    time_days = (censor_date - T0).astype("timedelta64[D]").astype(int)

    return pd.DataFrame(
        {
            "censor_date": censor_date.astype("datetime64[D]"),
            "time_days": time_days.astype(int),
            "status": target_death_in_followup.astype(int),
            "death_any_in_followup": death_any_in_followup.astype(int),
            "target_death_in_followup": target_death_in_followup.astype(int),
        },
        index=df.index,
    )


# ────────────────────────────────────────────────────────────────────────────
# Common-support harmonisation
# ────────────────────────────────────────────────────────────────────────────

def filter_to_common_support(ukb_df, pmr_df, vars_, max_iter=10):
    support_levels = {}
    ukb = ukb_df.copy()
    pmr = pmr_df.copy()

    for _ in range(max_iter):
        changed = False
        for var in vars_:
            ukb_vals = sorted(pd.Series(ukb[var].dropna().astype(str)).unique().tolist())
            pmr_vals = sorted(pd.Series(pmr[var].dropna().astype(str)).unique().tolist())
            common = sorted(set(ukb_vals).intersection(pmr_vals))
            if not common:
                raise RuntimeError(f"No common support for variable: {var}")
            ukb_keep = ukb[var].notna() & ukb[var].astype(str).isin(common)
            pmr_keep = pmr[var].notna() & pmr[var].astype(str).isin(common)
            if (~ukb_keep).any() or (~pmr_keep).any():
                changed = True
                ukb = ukb.loc[ukb_keep].copy()
                pmr = pmr.loc[pmr_keep].copy()
            support_levels[var] = common
        if not changed:
            return ukb, pmr, support_levels

    raise RuntimeError("Common-support filtering did not converge.")


def build_combined_levels(ukb_df, pmr_df):
    levels = {}
    for col in RHS_CATEGORICAL:
        present = set(pd.concat([ukb_df[col], pmr_df[col]], axis=0).dropna().astype(str).unique().tolist())
        ordered = [x for x in CATEGORY_LEVEL_ORDERS[col] if x in present]
        leftovers = sorted(present.difference(ordered))
        levels[col] = ordered + leftovers
    return levels


# ────────────────────────────────────────────────────────────────────────────
# Preprocessing (one-hot + standardisation)
# ────────────────────────────────────────────────────────────────────────────

@dataclass
class Preprocessor:
    numeric_features: list
    categorical_features: list
    categorical_levels: dict
    numeric_means: dict
    numeric_sds: dict
    feature_names: list

    def transform(self, df):
        parts = []
        feature_names = []
        for col in self.numeric_features:
            x = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
            mean = self.numeric_means[col]
            sd = self.numeric_sds[col]
            if not np.isfinite(sd) or sd <= 0:
                sd = 1.0
            z = (x - mean) / sd
            parts.append(z.reshape(-1, 1))
            feature_names.append(col)
        for col in self.categorical_features:
            vals = df[col].astype("string").to_numpy()
            for lev in self.categorical_levels[col]:
                parts.append((vals == lev).astype(float).reshape(-1, 1))
                feature_names.append(f"{col}__{lev}")
        X = np.hstack(parts).astype(np.float32) if parts else np.empty((len(df), 0), dtype=np.float32)
        return X, feature_names


def serialise_preprocessor(p: Preprocessor) -> dict:
    return {
        "numeric_features": list(p.numeric_features),
        "categorical_features": list(p.categorical_features),
        "categorical_levels": {k: list(v) for k, v in p.categorical_levels.items()},
        "numeric_means": {k: float(v) for k, v in p.numeric_means.items()},
        "numeric_sds": {k: float(v) for k, v in p.numeric_sds.items()},
        "feature_names": list(p.feature_names),
    }


def fit_preprocessor(train_df, combined_levels) -> Preprocessor:
    numeric_means, numeric_sds = {}, {}
    for col in RHS_NUMERIC:
        x = pd.to_numeric(train_df[col], errors="coerce").to_numpy(dtype=float)
        mean = float(np.nanmean(x))
        sd = float(np.nanstd(x))
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        numeric_means[col] = mean
        numeric_sds[col] = sd

    feature_names = list(RHS_NUMERIC)
    for col in RHS_CATEGORICAL:
        feature_names.extend([f"{col}__{lev}" for lev in combined_levels[col]])

    return Preprocessor(
        numeric_features=list(RHS_NUMERIC),
        categorical_features=list(RHS_CATEGORICAL),
        categorical_levels=combined_levels,
        numeric_means=numeric_means,
        numeric_sds=numeric_sds,
        feature_names=feature_names,
    )


# ────────────────────────────────────────────────────────────────────────────
# Network architecture (feed-forward; hidden_dims read from the locked spec)
# ────────────────────────────────────────────────────────────────────────────

class DeepSurvNet(nn.Module):
    """Multi-hidden-layer feed-forward network with optional input dropout,
    per-layer hidden dropout, and a scalar (linear-predictor) output.

    `hidden_dims` is a list (e.g., [256, 128, 64]) read from the
    cause-specific locked JSON spec.
    """

    def __init__(self, input_dim, hidden_dims, dropout=0.0, input_dropout=0.0):
        super().__init__()
        hidden_dims = [int(x) for x in hidden_dims if int(x) > 0]
        layers = []
        if input_dropout > 0:
            layers.append(nn.Dropout(float(input_dropout)))
        prev_dim = int(input_dim)
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(prev_dim, int(hidden_dim)))
            layers.append(nn.ReLU())
            if float(dropout) > 0:
                layers.append(nn.Dropout(float(dropout)))
            prev_dim = int(hidden_dim)
        self.backbone = nn.Sequential(*layers) if layers else nn.Identity()
        self.output = nn.Linear(prev_dim, 1)

    def forward(self, X):
        X = self.backbone(X)
        return self.output(X).squeeze(1)


def make_torch_model(input_dim, spec=None, seed=None):
    """Construct a DeepSurvNet from a (locked) spec dict.

    Required spec keys: hidden_dims, dropout. Optional: input_dropout, seed.
    """
    spec = spec or {}
    if seed is None:
        seed = int(spec.get("seed", 11))
    torch.manual_seed(int(seed))
    np.random.seed(int(seed))
    hidden_dims = spec.get("hidden_dims")
    if hidden_dims is None:
        raise RuntimeError("DeepSurv spec is missing 'hidden_dims'")
    return DeepSurvNet(
        input_dim=int(input_dim),
        hidden_dims=hidden_dims,
        dropout=float(spec.get("dropout", 0.0)),
        input_dropout=float(spec.get("input_dropout", 0.0)),
    )


def tensor_from_numpy(x, dtype=torch.float32):
    return torch.from_numpy(np.asarray(x)).to(dtype=dtype)


def make_optimizer(model, lr, weight_decay=0.0, name="adam"):
    name = str(name).lower()
    if name == "adamw":
        return torch.optim.AdamW(model.parameters(), lr=float(lr), weight_decay=float(weight_decay))
    if name == "adam":
        return torch.optim.Adam(model.parameters(), lr=float(lr), weight_decay=float(weight_decay))
    raise RuntimeError(f"Unsupported optimizer: {name}")


# ────────────────────────────────────────────────────────────────────────────
# Cox loss + design tensors + Breslow baseline
# ────────────────────────────────────────────────────────────────────────────

@dataclass
class CoxDesign:
    X_ord: torch.Tensor
    weight_ord: torch.Tensor
    event_weight_ord: torch.Tensor
    group_index: torch.Tensor
    group_end_idx: torch.Tensor
    group_event_weight: torch.Tensor
    event_group_mask: torch.Tensor
    total_event_weight: float
    event_times: np.ndarray


def build_cox_design(X, time_days, status, sample_weight) -> CoxDesign:
    time_days = np.asarray(time_days, dtype=float)
    status = np.asarray(status, dtype=int)
    sample_weight = np.asarray(sample_weight, dtype=np.float32)

    order = np.argsort(-time_days, kind="mergesort")
    time_ord = time_days[order]
    status_ord = status[order]
    weight_ord = sample_weight[order]
    X_ord = np.asarray(X, dtype=np.float32)[order]

    group_start = np.r_[True, time_ord[1:] != time_ord[:-1]]
    group_index = np.cumsum(group_start).astype(np.int64) - 1
    group_end_idx = np.where(np.r_[time_ord[1:] != time_ord[:-1], True])[0].astype(np.int64)
    event_weight_ord = (weight_ord * (status_ord == 1)).astype(np.float32)
    group_event_weight = np.bincount(group_index, weights=event_weight_ord, minlength=len(group_end_idx)).astype(np.float32)
    event_group_mask = group_event_weight > 0

    return CoxDesign(
        X_ord=tensor_from_numpy(X_ord, dtype=torch.float32),
        weight_ord=tensor_from_numpy(weight_ord, dtype=torch.float32),
        event_weight_ord=tensor_from_numpy(event_weight_ord, dtype=torch.float32),
        group_index=torch.from_numpy(group_index.astype(np.int64)),
        group_end_idx=torch.from_numpy(group_end_idx.astype(np.int64)),
        group_event_weight=tensor_from_numpy(group_event_weight, dtype=torch.float32),
        event_group_mask=torch.from_numpy(event_group_mask.astype(bool)),
        total_event_weight=float(np.sum(event_weight_ord)),
        event_times=time_ord[group_end_idx].astype(float),
    )


def weighted_cox_loss_torch_efron(model: nn.Module, design: CoxDesign) -> torch.Tensor:
    """Weighted Cox partial likelihood with **Efron** ties handling.

    For each event-time group g with d_g tied events:

      contribution_g =  Σ_{i in events of g} w_i * η_i
                      − (Σ_{i in events of g} w_i) / d_g
                          · Σ_{l=0..d_g-1} log( S_g − (l/d_g) · T_g )

    where
        S_g = Σ_{k in risk_set(g)} w_k * exp(η_k)
        T_g = Σ_{i in events(g)} w_i * exp(η_i)
        d_g = number of events in g (unweighted count)
        w_i = per-row sample weight (UKB unit, wUKB IPW, or PMR design 1/20)
    """
    lp = model(design.X_ord)
    weights = design.weight_ord

    # Stabilise exp(lp) by subtracting per-row log_weight before exponentiating.
    weighted_exp = weights * torch.exp(lp)

    # Cumulative sum in time-descending order yields the weighted risk-set
    # sum at each subject's time. Risk-set sum for each group = the cumulative
    # sum at group_end_idx (the last subject in that group).
    cum_weighted_exp = torch.cumsum(weighted_exp, dim=0)
    S = cum_weighted_exp[design.group_end_idx]  # (n_groups,)

    # Per-group weighted exp sum over tied events: T_g.
    is_event = design.event_weight_ord > 0
    event_weighted_exp = torch.where(is_event, weighted_exp, torch.zeros_like(weighted_exp))
    T = torch.zeros_like(S)
    T.index_add_(0, design.group_index, event_weighted_exp)

    # Per-group weighted sum of (event) lp's: numerator term.
    weighted_event_lp = design.event_weight_ord * lp
    group_event_eta = torch.zeros_like(S)
    group_event_eta.index_add_(0, design.group_index, weighted_event_lp)

    # Number of events per group (unweighted count).
    is_event_long = is_event.to(torch.long)
    d_g = torch.zeros(S.shape[0], dtype=torch.long, device=lp.device)
    d_g.scatter_add_(0, design.group_index, is_event_long)
    d_g_f = d_g.to(lp.dtype)

    sum_w_event = design.group_event_weight  # already the weighted sum per group

    mask = design.event_group_mask
    S_ev = S[mask]
    T_ev = T[mask]
    d_ev = d_g_f[mask]
    sum_w_ev = sum_w_event[mask]
    eta_w_ev = group_event_eta[mask]

    if d_ev.numel() == 0:
        return torch.zeros((), dtype=lp.dtype, device=lp.device)

    d_max = int(d_g_f.max().item())
    if d_max <= 0:
        return torch.zeros((), dtype=lp.dtype, device=lp.device)

    # Efron's adjustment: for the l-th tied event (l = 0..d_g-1),
    # subtract (l/d_g) * T_g from S_g before logging.
    l_idx = torch.arange(d_max, dtype=lp.dtype, device=lp.device)              # (d_max,)
    l_over_d = l_idx.unsqueeze(0) / d_ev.unsqueeze(1)                          # (n_g, d_max)
    denom = S_ev.unsqueeze(1) - l_over_d * T_ev.unsqueeze(1)                   # (n_g, d_max)
    valid = l_idx.unsqueeze(0) < d_ev.unsqueeze(1)                             # bool mask
    log_denom = torch.where(
        valid,
        torch.log(torch.clamp(denom, min=1e-12)),
        torch.zeros_like(denom),
    )
    sum_log_denom = log_denom.sum(dim=1)                                       # (n_g,)

    contrib = eta_w_ev - (sum_w_ev / d_ev) * sum_log_denom
    total = contrib.sum()
    denom_normaliser = max(design.total_event_weight, 1e-8)
    return -total / denom_normaliser


def weighted_cox_loss_torch(model: nn.Module, design: CoxDesign,
                            ties: str = "efron") -> torch.Tensor:
    """Weighted Cox partial likelihood with Efron ties handling."""
    ties = (ties or "efron").lower()
    if ties != "efron":
        raise RuntimeError(f"Unsupported ties handling: {ties!r}")
    return weighted_cox_loss_torch_efron(model, design)


def build_weighted_breslow_baseline(model, design: CoxDesign):
    model.eval()
    with torch.no_grad():
        lp = model(design.X_ord)
        log_weight = torch.log(torch.clamp(design.weight_ord, min=1e-8))
        log_risk_term = lp + log_weight
        log_cum_risk = torch.logcumsumexp(log_risk_term, dim=0)
        group_log_risk = log_cum_risk[design.group_end_idx]
        mask = design.event_group_mask
        event_times_desc = design.event_times[design.event_group_mask.cpu().numpy()].astype(float)
        delta_hazard_desc = design.group_event_weight[mask] / torch.exp(group_log_risk[mask])
        event_times = np.asarray(event_times_desc[::-1], dtype=float)
        delta_hazard = torch.flip(delta_hazard_desc, dims=[0])
        cum_hazard = torch.cumsum(delta_hazard, dim=0)
    return {"event_times": event_times,
            "cum_baseline_hazard": cum_hazard.cpu().numpy().astype(float)}


def predict_lp_torch(model, X, batch_size=32768):
    model.eval()
    X_t = tensor_from_numpy(X, dtype=torch.float32)
    out = []
    with torch.no_grad():
        for s in range(0, X_t.shape[0], batch_size):
            e = min(s + batch_size, X_t.shape[0])
            out.append(model(X_t[s:e]).cpu().numpy())
    return np.concatenate(out, axis=0) if out else np.empty((0,), dtype=np.float32)


def predict_risk_from_baseline(lp, baseline, horizons_days):
    lp = np.asarray(lp, dtype=float)
    event_times = np.asarray(baseline["event_times"], dtype=float)
    cum_hazard = np.asarray(baseline["cum_baseline_hazard"], dtype=float)
    out = {}
    for horizon in horizons_days:
        horizon = float(horizon)
        if len(event_times) == 0 or horizon <= 0:
            base_h = 0.0
        else:
            ix = np.searchsorted(event_times, horizon, side="right") - 1
            base_h = float(cum_hazard[ix]) if ix >= 0 else 0.0
        out[int(horizon)] = 1.0 - np.exp(-base_h * np.exp(lp))
    return out


# ────────────────────────────────────────────────────────────────────────────
# Stratified k-fold split + nested 80/10/10 within fold
# ────────────────────────────────────────────────────────────────────────────

def stratified_kfold_indices(status, k=5, seed=11):
    """Stratified k-fold splitter on a binary event indicator. Returns a list
    of length k, each element a numpy array of test-fold row indices."""
    rng = np.random.default_rng(int(seed))
    status = np.asarray(status, dtype=int)
    event_idx = np.where(status == 1)[0]
    nonevent_idx = np.where(status == 0)[0]
    rng.shuffle(event_idx)
    rng.shuffle(nonevent_idx)
    folds_event = np.array_split(event_idx, k)
    folds_nonevent = np.array_split(nonevent_idx, k)
    folds = [np.sort(np.concatenate([fe, fn])) for fe, fn in zip(folds_event, folds_nonevent)]
    return folds


def stratified_train_val_test_split_within_fold(status, fold_test_idx, seed=11):
    """In a 5-fold CV scheme, the per-fold "out-fold" partition is 20%
    of the cohort. This function splits that 20% into a 10% validation
    partition (early stopping / LR drops) and a 10% held-out test
    partition, each stratified on event status. Train is the remaining
    80% (the union of the other four folds).

    Returns (train_idx, val_idx, test_idx) in original-row coordinates.
    """
    rng = np.random.default_rng(int(seed))
    status = np.asarray(status, dtype=int)
    n = len(status)
    fold_idx = np.asarray(fold_test_idx, dtype=int)

    fold_status = status[fold_idx]
    event_pos = np.where(fold_status == 1)[0]
    nonevent_pos = np.where(fold_status == 0)[0]
    rng.shuffle(event_pos)
    rng.shuffle(nonevent_pos)

    # Half of each (event / non-event) goes to val, the other half to test.
    n_event_val = len(event_pos) // 2
    n_nonevent_val = len(nonevent_pos) // 2
    val_pos = np.concatenate([event_pos[:n_event_val], nonevent_pos[:n_nonevent_val]])
    test_pos = np.concatenate([event_pos[n_event_val:], nonevent_pos[n_nonevent_val:]])

    val_idx = np.sort(fold_idx[val_pos])
    test_idx = np.sort(fold_idx[test_pos])

    train_mask = np.ones(n, dtype=bool)
    train_mask[fold_idx] = False
    train_idx = np.where(train_mask)[0]

    return train_idx, val_idx, test_idx


# ────────────────────────────────────────────────────────────────────────────
# Wagner-style training with LR-decay schedule
# ────────────────────────────────────────────────────────────────────────────

def _resolve_training_schedule(spec):
    """Build the runtime training schedule: cause-specific locked-spec
    fields (`hidden_dims`, `dropout`, `optimizer`, `weight_decay`, `seed`)
    are read directly; LR schedule and per-stage epoch/patience values
    fall back to TRAINING_SCHEDULE_DEFAULTS."""
    out = dict(TRAINING_SCHEDULE_DEFAULTS)
    if spec:
        for k in ("optimizer", "weight_decay", "ties", "input_dropout",
                  "base_learning_rate", "lr_floor", "lr_decay_factor",
                  "epochs_per_stage", "early_stop_patience",
                  "grad_clip_max_norm", "seed"):
            if k in spec:
                out[k] = spec[k]
    return out


def train_with_lr_schedule(model, train_design, val_design, spec, seed=None):
    """Train `model` with the Wagner-style LR-decay schedule and Efron loss.

    Hyperparameters from `spec` (or TRAINING_SCHEDULE_DEFAULTS):
      - base_learning_rate, lr_floor, lr_decay_factor : LR schedule
      - epochs_per_stage, early_stop_patience          : per-stage stop rule
      - optimizer, weight_decay                        : optimiser
      - ties                                           : 'efron'
      - grad_clip_max_norm                             : per-step grad clip

    Tracks the lowest validation loss across **all stages**, restoring the
    corresponding model weights at the end.

    Returns (model, history dict, best_val_loss). The best-validation
    checkpoint weights are restored and deployed directly; the model is not
    refit on the full data.
    """
    sched = _resolve_training_schedule(spec)
    if seed is None:
        seed = int(sched.get("seed", 11))
    torch.manual_seed(int(seed))
    np.random.seed(int(seed))
    torch.set_num_threads(1)  # single-threaded CPU for reproducible training

    base_lr = float(sched["base_learning_rate"])
    lr_floor = float(sched["lr_floor"])
    lr_decay = float(sched["lr_decay_factor"])
    epochs_per_stage = int(sched["epochs_per_stage"])
    patience = int(sched["early_stop_patience"])
    grad_clip = float(sched.get("grad_clip_max_norm", 5.0))
    weight_decay = float(sched.get("weight_decay", 0.0))
    optimizer_name = str(sched.get("optimizer", "adamw"))
    ties = str(sched.get("ties", "efron"))

    history = {"stage": [], "lr": [], "epoch": [], "train_loss": [], "val_loss": []}
    best_val = float("inf")
    best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}


    lr = base_lr
    stage = 0
    while True:
        optimizer = make_optimizer(model, lr=lr, weight_decay=weight_decay, name=optimizer_name)
        no_improve = 0
        for epoch in range(1, epochs_per_stage + 1):
            model.train()
            optimizer.zero_grad(set_to_none=True)
            loss = weighted_cox_loss_torch(model, train_design, ties=ties)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
            optimizer.step()

            train_loss = float(loss.detach().cpu().item())
            if val_design is not None:
                model.eval()
                with torch.no_grad():
                    vloss = float(weighted_cox_loss_torch(model, val_design, ties=ties).detach().cpu().item())
            else:
                vloss = train_loss

            history["stage"].append(int(stage))
            history["lr"].append(float(lr))
            history["epoch"].append(int(epoch))
            history["train_loss"].append(train_loss)
            history["val_loss"].append(vloss)

            if vloss + 1e-6 < best_val:
                best_val = vloss
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                no_improve = 0
            else:
                no_improve += 1
                if no_improve >= patience:
                    break

        stage += 1
        if lr <= lr_floor + 1e-12:
            break
        lr = max(lr / lr_decay, lr_floor)

    model.load_state_dict(best_state)
    return model, history, best_val


# ────────────────────────────────────────────────────────────────────────────
# Path / config helpers
# ────────────────────────────────────────────────────────────────────────────

def find_framework_root(start_path: Path) -> Path:
    p = Path(start_path).resolve()
    for parent in [p, *p.parents]:
        if (parent / "framework_config.json").is_file():
            return parent
    raise RuntimeError(f"framework_config.json not found upward from {start_path}")


def read_framework_config(start_path: Path) -> dict:
    root = find_framework_root(start_path)
    with open(root / "framework_config.json") as f:
        cfg = json.load(f)
    cfg["framework_root"] = str(root)
    return cfg


def get_outcome_cfg_from_framework(cfg: dict, key: str) -> dict:
    outcomes = cfg["phase2"]["outcomes"]
    if key not in outcomes:
        raise RuntimeError(f"Outcome {key!r} not in framework_config.json")
    raw = outcomes[key]
    return {
        "key": key,
        "label": raw["label"],
        "type": raw["type"],
        "icd10": list(raw.get("icd10", []) or []),
    }


def get_weight_source_from_framework(cfg: dict, key: str) -> dict:
    sources = cfg["phase2"]["weight_sources"]
    if key not in sources:
        raise RuntimeError(f"Weight source {key!r} not in framework_config.json")
    return {"key": key, **sources[key]}


def get_evaluation_dir(framework_root: Path) -> Path:
    out = Path(framework_root) / "mortality_risk_prediction" / "results" / "evaluation" / "deep_surv"
    out.mkdir(parents=True, exist_ok=True)
    return out


def get_cause_dir(framework_root: Path, cause: str) -> Path:
    out = get_evaluation_dir(framework_root) / cause
    out.mkdir(parents=True, exist_ok=True)
    return out


def get_fold_dir(framework_root: Path, cause: str, source: str, fold: int) -> Path:
    out = get_cause_dir(framework_root, cause) / "folds" / source / f"fold_{fold}"
    out.mkdir(parents=True, exist_ok=True)
    return out


# ────────────────────────────────────────────────────────────────────────────
# Locked-spec loader
# ────────────────────────────────────────────────────────────────────────────
#
# A single architecture (selected once on all-cause mortality) is used for
# every outcome, read from deep_surv_specs/final_locked_v1.json.

LOCKED_SPEC_FILENAME = "final_locked_v1.json"


def load_locked_spec(framework_root: Path) -> dict:
    """Read the single locked spec (selected once on all-cause mortality)
    that is applied uniformly across all outcomes. Only `hidden_dims`,
    `dropout`, and (optionally) `seed` are used downstream; any other
    fields the spec records are superseded by TRAINING_SCHEDULE_DEFAULTS,
    which defines the LR schedule, optimizer, and ties handling."""
    path = Path(framework_root) / "mortality_risk_prediction" / "code" / "evaluate" / "deep_surv_specs" / LOCKED_SPEC_FILENAME
    if not path.exists():
        raise RuntimeError(f"Missing locked spec file: {path}")
    with open(path) as f:
        spec = json.load(f)
    if "hidden_dims" not in spec:
        raise RuntimeError(f"Locked spec at {path} is missing 'hidden_dims'")
    return spec
