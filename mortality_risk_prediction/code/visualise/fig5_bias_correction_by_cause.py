#!/usr/bin/env python3
"""Figure 5 - bias correction by cause (UKB, PMR, weighted UKB HSE Super Learner)."""

import gc
import matplotlib.ticker as mticker
from helpers import *

CAUSE_ORDER = [
    "all_cause_mortality",
    "cancer_mortality",
    "cardiovascular_mortality",
    "respiratory_mortality",
    "digestive_mortality",
]

DARK_BLUE = COLOURS["ukbw_hse_superlearner"]


def _apply_box_spine(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)
        spine.set_color("black")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=11)


def _load_bootstrap_overall(results_root, cause, hz):
    """Return (wukb_ci_lo, wukb_ci_hi) per-million from bootstrap overall row."""
    ci_path = results_root / "visualise" / "bootstrap_ci" / \
              f"bootstrap_ci__{cause}__{hz}y.csv"
    if not ci_path.exists():
        return np.nan, np.nan
    df = pd.read_csv(ci_path)
    row = df[(df["var"] == "overall") & (df["level"] == "all")]
    if row.empty:
        return np.nan, np.nan
    return (float(row["wukb_bias_ci_lo_per_mil"].iloc[0]),
            float(row["wukb_bias_ci_hi_per_mil"].iloc[0]))


def _correction_ci(ukb_bias, wukb_ci_lo, wukb_ci_hi):
    """Transform wUKB bias CI to bias correction % CI (UKB treated as fixed)."""
    abs_ukb = abs(ukb_bias)
    if abs_ukb == 0 or not np.isfinite(abs_ukb):
        return np.nan, np.nan
    # wukb values are negative; ci_lo is more negative (larger |wukb|) → less correction
    # wukb ci_hi is less negative (smaller |wukb|) → more correction
    corr_lo = (abs_ukb - abs(wukb_ci_hi)) / abs_ukb * 100
    corr_hi = (abs_ukb - abs(wukb_ci_lo)) / abs_ukb * 100
    return float(corr_lo), float(corr_hi)


def plot_one_horizon(records_df, hz, output_dir):
    sub = records_df[records_df["horizon_years"] == hz].copy()
    cause_rank = {c: i for i, c in enumerate(CAUSE_ORDER)}
    sub["_rank"] = sub["cause"].map(cause_rank)
    sub = sub.sort_values("_rank").reset_index(drop=True)
    y = np.arange(len(sub))

    sub["bias_correction_pct"] = np.where(
        np.abs(sub["ukb_bias_per_mil"]) > 0,
        (np.abs(sub["ukb_bias_per_mil"]) - np.abs(sub["wukb_bias_per_mil"])) /
        np.abs(sub["ukb_bias_per_mil"]) * 100,
        np.nan,
    )

    fig, (ax, ax_bc) = plt.subplots(1, 2, figsize=(8, 3),
                                     gridspec_kw={"width_ratios": [3, 1.2]}, sharey=True)

    for i, row in sub.iterrows():
        ax.scatter(row["pmr_bias_per_mil"], i, c=C_PMR, marker=MARKERS["pmr"],
                   s=50, zorder=5, edgecolors="none")
        ax.scatter(row["ukb_bias_per_mil"], i, c=C_UKB, marker=MARKERS["ukb"],
                   s=50, zorder=5, edgecolors="none")
        ukb_val  = row["ukb_bias_per_mil"]
        wukb_val = row["wukb_bias_per_mil"]
        if np.isfinite(ukb_val) and np.isfinite(wukb_val):
            ax.annotate("", xy=(wukb_val, i), xytext=(ukb_val, i),
                        arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8))
        ax.scatter(row["wukb_bias_per_mil"], i, c=DARK_BLUE,
                   marker=MARKERS["ukbw_hse_superlearner"], s=50, zorder=7, edgecolors="none")

    for i, row in sub.iterrows():
        bc = row["bias_correction_pct"]
        if not np.isfinite(bc):
            continue
        ax_bc.barh(i, bc, height=0.28, color=DARK_BLUE, zorder=4, left=0)


    _y_bot = len(sub) - 0.5 + 0.4
    _y_top = -0.5 - 0.4
    ax.plot([0, 0], [_y_bot, _y_top], color="grey", linewidth=1.6, alpha=0.7, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(sub["cause_label"].tolist(), fontsize=11)
    ax.invert_yaxis()
    ax.set_ylim(_y_bot, _y_top)

    _apply_box_spine(ax)
    _apply_box_spine(ax_bc)
    ax.tick_params(axis="y", length=4, width=1.2)
    ax_bc.tick_params(axis="y", left=False, labelleft=False)

    ax.set_xlabel("Prediction bias (deaths per 1,000,000)", fontsize=12)
    ax_bc.set_xlabel("Bias correction (%)", fontsize=12)
    ax_bc.set_xlim(0, 100)
    ax_bc.set_xticks([0, 25, 50, 75, 100])

    if hz == 10:
        ax.set_xticks([-20000, -15000, -10000, -5000, 0])
        ax.set_xlim(-23000, 3000)
    else:
        ax.set_xticks([-10000, -7500, -5000, -2500, 0])
        ax.set_xlim(-10500, 1300)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    fig.tight_layout()

    fname = f"fig5_bias_correction_by_cause__{hz}y"
    save_fig(fig, output_dir, fname)
    print(f"  Saved {fname}")


def main():
    results_root = get_results_root()
    output_dir   = get_output_dir("fig5_bias_correction_by_cause")

    pred_cols = {
        "pmr":  "pmr_pred_{hz}y",
        "ukb":  "ukb_pred_{hz}y",
        "wukb": "ukbw_hse_superlearner_pred_{hz}y",
    }
    all_records = []

    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            continue
        for hz in [5, 10]:
            event_col = f"event_{hz}y"
            pc = {k: v.replace("{hz}", str(hz)) for k, v in pred_cols.items()}
            rates = compute_rates(df, event_col, pc)
            rates["cause"]         = cause
            rates["cause_label"]   = CAUSE_LABELS[cause]
            rates["horizon_years"] = hz
            wukb_ci_lo, wukb_ci_hi = _load_bootstrap_overall(results_root, cause, hz)
            rates["wukb_ci_lo"] = wukb_ci_lo
            rates["wukb_ci_hi"] = wukb_ci_hi
            all_records.append(rates)
        del df; gc.collect()

    records_df = pd.DataFrame(all_records)
    records_df["pmr_bias_per_mil"]  = (records_df["pmr_pred_rate"]  - records_df["observed_rate"]) * 1_000_000
    records_df["ukb_bias_per_mil"]  = (records_df["ukb_pred_rate"]  - records_df["observed_rate"]) * 1_000_000
    records_df["wukb_bias_per_mil"] = (records_df["wukb_pred_rate"] - records_df["observed_rate"]) * 1_000_000

    records_df[["correction_ci_lo", "correction_ci_hi"]] = records_df.apply(
        lambda r: pd.Series(_correction_ci(r["ukb_bias_per_mil"],
                                           r["wukb_ci_lo"], r["wukb_ci_hi"])),
        axis=1
    )

    cause_rank = {c: i for i, c in enumerate(CAUSE_ORDER)}
    for hz in [5, 10]:
        sub = records_df[records_df["horizon_years"] == hz].copy()
        sub["_rank"] = sub["cause"].map(cause_rank)
        sub = sub.sort_values("_rank").reset_index(drop=True)
        ukb_abs  = np.abs(sub["ukb_bias_per_mil"])
        wukb_abs = np.abs(sub["wukb_bias_per_mil"])
        bc_pct   = np.where(ukb_abs > 0, (ukb_abs - wukb_abs) / ukb_abs * 100, np.nan)
        summary  = pd.DataFrame({
            "cause":               "    " + sub["cause_label"],
            "UKB predicted":       (sub["ukb_pred_rate"]  * 1_000_000).round(0),
            "wUKB predicted":      (sub["wukb_pred_rate"] * 1_000_000).round(0),
            "Observed":            (sub["observed_rate"]  * 1_000_000).round(0),
            "UKB relative bias":   (sub["ukb_relative_bias"]  * 100).round(0).astype(str) + "%",
            "wUKB relative bias":  (sub["wukb_relative_bias"] * 100).round(0).astype(str) + "%",
            "Bias correction":     pd.Series(bc_pct).apply(
                lambda x: f"{int(round(x))}%" if np.isfinite(x) else ""),
            "Correction CI":       sub.apply(
                lambda r: f"({int(round(r.correction_ci_lo))}%, {int(round(r.correction_ci_hi))}%)"
                          if np.isfinite(r.correction_ci_lo) else "", axis=1),
        })
        save_csv(summary, output_dir, f"fig5_bias_correction_by_cause__summary__{hz}y")

    plot_one_horizon(records_df, 5, output_dir)
    plot_one_horizon(records_df, 10, output_dir)

    fig_leg, ax_leg = plt.subplots(figsize=(2.2, 1.4))
    ax_leg.scatter([], [], c="#888888", marker=MARKERS["pmr"], s=30, label="PMR")
    ax_leg.scatter([], [], c="#8b0000", marker=MARKERS["ukb"], s=30, label="UKB")
    ax_leg.scatter([], [], c=DARK_BLUE, marker=MARKERS["ukbw_hse_superlearner"], s=30, label="wUKB")
    ax_leg.legend(fontsize=10, ncol=1, loc="center", frameon=False)
    ax_leg.axis("off")
    fig_leg.tight_layout()
    save_fig(fig_leg, output_dir, "fig5_bias_correction_by_cause__legend")

    print("Done: fig5_bias_correction_by_cause")


if __name__ == "__main__":
    main()
