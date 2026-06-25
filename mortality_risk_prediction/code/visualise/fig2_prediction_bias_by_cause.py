#!/usr/bin/env python3
"""Figure 2 - prediction bias by cause (UKB vs PMR, no bias correction)."""

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


def _fmt_ci(center, lo, hi):
    """Format as 'xxx (xxx, xxx)' with comma-separated integers."""
    def fi(x):
        return f"{int(round(x)):,}" if np.isfinite(x) else ""
    return f"{fi(center)} ({fi(lo)}, {fi(hi)})"


def _save_csv(df, output_dir, fname_base):
    df_r = df.copy()
    for col in df_r.columns:
        if pd.api.types.is_float_dtype(df_r[col]):
            df_r[col] = df_r[col].apply(
                lambda x: "" if pd.isna(x) else f"{int(round(x)):,}"
            )
    df_r.to_excel(output_dir / f"{fname_base}.xlsx", index=False)


def plot_one_horizon(records_df, hz, output_dir):
    sub = records_df[records_df["horizon_years"] == hz].copy()
    cause_rank = {c: i for i, c in enumerate(CAUSE_ORDER)}
    sub["_rank"] = sub["cause"].map(cause_rank)
    sub = sub.sort_values("_rank").reset_index(drop=True)
    y = np.arange(len(sub))

    fig, ax = plt.subplots(figsize=(4.5, 3))

    for i, row in sub.iterrows():
        ax.plot([0, row["ukb_bias_per_mil"]], [i, i], color=C_UKB, linewidth=1.6, zorder=2)
        ax.scatter(row["ukb_bias_per_mil"], i, c=C_UKB, marker=MARKERS["ukb"],
                   s=50, zorder=4, edgecolors="none")
        ax.plot([0, row["pmr_bias_per_mil"]], [i, i], color=C_PMR, linewidth=1.6, zorder=5)
        ax.scatter(row["pmr_bias_per_mil"], i, c=C_PMR, marker=MARKERS["pmr"],
                   s=50, zorder=6, edgecolors="none")

    ax.plot([0, 0], [len(sub), -1.7], color="black", linewidth=2.0, zorder=3, clip_on=False)

    ax.set_yticks(y)
    ax.set_yticklabels([""] * len(y))
    ax.invert_yaxis()
    ax.set_ylim(len(sub) - 0.5 + 0.5, -0.5 - 0.7)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.spines["bottom"].set_color("black")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", width=2.0, labelsize=11)

    if hz == 10:
        ax.set_xticks([-20000, -15000, -10000, -5000, 0])
        ax.set_xlim(-23000, 3000)
    else:
        ax.set_xticks([-10000, -7500, -5000, -2500, 0])
        ax.set_xlim(-10500, 1300)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_xlabel("Prediction bias (deaths per 1,000,000)", fontsize=12)

    fig.tight_layout()

    fname = f"fig2_prediction_bias_by_cause__{hz}y"
    save_fig(fig, output_dir, fname)
    print(f"  Saved {fname}")


def _load_bootstrap_ci(results_root, cause, hz):
    ci_path = results_root / "visualise" / "bootstrap_ci" / \
              f"bootstrap_ci__{cause}__{hz}y.csv"
    if not ci_path.exists():
        return np.nan, np.nan
    df = pd.read_csv(ci_path)
    row = df[(df["var"] == "overall") & (df["level"] == "all")]
    if row.empty:
        return np.nan, np.nan
    return float(row["ukb_bias_ci_lo_per_mil"].iloc[0]), \
           float(row["ukb_bias_ci_hi_per_mil"].iloc[0])


def main():
    results_root = get_results_root()
    output_dir   = get_output_dir("fig2_prediction_bias_by_cause")

    pred_cols = {"pmr": "pmr_pred_{hz}y", "ukb": "ukb_pred_{hz}y"}
    all_records = []

    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            print(f"[SKIP] {cause}")
            continue
        for hz in [5, 10]:
            event_col = f"event_{hz}y"
            pc = {k: v.replace("{hz}", str(hz)) for k, v in pred_cols.items()}
            rates = compute_rates(df, event_col, pc)
            rates["cause"]         = cause
            rates["cause_label"]   = CAUSE_LABELS[cause]
            rates["horizon_years"] = hz
            ci_lo, ci_hi = _load_bootstrap_ci(results_root, cause, hz)
            rates["ukb_bias_ci_lo"] = ci_lo
            rates["ukb_bias_ci_hi"] = ci_hi
            all_records.append(rates)
        del df; gc.collect()

    records_df = pd.DataFrame(all_records)
    records_df["pmr_bias_per_mil"] = (
        records_df["pmr_pred_rate"] - records_df["observed_rate"]
    ) * 1_000_000
    records_df["ukb_bias_per_mil"] = (
        records_df["ukb_pred_rate"] - records_df["observed_rate"]
    ) * 1_000_000

    cause_rank = {c: i for i, c in enumerate(CAUSE_ORDER)}
    for hz in [5, 10]:
        sub = records_df[records_df["horizon_years"] == hz].copy()
        sub["_rank"] = sub["cause"].map(cause_rank)
        sub = sub.sort_values("_rank").reset_index(drop=True)
        summary = pd.DataFrame({
            "cause": "    " + sub["cause_label"],
            "UKB bias per million (95% CI)": sub.apply(
                lambda r: _fmt_ci(r["ukb_bias_per_mil"],
                                  r["ukb_bias_ci_lo"],
                                  r["ukb_bias_ci_hi"]), axis=1
            ),
        })
        summary.to_excel(output_dir / f"fig2_prediction_bias_by_cause__summary__{hz}y.xlsx",
                         index=False)

    plot_one_horizon(records_df, 5, output_dir)
    plot_one_horizon(records_df, 10, output_dir)

    fig_leg, ax_leg = plt.subplots(figsize=(4, 0.8))
    ax_leg.scatter([], [], c=C_PMR, marker=MARKERS["pmr"], s=50, label="PMR")
    ax_leg.scatter([], [], c=C_UKB, marker=MARKERS["ukb"], s=50, label="UKB")
    ax_leg.legend(fontsize=14, ncol=2, loc="center", frameon=False)
    ax_leg.axis("off")
    fig_leg.tight_layout()
    save_fig(fig_leg, output_dir, "fig2_prediction_bias_by_cause__legend")

    print("Done: fig2_prediction_bias_by_cause")


if __name__ == "__main__":
    main()
