#!/usr/bin/env python3
"""Figure S13 - Cox vs deep survival prediction bias by cause (10y)."""

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

COL_UKB_LEFT  = C_SENS_UKB_ALT
COL_PMR_LEFT  = C_SENS_PMR_ALT
COL_UKB_RIGHT = C_UKB
COL_PMR_RIGHT = C_PMR

PAD = 0.45


def load_linearcox_records(results_root):
    """Both PMR and UKB from the linear_cox evaluation pipeline."""
    records = []
    for cause in CAUSES:
        if cause in ALL_CAUSE_KEYS:
            summary_path = (results_root / "evaluation" / "linear_cox" /
                            cause / "summaries")
        else:
            summary_path = (results_root / "evaluation" / "linear_cox" /
                            cause / "underlying" / "summaries")
        risk_files = list(summary_path.glob("pmr_reference_vs_transferred_risk__*.csv"))
        if not risk_files:
            continue
        risk_df = pd.read_csv(risk_files[0])
        row10 = risk_df[(risk_df["strata_variable"] == "Overall") & (risk_df["horizon_years"] == 10)]
        if row10.empty:
            continue
        row = row10.iloc[0]
        obs = float(row["pmr_obs"])
        records.append({
            "cause":        cause,
            "cause_label":  CAUSE_LABELS[cause],
            "pmr_bias_per_mil": (float(row["pmr_pred"]) - obs) * 1_000_000,
            "ukb_bias_per_mil": (float(row["ukb_pred"])  - obs) * 1_000_000,
        })
    return records


def load_deepsurv_records(results_root):
    pred_cols = {"pmr": "pmr_pred_10y", "ukb": "ukb_pred_10y"}
    records = []
    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            continue
        rates = compute_rates(df, "event_10y", pred_cols)
        records.append({
            "cause":        cause,
            "cause_label":  CAUSE_LABELS[cause],
            "pmr_bias_per_mil": (rates["pmr_pred_rate"] - rates["observed_rate"]) * 1_000_000,
            "ukb_bias_per_mil": (rates["ukb_pred_rate"] - rates["observed_rate"]) * 1_000_000,
        })
        del df; gc.collect()
    return records


def _draw_panel(ax, rdf, pmr_col, ukb_col, show_yticks=True):
    y = np.arange(len(rdf))
    for i, row in rdf.reset_index(drop=True).iterrows():
        pmr_val = row["pmr_bias_per_mil"]
        ukb_val = row["ukb_bias_per_mil"]
        pmr_closer = abs(pmr_val) < abs(ukb_val)
        ax.plot([0, ukb_val], [i, i], color=ukb_col, linewidth=1.6, zorder=2)
        ax.scatter(ukb_val, i, c=ukb_col, marker=MARKERS["ukb"],
                   s=50, zorder=4 if pmr_closer else 6, edgecolors="none")
        ax.plot([0, pmr_val], [i, i], color=pmr_col, linewidth=1.6, zorder=3)
        ax.scatter(pmr_val, i, c=pmr_col, marker=MARKERS["pmr"],
                   s=50, zorder=6 if pmr_closer else 4, edgecolors="none")

    _y_bot = len(rdf) - 0.5 + PAD
    _y_top = -0.5 - PAD
    ax.plot([0, 0], [_y_bot, _y_top], color="grey", linewidth=1.2, alpha=0.7)
    ax.set_yticks(y)
    if show_yticks:
        ax.set_yticklabels(rdf["cause_label"].tolist(), fontsize=11)
    else:
        ax.tick_params(axis="y", left=False, labelleft=False)
    ax.invert_yaxis()
    ax.set_ylim(_y_bot, _y_top)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_xlabel("Prediction bias (deaths per 1,000,000)", fontsize=10)


def main():
    results_root = get_results_root()
    output_dir   = get_output_dir("figS13_cox_vs_deepsurv")

    lc_records = load_linearcox_records(results_root)
    ds_records = load_deepsurv_records(results_root)

    if not lc_records:
        print("No linear_cox data found.")
        return

    cause_rank = {c: i for i, c in enumerate(CAUSE_ORDER)}
    lc_df = pd.DataFrame(lc_records)
    lc_df["_rank"] = lc_df["cause"].map(cause_rank)
    lc_df = lc_df.sort_values("_rank").reset_index(drop=True)

    ds_df = pd.DataFrame(ds_records)
    ds_df["_rank"] = ds_df["cause"].map(cause_rank)
    ds_df = ds_df.sort_values("_rank").reset_index(drop=True)

    n = len(lc_df)
    fig_height = max(3, n * 0.45 + 0.6)
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(8, fig_height),
                                             gridspec_kw={"width_ratios": [1, 1]},
                                             sharey=True)

    _draw_panel(ax_left,  lc_df, COL_PMR_LEFT,  COL_UKB_LEFT,  show_yticks=True)
    _draw_panel(ax_right, ds_df, COL_PMR_RIGHT, COL_UKB_RIGHT, show_yticks=False)

    fig.tight_layout()
    save_fig(fig, output_dir, "figS13_cox_vs_deepsurv")

    # Per-panel legends are saved separately so the two model families can be labelled independently
    fig_leg1, ax_leg1 = plt.subplots(figsize=(3.0, 1.0))
    ax_leg1.scatter([], [], c=COL_PMR_LEFT, marker=MARKERS["pmr"], s=50,
                    edgecolors="none", label="PMR (Cox PH)")
    ax_leg1.scatter([], [], c=COL_UKB_LEFT, marker=MARKERS["ukb"], s=50,
                    edgecolors="none", label="UKB (Cox PH)")
    ax_leg1.legend(fontsize=10, ncol=1, loc="center", frameon=False)
    ax_leg1.axis("off")
    fig_leg1.tight_layout()
    leg1_path = output_dir / "figS13_cox_vs_deepsurv__legend_left.svg"
    fig_leg1.savefig(leg1_path, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig_leg1)
    print(f"  Saved legend (left) to {leg1_path}")

    fig_leg2, ax_leg2 = plt.subplots(figsize=(3.0, 1.0))
    ax_leg2.scatter([], [], c=COL_PMR_RIGHT, marker=MARKERS["pmr"], s=50,
                    edgecolors="none", label="PMR (DeepSurv)")
    ax_leg2.scatter([], [], c=COL_UKB_RIGHT, marker=MARKERS["ukb"], s=50,
                    edgecolors="none", label="UKB (DeepSurv)")
    ax_leg2.legend(fontsize=10, ncol=1, loc="center", frameon=False)
    ax_leg2.axis("off")
    fig_leg2.tight_layout()
    leg2_path = output_dir / "figS13_cox_vs_deepsurv__legend_right.svg"
    fig_leg2.savefig(leg2_path, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig_leg2)
    print(f"  Saved legend (right) to {leg2_path}")

    print("Done: figS13_cox_vs_deepsurv")


if __name__ == "__main__":
    main()
