#!/usr/bin/env python3
"""Figure 7 sensitivity - alternative weighting model: HSE lasso-logit vs HSE Super Learner (10y)."""

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

PAD = 0.45


def load_hse_ll_records(results_root):
    pred_cols = {
        "pmr":  "pmr_pred_10y",
        "ukb":  "ukb_pred_10y",
        "wukb": "ukbw_hse_lassologit_pred_10y",
    }
    records = []
    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            continue
        rates = compute_rates(df, "event_10y", pred_cols)
        records.append({
            "cause":        cause,
            "cause_label":  CAUSE_LABELS[cause],
            "pmr_bias_per_mil":  (rates["pmr_pred_rate"]  - rates["observed_rate"]) * 1_000_000,
            "ukb_bias_per_mil":  (rates["ukb_pred_rate"]  - rates["observed_rate"]) * 1_000_000,
            "wukb_bias_per_mil": (rates["wukb_pred_rate"] - rates["observed_rate"]) * 1_000_000,
        })
        del df; gc.collect()
    return records


def load_hse_sl_records(results_root):
    pred_cols = {
        "pmr":  "pmr_pred_10y",
        "ukb":  "ukb_pred_10y",
        "wukb": "ukbw_hse_superlearner_pred_10y",
    }
    records = []
    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            continue
        rates = compute_rates(df, "event_10y", pred_cols)
        records.append({
            "cause":        cause,
            "cause_label":  CAUSE_LABELS[cause],
            "pmr_bias_per_mil":  (rates["pmr_pred_rate"]  - rates["observed_rate"]) * 1_000_000,
            "ukb_bias_per_mil":  (rates["ukb_pred_rate"]  - rates["observed_rate"]) * 1_000_000,
            "wukb_bias_per_mil": (rates["wukb_pred_rate"] - rates["observed_rate"]) * 1_000_000,
        })
        del df; gc.collect()
    return records


def _draw_panel(ax, rdf, has_wukb, wukb_col, wukb_marker, wukb_label,
                show_yticks=True, arrow_causes=None):
    """arrow_causes: set of cause strings to draw arrows for; None = all; empty set = none."""
    y = np.arange(len(rdf))
    for i, row in rdf.reset_index(drop=True).iterrows():
        ax.scatter(row["pmr_bias_per_mil"], i, c=C_PMR, marker=MARKERS["pmr"],
                   s=50, zorder=5, edgecolors="none")
        ax.scatter(row["ukb_bias_per_mil"], i, c=C_UKB, marker=MARKERS["ukb"],
                   s=50, zorder=5, edgecolors="none")
        if has_wukb and "wukb_bias_per_mil" in row:
            draw_arrow = (arrow_causes is None) or (row["cause"] in arrow_causes)
            if draw_arrow:
                ukb_val  = row["ukb_bias_per_mil"]
                wukb_val = row["wukb_bias_per_mil"]
                if np.isfinite(ukb_val) and np.isfinite(wukb_val):
                    ax.annotate("", xy=(wukb_val, i), xytext=(ukb_val, i),
                                arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8))
            if np.isfinite(row["wukb_bias_per_mil"]):
                ax.scatter(row["wukb_bias_per_mil"], i, c=wukb_col, marker=wukb_marker,
                           s=68, zorder=6, edgecolors="none")

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
    output_dir   = get_output_dir("fig7_sensitivity_hse_ll")

    cause_rank = {c: i for i, c in enumerate(CAUSE_ORDER)}

    ll_records = load_hse_ll_records(results_root)
    sl_records = load_hse_sl_records(results_root)

    if not ll_records:
        print("No HSE LL data found.")
        return

    ll_df = pd.DataFrame(ll_records)
    ll_df["_rank"] = ll_df["cause"].map(cause_rank)
    ll_df = ll_df.sort_values("_rank").reset_index(drop=True)

    sl_df = pd.DataFrame(sl_records)
    sl_df["_rank"] = sl_df["cause"].map(cause_rank)
    sl_df = sl_df.sort_values("_rank").reset_index(drop=True)

    n = len(ll_df)
    fig_height = max(3, n * 0.45 + 0.6)
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(8, fig_height),
                                             gridspec_kw={"width_ratios": [1, 1]},
                                             sharey=True)

    _draw_panel(ax_left,  ll_df, has_wukb=True,
                wukb_col=COLOURS["ukbw_hse_lassologit"],
                wukb_marker=MARKERS["ukbw_hse_lassologit"],
                wukb_label="wUKB (HSE LL)", show_yticks=True,
                arrow_causes={"all_cause_mortality", "cancer_mortality"})
    _draw_panel(ax_right, sl_df, has_wukb=True,
                wukb_col=COLOURS["ukbw_hse_superlearner"],
                wukb_marker=MARKERS["ukbw_hse_superlearner"],
                wukb_label="wUKB (HSE SL)", show_yticks=False,
                arrow_causes=None)

    fig.tight_layout()
    save_fig(fig, output_dir, "fig7_sensitivity_hse_ll")

    # Per-panel legends are saved separately so the two weighting models can be labelled independently
    fig_leg1, ax_leg1 = plt.subplots(figsize=(3.0, 1.0))
    ax_leg1.scatter([], [], c=C_PMR, marker=MARKERS["pmr"], s=50, edgecolors="none", label="PMR")
    ax_leg1.scatter([], [], c=C_UKB, marker=MARKERS["ukb"],  s=50, edgecolors="none", label="UKB")
    ax_leg1.scatter([], [], c=COLOURS["ukbw_hse_lassologit"],
                    marker=MARKERS["ukbw_hse_lassologit"], s=68,
                    edgecolors="none", label="wUKB (HSE LL)")
    ax_leg1.legend(fontsize=10, ncol=1, loc="center", frameon=False)
    ax_leg1.axis("off")
    fig_leg1.tight_layout()
    leg1_path = output_dir / "fig7_sensitivity_hse_ll__legend_left.svg"
    fig_leg1.savefig(leg1_path, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig_leg1)
    print(f"  Saved legend (left) to {leg1_path}")

    fig_leg2, ax_leg2 = plt.subplots(figsize=(3.0, 1.0))
    ax_leg2.scatter([], [], c=C_PMR, marker=MARKERS["pmr"], s=50, edgecolors="none", label="PMR")
    ax_leg2.scatter([], [], c=C_UKB, marker=MARKERS["ukb"],  s=50, edgecolors="none", label="UKB")
    ax_leg2.scatter([], [], c=COLOURS["ukbw_hse_superlearner"],
                    marker=MARKERS["ukbw_hse_superlearner"], s=68,
                    edgecolors="none", label="wUKB (HSE SL)")
    ax_leg2.legend(fontsize=10, ncol=1, loc="center", frameon=False)
    ax_leg2.axis("off")
    fig_leg2.tight_layout()
    leg2_path = output_dir / "fig7_sensitivity_hse_ll__legend_right.svg"
    fig_leg2.savefig(leg2_path, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig_leg2)
    print(f"  Saved legend (right) to {leg2_path}")

    print("Done: fig7_sensitivity_hse_ll")


if __name__ == "__main__":
    main()
