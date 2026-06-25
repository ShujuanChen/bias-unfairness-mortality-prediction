#!/usr/bin/env python3
"""Figure 6 and Figures S9-S12 - bias correction by strata (UKB, PMR, weighted UKB HSE Super Learner)."""

import gc
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
from helpers import *

DARK_BLUE = COLOURS["ukbw_hse_superlearner"]

VAR_GROUP_NAMES = {
    "disability": "Disability",
    "education":  "Education",
    "tenure":     "Tenure",
    "ruralurban": "Rural/Urban",
    "imd_decile": "Deprivation index",
}

LEVEL_ANNOTATIONS = {
    ("education", "Level 1"): "Level 1 (lowest)",
    ("education", "Level 4"): "Level 4 (highest)",
    ("imd_decile", "1"):      "1 (most deprived)",
    ("imd_decile", "10"):     "10 (least deprived)",
}

PAD = 0.45   # y-axis padding above first and below last row


def _col_to_str(series):
    """Convert column to string, coercing whole-number floats (e.g. 1.0 -> '1')."""
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


def _get_display(var, level):
    if (var, level) in LEVEL_ANNOTATIONS:
        return LEVEL_ANNOTATIONS[(var, level)]
    return LEVEL_DISPLAY_NAMES.get(level, level)


def main():
    results_root = get_results_root()
    output_dir = get_output_dir("fig6_bias_correction_by_strata")

    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            print(f"[SKIP] {cause}")
            continue

        for hz in [5, 10]:
            pred_cols = {
                "pmr":  f"pmr_pred_{hz}y",
                "ukb":  f"ukb_pred_{hz}y",
                "wukb": f"ukbw_hse_superlearner_pred_{hz}y",
            }
            event_col = f"event_{hz}y"

            strata_labels = build_strata_labels(STRATA_SUBSET)
            records = []
            prev_var = None

            for var, level in strata_labels:
                if var == "__spacer__":
                    records.append({"var": "__spacer__", "level": "", "label": "",
                                    "pmr_bias_per_mil": np.nan, "ukb_bias_per_mil": np.nan,
                                    "wukb_bias_per_mil": np.nan, "bias_correction_pct": np.nan})
                    continue
                if prev_var is None:
                    # Header spacer before the very first group
                    records.append({"var": "__spacer__", "level": "", "label": "",
                                    "pmr_bias_per_mil": np.nan, "ukb_bias_per_mil": np.nan,
                                    "wukb_bias_per_mil": np.nan, "bias_correction_pct": np.nan})
                prev_var = var

                # Coerce to string so float-coded deprivation deciles (1.0) match level keys ("1")
                mask = _col_to_str(df[var]) == level
                if mask.sum() < 50:
                    # Keep the row as NaN so every cause has identical strata structure
                    # (small subgroups are suppressed but the layout stays aligned)
                    records.append({
                        "var": var, "level": level, "label": level,
                        "pmr_bias_per_mil": np.nan, "ukb_bias_per_mil": np.nan,
                        "wukb_bias_per_mil": np.nan, "bias_correction_pct": np.nan,
                        "observed_rate": np.nan, "pmr_pred_rate": np.nan,
                        "ukb_pred_rate": np.nan, "wukb_pred_rate": np.nan,
                        "n": int(mask.sum()), "n_events": 0,
                    })
                    continue

                rates = compute_rates(df[mask], event_col, pred_cols)
                ukb_bias  = (rates["ukb_pred_rate"]  - rates["observed_rate"]) * 1_000_000
                wukb_bias = (rates["wukb_pred_rate"] - rates["observed_rate"]) * 1_000_000
                pmr_bias  = (rates["pmr_pred_rate"]  - rates["observed_rate"]) * 1_000_000
                ukb_abs   = abs(ukb_bias)  if np.isfinite(ukb_bias)  else 0
                wukb_abs  = abs(wukb_bias) if np.isfinite(wukb_bias) else 0
                bc_pct    = (ukb_abs - wukb_abs) / ukb_abs * 100 if ukb_abs > 0 else np.nan
                records.append({
                    "var": var, "level": level, "label": level,
                    "pmr_bias_per_mil":    pmr_bias,
                    "ukb_bias_per_mil":    ukb_bias,
                    "wukb_bias_per_mil":   wukb_bias,
                    "bias_correction_pct": bc_pct,
                    "observed_rate":       rates["observed_rate"],
                    "pmr_pred_rate":       rates["pmr_pred_rate"],
                    "ukb_pred_rate":       rates["ukb_pred_rate"],
                    "wukb_pred_rate":      rates["wukb_pred_rate"],
                    "n": rates["n"], "n_events": rates["n_events"],
                })

            # Label each spacer row with the name of the group that follows it
            for i, rec in enumerate(records):
                if rec["var"] == "__spacer__":
                    for j in range(i + 1, len(records)):
                        if records[j]["var"] != "__spacer__":
                            rec["_group_header"] = VAR_GROUP_NAMES.get(
                                records[j]["var"], records[j]["var"])
                            break
                    else:
                        rec["_group_header"] = ""

            rdf = pd.DataFrame(records).reset_index(drop=True)
            n = len(rdf)
            fig_height = max(3, n * 0.22 + 0.5)

            fig, (ax, ax_bc) = plt.subplots(
                1, 2, figsize=(8.0, fig_height),
                gridspec_kw={"width_ratios": [3, 1.2], "wspace": 0.08},
                sharey=True,
            )

            show_pmr = (cause == "all_cause_mortality")

            for i, row in rdf.iterrows():
                if row["var"] == "__spacer__":
                    continue
                pmr_val  = row["pmr_bias_per_mil"]
                ukb_val  = row["ukb_bias_per_mil"]
                wukb_val = row["wukb_bias_per_mil"]

                if show_pmr and np.isfinite(pmr_val):
                    ax.scatter(pmr_val, i, c=C_PMR, marker=MARKERS["pmr"],
                               s=45, zorder=5, edgecolors="none")
                if np.isfinite(ukb_val):
                    ax.scatter(ukb_val, i, c=C_UKB, marker=MARKERS["ukb"],
                               s=45, zorder=5, edgecolors="none")
                if np.isfinite(ukb_val) and np.isfinite(wukb_val):
                    ax.annotate("", xy=(wukb_val, i), xytext=(ukb_val, i),
                                arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8))
                if np.isfinite(wukb_val):
                    ax.scatter(wukb_val, i, c=DARK_BLUE,
                               marker=MARKERS["ukbw_hse_superlearner"],
                               s=45, zorder=7, edgecolors="none")

                bc = row["bias_correction_pct"]
                if np.isfinite(bc):
                    ax_bc.barh(i, bc, height=0.45, color=DARK_BLUE, zorder=5, left=0)

            _y_bot = n - 0.5 + PAD
            _y_top = -0.5 - PAD
            ax.plot([0, 0], [_y_bot, _y_top], color="grey", linewidth=1.6, alpha=0.7, zorder=3)
            ax_bc.plot([0, 0], [_y_bot, _y_top], color="grey", linewidth=1.6, alpha=0.7, zorder=3)

            data_rows = [i for i, row in rdf.iterrows() if row["var"] != "__spacer__"]
            ax.set_yticks(data_rows)
            ax.set_yticklabels([""] * len(data_rows))
            ax.tick_params(axis="y", length=3.5, width=0.8)
            ax.invert_yaxis()
            ax.set_ylim(_y_bot, _y_top)
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

            # Labels: group headers and level names right-aligned outside left axis edge
            label_trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)

            for i, row in rdf.iterrows():
                if row["var"] == "__spacer__":
                    grp = row.get("_group_header", "")
                    if grp:
                        ax.text(-0.02, i, grp, transform=label_trans,
                                ha="right", va="center", fontweight="bold",
                                fontsize=8.5, clip_on=False)
                else:
                    display = _get_display(row["var"], row["level"])
                    ax.text(-0.02, i, "  " + display, transform=label_trans,
                            ha="right", va="center", fontsize=8.5, clip_on=False)

            ax.set_xlabel("Prediction bias (deaths per 1,000,000)", fontsize=9)
            ax_bc.set_xlabel("Bias correction (%)", fontsize=9)
            ax_bc.set_xlim(0, 100)
            ax_bc.set_xticks([0, 25, 50, 75, 100])
            ax_bc.tick_params(axis="y", left=False, labelleft=False)

            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(1.0)
            for spine in ax_bc.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(1.0)

            fig.subplots_adjust(left=0.22, right=0.97, top=0.97, bottom=0.12, wspace=0.08)

            fname = f"fig6_bias_correction_by_strata__{cause}__{hz}y"
            save_fig(fig, output_dir, fname)
            plt.close(fig)
            print(f"  Saved {fname}")

        del df; gc.collect()

    print("Done: fig6_bias_correction_by_strata")


if __name__ == "__main__":
    main()
