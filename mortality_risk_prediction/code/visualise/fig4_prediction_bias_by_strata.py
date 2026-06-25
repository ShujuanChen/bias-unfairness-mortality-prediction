#!/usr/bin/env python3
"""Figure 4 and Figures S5-S8 - prediction bias by strata (UKB vs PMR, no bias correction)."""

import gc
import math
import matplotlib.ticker as mticker
from helpers import *

TICK_STEP = {
    ("all_cause_mortality",       5):  10000,
    ("all_cause_mortality",       10): 20000,
    ("cancer_mortality",          5):  2500,
    ("cancer_mortality",          10): 5000,
    ("cardiovascular_mortality",  5):  2000,
    ("cardiovascular_mortality",  10): 5000,
    ("digestive_mortality",       5):  1000,
    ("digestive_mortality",       10): 2000,
    ("respiratory_mortality",     5):  2000,
    ("respiratory_mortality",     10): 2500,
}

LEVEL_ANNOTATIONS_CSV = {
    "imd_decile": {"1": "1 (most deprived)", "10": "10 (least deprived)"},
    "education":  {"Level 1": "Level 1 (lowest)", "Level 4": "Level 4 (highest)"},
}

VAR_GROUP_NAMES = {
    "disability": "Disability",
    "education":  "Education",
    "tenure":     "Tenure",
    "ruralurban": "Rural/Urban",
    "imd_decile": "Deprivation index",
}


def _fmt_ci(center, lo, hi):
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


def _load_bootstrap_index(results_root, cause, hz):
    ci_path = results_root / "visualise" / "bootstrap_ci" / \
              f"bootstrap_ci__{cause}__{hz}y.csv"
    if not ci_path.exists():
        return {}
    df = pd.read_csv(ci_path)
    df = df[df["var"] != "overall"]
    out = {}
    for _, row in df.iterrows():
        out[(str(row["var"]), str(row["level"]))] = (
            float(row["ukb_bias_ci_lo_per_mil"]),
            float(row["ukb_bias_ci_hi_per_mil"]),
        )
    return out


def plot_one_horizon(rdf, cause, hz, output_dir):
    show_pmr = (cause == "all_cause_mortality")
    n = len(rdf)
    fig_height = max(4, n * 0.28 + 0.8)
    fig, ax = plt.subplots(figsize=(4, fig_height))
    y = np.arange(n)

    for i, row in rdf.iterrows():
        if row["var"] == "__spacer__":
            continue
        ax.plot([0, row["ukb_bias_per_mil"]], [i, i], color=C_UKB, linewidth=1.6, zorder=2)
        ax.scatter(row["ukb_bias_per_mil"], i, c=C_UKB, marker=MARKERS["ukb"],
                   s=50, zorder=4, edgecolors="none")
        if show_pmr:
            ax.plot([0, row["pmr_bias_per_mil"]], [i, i], color=C_PMR, linewidth=1.6, zorder=5)
            ax.scatter(row["pmr_bias_per_mil"], i, c=C_PMR, marker=MARKERS["pmr"],
                       s=50, zorder=6, edgecolors="none")

    ax.plot([0, 0], [n, -1.7], color="black", linewidth=2.0, zorder=7, clip_on=False)

    ax.set_yticks(y)
    ax.set_yticklabels([""] * n)
    ax.invert_yaxis()
    ax.set_ylim(n - 0.5 + 0.5, -0.5 - 0.7)

    step = TICK_STEP.get((cause, hz), 5000)
    cols = ["pmr_bias_per_mil", "ukb_bias_per_mil"] if show_pmr else ["ukb_bias_per_mil"]
    data_vals = rdf[cols].values.flatten()
    data_vals = data_vals[np.isfinite(data_vals)]
    if len(data_vals):
        d_min = min(float(np.min(data_vals)), 0)
        d_max = max(float(np.max(data_vals)), 0)
    else:
        d_min, d_max = -step * 4, 0
    pad   = (d_max - d_min) * 0.05
    x_min = d_min - pad
    x_max = d_max + pad

    first_tick = math.floor(x_min / step) * step
    last_tick  = math.ceil(x_max / step) * step
    ticks = [int(first_tick + step * i)
             for i in range(int(round((last_tick - first_tick) / step)) + 1)]
    ax.set_xticks(ticks)
    ax.set_xlim(x_min, x_max)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.spines["bottom"].set_color("black")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", width=2.0, labelsize=11)
    ax.set_xlabel("Prediction bias (deaths per 1,000,000)", fontsize=12)

    fig.tight_layout()

    fname = f"fig4_prediction_bias_by_strata__{cause}__{hz}y"
    save_fig(fig, output_dir, fname)
    print(f"  Saved {fname}")


def build_summary_csv(rdf, cause, hz, output_dir):
    csv_rows = []
    current_var = None
    for _, row in rdf.iterrows():
        if row["var"] == "__spacer__":
            continue
        if row["var"] != current_var:
            current_var = row["var"]
            csv_rows.append({
                "strata": VAR_GROUP_NAMES.get(current_var, current_var),
                "UKB bias per million (95% CI)": "",
            })
        annot = LEVEL_ANNOTATIONS_CSV.get(row["var"], {})
        display_level = annot.get(row["level"], row["level"])
        display_level = LEVEL_DISPLAY_NAMES.get(display_level, display_level)
        csv_rows.append({
            "strata": "    " + display_level,
            "UKB bias per million (95% CI)": _fmt_ci(
                row.get("ukb_bias_per_mil", np.nan),
                row.get("ukb_bias_ci_lo", np.nan),
                row.get("ukb_bias_ci_hi", np.nan),
            ),
        })
    out = pd.DataFrame(csv_rows)
    out.to_excel(output_dir / f"fig4_prediction_bias_by_strata__{cause}__{hz}y__summary.xlsx",
                 index=False)


def main():
    results_root    = get_results_root()
    base_output_dir = get_output_dir("fig4_prediction_bias_by_strata")

    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            print(f"[SKIP] {cause}")
            continue

        cause_dir = base_output_dir / cause
        cause_dir.mkdir(parents=True, exist_ok=True)

        for hz in [5, 10]:
            event_col  = f"event_{hz}y"
            pred_cols  = {"pmr": f"pmr_pred_{hz}y", "ukb": f"ukb_pred_{hz}y"}
            ci_index   = _load_bootstrap_index(results_root, cause, hz)

            strata_labels = build_strata_labels(STRATA_SUBSET)
            records = []

            for var, level in strata_labels:
                if var == "__spacer__":
                    records.append({"var": "__spacer__", "level": "", "label": "",
                                    "pmr_bias_per_mil": np.nan,
                                    "ukb_bias_per_mil": np.nan,
                                    "ukb_bias_ci_lo":   np.nan,
                                    "ukb_bias_ci_hi":   np.nan})
                    continue
                mask = df[var].astype(str) == level
                # Suppress small subgroups (< 50 members) to avoid over-interpreting
                # unstable stratum estimates.
                if mask.sum() < 50:
                    continue
                rates = compute_rates(df[mask], event_col, pred_cols)
                ci_lo, ci_hi = ci_index.get((var, level), (np.nan, np.nan))
                records.append({
                    "var":               var,
                    "level":             level,
                    "label":             level,
                    "pmr_bias_per_mil":  (rates["pmr_pred_rate"] - rates["observed_rate"]) * 1_000_000,
                    "ukb_bias_per_mil":  (rates["ukb_pred_rate"] - rates["observed_rate"]) * 1_000_000,
                    "ukb_bias_ci_lo":    ci_lo,
                    "ukb_bias_ci_hi":    ci_hi,
                    "pmr_relative_bias": rates["pmr_relative_bias"],
                    "ukb_relative_bias": rates["ukb_relative_bias"],
                    "observed_rate":     rates["observed_rate"],
                    "pmr_pred_rate":     rates["pmr_pred_rate"],
                    "ukb_pred_rate":     rates["ukb_pred_rate"],
                    "n":                 rates["n"],
                    "n_events":          rates["n_events"],
                })

            rdf = pd.DataFrame(records)
            plot_one_horizon(rdf, cause, hz, cause_dir)
            build_summary_csv(rdf, cause, hz, cause_dir)

        del df; gc.collect()

    fig_leg, ax_leg = plt.subplots(figsize=(4, 0.8))
    ax_leg.scatter([], [], c=C_PMR, marker=MARKERS["pmr"], s=50, label="PMR")
    ax_leg.scatter([], [], c=C_UKB, marker=MARKERS["ukb"], s=50, label="UKB")
    ax_leg.legend(fontsize=14, ncol=2, loc="center", frameon=False)
    ax_leg.axis("off")
    fig_leg.tight_layout()
    save_fig(fig_leg, base_output_dir, "fig4_prediction_bias_by_strata__legend")

    print("Done: fig4_prediction_bias_by_strata")


if __name__ == "__main__":
    main()
