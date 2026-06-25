#!/usr/bin/env python3
"""Figure 3 and Figure S4 - calibration of UKB vs PMR predictions, overall."""

import gc
from helpers import *

N_BINS_10 = 10
N_BINS_1 = 100

# PMR shades override the helpers defaults: lighter grey for the connecting line, darker grey for the fine background scatter.
C_PMR = "#b8b8b8"
C_PMR_SCATTER = "#3a3a3a"


def compute_calibration(pred, obs, w, n_bins):
    valid = np.isfinite(pred) & np.isfinite(obs) & np.isfinite(w)
    pred, obs, w = pred[valid], obs[valid], w[valid]
    if len(pred) < n_bins * 5:
        return None
    edges = weighted_percentile_edges(pred, w, n_bins)
    bin_idx = np.digitize(pred, edges) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    rows = []
    for b in range(n_bins):
        mask = bin_idx == b
        sw = w[mask].sum()
        if sw < 5:
            rows.append({"bin": b, "mean_pred": np.nan, "obs_rate": np.nan, "n_eff": 0})
            continue
        mp = np.sum(pred[mask] * w[mask]) / sw
        mr = np.sum(obs[mask] * w[mask]) / sw
        n_eff = sw ** 2 / np.sum(w[mask] ** 2)
        rows.append({"bin": b, "mean_pred": mp, "obs_rate": mr, "n_eff": n_eff})
    return pd.DataFrame(rows)


def plot_calibration(df, cause, output_dir, horizon):
    event_col = f"event_{horizon}"
    if event_col not in df.columns:
        return

    model_sets = {
        "pmr": (f"pmr_pred_{horizon}", C_PMR, MARKERS["pmr"], "PMR"),
        "ukb": (f"ukb_pred_{horizon}", C_UKB, MARKERS["ukb"], "UKB"),
    }

    if len(df) < 200:
        return

    obs = df[event_col].to_numpy(dtype=float)
    w = df["w"].to_numpy(dtype=float) if "w" in df.columns else np.ones(len(df))

    fig, ax = plt.subplots(figsize=(5, 5))
    main_curve_max = 0.0

    for model_key, (pred_col, colour, marker, label) in model_sets.items():
        if pred_col not in df.columns:
            continue
        pred = df[pred_col].to_numpy(dtype=float)

        cal_1p = compute_calibration(pred, obs, w, N_BINS_1)
        if cal_1p is not None:
            if model_key == "pmr":
                ax.scatter(cal_1p["mean_pred"], cal_1p["obs_rate"], c=C_PMR_SCATTER,
                           marker="x", s=20, alpha=0.7, zorder=3, linewidths=0.7)
            else:
                ax.scatter(cal_1p["mean_pred"], cal_1p["obs_rate"], c=colour, marker="x",
                           s=18, alpha=0.65, zorder=2, linewidths=0.6)

        cal_10p = compute_calibration(pred, obs, w, N_BINS_10)
        if cal_10p is not None:
            valid_main = cal_10p[np.isfinite(cal_10p["mean_pred"]) & np.isfinite(cal_10p["obs_rate"])]
            if not valid_main.empty:
                main_curve_max = max(
                    main_curve_max,
                    float(valid_main["mean_pred"].max()),
                    float(valid_main["obs_rate"].max()),
                )
            linestyle = "-." if model_key == "pmr" else "-"
            ax.plot(cal_10p["mean_pred"], cal_10p["obs_rate"], color=colour, marker=marker,
                    markersize=4, linewidth=1.2, linestyle=linestyle, label=label, zorder=5)
            if not valid_main.empty:
                se = np.sqrt(
                    valid_main["obs_rate"] * (1 - valid_main["obs_rate"]) / valid_main["n_eff"].clip(lower=1)
                )
                ax.errorbar(
                    valid_main["mean_pred"], valid_main["obs_rate"],
                    yerr=1.96 * se, fmt="none", ecolor=colour,
                    capsize=0, linewidth=1.8, zorder=6,
                )

    if main_curve_max > 0:
        axis_upper = main_curve_max * 1.10
    else:
        axis_upper = max(ax.get_xlim()[1], ax.get_ylim()[1])
        if not np.isfinite(axis_upper) or axis_upper <= 0:
            axis_upper = 1.0

    ax.plot([0, axis_upper], [0, axis_upper], color="black", linestyle="-", linewidth=0.7, alpha=0.7)
    ax.set_xlim(0, axis_upper)
    ax.set_ylim(0, axis_upper)
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=5))
    ax.set_xlabel("Predicted risk", fontsize=9)
    ax.set_ylabel("Observed risk", fontsize=9)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()

    fname = f"fig3_calibration__{cause}__{horizon}__Overall__All"
    save_fig(fig, output_dir, fname)


def main():
    results_root = get_results_root()

    for cause in CAUSES:
        df = load_cause(results_root, cause)
        if df is None:
            print(f"[SKIP] {cause}")
            continue
        print(f"Processing {cause}...")

        for horizon in ["5y", "10y"]:
            output_dir = get_output_dir("fig3_calibration") / horizon
            output_dir.mkdir(parents=True, exist_ok=True)
            plot_calibration(df, cause, output_dir, horizon=horizon)

        del df; gc.collect()

    print("Done: fig3_calibration")


if __name__ == "__main__":
    main()
