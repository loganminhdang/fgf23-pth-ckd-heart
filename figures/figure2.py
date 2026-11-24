"""Figure 2: FGF23 and LV structure (panel functions and source data).

  a  Adjusted LV mass index against FGF23 (partial residuals), with the fitted
     line and its 95% confidence band.
  b  Per-SD associations of log10 FGF23 with LV mass index, wall thickness,
     relative wall thickness, LV internal dimensions and fractional shortening.
  c  The FGF23 association at each echocardiographic visit, and its
     between-patient and within-patient parts.
  d  Sensitivity to assay handling (from outputs/07_figure2d.csv).

submission_figures.py draws the submitted figure with these panel functions.
Running this script writes the source data of panels a to c (and a draft
layout of the figure). The points in panel a are individual echocardiograms
and are not written to the source data.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from kipanel import (ACCENT, HAIR, LW_DATA, NEUTRAL, OUT, PRIMARY, TABLES, forest, frame,
                     letters, mm_figure, plt, save)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
import prep  # noqa: E402

COVS = prep.PRIMARY_COVS
BAND_FC, DOT = "#E7C9C7", "#B4BAC2"
XLAB_G = "Change in LV mass index (g/m\u00b2) per 10-fold FGF23"
XLAB_Z = "Standardized change per SD log\u2081\u2080 FGF23"
PHENOTYPES = [("LVMI", "LV mass index"), ("IVSd", "Septum, diastole"),
              ("IVSs", "Septum, systole"), ("PWTd", "Posterior wall, diastole"),
              ("RWT", "Relative wall thickness"), ("LVIDd", "LV diameter, diastole"),
              ("LVIDs", "LV diameter, systole"), ("FS", "Fractional shortening")]


def term(frame_, outcome, exposures, index):
    """Estimate and 95% CI of exposures[index], adjusted for the primary
    covariates, with patient-clustered standard errors."""
    y = frame_[outcome].to_numpy(float)
    columns = [frame_[c] for c in COVS] + [frame_[e] for e in exposures]
    beta, V = prep.cluster_ols(y, columns, frame_["Patient_ID"].to_numpy())
    i = 1 + len(COVS) + index
    se = float(np.sqrt(V[i, i]))
    return beta, V, float(beta[i]), float(beta[i] - 1.96 * se), float(beta[i] + 1.96 * se)


def panel_a(ax, qc):
    """Partial-residual plot of LV mass index against FGF23."""
    s = qc.dropna(subset=["LVMI", "logFGF23"] + COVS)
    beta, V, b, lo, hi = term(s, "LVMI", ["logFGF23"], 0)
    i = 1 + len(COVS)
    se = float(np.sqrt(V[i, i]))
    X = np.column_stack([np.ones(len(s))] + [s[c] for c in COVS] + [s["logFGF23"]])
    y, x = s["LVMI"].to_numpy(float), s["logFGF23"].to_numpy(float)
    adjusted = (y - X @ beta) + b * (x - x.mean()) + y.mean()
    grid = np.linspace(x.min(), x.max(), 120)
    fit = y.mean() + b * (grid - x.mean())
    half = 1.96 * se * np.abs(grid - x.mean())
    ax.scatter(10 ** x, adjusted, s=2.4, color=DOT, lw=0, zorder=2)
    ax.fill_between(10 ** grid, fit - half, fit + half, color=BAND_FC, alpha=0.55, lw=0, zorder=3)
    ax.plot(10 ** grid, fit, color=ACCENT, lw=LW_DATA, zorder=5)
    ax.set_xscale("log")
    ax.set_ylabel("Adjusted LV mass index (g/m\u00b2)")
    ax.grid(color=HAIR, lw=0.35, alpha=0.9)
    frame(ax, "FGF23 (pg/ml, log scale)", "No association with LV mass")
    return pd.DataFrame({"fgf23_pg_ml": (10 ** grid).round(3), "fitted": fit.round(4),
                         "ci_low": (fit - half).round(4), "ci_high": (fit + half).round(4)})


def panel_b(ax, qc):
    """Forest plot: FGF23 against each cardiac measure, in SD units."""
    rows, records = [], []
    for column, label in PHENOTYPES:
        s = qc.dropna(subset=[column, "logFGF23"] + COVS).copy()
        s["z_y"] = (s[column] - s[column].mean()) / s[column].std()
        s["z_x"] = (s["logFGF23"] - s["logFGF23"].mean()) / s["logFGF23"].std()
        *_, b, lo, hi = term(s, "z_y", ["z_x"], 0)
        rows.append((label, b, lo, hi, NEUTRAL))
        records.append({"phenotype": label, "estimate": b, "ci_low": lo, "ci_high": hi,
                        "n": len(s), "patients": s["Patient_ID"].nunique()})
    forest(ax, rows, XLAB_Z, "Across cardiac phenotypes", sep_after=(0, 4, 6))
    return pd.DataFrame(records)


def panel_c(ax, qc):
    """Forest plot: FGF23 by visit, and between and within patients."""
    rows, records = [], []
    for visit in (1, 3, 5):
        s = qc[qc["visit"] == visit].dropna(subset=["LVMI", "logFGF23"] + COVS)
        *_, b, lo, hi = term(s, "LVMI", ["logFGF23"], 0)
        rows.append((f"Visit {visit}", b, lo, hi, NEUTRAL))
        records.append({"specification": f"Visit {visit}", "estimate": b, "ci_low": lo,
                        "ci_high": hi, "n": len(s)})
    s = qc.dropna(subset=["LVMI", "logFGF23"] + COVS).copy()
    s["between"] = s.groupby("Patient_ID")["logFGF23"].transform("mean")
    s["within"] = s["logFGF23"] - s["between"]
    for index, label in ((0, "Between patients"), (1, "Within patients")):
        *_, b, lo, hi = term(s, "LVMI", ["between", "within"], index)
        rows.append((label, b, lo, hi, PRIMARY))
        records.append({"specification": label, "estimate": b, "ci_low": lo, "ci_high": hi,
                        "n": len(s)})
    forest(ax, rows, XLAB_G, "By occasion and by contrast", sep_after=(2,))
    return pd.DataFrame(records)


def panel_d(ax):
    """Forest plot: assay sensitivity analyses."""
    t = pd.read_csv(TABLES / "07_figure2d.csv")
    rows = [(r.Specification, r.estimate, r.ci_low, r.ci_high, NEUTRAL) for r in t.itertuples()]
    forest(ax, rows, XLAB_G, "Assay-aware sensitivity")


def main():
    plt.rcParams.update({"savefig.bbox": "tight", "savefig.pad_inches": 0.015})
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    fig = mm_figure(180.0, 96.0)
    gs = fig.add_gridspec(nrows=2, ncols=2, left=0.085, right=0.985, top=0.89, bottom=0.10,
                          wspace=0.50, hspace=0.66)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    fit = panel_a(axes[0], qc)
    phenotypes = panel_b(axes[1], qc)
    repeated = panel_c(axes[2], qc)
    panel_d(axes[3])
    letters(fig, axes, "abcd", dx=0.052, dy=0.030)
    save(fig, "Figure2")
    fit.to_csv(TABLES / "figure2a_fit.csv", index=False)
    phenotypes.round(4).to_csv(TABLES / "figure2b_phenotypes.csv", index=False)
    repeated.round(4).to_csv(TABLES / "figure2c_repeated.csv", index=False)


if __name__ == "__main__":
    main()
