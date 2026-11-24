"""Shape of the PTH-LV mass association (Figure 4a, b and e; Supplementary
Table S4e; Supplementary Figure S2b).

Spline model: LV mass index is modelled on log10 PTH with a restricted cubic
spline (5 knots, at the 5th, 27.5th, 50th, 72.5th and 95th percentiles of
log10 PTH in the echocardiograms with normal LV geometry), adjusted for the
primary covariates, with patient-clustered standard errors. The same knots are
used for all echocardiograms (Figure 4a) and for normal geometry only
(Figure 4b). Departure from a straight line is tested with a Wald test of the
non-linear spline terms.

From the Figure 4a model:
  Table S4e      LV mass index at PTH 100, 150, 200 and 300 pg/mL compared
                 with 65 pg/mL (300 pg/mL lies beyond the last knot, where the
                 spline is a straight line)
  Figure S2b     the same comparison as a curve
  PTH threshold  the PTH concentration above which the 95% CI of the slope
                 excludes zero (descriptive)

Figure 4e: echocardiograms grouped by PTH relative to the upper reference
limit of the assay (68.3 pg/mL): up to 1, >1 to 2, >2 to 4 and >4 times the
limit; each group compared with the first, adjusted for the primary
covariates.

The script stops if a P value or estimate differs from the reported value.

Writes:
  outputs/05_spline_knots.csv          knots (Supplementary Methods S4)
  outputs/05_nonlinearity_tests.csv    P for non-linearity (Figure 4a, b)
  outputs/05_figure4ab_curves.csv      adjusted LV mass index across PTH (Figure 4a, b)
  outputs/05_figure4ab_strips.csv      number of echocardiograms across PTH (Figure 4a, b)
  outputs/05_spline_contrasts.csv      Table S4e
  outputs/05_spline_curve.csv          curve for Supplementary Figure S2b
  outputs/05_pth_histogram.csv         number of echocardiograms across PTH (Figure S2b)
  outputs/05_pth_uln_categories.csv    Figure 4e
  and Table S4e and Figure 4e unrounded in outputs/unrounded/.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"
KNOT_QUANTILES = [0.05, 0.275, 0.50, 0.725, 0.95]
REFERENCE_PTH = 65.0
CONTRAST_PTH = [100, 150, 200, 300]

ULN_GROUPS = [("Up to 1x", 0.0, 1.0), ("1-2x", 1.0, 2.0), ("2-4x", 2.0, 4.0), ("Above 4x", 4.0, np.inf)]
REPORTED_ULN = {"1-2x": (-0.18, -4.26, 3.9), "2-4x": (3.87, -0.98, 8.72),   # vs up to 1x, g/m2
                 "Above 4x": (15.65, 7.29, 24.02)}
REPORTED_NONLINEARITY = (0.032, 0.158)               # Figure 4a, Figure 4b (P rounded to 3 decimals)
REPORTED_S4E = {100: (-0.14, -3.44, 3.16),
                 150: (1.35, -3.49, 6.18),
                 200: (5.15, 0.14, 10.15),
                 300: (12.0, 5.95, 18.06)}


def rcs_basis(x, knots):
    """Non-linear terms of a restricted cubic spline (Harrell's formulation)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    scale = (knots[-1] - knots[0]) ** 2
    cube = lambda u: np.clip(u, 0, None) ** 3
    return np.column_stack([
        (cube(x - knots[j]) - cube(x - knots[-2]) * (knots[-1] - knots[j]) / (knots[-1] - knots[-2])
         + cube(x - knots[-1]) * (knots[-2] - knots[j]) / (knots[-1] - knots[-2])) / scale
        for j in range(len(knots) - 2)])


def rcs_derivative(x, knots):
    """Slope of each non-linear spline term with respect to log10 PTH."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    scale = (knots[-1] - knots[0]) ** 2
    square = lambda u: 3 * np.clip(u, 0, None) ** 2
    return np.column_stack([
        (square(x - knots[j]) - square(x - knots[-2]) * (knots[-1] - knots[j]) / (knots[-1] - knots[-2])
         + square(x - knots[-1]) * (knots[-2] - knots[j]) / (knots[-1] - knots[-2])) / scale
        for j in range(len(knots) - 2)])


def cluster_fit(X, y, groups):
    """Linear regression with patient-clustered standard errors (the same
    calculation as prep.cluster_ols, with the design matrix given directly)."""
    n, k = X.shape
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y - X @ beta
    bread = np.linalg.solve(X.T @ X, np.eye(k))
    bread = (bread + bread.T) / 2
    meat = np.zeros((k, k))
    for g in np.unique(groups):
        u = X[groups == g].T @ residual[groups == g]
        meat += np.outer(u, u)
    G = len(np.unique(groups))
    V = bread @ meat @ bread * (G / (G - 1)) * ((n - 1) / (n - k))
    return beta, V


def fit_spline(frame, knots):
    """Fit the spline model; return the coefficients, their covariance, the
    positions of the PTH terms and the P value for non-linearity."""
    covariates = frame[prep.PRIMARY_COVS].to_numpy(float)
    X = np.column_stack([np.ones(len(frame)), covariates, frame["logPTH"].to_numpy(float),
                         rcs_basis(frame["logPTH"], knots)])
    beta, V = cluster_fit(X, frame["LVMI"].to_numpy(float), frame["Patient_ID"].to_numpy())
    first = 1 + len(prep.PRIMARY_COVS)                 # index of the linear PTH term
    exposure = np.arange(first, X.shape[1])
    nonlinear = exposure[1:]
    b = beta[nonlinear]
    wald = float(b @ np.linalg.solve(V[np.ix_(nonlinear, nonlinear)], b))
    p_nonlinear = 1 - stats.chi2.cdf(wald, len(nonlinear))
    return beta, V, exposure, p_nonlinear


def uln_categories(d):
    """Figure 4e: LV mass index by PTH relative to the upper reference limit,
    each group compared with PTH up to the limit."""
    d = d.copy()
    ratio = d["PTH"] / prep.PTH_UPPER_REFERENCE
    d["group"] = pd.Series(pd.NA, index=d.index, dtype="object")
    for i, (_, lo, hi) in enumerate(ULN_GROUPS):
        d.loc[(ratio > lo) & (ratio <= hi), "group"] = f"g{i}"
    assert d["group"].notna().all()
    formula = "LVMI ~ C(group, Treatment('g0')) + " + " + ".join(prep.PRIMARY_COVS)
    model = smf.ols(formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["Patient_ID"]})
    rows, raw = [], []
    for i, (label, lo, hi) in enumerate(ULN_GROUPS):
        g = d[d["group"] == f"g{i}"]
        upper = "" if np.isinf(hi) else f"{hi * prep.PTH_UPPER_REFERENCE:.1f}"
        row = {"group": label, "pth_from_pg_ml": round(lo * prep.PTH_UPPER_REFERENCE, 1),
               "pth_to_pg_ml": upper, "n": len(g), "patients": g["Patient_ID"].nunique(),
               "mean_lvmi": round(g["LVMI"].mean(), 1)}
        if i == 0:
            row.update(estimate=0.0, ci_low=np.nan, ci_high=np.nan, p=np.nan)
        else:
            term = f"C(group, Treatment('g0'))[T.g{i}]"
            lo_ci, hi_ci = model.conf_int().loc[term]
            raw.append({"group": label, "n": len(g), "patients": g["Patient_ID"].nunique(),
                        "estimate": model.params[term], "ci_low": lo_ci, "ci_high": hi_ci,
                        "p": model.pvalues[term]})
            row.update(estimate=round(model.params[term], 2), ci_low=round(lo_ci, 2),
                       ci_high=round(hi_ci, 2), p=float(f"{model.pvalues[term]:.3g}"))
            reported = REPORTED_ULN.get(label)
            if reported is not None:
                assert np.allclose([row["estimate"], row["ci_low"], row["ci_high"]], reported,
                                   atol=0.006), (label, row, reported)
        rows.append(row)
    assert len(REPORTED_ULN) == 3, "PTH upper-reference-limit groups are not pinned"
    return pd.DataFrame(rows), pd.DataFrame(raw)


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    d = qc.dropna(subset=["LVMI", "logPTH"] + prep.PRIMARY_COVS).reset_index(drop=True)
    normal = d[d["geometry"] == "Normal geometry"].reset_index(drop=True)
    # Knots: percentiles of log10 PTH in the normal-geometry echocardiograms.
    figure4_knots = np.quantile(normal["logPTH"], KNOT_QUANTILES)
    OUT.mkdir(exist_ok=True)
    pd.DataFrame({"percentile": [round(100 * q, 1) for q in KNOT_QUANTILES],
                  "log10_pth": figure4_knots.round(6),
                  "pth_pg_ml": (10 ** figure4_knots).round(3)}).to_csv(
        OUT / "05_spline_knots.csv", index=False)   # knots quoted in Supplementary Methods S4
    _, _, _, p_all = fit_spline(d, figure4_knots)
    _, _, _, p_normal = fit_spline(normal, figure4_knots)
    assert (round(p_all, 3), round(p_normal, 3)) == REPORTED_NONLINEARITY, (p_all, p_normal)

    # Figure 4a, b: adjusted LV mass index (covariates at their means) from the
    # 1st to the 99th percentile of PTH, and the number of echocardiograms
    # across PTH, drawn below the curve.
    curves, strips = [], []
    for panel, frame in (("4a", d), ("4b", normal)):
        b4, V4, e4, _ = fit_spline(frame, figure4_knots)
        x = frame["logPTH"].to_numpy(float)
        grid = np.linspace(np.quantile(x, 0.01), np.quantile(x, 0.99), 160)
        rows4 = np.column_stack([np.ones(len(grid))]
                                + [np.full(len(grid), frame[c].mean()) for c in prep.PRIMARY_COVS]
                                + [grid, rcs_basis(grid, figure4_knots)])
        pred = rows4 @ b4
        se = np.sqrt(np.einsum("ij,jk,ik->i", rows4, V4, rows4))
        curves.append(pd.DataFrame({"figure": panel, "pth_pg_ml": (10 ** grid).round(3),
                                    "adjusted_lvmi": pred.round(4),
                                    "ci_low": (pred - 1.96 * se).round(4),
                                    "ci_high": (pred + 1.96 * se).round(4)}))
        counts, edges = np.histogram(x, bins=np.linspace(x.min(), x.max(), 26))
        strips.append(pd.DataFrame({"figure": panel, "log10_pth_left": edges[:-1].round(5),
                                    "log10_pth_right": edges[1:].round(5), "count": counts}))

    # Table S4e and Supplementary Figure S2b use the Figure 4a model.
    knots = figure4_knots
    beta, V, exposure, _ = fit_spline(d, knots)
    tests = pd.DataFrame([
        {"figure": "4a", "sample": "All observations", "n": len(d),
         "patients": d.Patient_ID.nunique(), "nonlinearity_p": round(p_all, 4)},
        {"figure": "4b", "sample": "Normal LV geometry", "n": len(normal),
         "patients": normal.Patient_ID.nunique(), "nonlinearity_p": round(p_normal, 4)}])

    def exposure_row(log_pth):
        """Values of the PTH terms of the model at one PTH concentration."""
        row = np.zeros(len(beta))
        row[exposure[0]] = log_pth
        row[exposure[1:]] = rcs_basis(log_pth, knots)[0]
        return row

    # Table S4e: each PTH concentration against 65 pg/mL.
    reference = exposure_row(np.log10(REFERENCE_PTH))
    pth = 10 ** d["logPTH"]
    rows, raw = [], []
    for value in CONTRAST_PTH:
        delta = exposure_row(np.log10(value)) - reference
        estimate, se = float(delta @ beta), float(np.sqrt(delta @ V @ delta))
        percentile = float(stats.percentileofscore(pth, value))
        rows.append({"contrast": f"PTH {value} vs {REFERENCE_PTH:.0f} pg/mL",
                     "pth_pg_ml": value,
                     "cohort_percentile": round(percentile, 1),
                     "estimate": round(estimate, 2),
                     "ci_low": round(estimate - 1.96 * se, 2),
                     "ci_high": round(estimate + 1.96 * se, 2),
                     "beyond_outer_knot": bool(np.log10(value) > knots[-1])})
        raw.append({**rows[-1], "cohort_percentile": percentile, "estimate": estimate,
                    "ci_low": estimate - 1.96 * se, "ci_high": estimate + 1.96 * se})
        reported = REPORTED_S4E[value]
        assert np.allclose([rows[-1]["estimate"], rows[-1]["ci_low"], rows[-1]["ci_high"]],
                           reported, atol=0.006), (rows[-1], reported)
    contrasts, contrasts_raw = pd.DataFrame(rows), pd.DataFrame(raw)

    # PTH concentration above which the 95% CI of the slope excludes zero.
    grid = np.linspace(knots[0], knots[-1], 200)
    D = np.zeros((len(grid), len(beta)))
    D[:, exposure[0]] = 1.0
    D[:, exposure[1:]] = rcs_derivative(grid, knots)
    slope = D @ beta
    slope_se = np.sqrt(np.einsum("ij,jk,ik->i", D, V, D))
    positive = slope - 1.96 * slope_se > 0
    onset = float(10 ** grid[np.argmax(positive)]) if positive.any() else np.nan
    contrasts["slope_ci_excludes_zero_above_pg_ml"] = round(onset)
    contrasts_raw["slope_ci_excludes_zero_above_pg_ml"] = onset

    # Supplementary Figure S2b: difference from PTH 65 pg/mL across PTH.
    x = np.linspace(np.log10(REFERENCE_PTH * 0.6),
                    max(knots[-1], np.log10(max(CONTRAST_PTH))) + 0.02, 160)
    C = np.array([exposure_row(v) - reference for v in x])
    estimate = C @ beta
    se = np.sqrt(np.einsum("ij,jk,ik->i", C, V, C))
    curve = pd.DataFrame({"log10_pth": x.round(5), "difference": estimate.round(4),
                          "ci_low": (estimate - 1.96 * se).round(4),
                          "ci_high": (estimate + 1.96 * se).round(4)})

    counts, edges = np.histogram(d["logPTH"], bins=40)
    histogram = pd.DataFrame({"log10_pth_left": edges[:-1].round(5),
                              "log10_pth_right": edges[1:].round(5), "count": counts})

    uln, uln_raw = uln_categories(d)

    # Write the tables.
    OUT.mkdir(exist_ok=True)
    uln.to_csv(OUT / "05_pth_uln_categories.csv", index=False)
    (OUT / "unrounded").mkdir(exist_ok=True)
    uln_raw.to_csv(OUT / "unrounded" / "05_pth_uln_categories.csv", index=False)
    tests.to_csv(OUT / "05_nonlinearity_tests.csv", index=False)
    contrasts.to_csv(OUT / "05_spline_contrasts.csv", index=False)
    (OUT / "unrounded").mkdir(exist_ok=True)
    contrasts_raw.to_csv(OUT / "unrounded" / "05_spline_contrasts.csv", index=False)
    curve.to_csv(OUT / "05_spline_curve.csv", index=False)
    histogram.to_csv(OUT / "05_pth_histogram.csv", index=False)
    pd.concat(curves).to_csv(OUT / "05_figure4ab_curves.csv", index=False)
    pd.concat(strips).to_csv(OUT / "05_figure4ab_strips.csv", index=False)
    print(tests.to_string(index=False), "\n")
    print(contrasts.to_string(index=False))
    print("\nKnots, Figure 4a-b, Table S4e and Supplementary Figure S2b (pg/mL):",
          [round(float(10 ** k), 1) for k in figure4_knots])
    print("Nonlinearity tests and Table S4e reproduced.")


if __name__ == "__main__":
    main()
