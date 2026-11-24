"""Mineral-metabolism markers and LV structure (Figure 3, Supplementary Table S3).

  S3a  Each marker against LV mass index, per standard deviation (SD) of the
       marker.
  S3b  PTH and FGF23 in the same model: the difference between their per-SD
       associations, with a 95% CI from 2,000 bootstrap samples of patients.
  S3c  Each marker against each wall and cavity dimension, in mm (or RWT units)
       per SD of the marker.
  S3d  PTH and phosphate in the same model.

PTH and FGF23 are analysed as log10. All models are adjusted for the primary
covariates (a covariate that is itself the marker is left out) with
patient-clustered standard errors. Figure 3b adds a fourth set, baseline
echocardiograms only; it is run after the three sets of Table S3b so that
those three use the same bootstrap samples.

The script stops if any estimate differs from the value reported in the
Supplementary Tables (the REPORTED_ values below).

Writes outputs/04_table_s3a.csv to 04_table_s3d.csv and 04_figure3b.csv
(rounded), and the same tables unrounded to outputs/unrounded/.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"
UNROUNDED = OUT / "unrounded"
SEED = 20261021
N_BOOT = 2000
TOLERANCE = 0.006

REPORTED_S3A = {'PTH': (3.44, 1.33, 5.54),
                 'Phosphate': (2.44, 0.49, 4.38),
                 'Calcium': (-1.55, -3.37, 0.27),
                 'eGFR': (-3.25, -5.95, -0.55),
                 'Albumin': (-1.93, -4.57, 0.71),
                 'FGF23': (-0.1, -1.78, 1.57)}
REPORTED_S3B = {'All observations': (3.62, 0.76, 6.33),
                 'Without LVH': (3.15, 1.02, 5.26),
                 'Normal LV geometry': (2.41, 0.23, 4.65)}

# Tables S3c and S3d are reported to three decimals (mm, or RWT units, per SD).
TOLERANCE_3DP = 0.0006
REPORTED_S3C = {
    'PTH': {'Septum, diastole': (0.203, 0.052, 0.354),
        'Septum, systole': (0.3, 0.115, 0.485),
        'Posterior wall, diastole': (0.163, 0.02, 0.307),
        'LV diameter, diastole': (0.225, -0.191, 0.641),
        'LV diameter, systole': (0.049, -0.315, 0.412),
        'Relative wall thickness': (0.006, -0.001, 0.012)},
    'Phosphate': {'Septum, diastole': (0.112, -0.025, 0.248),
        'Septum, systole': (0.016, -0.143, 0.175),
        'Posterior wall, diastole': (0.102, -0.029, 0.233),
        'LV diameter, diastole': (0.315, -0.039, 0.668),
        'LV diameter, systole': (0.339, 0.063, 0.616),
        'Relative wall thickness': (0.002, -0.004, 0.008)},
    'Calcium': {'Septum, diastole': (-0.025, -0.14, 0.091),
        'Septum, systole': (-0.007, -0.164, 0.15),
        'Posterior wall, diastole': (0.016, -0.091, 0.124),
        'LV diameter, diastole': (-0.359, -0.701, -0.017),
        'LV diameter, systole': (-0.251, -0.536, 0.034),
        'Relative wall thickness': (0.004, -0.001, 0.008)},
    'Albumin': {'Septum, diastole': (-0.017, -0.166, 0.132),
        'Septum, systole': (-0.02, -0.192, 0.152),
        'Posterior wall, diastole': (-0.011, -0.129, 0.107),
        'LV diameter, diastole': (-0.513, -0.931, -0.096),
        'LV diameter, systole': (-0.4, -0.712, -0.087),
        'Relative wall thickness': (0.003, -0.002, 0.008)},
}
REPORTED_S3D = {  # (outcome, exposure): estimate, CI; LV mass index to two decimals
    ('Septum, diastole', 'PTH'): (0.193, 0.046, 0.339),
    ('Septum, diastole', 'Phosphate'): (0.094, -0.032, 0.221),
    ('Septum, systole', 'PTH'): (0.302, 0.117, 0.487),
    ('Septum, systole', 'Phosphate'): (-0.011, -0.163, 0.141),
    ('Posterior wall, diastole', 'PTH'): (0.154, 0.015, 0.293),
    ('Posterior wall, diastole', 'Phosphate'): (0.088, -0.034, 0.211),
    ('LV diameter, diastole', 'PTH'): (0.193, -0.229, 0.615),
    ('LV diameter, diastole', 'Phosphate'): (0.297, -0.062, 0.656),
    ('LV diameter, systole', 'PTH'): (0.012, -0.354, 0.378),
    ('LV diameter, systole', 'Phosphate'): (0.338, 0.057, 0.619),
    ('Relative wall thickness', 'PTH'): (0.005, -0.002, 0.012),
    ('Relative wall thickness', 'Phosphate'): (0.001, -0.005, 0.007),
    ('LV mass index', 'PTH'): (3.21, 1.2, 5.21),
    ('LV mass index', 'Phosphate'): (2.15, 0.37, 3.93),
}


def check(got, reported, label, tolerance=TOLERANCE):
    """Stop if an estimate or CI limit differs from the reported value."""
    assert all(abs(g - p) <= tolerance for g, p in zip(got, reported)), \
        f"{label}: {tuple(round(g, 2) for g in got)} differs from reported {reported}"


MARKERS = {"PTH": "logPTH", "Phosphate": "Phosphate", "Calcium": "Calcium",
           "eGFR": "eGFR", "Albumin": "Albumin", "FGF23": "logFGF23"}
DIMENSIONS = {"Septum, diastole": "IVSd", "Septum, systole": "IVSs",
              "Posterior wall, diastole": "PWTd",
              "LV diameter, diastole": "LVIDd",
              "LV diameter, systole": "LVIDs",
              "Relative wall thickness": "RWT"}


def standardise(frame, columns):
    """Add z_<column>: the column in SD units."""
    out = frame.copy()
    for c in columns:
        out[f"z_{c}"] = (out[c] - out[c].mean()) / out[c].std()
    return out


def adjustment_set(exposure_column):
    """Primary covariates, excluding the exposure if it is one of them."""
    return [c for c in prep.PRIMARY_COVS if c != exposure_column]


def fit(outcome, exposure, data, covariates):
    """Linear regression with patient-clustered standard errors. Rows missing
    the outcome are dropped first, so that the patient clusters match the rows
    used."""
    data = data.dropna(subset=[outcome])
    formula = f"{outcome} ~ {exposure}" + ("".join(f" + {c}" for c in covariates))
    model = smf.ols(formula, data=data).fit(
        cov_type="cluster", cov_kwds={"groups": data["Patient_ID"]})
    lo, hi = model.conf_int().loc[exposure]
    return model.params[exposure], lo, hi, model.pvalues[exposure], int(model.nobs)


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    d = standardise(qc, [v for v in MARKERS.values()])
    d = d.dropna(subset=[f"z_{v}" for v in MARKERS.values()] + prep.PRIMARY_COVS)

    # S3a: each marker against LV mass index.
    rows, raw = [], []
    for name, col in MARKERS.items():
        est, lo, hi, p, n = fit("LVMI", f"z_{col}", d, adjustment_set(col))
        check((est, lo, hi), REPORTED_S3A[name], f"Table S3a, {name}")
        rows.append({"Marker": name, "n": n, "patients": d.Patient_ID.nunique(),
                     "estimate": round(est, 2), "ci_low": round(lo, 2),
                     "ci_high": round(hi, 2), "p": float(f"{p:.3g}")})
        raw.append({"Marker": name, "n": n, "patients": d.Patient_ID.nunique(),
                    "estimate": est, "ci_low": lo, "ci_high": hi, "p": p})
    a, a_raw = pd.DataFrame(rows), pd.DataFrame(raw)

    # S3b: PTH and FGF23 in one model; difference of their per-SD estimates.
    # PTH and FGF23 are put in SD units within each set of echocardiograms and
    # within each bootstrap sample.
    columns = prep.PRIMARY_COVS + ["logPTH", "logFGF23"]

    def difference(frame):
        """PTH estimate minus FGF23 estimate, both per SD."""
        x = frame[columns].to_numpy(float).copy()
        x[:, -2:] = (x[:, -2:] - x[:, -2:].mean(axis=0)) / x[:, -2:].std(axis=0, ddof=1)
        X = np.column_stack([np.ones(len(x)), x])
        beta = np.linalg.lstsq(X, frame["LVMI"].to_numpy(float), rcond=None)[0]
        return beta[-2] - beta[-1]

    rng = np.random.default_rng(SEED)
    complete = qc.dropna(subset=["LVMI"] + columns)
    subsets = {"All observations": complete,
               "Without LVH": complete[complete["LVH"].eq(0)],
               "Normal LV geometry": complete[complete["geometry"] == "Normal geometry"],
               "Baseline only": complete[complete["visit"] == 1]}
    rows, raw = [], []
    for label, frame in subsets.items():
        frame = frame.reset_index(drop=True)
        point = difference(frame)
        patients = frame["Patient_ID"].unique()
        rows_of = {pid: np.flatnonzero(frame["Patient_ID"].to_numpy() == pid)
                   for pid in patients}
        # Bootstrap: resample whole patients (all their echocardiograms).
        draws = np.empty(N_BOOT)
        for i in range(N_BOOT):
            picked = rng.choice(patients, size=len(patients), replace=True)
            draws[i] = difference(frame.iloc[np.concatenate([rows_of[p] for p in picked])])
        lo, hi = np.percentile(draws, [2.5, 97.5])
        if label in REPORTED_S3B:
            check((point, lo, hi), REPORTED_S3B[label], f"Table S3b, {label}")
        rows.append({"Observations": label, "n": len(frame),
                     "patients": frame["Patient_ID"].nunique(),
                     "difference": round(point, 2), "ci_low": round(lo, 2),
                     "ci_high": round(hi, 2)})
        raw.append({"Observations": label, "n": len(frame),
                    "patients": frame["Patient_ID"].nunique(),
                    "difference": point, "ci_low": lo, "ci_high": hi})
    figure3b, figure3b_raw = pd.DataFrame(rows), pd.DataFrame(raw)
    b = figure3b[figure3b["Observations"].isin(REPORTED_S3B)]
    b_raw = figure3b_raw[figure3b_raw["Observations"].isin(REPORTED_S3B)]

    # S3c: each marker against each wall and cavity dimension.
    rows, raw = [], []
    for marker in ["PTH", "Phosphate", "Calcium", "Albumin"]:
        col = MARKERS[marker]
        for label, outcome in DIMENSIONS.items():
            est, lo, hi, p, n = fit(outcome, f"z_{col}", d, adjustment_set(col))
            check((est, lo, hi), REPORTED_S3C[marker][label],
                  f"Table S3c, {marker}, {label}", TOLERANCE_3DP)
            rows.append({"Marker": marker, "Outcome": label, "n": n,
                         "estimate": round(est, 3), "ci_low": round(lo, 3),
                         "ci_high": round(hi, 3), "p": float(f"{p:.3g}")})
            raw.append({"Marker": marker, "Outcome": label, "n": n,
                        "patients": d.dropna(subset=[outcome])["Patient_ID"].nunique(),
                        "estimate": est, "ci_low": lo, "ci_high": hi, "p": p})
    c, c_raw = pd.DataFrame(rows), pd.DataFrame(raw)

    # S3d: PTH and phosphate in the same model.
    rows, raw = [], []
    for label, outcome in {**DIMENSIONS, "LV mass index": "LVMI"}.items():
        for exposure, other in (("z_logPTH", "z_Phosphate"),
                                ("z_Phosphate", "z_logPTH")):
            covs = adjustment_set("logPTH" if exposure == "z_logPTH" else "Phosphate")
            est, lo, hi, p, n = fit(outcome, exposure, d, covs + [other])
            name = "PTH" if exposure == "z_logPTH" else "Phosphate"
            check((est, lo, hi), REPORTED_S3D[(label, name)], f"Table S3d, {label}, {name}",
                  TOLERANCE if outcome == "LVMI" else TOLERANCE_3DP)
            rows.append({"Outcome": label,
                         "Exposure": name,
                         "n": n, "estimate": round(est, 3), "ci_low": round(lo, 3),
                         "ci_high": round(hi, 3), "p": float(f"{p:.3g}")})
            raw.append({"Outcome": label, "Exposure": name, "n": n,
                        "patients": d.dropna(subset=[outcome])["Patient_ID"].nunique(),
                        "estimate": est, "ci_low": lo, "ci_high": hi, "p": p})
    dd, dd_raw = pd.DataFrame(rows), pd.DataFrame(raw)

    # Write the rounded and unrounded tables.
    OUT.mkdir(exist_ok=True)
    UNROUNDED.mkdir(exist_ok=True)
    figure3b.to_csv(OUT / "04_figure3b.csv", index=False)
    figure3b_raw.to_csv(UNROUNDED / "04_figure3b.csv", index=False)
    for name, table, unrounded in (("a", a, a_raw), ("b", b, b_raw), ("c", c, c_raw),
                                   ("d", dd, dd_raw)):
        table.to_csv(OUT / f"04_table_s3{name}.csv", index=False)
        unrounded.to_csv(UNROUNDED / f"04_table_s3{name}.csv", index=False)
        print(f"--- Table S3{name} ---")
        print(table.to_string(index=False), "\n")
    print("Tables S3a to S3d reproduced.")


if __name__ == "__main__":
    main()
