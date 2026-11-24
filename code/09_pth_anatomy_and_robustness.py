"""PTH, wall and cavity dimensions, and robustness of the PTH association
(Supplementary Tables S4a, S4b, S4d; Figure 4c; Supplementary Figure S2a).

  S4a  log10 PTH against each wall and cavity dimension, all echocardiograms.
  S4b  The same in echocardiograms with normal LV geometry.
  S4d  log10 PTH against LV mass index under alternative model specifications.

Estimates are per 10-fold higher PTH, in mm for dimensions and g/m2 for LV
mass index. Models are adjusted for the primary covariates unless stated, with
patient-clustered standard errors.

Specifications in S4d:
  Prespecified covariates   age, sex, BMI, systolic blood pressure, eGFR
  Visit 1, 3, 5             one echocardiographic visit at a time
  Non-hypertrophic only     echocardiograms without LVH
  Normal LV geometry        no LVH and relative wall thickness 0.42 or less
  Volume proxies adjusted   primary covariates plus weight, albumin and calcium
  Patient means             one row per patient, all variables averaged
  Between, within patients  each patient's mean log10 PTH, and each visit's
                            difference from that mean, in the same model
  Visit-adjusted            primary covariates plus study visit
  Adjusted for FGF23, phosphate and calcium
  LV mass indexed to height^2.7
                            LV mass (LV mass index x body surface area)
                            divided by height^2.7, in g/m^2.7

Also written:
  09_within_visit_adjusted.csv  within-patient estimate adjusted for study
                                visit (Supplementary Methods S4)
  09_kidney_function.csv        eGFR in the primary model with and without
                                PTH (per 10 ml/min/1.73 m2), the PTH x eGFR
                                interaction, and the PTH association at eGFR
                                60, 45, 30 and 15 ml/min/1.73 m2

The script stops if any estimate differs from the reported value.

Writes outputs/09_table_s4a.csv, 09_table_s4b.csv, 09_table_s4d.csv,
09_within_visit_adjusted.csv and 09_kidney_function.csv (rounded), and the
same tables unrounded to outputs/unrounded/.
"""

from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"
TOLERANCE = 0.006

DIMENSIONS = [("IVSd", "Septum, diastole"), ("IVSs", "Septum, systole"),
              ("PWTd", "Posterior wall, diastole"), ("LVIDd", "LV diameter, diastole"),
              ("LVIDs", "LV diameter, systole")]

REPORTED_S4A = {"IVSd": (0.73, 0.19, 1.28),
                 "IVSs": (1.09, 0.42, 1.76),
                 "PWTd": (0.59, 0.07, 1.11),
                 "LVIDd": (0.82, -0.69, 2.32),
                 "LVIDs": (0.18, -1.14, 1.49)}
REPORTED_S4B = {"IVSd": (0.4, -0.03, 0.82),
                 "IVSs": (0.78, 0.08, 1.47),
                 "PWTd": (0.25, -0.11, 0.61),
                 "LVIDd": (1.02, -0.53, 2.57),
                 "LVIDs": (0.4, -0.96, 1.76)}
REPORTED_S4D = {
    "Primary model": (12.44, 4.81, 20.07),
    "Prespecified covariates": (12.95, 5.51, 20.4),
    "Visit 1": (12.75, 2.56, 22.94),
    "Visit 3": (16.0, 5.31, 26.7),
    "Visit 5": (13.16, -1.22, 27.54),
    "Non-hypertrophic only": (9.73, 4.15, 15.31),
    "Normal LV geometry": (8.08, 2.17, 13.99),
    "Volume proxies adjusted": (12.13, 4.02, 20.25),
    "Patient means": (13.69, 4.79, 22.58),
    "Between patients": (14.1, 4.79, 23.41),
    "Within patients": (7.58, 0.5, 14.67),
    "Visit-adjusted": (13.15, 5.05, 21.24),
    "Adjusted for FGF23, phosphate and calcium": (10.2, 2.11, 18.29),
    "LV mass indexed to height^2.7": (5.96, 2.58, 9.34)}   # g/m^2.7
REPORTED_WITHIN_VISIT = (9.4, 1.27, 17.52)
REPORTED_KIDNEY = {
    "eGFR per 10 ml/min/1.73 m2 higher, without PTH": (-1.55, -2.84, -0.26),
    "eGFR per 10 ml/min/1.73 m2 higher, with PTH": (-0.63, -1.88, 0.62),
    "Interaction log10 PTH x eGFR (per 10 ml/min/1.73 m2)": (-5.18, -8.61, -1.75),
    "PTH at eGFR 60 ml/min/1.73 m2": (1.19, -10.2, 12.59),
    "PTH at eGFR 45 ml/min/1.73 m2": (8.97, 0.67, 17.26),
    "PTH at eGFR 30 ml/min/1.73 m2": (16.74, 8.95, 24.53),
    "PTH at eGFR 15 ml/min/1.73 m2": (24.51, 14.23, 34.79)}
EGFR_POINTS = (60, 45, 30, 15)


def fit(data, outcome, exposure, covariates):
    """Linear regression with patient-clustered standard errors; returns n,
    patients, estimate, 95% CI and P of the exposure."""
    data = data.dropna(subset=[outcome, exposure] + [c for c in covariates if not c.startswith("C(")])
    formula = f"{outcome} ~ {exposure} + " + " + ".join(covariates)
    model = smf.ols(formula, data=data).fit(
        cov_type="cluster", cov_kwds={"groups": data["Patient_ID"]})
    lo, hi = model.conf_int().loc[exposure]
    return {"n": int(model.nobs), "patients": data["Patient_ID"].nunique(),
            "estimate": model.params[exposure], "ci_low": lo, "ci_high": hi,
            "p": model.pvalues[exposure]}


def check(row, reported, label):
    """Stop if an estimate or CI limit differs from the reported value."""
    got = (row["estimate"], row["ci_low"], row["ci_high"])
    assert all(abs(g - p) <= TOLERANCE for g, p in zip(got, reported)), \
        f"{label}: {tuple(round(g, 2) for g in got)} differs from reported {reported}"


def table(rows):
    """Round estimates to 2 decimals and P to 3 significant figures."""
    out = pd.DataFrame(rows)
    out[["estimate", "ci_low", "ci_high"]] = out[["estimate", "ci_low", "ci_high"]].round(2)
    out["p"] = out["p"].map(lambda v: float(f"{v:.3g}"))
    return out


def kidney_function(d, covs):
    """eGFR with and without PTH in the model, and the PTH x eGFR interaction.
    The interaction model is refitted with eGFR centred at 60, 45, 30 and 15,
    so that the PTH term gives the PTH association at that eGFR."""
    others = [c for c in covs if c != "eGFR"]
    d = d.assign(eGFR10=d["eGFR"] / 10)
    rows = [{"Analysis": "eGFR per 10 ml/min/1.73 m2 higher, without PTH",
             **fit(d, "LVMI", "eGFR10", others)},
            {"Analysis": "eGFR per 10 ml/min/1.73 m2 higher, with PTH",
             **fit(d, "LVMI", "eGFR10", others + ["logPTH"])}]
    for point in EGFR_POINTS:
        e = d.assign(eGFR_c=(d["eGFR"] - point) / 10)
        formula = "LVMI ~ logPTH * eGFR_c + " + " + ".join(others)
        model = smf.ols(formula, data=e).fit(cov_type="cluster", cov_kwds={"groups": e["Patient_ID"]})
        if point == EGFR_POINTS[0]:
            lo, hi = model.conf_int().loc["logPTH:eGFR_c"]
            rows.append({"Analysis": "Interaction log10 PTH x eGFR (per 10 ml/min/1.73 m2)",
                         "n": int(model.nobs), "patients": e["Patient_ID"].nunique(),
                         "estimate": model.params["logPTH:eGFR_c"], "ci_low": lo, "ci_high": hi,
                         "p": model.pvalues["logPTH:eGFR_c"]})
        lo, hi = model.conf_int().loc["logPTH"]
        rows.append({"Analysis": f"PTH at eGFR {point} ml/min/1.73 m2", "n": int(model.nobs),
                     "patients": e["Patient_ID"].nunique(), "estimate": model.params["logPTH"],
                     "ci_low": lo, "ci_high": hi, "p": model.pvalues["logPTH"]})
    return rows


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    covs = prep.PRIMARY_COVS
    normal = qc[qc["geometry"] == "Normal geometry"]

    # S4a and S4b: PTH against each dimension.
    a, b = [], []
    for column, label in DIMENSIONS:
        for rows, data, reported in ((a, qc, REPORTED_S4A), (b, normal, REPORTED_S4B)):
            row = fit(data, column, "logPTH", covs)
            check(row, reported[column], label)
            rows.append({"Outcome": label, **row})

    # S4d: the primary model and its alternative specifications.
    d = qc.dropna(subset=["LVMI", "logPTH"] + covs).copy()
    d["mean_lp"] = d.groupby("Patient_ID")["logPTH"].transform("mean")
    d["dev_lp"] = d["logPTH"] - d["mean_lp"]
    means = d.groupby("Patient_ID")[["LVMI", "logPTH"] + covs].mean().reset_index()

    specifications = [
        ("Primary model", d, "logPTH", covs),
        ("Prespecified covariates", d, "logPTH", prep.PRESPECIFIED_COVS),
        ("Visit 1", d[d["visit"] == 1], "logPTH", covs),
        ("Visit 3", d[d["visit"] == 3], "logPTH", covs),
        ("Visit 5", d[d["visit"] == 5], "logPTH", covs),
        ("Non-hypertrophic only", d[d["LVH"].eq(0)], "logPTH", covs),
        ("Normal LV geometry", d[d["geometry"] == "Normal geometry"], "logPTH", covs),
        ("Volume proxies adjusted", d, "logPTH", covs + prep.VOLUME_PROXIES),
        ("Patient means", means, "logPTH", covs),
        ("Between patients", d, "mean_lp", covs + ["dev_lp"]),
        ("Within patients", d, "dev_lp", covs + ["mean_lp"]),
        ("Visit-adjusted", d, "logPTH", covs + ["C(visit)"]),
        ("Adjusted for FGF23, phosphate and calcium", d, "logPTH",
         covs + ["logFGF23", "Phosphate", "Calcium"]),
        ("LV mass indexed to height^2.7", d, "logPTH", covs),
    ]
    d["LVMI_ht27"] = d["LVMI"] * d["BSA"] / (d["Height"] / 100) ** 2.7
    s4d = []
    for label, data, exposure, covariates in specifications:
        outcome = "LVMI_ht27" if label == "LV mass indexed to height^2.7" else "LVMI"
        row = fit(data, outcome, exposure, covariates)
        check(row, REPORTED_S4D[label], label)
        s4d.append({"Specification": label, **row})

    # Within-patient estimate adjusted for study visit, and kidney function.
    within_visit = fit(d, "LVMI", "dev_lp", covs + ["mean_lp", "C(visit)"])
    check(within_visit, REPORTED_WITHIN_VISIT, "Within patients, visit-adjusted")

    kidney = kidney_function(d, covs)
    assert len(REPORTED_KIDNEY) == len(kidney), "kidney-function results are not pinned"
    for row in kidney:
        check(row, REPORTED_KIDNEY[row["Analysis"]], row["Analysis"])

    # Write the tables.
    OUT.mkdir(exist_ok=True)
    table([{"Specification": "Within patients, visit-adjusted", **within_visit}]).to_csv(
        OUT / "09_within_visit_adjusted.csv", index=False)
    table(kidney).to_csv(OUT / "09_kidney_function.csv", index=False)
    (OUT / "unrounded").mkdir(exist_ok=True)
    pd.DataFrame([{"Specification": "Within patients, visit-adjusted", **within_visit}]).to_csv(
        OUT / "unrounded" / "09_within_visit_adjusted.csv", index=False)
    pd.DataFrame(kidney).to_csv(OUT / "unrounded" / "09_kidney_function.csv", index=False)
    for name, rows in (("a", a), ("b", b), ("d", s4d)):
        out = table(rows)
        out.to_csv(OUT / f"09_table_s4{name}.csv", index=False)
        pd.DataFrame(rows).to_csv(OUT / "unrounded" / f"09_table_s4{name}.csv", index=False)
        print(f"--- Table S4{name} ---\n{out.to_string(index=False)}\n")
    print("Tables S4a, S4b and S4d reproduced.")


if __name__ == "__main__":
    main()
