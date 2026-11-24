"""Primary models: LV mass index against PTH and against FGF23.

Fits the primary model (LV mass index on log10 PTH or log10 FGF23, adjusted
for age, sex, BMI, blood pressure, heart rate and eGFR, with standard errors
allowing for repeated echocardiograms of the same patient), the PTH model at
each echocardiographic visit, and the PTH association split into its
between-patient and within-patient parts.

The script first checks the analysis set (200 participants, 865 visits, 524
echocardiograms, 523 passing quality control) and stops if any count or any
estimate below differs from the value reported in the manuscript:
  PTH primary model     12.44 (4.81 to 20.07) g/m2 per 10-fold higher PTH
  FGF23 primary model   -0.12 (-2.07 to 1.83) g/m2 per 10-fold higher FGF23
  PTH within / between patients   7.58 / 14.10

Writes outputs/01_primary_models.csv.
"""

from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"

COVARIATES = "Age + male + BMI + SBP + DBP + Heart_rate + eGFR"


def fit(formula, data, term):
    """Linear regression with patient-clustered standard errors; returns the
    estimate, 95% CI and P value of one term."""
    model = smf.ols(formula, data=data).fit(
        cov_type="cluster", cov_kwds={"groups": data["Patient_ID"]})
    lo, hi = model.conf_int().loc[term]
    return {"term": term, "n": int(model.nobs),
            "patients": data["Patient_ID"].nunique(),
            "estimate": model.params[term], "ci_low": lo, "ci_high": hi,
            "p": model.pvalues[term]}


def main():
    # Analysis set: participants, visits and echocardiograms.
    counts = prep.enrolment_counts()
    assert counts["analysed"] == counts["enrolled"] == 200, counts
    df = prep.load_visits()
    assert len(df) == 865 and df["Patient_ID"].nunique() == 200, (len(df), df["Patient_ID"].nunique())
    echo, qc = prep.echo_qc(df)
    assert len(echo) == 524, f"expected 524 echo rows, got {len(echo)}"
    assert echo["visit"].value_counts().to_dict() == {1: 200, 3: 169, 5: 155}
    assert len(qc) == 523, f"expected 523 QC-passing rows, got {len(qc)}"
    assert qc.groupby("Patient_ID").size().value_counts().to_dict() == {3: 154, 2: 15, 1: 31}
    assert echo["IVSs"].isna().sum() == 1, "expected 1 echocardiogram without IVSs"
    assert qc["IVSs"].isna().sum() == 0, "expected no QC-passing echocardiogram without IVSs"
    d = prep.primary_frame(qc)
    assert d["Patient_ID"].nunique() == 200

    rows = []

    # PTH: primary model, then each echocardiographic visit separately.
    r = fit(f"LVMI ~ logPTH + {COVARIATES}", d, "logPTH")
    r["model"] = "PTH primary"
    assert abs(r["estimate"] - 12.44) < 0.01, r
    assert abs(r["ci_low"] - 4.81) < 0.01 and abs(r["ci_high"] - 20.07) < 0.01, r
    rows.append(r)

    for v in (1, 3, 5):
        rv = fit(f"LVMI ~ logPTH + {COVARIATES}", d[d["visit"] == v], "logPTH")
        rv["model"] = f"PTH visit {v}"
        rows.append(rv)

    # PTH: within-patient change and between-patient difference, fitted together.
    r = fit(f"LVMI ~ mean_lp + dev_lp + {COVARIATES}", d, "dev_lp")
    r["model"] = "PTH within-participant"
    assert abs(r["estimate"] - 7.58) < 0.01, r
    rows.append(r)
    r = fit(f"LVMI ~ mean_lp + dev_lp + {COVARIATES}", d, "mean_lp")
    r["model"] = "PTH between-participant"
    assert abs(r["estimate"] - 14.10) < 0.01, r
    rows.append(r)

    # FGF23: primary model.
    dg = qc.dropna(subset=["logFGF23", "Age", "BMI", "SBP", "DBP",
                           "Heart_rate", "eGFR"])
    r = fit(f"LVMI ~ logFGF23 + {COVARIATES}", dg, "logFGF23")
    r["model"] = "FGF23 primary"
    assert abs(r["estimate"] - (-0.12)) < 0.006, r
    assert abs(r["ci_low"] - (-2.07)) < 0.006 and abs(r["ci_high"] - 1.83) < 0.006, r
    rows.append(r)

    out = pd.DataFrame(rows)[["model", "n", "patients", "estimate",
                              "ci_low", "ci_high", "p"]].round(4)
    OUT.mkdir(exist_ok=True)
    out.to_csv(OUT / "01_primary_models.csv", index=False)
    print(out.to_string(index=False))
    print("\nAll reported values reproduced.")


if __name__ == "__main__":
    main()
