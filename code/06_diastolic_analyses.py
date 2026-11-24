"""Diastolic function and mineral metabolism (Figure 4d, Supplementary Table S4c).

Baseline (visit 1) echocardiograms of participants with normal LV geometry.
Outcomes: average E/e' (primary), septal e', lateral e', left atrial volume
index and E/A (calculated from the E and A velocities).

Exposures: log10 PTH and log10 FGF23. Outcomes are in SD units, so each
estimate is the change in SD per 10-fold higher PTH or FGF23. Age is included
as a positive control (a known determinant of diastolic function), in native
units per year; Table S4c shows age for E/e', septal e' and lateral e' only.

Each model is run with the primary covariates and with the prespecified
(smaller) covariate set: age, sex, BMI, systolic blood pressure and eGFR.

Writes outputs/06_table_s4c.csv (rounded) and outputs/unrounded/06_table_s4c.csv.
"""

from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"

OUTCOMES = {"E/e' (average)": "E_e_avg", "e' septal": "e_septal",
            "e' lateral": "e_lateral", "LA volume index": "LAVI",
            "E/A (derived)": "EA"}
EXPOSURES = {"PTH": "logPTH", "FGF23": "logFGF23", "Age (positive control)": "Age"}
COVARIATE_SETS = {"primary": prep.PRIMARY_COVS,
                  "prespecified": prep.PRESPECIFIED_COVS}


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    baseline = qc[(qc["visit"] == 1) & (qc["geometry"] == "Normal geometry")].copy()

    # One model per exposure, covariate set and outcome.
    rows, raw = [], []
    for exposure_label, exposure in EXPOSURES.items():
        for set_label, covariates in COVARIATE_SETS.items():
            covs = [c for c in covariates if c != exposure]
            for outcome_label, outcome in OUTCOMES.items():
                data = baseline.dropna(subset=[outcome, exposure] + covs).copy()
                if exposure == "Age":
                    target = outcome          # positive control: native units
                else:
                    target = f"z_{outcome}"   # outcome in SD units
                    data[target] = (data[outcome] - data[outcome].mean()) / data[outcome].std()
                formula = f"{target} ~ {exposure}" + "".join(f" + {c}" for c in covs)
                model = smf.ols(formula, data=data).fit(
                    cov_type="cluster", cov_kwds={"groups": data["Patient_ID"]})
                lo, hi = model.conf_int().loc[exposure]
                decimals = 3 if exposure == "Age" else 2
                rows.append({"Exposure": exposure_label, "Outcome": outcome_label,
                             "Covariates": set_label, "n": int(model.nobs),
                             "estimate": round(model.params[exposure], decimals),
                             "ci_low": round(lo, decimals), "ci_high": round(hi, decimals),
                             "p": float(f"{model.pvalues[exposure]:.3g}")})
                raw.append({**rows[-1], "estimate": model.params[exposure],
                            "ci_low": lo, "ci_high": hi, "p": model.pvalues[exposure]})

    def reported(table):
        # Age is reported under the primary covariate set only.
        return table[~((table["Exposure"].str.startswith("Age"))
                       & (table["Covariates"] == "prespecified"))]

    table = reported(pd.DataFrame(rows))

    OUT.mkdir(exist_ok=True)
    table.to_csv(OUT / "06_table_s4c.csv", index=False)
    (OUT / "unrounded").mkdir(exist_ok=True)
    reported(pd.DataFrame(raw)).to_csv(OUT / "unrounded" / "06_table_s4c.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
