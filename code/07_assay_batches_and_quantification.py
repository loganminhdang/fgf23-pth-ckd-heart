"""FGF23 assay batches and limits of quantification (Supplementary Table S2,
Figure 2d, Supplementary Figure S1).

Each FGF23 measurement carries its assay batch and, for values the laboratory
reported below the limit of quantification, that limit (78 pg/mL at standard
sample volume, 39 pg/mL at doubled volume, 20 pg/mL in the earliest batches).
For those values the dataset holds the concentration extrapolated from the
standard curve and corrected for dilution.

  Table S2    FGF23 and LV mass index in the primary model, at each visit,
              using quantifiable values only (values reported below the limit
              excluded), within the calibrated range (78.1 to 5,000 pg/mL) and
              in normal LV geometry.
  Figure 2d   primary model; adjusted for assay batch; calibrated range;
              quantifiable values only.
  Figure S1a  number of measurements by assay batch and study visit. Batches
              are numbered B01 to B15 in the order assayed and labelled with
              the year of their first plate.
  Figure S1b  FGF23 by assay batch (box plots, all 865 measurements).

The script stops if a Table S2 estimate or n differs from the reported value.

Writes outputs/07_table_s2.csv, 07_figure2d.csv, 07_batch_by_visit.csv,
07_batch_distribution.csv and 07_linkage_summary.csv (counts quoted in
Supplementary Methods S3), and Table S2 and Figure 2d unrounded to
outputs/unrounded/.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from matplotlib import cbook

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"
TOLERANCE = 0.006
LOW, HIGH = prep.CALIBRATED_RANGE

REPORTED_S2 = {
    "Primary model": (-0.12, -2.07, 1.83, 523),
    "Visit 1": (1.83, -1.23, 4.9, 199),
    "Visit 3": (-1.61, -3.63, 0.41, 169),
    "Visit 5": (-0.55, -3.28, 2.18, 155),
    "Quantifiable values only": (-1.11, -3.86, 1.64, 398),
    "Within calibrated range": (-0.08, -4.86, 4.7, 287),
    "Normal LV geometry": (-0.16, -1.9, 1.58, 450)}


def fit(frame, covariates):
    """LV mass index on log10 FGF23 with patient-clustered standard errors."""
    data = frame.dropna(subset=["LVMI", "logFGF23"] + prep.PRIMARY_COVS)
    model = smf.ols("LVMI ~ logFGF23 + " + " + ".join(covariates), data=data).fit(
        cov_type="cluster", cov_kwds={"groups": data["Patient_ID"]})
    lo, hi = model.conf_int().loc["logFGF23"]
    return {"n": int(model.nobs), "patients": data["Patient_ID"].nunique(),
            "estimate": model.params["logFGF23"], "ci_low": lo, "ci_high": hi,
            "p": model.pvalues["logFGF23"]}


def table(rows):
    """Round estimates to 2 decimals and P to 3 significant figures."""
    out = pd.DataFrame(rows)
    out[["estimate", "ci_low", "ci_high"]] = out[["estimate", "ci_low", "ci_high"]].round(2)
    out["p"] = out["p"].map(lambda v: float(f"{v:.3g}"))
    return out


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    # Assay batch of each measurement, matched to the echocardiograms.
    link = prep.load_linkage()
    assert len(link) == 865 and link["run"].nunique() == 15, \
        f"linkage has {len(link)} measurements in {link['run'].nunique()} batches"
    assert link["plate_year"].notna().all(), "measurements without a plate year"

    keys = ["Patient_ID", "visit"]
    analysed = link[link["Patient_ID"].isin(df["Patient_ID"])]
    assert len(analysed) == len(link) == len(df), "linkage and visits do not match one to one"
    qc = qc.merge(link[keys + ["run", "limit", "bound"]], on=keys, how="left")
    assert qc["run"].notna().all(), "echocardiographic observations without an assay batch"
    qc["below_limit"] = qc["limit"] == "below"
    covs = prep.PRIMARY_COVS

    # Table S2: the primary FGF23 model under each specification.
    specifications = {
        "Primary model": qc,
        "Visit 1": qc[qc["visit"] == 1],
        "Visit 3": qc[qc["visit"] == 3],
        "Visit 5": qc[qc["visit"] == 5],
        "Quantifiable values only": qc[~qc["below_limit"]],
        "Within calibrated range": qc[qc["FGF23"].between(LOW, HIGH)],
        "Normal LV geometry": qc[qc["geometry"] == "Normal geometry"],
    }
    s2 = []
    for label, frame in specifications.items():
        row = fit(frame, covs)
        est, lo, hi, n = REPORTED_S2[label]
        assert row["n"] == n and all(abs(g - p) <= TOLERANCE for g, p in
                                     zip((row["estimate"], row["ci_low"], row["ci_high"]),
                                         (est, lo, hi))), f"Table S2, {label}: {row}"
        s2.append({"Specification": label, **row})

    # Figure 2d: sensitivity to assay batch and quantification.
    batch = fit(qc, covs + ["C(run)"])
    figure2d = [{"Specification": "Primary model", **fit(qc, covs)},
                {"Specification": "Adjusted for assay batch", **batch},
                {"Specification": "Within calibrated range",
                 **fit(specifications["Within calibrated range"], covs)},
                {"Specification": "Quantifiable values only",
                 **fit(specifications["Quantifiable values only"], covs)}]

    # Figure S1a: measurements by batch and visit.
    first_year = link.groupby("run")["plate_year"].min()
    labels = {run: f"{run} ({int(year)})" for run, year in first_year.items()}
    by_visit = (pd.crosstab(analysed["run"], analysed["visit"])
                  .reindex(columns=range(1, 6), fill_value=0))
    by_visit.insert(0, "batch", by_visit.index.map(labels))
    spanning = int((by_visit[list(range(1, 6))] > 0).sum(axis=1).gt(1).sum())

    # Figure S1b: FGF23 by batch, all 865 measurements.
    values = analysed[keys + ["run"]].merge(df[keys + ["FGF23"]], on=keys, how="left")
    assert len(values) == len(analysed) and values["FGF23"].notna().all()
    rows = []
    for run, group in values.groupby("run"):
        logs = np.log10(group["FGF23"].dropna().to_numpy())
        box = cbook.boxplot_stats(logs, whis=1.5)[0]
        rows.append({"run": run, "batch": labels[run], "n": len(logs),
                     "median_pg_ml": round(float(10 ** box["med"]), 1),
                     **{f"log10_{k}": round(float(box[k]), 5)
                        for k in ("whislo", "q1", "med", "q3", "whishi")}})
    distribution = pd.DataFrame(rows)

    # Counts quoted in Supplementary Methods S3.
    limits = qc.loc[qc["below_limit"], "bound"].value_counts().sort_index()
    summary = pd.DataFrame([{
        "measurements": len(link), "patients": link["Patient_ID"].nunique(),
        "batches": link["run"].nunique(), "batches_spanning_visits": spanning,
        "first_plate_year": int(link["plate_year"].min()),
        "last_plate_year": int(link["plate_year"].max()),
        "below_limit": int((link["limit"] == "below").sum()),
        "below_limit_analysed": int(qc["below_limit"].sum()),
        "limits_applied_analysed": "; ".join(f"{k:g} pg/mL: {v}" for k, v in limits.items())}])

    OUT.mkdir(exist_ok=True)
    table(s2).to_csv(OUT / "07_table_s2.csv", index=False)
    table(figure2d).to_csv(OUT / "07_figure2d.csv", index=False)
    (OUT / "unrounded").mkdir(exist_ok=True)
    pd.DataFrame(s2).to_csv(OUT / "unrounded" / "07_table_s2.csv", index=False)
    pd.DataFrame(figure2d).to_csv(OUT / "unrounded" / "07_figure2d.csv", index=False)
    by_visit.to_csv(OUT / "07_batch_by_visit.csv")
    distribution.to_csv(OUT / "07_batch_distribution.csv", index=False)
    summary.to_csv(OUT / "07_linkage_summary.csv", index=False)
    print(table(s2).to_string(index=False), "\n")
    print(table(figure2d).to_string(index=False), "\n")
    print(summary.T.to_string(header=False), "\n")
    print("Table S2 reproduced.")


if __name__ == "__main__":
    main()
