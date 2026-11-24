"""Baseline characteristics (Table 1), LVH at baseline (Supplementary Table S1b)
and LV mass index over follow-up.

Table 1 describes the 200 participants at visit 1: clinical and laboratory
rows from all baseline visits, echocardiographic rows from the baseline
echocardiograms passing quality control. Values are mean +/- SD, or median
(interquartile range) for PTH, FGF23 and time since CKD diagnosis; n is the
number of participants with the value.

Primary kidney disease and the year of CKD diagnosis are read from
data/study_workbook.xlsx. Time since diagnosis = 2023 (first year of
enrolment) minus the year of diagnosis. Without the workbook these rows are
left out, with a message.

Table S1b lists the participants whose baseline LV mass index is above the
sex-specific LVH threshold, and by how much (%). Its output files list
individual participants and are not included in the public repository.

LV mass index over follow-up: mean at visits 1, 3 and 5 in participants with
an analysable echocardiogram at all three.

The script stops if the participant count, time since CKD diagnosis, the
number of participants in Table S1b or the mean LV mass index at visits 1, 3
and 5 differ from the reported values.

Writes outputs/08_table1.csv, 08_table_s1b.csv and 08_lvmi_by_visit.csv, and
Table 1 and Table S1b unrounded to outputs/unrounded/.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"
WORKBOOK = prep.DATA / "study_workbook.xlsx"
REFERENCE_YEAR = 2023
REPORTED_DURATION = (197, 3.0, 1.0, 7.0)   # participants, median, quartiles (years)

# Primary kidney disease: workbook labels (Vietnamese) and their English names.
AETIOLOGY = {"Viêm cầu thận mạn": "Chronic glomerulonephritis",
             "Bệnh kẽ thận": "Interstitial kidney disease",
             "Khác (ghi rõ)": "Other",
             "Viêm thận bể thận": "Pyelonephritis",
             "Bệnh thận đái tháo đường": "Diabetic nephropathy"}
AETIOLOGY_ORDER = ["Chronic glomerulonephritis", "Interstitial kidney disease",
                   "Pyelonephritis", "Diabetic nephropathy", "Other"]
# KDIGO eGFR categories (ml/min/1.73 m2).
STAGES = [(90, np.inf, "G1 (\u226590)"), (60, 90, "G2 (60\u201389)"), (45, 60, "G3a (45\u201359)"), (30, 45, "G3b (30\u201344)"),
          (15, 30, "G4 (15\u201329)"), (-np.inf, 15, "G5 (<15)")]


def workbook_rows(patients):
    """Primary kidney disease and time since CKD diagnosis for the analysed
    participants. The workbook has column labels in row 2 and one row per
    participant from row 3; columns are found by their Vietnamese labels."""
    if not WORKBOOK.exists():
        print(f"{WORKBOOK.name} not found: CKD aetiology and time since diagnosis omitted",
              file=sys.stderr)
        return None, None, 0
    from openpyxl import load_workbook
    book = load_workbook(WORKBOOK, data_only=True)
    sheet = next((s for s in book.sheetnames if "Lần 5" in s or "Lan 5" in s), book.sheetnames[0])
    ws = book[sheet]
    header = [(ws.cell(row=2, column=c).value or f"c{c}") for c in range(1, ws.max_column + 1)]
    w = pd.DataFrame(list(ws.iter_rows(min_row=3, values_only=True)), columns=header).dropna(how="all")
    pid = next(c for c in w.columns if "PID" in str(c))
    w = w[w[pid].notna()].copy()
    w["pid"] = w[pid].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    matched = w[w["pid"].isin(patients)]
    aetiology_col = next((c for c in w.columns if "Chẩn đoán bệnh thận" in str(c)), None)
    year_col = next((c for c in w.columns if "Thời gian mắc bệnh" in str(c)), None)
    aetiology = (matched[aetiology_col].map(lambda x: AETIOLOGY.get(str(x).strip(), str(x).strip()))
                 if aetiology_col else None)
    year = (matched[year_col].map(lambda v: v.year if hasattr(v, "year") else v)
            if year_col else None)
    duration = ((REFERENCE_YEAR - pd.to_numeric(year, errors="coerce")).dropna()
                if year_col else None)
    return aetiology, duration, len(matched)


def full(x):
    """A number written in full, for the unrounded tables."""
    return repr(float(x))


def table1(df, qc):
    """Table 1, rounded and unrounded (the same rows and columns)."""
    # Baseline visits (b) and baseline echocardiograms passing quality control (e).
    b = df[df["visit"] == 1]
    e = qc[qc["visit"] == 1].copy()
    e["LVMI_ht27"] = e["LVMI"] * e["BSA"] / (e["Height"] / 100) ** 2.7
    n = len(b)
    assert n == 200, f"expected 200 baseline visits, got {n}"
    rows = [("Participants", f"{n}", "")]
    raw = [("Participants", f"{n}", "")]

    def section(title):
        """A heading row."""
        rows.append((title, "", ""))
        raw.append((title, "", ""))

    def num(label, s, dp=1, median=False):
        """A row "mean +/- SD", or "median (Q1-Q3)" when median=True."""
        v = s.dropna()
        value = (f"{v.median():.{dp}f} ({v.quantile(.25):.{dp}f}\u2013{v.quantile(.75):.{dp}f})"
                 if median else f"{v.mean():.{dp}f} \u00b1 {v.std():.{dp}f}")
        exact = (f"{full(v.median())} ({full(v.quantile(.25))}\u2013{full(v.quantile(.75))})"
                 if median else f"{full(v.mean())} \u00b1 {full(v.std())}")
        rows.append((label, value, str(len(v))))
        raw.append((label, exact, str(len(v))))

    def pct(label, k, denominator):
        """A row "count (%)"."""
        rows.append((label, f"{k} ({100 * k / denominator:.1f})", str(denominator)))
        raw.append((label, f"{k} ({full(100 * k / denominator)})", str(denominator)))

    section("Demographic and anthropometric")
    num("Age, years", b["Age"])
    pct("Male sex, n (%)", int(b["male"].sum()), n)
    num("Body-mass index, kg m\u207b\u00b2", b["BMI"])
    num("Body surface area, m\u00b2", b["BSA"], 2)

    section("Haemodynamic")
    num("Systolic blood pressure, mm Hg", b["SBP"], 0)
    num("Diastolic blood pressure, mm Hg", b["DBP"], 0)
    num("Heart rate, beats min\u207b\u00b9", b["Heart_rate"], 0)

    section("Kidney disease")
    num("eGFR, ml min\u207b\u00b9 1.73 m\u207b\u00b2", b["eGFR"])
    num("Serum creatinine, \u00b5mol l\u207b\u00b9", b["Creatinine"], 0)
    for lo, hi, label in STAGES:
        pct(f"   {label}, n (%)", int(((b["eGFR"] >= lo) & (b["eGFR"] < hi)).sum()), n)
    aetiology, duration, matched = workbook_rows(set(b["Patient_ID"].astype(str)))
    if aetiology is not None:
        counts = aetiology.value_counts()
        for label in AETIOLOGY_ORDER + [k for k in counts.index if k not in AETIOLOGY_ORDER]:
            if label in counts:
                pct(f"   {label}, n (%)", int(counts[label]), matched)
    if duration is not None and len(duration):
        got = (len(duration), duration.median(), duration.quantile(.25), duration.quantile(.75))
        assert got == REPORTED_DURATION, f"time since CKD diagnosis: {got}"
        rows.append(("Time since CKD diagnosis, years",
                     f"{duration.median():.0f} ({duration.quantile(.25):.0f}\u2013"
                     f"{duration.quantile(.75):.0f})", str(len(duration))))
        raw.append(("Time since CKD diagnosis, years",
                    f"{full(duration.median())} ({full(duration.quantile(.25))}\u2013"
                    f"{full(duration.quantile(.75))})", str(len(duration))))

    section("Mineral metabolism")
    num("Parathyroid hormone, pg ml\u207b\u00b9", b["PTH"], 0, median=True)
    num("FGF23, pg ml\u207b\u00b9", b["FGF23"], 1, median=True)
    num("Calcium, mmol l\u207b\u00b9", b["Calcium"], 2)
    num("Phosphate, mmol l\u207b\u00b9", b["Phosphate"], 2)
    num("Albumin, g l\u207b\u00b9", b["Albumin"], 1)

    section("Cardiac structure")
    num("LV mass index, g m\u207b\u00b2", e["LVMI"])
    num("LV mass index, g m\u207b\u00b2\u02d9\u2077", e["LVMI_ht27"])
    num("Relative wall thickness", e["RWT"], 3)
    num("Septal thickness in diastole, mm", e["IVSd"])
    num("Posterior wall thickness in diastole, mm", e["PWTd"])
    num("LV internal dimension in diastole, mm", e["LVIDd"])
    pct("Left ventricular hypertrophy, n (%)", int(e["LVH"].sum()), len(e))
    pct("Increased relative wall thickness, n (%)", int((e["RWT"] > prep.RWT_INCREASED).sum()), len(e))

    section("Diastolic function")
    num("E/e\u2032, average", b["E_e_avg"], 2)
    num("e\u2032 septal, cm s\u207b\u00b9", b["e_septal"], 1)
    num("e\u2032 lateral, cm s\u207b\u00b9", b["e_lateral"], 1)
    num("Left atrial volume index, ml m\u207b\u00b2", b["LAVI"], 1)
    columns = ["Characteristic", "Value", "n"]
    return pd.DataFrame(rows, columns=columns), pd.DataFrame(raw, columns=columns)


def table_s1b(df):
    """Table S1b: baseline LV mass index above the LVH threshold, rounded and
    unrounded. Participants are numbered 1 to 11 in order of LV mass index."""
    baseline = df[df["visit"] == 1]
    lvh = baseline[baseline["LVH"].eq(1)].copy()
    lvh["pct_above"] = 100 * (lvh["LVMI"] - lvh["lvh_threshold"]) / lvh["lvh_threshold"]
    lvh = lvh.sort_values("LVMI").reset_index(drop=True)
    assert len(lvh) == 11, f"expected 11 participants meeting LVH criteria, got {len(lvh)}"

    def table(decimals):
        def value(series, dp):
            return series.round(dp) if decimals else series
        return pd.DataFrame({
            "Participant": range(1, len(lvh) + 1),
            "Sex": np.where(lvh["male"].eq(1), "M", "F"),
            "LVMI (g/m2)": value(lvh["LVMI"], 1),
            "Threshold (g/m2)": value(lvh["lvh_threshold"], 0),
            "% above threshold": value(lvh["pct_above"], 1),
            "RWT": value(lvh["RWT"], 2),
            "Geometry": lvh["geometry"],
        })
    return table(True), table(False)


REPORTED_LVMI_BY_VISIT = (154, 78.0, 78.0, 78.2)   # participants, mean at visits 1, 3, 5


def lvmi_by_visit(qc):
    """Mean LV mass index at visits 1, 3 and 5, in participants with all three."""
    wide = qc.pivot_table(index="Patient_ID", columns="visit", values="LVMI").dropna()
    out = pd.DataFrame({"visit": [1, 3, 5], "participants": len(wide),
                        "mean_lvmi": [round(wide[v].mean(), 1) for v in (1, 3, 5)],
                        "sd_lvmi": [round(wide[v].std(), 1) for v in (1, 3, 5)]})
    got = (len(wide), *out["mean_lvmi"])
    assert REPORTED_LVMI_BY_VISIT is not None and got == REPORTED_LVMI_BY_VISIT, got
    return out


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    (t1, t1_raw), (s1b, s1b_raw), by_visit = table1(df, qc), table_s1b(df), lvmi_by_visit(qc)
    OUT.mkdir(exist_ok=True)
    (OUT / "unrounded").mkdir(exist_ok=True)
    t1.to_csv(OUT / "08_table1.csv", index=False)
    t1_raw.to_csv(OUT / "unrounded" / "08_table1.csv", index=False)
    s1b.to_csv(OUT / "08_table_s1b.csv", index=False)
    s1b_raw.to_csv(OUT / "unrounded" / "08_table_s1b.csv", index=False)
    by_visit.to_csv(OUT / "08_lvmi_by_visit.csv", index=False)
    print(t1.to_string(index=False), "\n")
    print(s1b.to_string(index=False))


if __name__ == "__main__":
    main()
