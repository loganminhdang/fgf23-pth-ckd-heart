"""Supplementary Table S1c: characteristics at each study visit.

Clinical and laboratory rows: median (n) of all participants attending each
of visits 1 to 5. Echocardiographic rows: mean (n) of the echocardiograms
passing quality control at visits 1, 3 and 5 (no echocardiography at visits
2 and 4, shown as a dash).

Before writing, the script checks that no echocardiographic n exceeds the
number attending that visit and that the dashes are at visits 2 and 4 only.

Writes outputs/03_table_s1c.csv and outputs/03_table_s1c.md.
"""

import re
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import pandas as pd

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"
DASH = "\u2014"

CLINICAL = [("eGFR (ml/min/1.73 m\u00b2)", "eGFR", 1),
            ("Creatinine (\u00b5mol/L)", "Creatinine", 1),
            ("PTH (pg/mL)", "PTH", 1),
            ("FGF23 (pg/mL)", "FGF23", 1),
            ("Phosphate (mmol/L)", "Phosphate", 1),
            ("Calcium (mmol/L)", "Calcium", 1),
            ("Albumin (g/L)", "Albumin", 1),
            ("Systolic BP (mm Hg)", "SBP", 0),
            ("Diastolic BP (mm Hg)", "DBP", 0)]

LV = [("LV mass index (g/m\u00b2)", "LVMI", 1),
      ("Relative wall thickness", "RWT", 2),
      ("Septal thickness, diastole (mm)", "IVSd", 1),
      ("LV diameter, diastole (mm)", "LVIDd", 1)]

FOOTNOTE = ("Values are median (n) except LV measures, mean (n). "
             "Echocardiography was scheduled at visits 1, 3 and 5 only "
             f"({DASH}). LV measures are from the {{n}} echocardiograms "
             "passing prespecified quality control. Assay batches "
             "contributing to each visit are shown in Supplementary Figure S1.")


def half_up(value, decimals):
    """Round half away from zero, as in all other tables (Python's format rounds half to even)."""
    return str(Decimal(repr(float(value))).quantize(Decimal(1).scaleb(-decimals), rounding=ROUND_HALF_UP))


def cell(series, decimals, statistic):
    """One table cell, "value (n)"."""
    s = series.dropna()
    value = s.median() if statistic == "median" else s.mean()
    return f"{half_up(value, decimals)} ({len(s)})"


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    attendance = df.groupby("visit")["Patient_ID"].nunique()

    # Attendance, then the clinical rows (medians) and the LV rows (means).
    table = {"Attended, n": {v: str(attendance[v]) for v in range(1, 6)}}
    for label, var, dec in CLINICAL:
        table[label] = {v: cell(df.loc[df["visit"] == v, var], dec, "median")
                        for v in range(1, 6)}
    for label, var, dec in LV:
        table[label] = {v: (cell(qc.loc[qc["visit"] == v, var], dec, "mean")
                            if v in (1, 3, 5) else DASH)
                        for v in range(1, 6)}

    # Check 1: LV n does not exceed attendance.
    for label, _, _ in LV:
        for v in (1, 3, 5):
            n = int(re.search(r"\((\d+)\)", table[label][v]).group(1))
            assert n <= attendance[v], \
                f"{label} visit {v}: n={n} exceeds attendance {attendance[v]}"
    # Check 2: dashes at visits 2 and 4 only.
    for label, _, _ in LV:
        for v in range(1, 6):
            is_dash = table[label][v] == DASH
            assert is_dash == (v in (2, 4)), \
                f"{label} visit {v}: placeholder misplaced"

    out = pd.DataFrame(table).T
    out.columns = [f"Visit {v}" for v in range(1, 6)]
    out.index.name = "Variable"

    OUT.mkdir(exist_ok=True)
    out.to_csv(OUT / "03_table_s1c.csv")
    md = ["**(c) Characteristics at each study visit**", "",
          out.reset_index().to_markdown(index=False), "",
          f"(c) {FOOTNOTE.format(n=len(qc))}"]
    (OUT / "03_table_s1c.md").write_text("\n".join(md), encoding="utf-8")
    print(out.to_string())
    print("\nChecks passed; table written.")


if __name__ == "__main__":
    main()
