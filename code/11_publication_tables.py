"""Tables of the manuscript and the Supplementary Information, as printed.

Writes Table 1 and Supplementary Tables S1b to S5 to outputs/tables/, one CSV
file per printed table, with the printed headers, row labels, row order and
cell text. The manuscript tables are filled from these files.

Numbers are rounded half-up (0.125 -> 0.13) from the unrounded results in
outputs/unrounded/. Estimates and CI limits have 2 decimals, or 3 where printed
so (Tables S3c and S3d except the LV mass index rows, and the age rows of
Table S4c); P values have 3 decimals, and are shown as <0.001 below 0.001.
Negative numbers use the typographic minus sign.

  table_1.csv          Table 1          outputs/unrounded/08_table1.csv
  table_s1b.csv        Table S1b        outputs/unrounded/08_table_s1b.csv
                                        (lists participants; not in the public repository)
  table_s1c.csv        Table S1c        outputs/03_table_s1c.csv
  table_s2.csv         Table S2         outputs/unrounded/07_table_s2.csv, and the
                                        assay-batch row of outputs/unrounded/07_figure2d.csv
  table_s3a.csv to     Tables S3a-S3d   outputs/unrounded/04_table_s3a.csv to
  table_s3d.csv                         04_table_s3d.csv, and the baseline row of
                                        Table S3b from outputs/unrounded/04_figure3b.csv
  table_s4a.csv,       Tables S4a, S4b, outputs/unrounded/09_table_s4a.csv,
  table_s4b.csv,       S4d              09_table_s4b.csv, 09_table_s4d.csv
  table_s4d.csv
  table_s4c.csv        Table S4c        outputs/unrounded/06_table_s4c.csv
  table_s4e.csv        Table S4e        outputs/unrounded/05_spline_contrasts.csv
  table_s5.csv         Table S5         outputs/single_cell/tables/table_s5.csv and
                                        by_celltype.csv (single-cell analysis)

Table S1a (patients excluded at screening) comes from the screening log and is
not produced by this public code. Table S5 is written only after
single_cell/run_single_cell.sh has been run; run this script again after it.

Before formatting, each unrounded table is checked against the rounded table
of the same name in outputs/ (same rows, values agreeing to rounding), and
every row of each input must be either printed or listed here as not printed.
"""

import csv
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import pandas as pd

import prep

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
UNROUNDED = OUTPUTS / "unrounded"
TABLES = OUTPUTS / "tables"
SINGLE_CELL = OUTPUTS / "single_cell" / "tables"

MINUS = "−"
DASH = "—"
CI = "Estimate (95% CI)"
N_PATIENTS = "n (patients)"


# Formatting ----------------------------------------------------------------

def fixed(value, decimals):
    """Round half-up to a number of decimals, with the minus sign U+2212."""
    exact = Decimal(str(value).strip())
    rounded = exact.quantize(Decimal(1).scaleb(-decimals), rounding=ROUND_HALF_UP)
    if rounded == 0:
        rounded = abs(rounded)
    return f"{rounded:f}".replace("-", MINUS)


def p_value(value):
    """P to 3 decimals, <0.001 below 0.001."""
    return "<0.001" if Decimal(str(value).strip()) < Decimal("0.001") else fixed(value, 3)


def estimate_ci(row, decimals, estimate="estimate"):
    """Estimate and 95% CI as printed: estimate (low to high)."""
    return (f"{fixed(row[estimate], decimals)} ({fixed(row['ci_low'], decimals)} to "
            f"{fixed(row['ci_high'], decimals)})")


def n_patients(row):
    """n (patients): echocardiograms, and the patients they come from."""
    return f"{row['n']} ({row['patients']})"


# Inputs and checks ---------------------------------------------------------

def read(path):
    """Read a CSV file as text, so that values are not changed on reading."""
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")


def unrounded(name):
    """An unrounded table, checked against the rounded table of the same name.

    The two must have the same rows, and each number of the unrounded table
    must lie within rounding of the rounded one (half a unit in the last
    decimal written, or equal at 3 significant figures for P values written in
    exponent notation). This catches unrounded and rounded tables that come
    from different runs.
    """
    exact, rounded = read(UNROUNDED / name), read(OUTPUTS / name)
    if len(exact) != len(rounded):
        raise SystemExit(f"outputs/unrounded/{name} and outputs/{name} differ in rows; "
                         "run the analysis scripts again")
    for column in rounded.columns:
        for i, (u, r) in enumerate(zip(exact[column], rounded[column])):
            if u == r:
                continue
            try:
                u_value, r_value = float(u), float(r)
            except ValueError:
                raise SystemExit(f"{name}, row {i + 1}, {column}: {u!r} and {r!r} differ")
            if "e" in r.lower():
                close = float(f"{u_value:.3g}") == r_value
            else:
                decimals = len(r.split(".")[1]) if "." in r else 0
                close = abs(u_value - r_value) <= 0.5 * 10 ** -decimals + 1e-9
            if not close:
                raise SystemExit(f"{name}, row {i + 1}, {column}: unrounded {u} does not round "
                                 f"to {r}; run the analysis scripts again")
    return exact


def select(table, column, printed, omitted=(), name=""):
    """Rows of `table` in the order of `printed` (a mapping from source label
    to printed label), checking that every source row is printed or omitted."""
    labels = list(table[column])
    unexpected = [x for x in labels if x not in printed and x not in omitted]
    missing = [x for x in printed if x not in labels]
    if unexpected or missing:
        raise SystemExit(f"{name}: rows not in the printed layout {unexpected}, "
                         f"printed rows missing {missing}")
    return [(printed[label], table[table[column] == label].iloc[0]) for label in printed]


def write(name, header, rows):
    """Write one printed table to outputs/tables/."""
    TABLES.mkdir(parents=True, exist_ok=True)
    width = len(header)
    assert all(len(r) == width for r in rows), f"{name}: rows of unequal length"
    with open(TABLES / name, "w", encoding="utf-8", newline="") as handle:
        out = csv.writer(handle, lineterminator="\n")
        out.writerow(header)
        out.writerows(rows)
    print(f"outputs/tables/{name}: {len(rows)} rows x {width} columns")


def absent(name, reason):
    """Say that a table cannot be written (its input is missing) and remove any
    old copy of it."""
    stale = TABLES / name
    if stale.exists():
        stale.unlink()
        reason += f"; outputs/tables/{name} from an earlier run removed"
    print(f"outputs/tables/{name} not written: {reason}")


# Table 1 -------------------------------------------------------------------

# Printed label, label in outputs/08_table1.csv (None: a heading added in
# print), and decimals (None: heading or count).
TABLE_1 = [
    ("Participants", "Participants", None),
    ("Demographic and anthropometric", "Demographic and anthropometric", None),
    ("Age, years", "Age, years", 1),
    ("Male sex, n (%)", "Male sex, n (%)", 1),
    ("Body-mass index, kg/m²", "Body-mass index, kg m⁻²", 1),
    ("Body surface area, m²", "Body surface area, m²", 2),
    ("Hemodynamic", "Haemodynamic", None),
    ("Systolic blood pressure, mm Hg", "Systolic blood pressure, mm Hg", 0),
    ("Diastolic blood pressure, mm Hg", "Diastolic blood pressure, mm Hg", 0),
    ("Heart rate, beats/min", "Heart rate, beats min⁻¹", 0),
    ("Kidney disease", "Kidney disease", None),
    ("eGFR, ml/min/1.73 m²", "eGFR, ml min⁻¹ 1.73 m⁻²", 1),
    ("Serum creatinine, µmol/l", "Serum creatinine, µmol l⁻¹", 0),
    ("eGFR category, ml/min/1.73 m²", None, None),
    ("G1 (≥90), n (%)", "   G1 (≥90), n (%)", 1),
    ("G2 (60–89), n (%)", "   G2 (60–89), n (%)", 1),
    ("G3a (45–59), n (%)", "   G3a (45–59), n (%)", 1),
    ("G3b (30–44), n (%)", "   G3b (30–44), n (%)", 1),
    ("G4 (15–29), n (%)", "   G4 (15–29), n (%)", 1),
    ("G5 (<15), n (%)", "   G5 (<15), n (%)", 1),
    ("Primary kidney disease", None, None),
    ("Chronic glomerulonephritis, n (%)", "   Chronic glomerulonephritis, n (%)", 1),
    ("Interstitial kidney disease, n (%)", "   Interstitial kidney disease, n (%)", 1),
    ("Other, n (%)", "   Other, n (%)", 1),
    ("Pyelonephritis, n (%)", "   Pyelonephritis, n (%)", 1),
    ("Diabetic nephropathy, n (%)", "   Diabetic nephropathy, n (%)", 1),
    ("Time since CKD diagnosis, years", "Time since CKD diagnosis, years", 0),
    ("Mineral metabolism", "Mineral metabolism", None),
    ("Parathyroid hormone, pg/ml", "Parathyroid hormone, pg ml⁻¹", 0),
    ("FGF23, pg/ml", "FGF23, pg ml⁻¹", 1),
    ("Calcium, mmol/l", "Calcium, mmol l⁻¹", 2),
    ("Phosphate, mmol/l", "Phosphate, mmol l⁻¹", 2),
    ("Albumin, g/l", "Albumin, g l⁻¹", 1),
    ("Cardiac structure", "Cardiac structure", None),
    ("LV mass index, g/m²", "LV mass index, g m⁻²", 1),
    ("Relative wall thickness", "Relative wall thickness", 3),
    ("Septal thickness in diastole, mm", "Septal thickness in diastole, mm", 1),
    ("Posterior wall thickness in diastole, mm", "Posterior wall thickness in diastole, mm", 1),
    ("LV internal dimension in diastole, mm", "LV internal dimension in diastole, mm", 1),
    ("Left ventricular hypertrophy, n (%)", "Left ventricular hypertrophy, n (%)", 1),
    ("Increased relative wall thickness, n (%)", "Increased relative wall thickness, n (%)", 1),
    ("Diastolic function", "Diastolic function", None),
    ("E/e′, average", "E/e′, average", 2),
    ("e′ septal, cm/s", "e′ septal, cm s⁻¹", 1),
    ("e′ lateral, cm/s", "e′ lateral, cm s⁻¹", 1),
    ("Left atrial volume index, ml/m²", "Left atrial volume index, ml m⁻²", 1),
]
TABLE_1_NOT_PRINTED = ["LV mass index, g m⁻²˙⁷"]   # LV mass / height^2.7


def table1_value(text, decimals):
    """Round one unrounded Table 1 value: "mean ± SD", "median (Q1–Q3)" or
    "count (percentage)"."""
    if decimals is None:
        return text
    if " ± " in text:
        mean, sd = text.split(" ± ")
        return f"{fixed(mean, decimals)} ± {fixed(sd, decimals)}"
    centre, spread = text.rstrip(")").split(" (")
    if "–" in spread:
        q1, q3 = spread.split("–")
        return f"{fixed(centre, decimals)} ({fixed(q1, decimals)}–{fixed(q3, decimals)})"
    return f"{int(centre)} ({fixed(spread, decimals)})"


def table_1():
    """Table 1, rounded half-up from the unrounded values; any value that
    differs from outputs/08_table1.csv is reported on screen."""
    exact, rounded = read(UNROUNDED / "08_table1.csv"), read(OUTPUTS / "08_table1.csv")
    assert list(exact["Characteristic"]) == list(rounded["Characteristic"]), \
        "outputs/unrounded/08_table1.csv and outputs/08_table1.csv differ in rows"
    printed = {source: label for label, source, _ in TABLE_1 if source}
    select(exact, "Characteristic", printed, TABLE_1_NOT_PRINTED, "Table 1")
    exact, rounded = exact.set_index("Characteristic"), rounded.set_index("Characteristic")
    rows = []
    for label, source, decimals in TABLE_1:
        if source is None:
            rows.append([label, "", ""])
            continue
        value = table1_value(exact.loc[source, "Value"], decimals)
        if value != rounded.loc[source, "Value"]:
            print(f"  Table 1, {label}: rounded half-up {value}; "
                  f"outputs/08_table1.csv has {rounded.loc[source, 'Value']}")
        rows.append([label, value, exact.loc[source, "n"]])
    write("table_1.csv", ["Characteristic", "Value", "n"], rows)


# Table S1 ------------------------------------------------------------------

def table_s1b():
    table = unrounded("08_table_s1b.csv")
    header = ["Participant", "Sex", "LVMI (g/m²)", "Threshold (g/m²)",
              "% above threshold", "RWT", "Geometry"]
    rows = [[r["Participant"], r["Sex"], fixed(r["LVMI (g/m2)"], 1),
             fixed(r["Threshold (g/m2)"], 0), fixed(r["% above threshold"], 1),
             fixed(r["RWT"], 2), r["Geometry"]] for _, r in table.iterrows()]
    write("table_s1b.csv", header, rows)


S1C_ROWS = {   # label in outputs/03_table_s1c.csv: printed label
    "Attended, n": "Attended patients",
    "eGFR (ml/min/1.73 m²)": "eGFR (ml/min/1.73 m²)",
    "Creatinine (µmol/L)": "Creatinine (µmol/l)",
    "PTH (pg/mL)": "PTH (pg/ml)",
    "FGF23 (pg/mL)": "FGF23 (pg/ml)",
    "Phosphate (mmol/L)": "Phosphate (mmol/l)",
    "Calcium (mmol/L)": "Calcium (mmol/l)",
    "Albumin (g/L)": "Albumin (g/l)",
    "Systolic BP (mm Hg)": "Systolic BP (mm Hg)",
    "Diastolic BP (mm Hg)": "Diastolic BP (mm Hg)",
    "LV mass index (g/m²)": "LV mass index (g/m²)",
    "Relative wall thickness": "Relative wall thickness",
    "Septal thickness, diastole (mm)": "Septal thickness, diastole (mm)",
    "LV diameter, diastole (mm)": "LV diameter, diastole (mm)",
}


def table_s1c():
    table = read(OUTPUTS / "03_table_s1c.csv")
    visits = [f"Visit {v}" for v in range(1, 6)]
    assert list(table.columns) == ["Variable"] + visits, list(table.columns)
    rows = [[label] + [r[v] for v in visits]
            for label, r in select(table, "Variable", S1C_ROWS, name="Table S1c")]
    write("table_s1c.csv", ["Variable"] + visits, rows)


# Table S2 ------------------------------------------------------------------

def table_s2():
    low, high = prep.CALIBRATED_RANGE
    table = unrounded("07_table_s2.csv")
    batch = unrounded("07_figure2d.csv")
    batch = batch[batch["Specification"] == "Adjusted for assay batch"]
    table = pd.concat([table.iloc[:4], batch, table.iloc[4:]], ignore_index=True)
    printed = {
        "Primary model": "Primary model",
        "Visit 1": "Visit 1",
        "Visit 3": "Visit 3",
        "Visit 5": "Visit 5",
        "Adjusted for assay batch": "Adjusted for assay batch",
        "Quantifiable values only": "Quantifiable values only (excluding values reported "
                                    "below the limit of quantification)",
        "Within calibrated range": f"Within calibrated range ({low:g}–{high:g} pg/ml)",
        "Normal LV geometry": "Normal LV geometry",
    }
    rows = [[label, n_patients(r), estimate_ci(r, 2), p_value(r["p"])]
            for label, r in select(table, "Specification", printed, name="Table S2")]
    write("table_s2.csv", ["Specification", N_PATIENTS, CI, "P"], rows)


# Table S3 ------------------------------------------------------------------

def table_s3a():
    table = unrounded("04_table_s3a.csv")
    markers = ["PTH", "Phosphate", "Calcium", "eGFR", "Albumin", "FGF23"]
    rows = [[label, n_patients(r), estimate_ci(r, 2), p_value(r["p"])]
            for label, r in select(table, "Marker", {m: m for m in markers}, name="Table S3a")]
    write("table_s3a.csv", ["Marker", N_PATIENTS, "Estimate (95% CI), g/m² per SD", "P"], rows)


def table_s3b():
    table = unrounded("04_table_s3b.csv")
    baseline = unrounded("04_figure3b.csv")
    table = pd.concat([table, baseline[baseline["Observations"] == "Baseline only"]],
                      ignore_index=True)
    sets = ["All observations", "Without LVH", "Normal LV geometry", "Baseline only"]
    rows = [[label, n_patients(r), estimate_ci(r, 2, "difference"), DASH]
            for label, r in select(table, "Observations", {s: s for s in sets},
                                   name="Table S3b")]
    write("table_s3b.csv", ["Observations", N_PATIENTS,
                            "Difference, PTH − FGF23 (95% CI), g/m² per SD", "P"], rows)


S3_OUTCOMES = ["Septum, diastole", "Septum, systole", "Posterior wall, diastole",
               "LV diameter, diastole", "LV diameter, systole", "Relative wall thickness"]


def table_s3c():
    table = unrounded("04_table_s3c.csv")
    table["key"] = table["Marker"] + " | " + table["Outcome"]
    keys = [f"{m} | {o}" for m in ["PTH", "Phosphate", "Calcium", "Albumin"] for o in S3_OUTCOMES]
    rows = [[r["Marker"], r["Outcome"], n_patients(r), estimate_ci(r, 3), p_value(r["p"])]
            for _, r in select(table, "key", {k: k for k in keys}, name="Table S3c")]
    write("table_s3c.csv", ["Marker", "Outcome", N_PATIENTS,
                            "Estimate (95% CI), native units per SD", "P"], rows)


def table_s3d():
    table = unrounded("04_table_s3d.csv")
    table["key"] = table["Outcome"] + " | " + table["Exposure"]
    keys = [f"{o} | {e}" for o in S3_OUTCOMES + ["LV mass index"] for e in ["PTH", "Phosphate"]]
    rows = [[r["Outcome"], r["Exposure"],
             estimate_ci(r, 2 if r["Outcome"] == "LV mass index" else 3), p_value(r["p"])]
            for _, r in select(table, "key", {k: k for k in keys}, name="Table S3d")]
    write("table_s3d.csv", ["Outcome", "Exposure", CI, "P"], rows)


# Table S4 ------------------------------------------------------------------

def table_s4ab(name, printed):
    table = unrounded(f"09_table_s4{name}.csv")
    rows = [[label, n_patients(r), estimate_ci(r, 2), p_value(r["p"])]
            for label, r in select(table, "Outcome", printed, name=f"Table S4{name}")]
    write(f"table_s4{name}.csv", ["Outcome", N_PATIENTS,
                                  "Estimate (95% CI), mm per 10-fold PTH", "P"], rows)


def table_s4a():
    table_s4ab("a", {"Septum, diastole": "Septal thickness, diastole",
                     "Septum, systole": "Septal thickness, systole",
                     "Posterior wall, diastole": "Posterior wall, diastole",
                     "LV diameter, diastole": "LV diameter, diastole",
                     "LV diameter, systole": "LV diameter, systole"})


def table_s4b():
    table_s4ab("b", {"Septum, diastole": "Septum, diastole",
                     "Septum, systole": "Septum, systole",
                     "Posterior wall, diastole": "Posterior wall, diastole",
                     "LV diameter, diastole": "Diameter, diastole",
                     "LV diameter, systole": "Diameter, systole"})


S4C_OUTCOMES = {"E/e' (average)": "E/e′ (average)", "e' septal": "e′ septal",
                "e' lateral": "e′ lateral", "LA volume index": "LA volume index",
                "E/A (derived)": "E/A (derived)"}
AGE = "Age (positive control)"


def table_s4c():
    table = unrounded("06_table_s4c.csv")
    table["key"] = table["Exposure"] + " | " + table["Outcome"] + " | " + table["Covariates"]
    keys = [f"{e} | {o} | {c}" for e in ["PTH", "FGF23"] for c in ["primary", "prespecified"]
            for o in S4C_OUTCOMES]
    keys += [f"{AGE} | {o} | primary" for o in ["E/e' (average)", "e' septal", "e' lateral"]]
    omitted = [f"{AGE} | {o} | primary" for o in ["LA volume index", "E/A (derived)"]]
    rows = [[r["Exposure"], S4C_OUTCOMES[r["Outcome"]], r["Covariates"], r["n"],
             estimate_ci(r, 3 if r["Exposure"] == AGE else 2), p_value(r["p"])]
            for _, r in select(table, "key", {k: k for k in keys}, omitted, "Table S4c")]
    write("table_s4c.csv", ["Exposure", "Outcome", "Covariates", "n", CI, "P"], rows)


def table_s4d():
    table = unrounded("09_table_s4d.csv")
    same = ["Primary model", "Prespecified covariates", "Visit 1", "Visit 3", "Visit 5",
            "Non-hypertrophic only", "Normal LV geometry", "Volume proxies adjusted",
            "Patient means", "Between patients", "Within patients", "Visit-adjusted"]
    printed = {**{s: s for s in same},
               "Adjusted for FGF23, phosphate and calcium": "FGF23, phosphate and calcium adjusted",
               "LV mass indexed to height^2.7": "Indexed to height²·⁷ (g/m²·⁷)"}
    rows = [[label, n_patients(r), estimate_ci(r, 2), p_value(r["p"])]
            for label, r in select(table, "Specification", printed, name="Table S4d")]
    write("table_s4d.csv", ["Specification", N_PATIENTS,
                            "Estimate (95% CI), g/m² per 10-fold PTH", "P"], rows)


def table_s4e():
    table = unrounded("05_spline_contrasts.csv")
    printed = {c: c.replace("pg/mL", "pg/ml") for c in table["contrast"]}
    assert list(printed) == [f"PTH {v} vs 65 pg/mL" for v in (100, 150, 200, 300)], list(printed)
    rows = [[label, f"p{fixed(r['cohort_percentile'], 0)}", estimate_ci(r, 2)]
            for label, r in select(table, "contrast", printed, name="Table S4e")]
    write("table_s4e.csv", ["Contrast", "Cohort percentile", "Estimate (95% CI), g/m²"], rows)


# Table S5 ------------------------------------------------------------------

S5_GENES = ["FGFR1", "FGFR2", "FGFR3", "KLB", "KL", "FGFR4", "FGF23"]
S5_NOT_PRINTED = ["cell"]   # atlas label of nuclei without a resolved cell type


def table_s5():
    source = SINGLE_CELL / "table_s5.csv"
    if not source.exists():
        absent("table_s5.csv", "outputs/single_cell/tables/table_s5.csv not found "
               "(run single_cell/run_single_cell.sh, then this script)")
        return
    table = read(source)
    columns = [f"{g}_detect%" for g in S5_GENES] + [f"{g}_mean" for g in S5_GENES]
    missing = [c for c in ["cell_type"] + columns if c not in table.columns]
    if missing:
        raise SystemExit(f"outputs/single_cell/tables/table_s5.csv lacks columns {missing}")
    cell_types = [c for c in table["cell_type"] if c not in S5_NOT_PRINTED]
    values = table.set_index("cell_type")
    long = SINGLE_CELL / "by_celltype.csv"
    if long.exists():
        # Unrounded values; table_s5.csv holds them to 4 decimals.
        exact = read(long)
        exact = exact[exact["gene"].isin(S5_GENES)].set_index(["cell_type", "gene"])
        for cell_type in values.index:
            for g in S5_GENES:
                for column, field in ((f"{g}_detect%", "detection_pct"), (f"{g}_mean", "mean_expr")):
                    u = exact.loc[(cell_type, g), field]
                    if abs(float(u) - float(values.loc[cell_type, column])) > 0.00005 + 1e-9:
                        raise SystemExit(f"Table S5, {cell_type}, {column}: by_celltype.csv "
                                         f"has {u}, table_s5.csv {values.loc[cell_type, column]}")
                    values.loc[cell_type, column] = u
    else:
        print("  Table S5: by_celltype.csv not found; detection percentages rounded to 3 "
              "decimals from the 4 decimals of table_s5.csv")
    header = (["Cell type"] + [f"{g} detected (%)" for g in S5_GENES]
              + [f"{g} mean" for g in S5_GENES])
    rows = [[c[0].upper() + c[1:]]
            + [fixed(values.loc[c, f"{g}_detect%"], 3) for g in S5_GENES]
            + [fixed(values.loc[c, f"{g}_mean"], 4) for g in S5_GENES] for c in cell_types]
    write("table_s5.csv", header, rows)


def main():
    if not UNROUNDED.is_dir():
        sys.exit("outputs/unrounded/ not found: run the analysis scripts first")
    for make in (table_1, table_s1b, table_s1c, table_s2,
                 table_s3a, table_s3b, table_s3c, table_s3d,
                 table_s4a, table_s4b, table_s4c, table_s4d, table_s4e, table_s5):
        make()
    print("Publication tables written to outputs/tables/.")


if __name__ == "__main__":
    main()
