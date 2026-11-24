"""Shared steps used by every analysis and figure script.

Reads the study dataset (data/visits.csv, one row per participant visit),
applies the echocardiographic quality control and adds the derived variables
(log PTH, log FGF23, body surface area, LVH, LV geometry). Nothing is written
to disk.

Age, eGFR, BMI, relative wall thickness and LV mass index are taken from the
dataset as they are:
  Age    completed years at the visit
  eGFR   CKD-EPI 2021 creatinine equation without race
  BMI    weight / height squared
  RWT    2 x posterior wall thickness / LV internal dimension, in diastole
  LVMI   Devereux LV mass divided by Du Bois body surface area; blank where
         the echocardiogram fails quality control
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

# Covariates of the primary model, of the prespecified (smaller) model, and the
# extra covariates of the "volume proxies" sensitivity analysis.
PRIMARY_COVS = ["Age", "male", "BMI", "SBP", "DBP", "Heart_rate", "eGFR"]
PRESPECIFIED_COVS = ["Age", "male", "BMI", "SBP", "eGFR"]
VOLUME_PROXIES = ["Weight", "Albumin", "Calcium"]

# LVH: LV mass index above 115 g/m2 in men and 95 g/m2 in women.
# Increased relative wall thickness: above 0.42.
LVH_THRESHOLD = {1: 115.0, 0: 95.0}   # key: 1 = male, 0 = female
RWT_INCREASED = 0.42

# Numeric columns of data/visits.csv. E/A is recalculated from the E and A
# velocities rather than taken from the reported ratio (EA_reported).
NUMERIC = ["Age", "Height", "Weight", "BMI", "SBP", "DBP", "Heart_rate",
           "Creatinine", "Calcium", "Phosphate", "Albumin", "eGFR", "PTH", "FGF23",
           "LVIDd", "LVIDs", "IVSd", "IVSs", "PWd_used", "LVMI", "RWT",
           "E", "A", "EA_reported", "e_septal", "e_lateral", "E_e_avg", "TRV", "LAVI", "GLS"]
REQUIRED = ["Patient_ID", "Visit", "Sex"] + NUMERIC

# Participant flow for Figure 1a: patients screened, excluded at screening,
# enrolled and analysed. These counts come from the screening log, which is not
# public; the file holds the four counts only.
FLOW_FILE = ROOT / "outputs" / "02_participant_flow.csv"
FLOW_STAGES = ["screened", "excluded_at_screening", "enrolled", "analysed"]

# FGF23 assay: calibrated range of the standard curve (pg/mL).
CALIBRATED_RANGE = (78.1, 5000.0)

# PTH assay (Architect i2000, Abbott): upper limit of the reference interval
# (15.00-68.30 pg/mL).
PTH_UPPER_REFERENCE = 68.3

# Echocardiographic quality control: plausible range of each measurement (mm),
# set before the analyses were run.
QC_RANGES = {"IVSd": (5, 20), "IVSs": (6, 25), "LVIDd": (30, 60), "LVIDs": (15, 45)}


def read_visits(path=None):
    """Read data/visits.csv and check that every required column is present."""
    df = pd.read_csv(path or DATA / "visits.csv", encoding="utf-8")
    df.columns = [c.strip() for c in df.columns]
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"data/visits.csv lacks columns {missing}")
    return df


def read_flow(path=None):
    """Participant flow counts for Figure 1a, checked to add up."""
    flow = pd.read_csv(path or FLOW_FILE, dtype={"stage": str, "n": int})
    if list(flow.columns) != ["stage", "n"] or list(flow["stage"]) != FLOW_STAGES:
        raise ValueError("outputs/02_participant_flow.csv: unexpected layout")
    n = {stage: int(count) for stage, count in zip(flow["stage"], flow["n"])}
    if n["screened"] != n["enrolled"] + n["excluded_at_screening"] or \
            n["analysed"] != n["enrolled"]:
        raise ValueError("outputs/02_participant_flow.csv: counts do not add up")
    return n


def enrolment_counts(path=None, flow_path=None):
    """Participants enrolled (from the flow counts) and analysed (counted in
    data/visits.csv); the two must agree."""
    flow = read_flow(flow_path)
    analysed = read_visits(path)["Patient_ID"].nunique()
    if analysed != flow["analysed"]:
        raise ValueError(f"data/visits.csv holds {analysed} participants; "
                         f"outputs/02_participant_flow.csv gives {flow['analysed']} analysed")
    return {"enrolled": flow["enrolled"], "analysed": analysed}


def load_visits(path=None):
    """Load all visits (one row per participant visit) and add the derived
    variables used in the analyses."""
    df = read_visits(path).copy()
    df["visit"] = pd.to_numeric(df["Visit"], errors="raise").astype(int)
    for c in NUMERIC:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Sex: "Male" or "Female" in every row.
    sex = df["Sex"].str.strip().str.lower()
    if not sex.isin(["male", "female"]).all():
        raise ValueError("Sex must be Male or Female in every row")
    df["male"] = (sex == "male").astype(int)

    # Biomarkers on the log10 scale, and simple derived measures.
    df["logPTH"] = np.log10(df["PTH"])
    df["logFGF23"] = np.log10(df["FGF23"])
    df["PWTd"] = df["PWd_used"]
    df["FS"] = 100 * (df["LVIDd"] - df["LVIDs"]) / df["LVIDd"]   # fractional shortening, %
    df["EA"] = df["E"] / df["A"]

    # Body surface area (Du Bois), LVH by the sex-specific threshold, and
    # LV geometry from LVH and relative wall thickness.
    df["BSA"] = 0.007184 * df["Height"] ** 0.725 * df["Weight"] ** 0.425
    df["lvh_threshold"] = df["male"].map(LVH_THRESHOLD)
    df["LVH"] = (df["LVMI"] > df["lvh_threshold"]).astype("Int64").where(df["LVMI"].notna())
    df["RWT_high"] = (df["RWT"] > RWT_INCREASED).astype("Int64").where(df["RWT"].notna())
    lvh = df["LVH"].fillna(-1).astype(int)
    rwt_high = df["RWT_high"].fillna(-1).astype(int)
    df["geometry"] = np.select(
        [(lvh == 0) & (rwt_high == 0), (lvh == 0) & (rwt_high == 1),
         (lvh == 1) & (rwt_high == 0), (lvh == 1) & (rwt_high == 1)],
        ["Normal geometry", "Concentric remodelling", "Eccentric hypertrophy",
         "Concentric hypertrophy"], default=None)
    return df


def echo_qc(df):
    """Return all echocardiograms and those that pass quality control.

    An echocardiogram is a visit with LV internal dimension in diastole
    recorded (visits 1, 3 and 5). It passes when each measurement is within its
    range in QC_RANGES (a missing systolic septal thickness is allowed) and LV
    mass index is recorded. The script stops if these two conditions disagree.
    """
    echo = df[df["LVIDd"].notna()].copy()
    ok = pd.Series(True, index=echo.index)
    for var, (lo, hi) in QC_RANGES.items():
        ok &= echo[var].between(lo, hi) | (echo[var].isna() & (var == "IVSs"))
    has_lvmi = echo["LVMI"].notna()
    assert (ok == has_lvmi).all(), \
        f"quality control and recorded LV mass index disagree in {int((ok != has_lvmi).sum())} echocardiograms"
    return echo, echo[ok & has_lvmi].copy()


def load_linkage(path=None):
    """Load the FGF23 assay batch of each measurement (data/lab_linked_v2.csv).

    One row per FGF23 measurement (865). Columns used:
      pid, tp          participant and study visit
      batch            assay batch, B01 to B15 in the order assayed
      plate_year       year in which the sample was assayed
      limit, bound     "below" and the limit of quantification, for values
                       reported below that limit
    """
    link = pd.read_csv(path or DATA / "lab_linked_v2.csv",
                       dtype={"pid": str, "batch": str, "limit": str})
    required = {"pid", "tp", "batch", "plate_year", "limit", "bound"}
    missing = required - set(link.columns)
    if missing:
        raise ValueError(f"data/lab_linked_v2.csv lacks columns {sorted(missing)}")
    link["Patient_ID"] = pd.to_numeric(link["pid"].str.strip())
    link["visit"] = link["tp"].astype(int)
    link["run"] = link["batch"].astype(str).str.strip()
    link["limit"] = link["limit"].fillna("")
    link["plate_year"] = pd.to_numeric(link["plate_year"])
    return link.drop_duplicates(["Patient_ID", "visit"])


def cluster_ols(y, columns, groups):
    """Linear regression with standard errors that allow for repeated
    echocardiograms of the same patient (patient-clustered, small-sample
    corrected).

    Returns the coefficients (intercept first) and their covariance matrix.
    """
    X = np.column_stack([np.ones(len(y))] + [np.asarray(c, dtype=float) for c in columns])
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


def primary_frame(qc):
    """Echocardiograms with complete primary-model covariates, with each
    patient's mean log10 PTH and each visit's difference from that mean
    (used to separate between-patient from within-patient associations)."""
    covs = ["logPTH", "Age", "BMI", "SBP", "DBP", "Heart_rate", "eGFR"]
    d = qc.dropna(subset=covs).copy()
    d["mean_lp"] = d.groupby("Patient_ID")["logPTH"].transform("mean")
    d["dev_lp"] = d["logPTH"] - d["mean_lp"]
    return d
