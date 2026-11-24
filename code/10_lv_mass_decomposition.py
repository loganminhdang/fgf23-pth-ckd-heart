"""How much of each marker's LV mass association comes from wall thickness and
how much from cavity size (Figure 3 legend).

LV mass (Devereux) = 0.8 x 1.04 x [(LVIDd + IVSd + PWTd)^3 - LVIDd^3] + 0.6 g,
with dimensions in cm. Around the mean dimensions of the echocardiograms, a
1 mm change in each dimension changes LV mass by a fixed number of grams.
Multiplying the change in each dimension associated with a marker by that
number gives the grams contributed by the septum, the posterior wall and the
cavity. Their sum is compared with the marker's association with LV mass
itself (its LV mass index association x mean body surface area).

Associations are per SD of each marker (log10 PTH for PTH), adjusted for the
primary covariates. The script stops if the reported percentages (PTH 84%,
calcium 70%, albumin 85%) are not reproduced.

Writes outputs/10_lv_mass_decomposition.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import prep

OUT = Path(__file__).resolve().parents[1] / "outputs"
DEVEREUX = 0.8 * 1.04
MARKERS = [("PTH", "logPTH"), ("Phosphate", "Phosphate"),
           ("Calcium", "Calcium"), ("Albumin", "Albumin")]
REPORTED_RECONSTRUCTED = {"PTH": 84, "Calcium": 70, "Albumin": 85}   # per cent


def per_sd(frame, marker, outcome):
    """Change in the outcome per SD of the marker, adjusted for the primary
    covariates."""
    s = frame.dropna(subset=[outcome, marker] + prep.PRIMARY_COVS).copy()
    s["z"] = (s[marker] - s[marker].mean()) / s[marker].std()
    beta, _ = prep.cluster_ols(s[outcome].to_numpy(float),
                               [s[c] for c in prep.PRIMARY_COVS] + [s["z"]],
                               s["Patient_ID"].to_numpy())
    return float(beta[-1])


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    # Grams of LV mass per 1 mm change in each dimension, at the mean dimensions.
    s = qc.dropna(subset=["LVIDd", "IVSd", "PWTd", "BSA"])
    outer = ((s["LVIDd"] + s["IVSd"] + s["PWTd"]) / 10).mean()
    cavity = (s["LVIDd"] / 10).mean()
    sensitivity = {"IVSd": 3 * DEVEREUX * outer ** 2 / 10,                 # g per mm
                   "PWTd": 3 * DEVEREUX * outer ** 2 / 10,
                   "LVIDd": 3 * DEVEREUX * (outer ** 2 - cavity ** 2) / 10}

    # Contribution of each dimension, and the total, for each marker.
    rows = []
    for name, column in MARKERS:
        observed = per_sd(qc, column, "LVMI") * s["BSA"].mean()
        parts = {k: per_sd(qc, column, k) * v for k, v in sensitivity.items()}
        total = sum(parts.values())
        rows.append({"marker": name, "observed_g": observed,
                     "septum_g": parts["IVSd"], "posterior_wall_g": parts["PWTd"],
                     "cavity_g": parts["LVIDd"], "reconstructed_g": total,
                     "reconstructed_pct": 100 * total / observed,
                     "walls_pct_of_reconstruction": 100 * (parts["IVSd"] + parts["PWTd"]) / total})
        if name in REPORTED_RECONSTRUCTED:
            assert round(rows[-1]["reconstructed_pct"]) == REPORTED_RECONSTRUCTED[name], rows[-1]
    out = pd.DataFrame(rows).round(2)
    OUT.mkdir(exist_ok=True)
    out.to_csv(OUT / "10_lv_mass_decomposition.csv", index=False)
    print(out.to_string(index=False))
    print("\nReconstruction percentages reproduced.")


if __name__ == "__main__":
    main()
