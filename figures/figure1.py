"""Figure 1: the study cohort (panels and their source data).

  a  Participant flow. The screening and enrolment counts are read from
     outputs/02_participant_flow.csv; all other counts are computed from the
     dataset. The submitted panel a is drawn by submission_figures.py from the
     counts written here (figure1a_flow.csv).
  b  eGFR at visit 1, with KDIGO categories.
  c  Baseline LV mass index minus the sex-specific LVH threshold (115 g/m2 in
     men, 95 g/m2 in women), stacked by sex; the shaded band is +/-10 g/m2.
  d  LV geometry of all echocardiograms passing quality control.
  e  E/e' and left atrial volume index at baseline, with the ASE/EACVI cut-offs
     (14 and 34 ml/m2).
  f  Change in LV geometry over 24 months in participants with normal geometry
     at baseline and an analysable echocardiogram at visits 1, 3 and 5.

Each panel is saved at its size in the submitted figure; submission_figures.py
places panels b to f. The counts and histograms behind each panel are written
to outputs/figure1*.csv.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

from kipanel import (FS_ANNOT, FS_TINY, FS_TITLE, INK, LW_AXIS, SUBTLE, HAIR, TABLES,
                     mm_figure, plt, save)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
import prep  # noqa: E402

# Reasons for exclusion at screening, as shown in panel a.
EXCLUSION_REASONS = ["Prevalent LVH", "Established cardiac disease", "Hypertension predating CKD",
                     "eGFR outside 15\u201390"]
BLUE, GOLD, ZONE = "#2C5F8A", "#C89B4A", "#F0F2F5"
STAGES = [(0, 15, "G5", "#D6E1ED"), (15, 30, "G4", "#DFE8F1"), (30, 45, "G3b", "#E7EDF4"),
          (45, 60, "G3a", "#EFF3F7"), (60, 90, "G2", "#F5F7FA"), (90, 1e3, "G1", "#FAFBFC")]
LVH_CUT = {1: 115.0, 0: 95.0}
NEAR_THRESHOLD = 10.0
PRIME = "\u2032"


def style_axis(ax):
    """Grid, ticks and axis lines shared by panels b to f."""
    ax.grid(axis="y", color=HAIR, lw=0.35, alpha=0.9)
    ax.tick_params(length=2.0, width=LW_AXIS, pad=1.8, colors=INK)
    for spine in ax.spines.values():
        spine.set_linewidth(LW_AXIS)
        spine.set_color("#4A4A4A")


# ─────────────────────────────────────────────────────────────── panel a
def panel_a(df):
    """Participant flow chart, and its counts (figure1a_flow.csv)."""
    echo, qc = prep.echo_qc(df)
    counts = prep.enrolment_counts()
    n_enrol = counts["enrolled"]
    n_analysed = df["Patient_ID"].nunique()
    if not n_analysed == counts["analysed"] == n_enrol:
        raise ValueError("outputs/02_participant_flow.csv differs from data/visits.csv")
    recorded = prep.read_flow()
    n_screened, n_excluded = recorded["screened"], recorded["excluded_at_screening"]
    n_visits = df["visit"].nunique()
    n_all = int((df.groupby("Patient_ID")["visit"].nunique() == n_visits).sum())
    contrib = qc.groupby("Patient_ID").size().value_counts().sort_index().to_dict()
    n_rep = sum(v for k, v in contrib.items() if k >= 2)
    n_pat = qc["Patient_ID"].nunique()
    assert sum(contrib.values()) == n_pat and sum(k * v for k, v in contrib.items()) == len(qc)

    # Layout of the flow chart (font sizes, line heights, box sizes in mm).
    FS_MAIN, FS_SIDE, FS_NOTE, FS_FINAL = 6.8, 6.3, 6.0, 7.8
    LH_MAIN, LH_SIDE, LH_NOTE = 3.3, 3.1, 2.9
    PADX, PADY, INDENT = 2.8, 2.7, 3.0
    NUMW_MAIN, NUMW_SIDE, NUMW_NOTE, NUMGAP = 8.0, 7.0, 8.0, 2.2
    GREY = "#6B7280"
    EDGE_MAIN, EDGE_SIDE, TINT, RULE = "#2F3B47", "#C3CBD4", "#EEF3F8", "#DCE1E6"

    A = [(n_screened, "Patients assessed for"), (None, "eligibility")]
    SA = [(n_excluded, "Were excluded")] + [(None, r, 1) for r in EXCLUSION_REASONS]
    B = [(n_enrol, "Patients with predialysis"), (None, "CKD were enrolled")]
    SB = [(n_all, f"Completed all {n_visits} visits"), (n_analysed - n_all, "Had partial follow-up")]
    C = [(len(echo), "Echocardiograms obtained at"), (None, "visits 1, 3, and 5")]
    SC = [(len(echo) - len(qc), "Excluded owing to"), (None, "internally inconsistent wall", 1),
          (None, "measurements", 1)]
    D = [(len(qc), "Echocardiograms analyzed,"), (None, f"from {n_pat} patients")]

    MX, MW, SX, SW = 0.0, 60.0, 66.0, 46.0
    CX, W = MX + MW / 2, SX + SW
    h_of = lambda rows, lh: len(rows) * lh + 2 * PADY                    # noqa: E731
    gap = lambda hs: max(14.0, hs + 5.0)                                  # noqa: E731
    mains, sides = [A, B, C, D], [SA, SB, SC]
    hm = [h_of(r, LH_MAIN) for r in mains]
    hs = [h_of(r, LH_SIDE) for r in sides]
    ys, mids = [3.0], []
    for h_main, h_side in zip(hm[:-1], hs):
        g = gap(h_side)
        mids.append(ys[-1] + h_main + g / 2)
        ys.append(ys[-1] + h_main + g)
    yD, hD = ys[-1], hm[-1]
    y_rule = yD + hD + 8.0
    y_note = y_rule + 3.4
    H = y_note + (len(contrib) + 1) * LH_NOTE + 5.0

    fig = mm_figure(W, H)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    def box(x, y, w, rows, kind):
        side = kind == "side"
        lh, fs = (LH_SIDE, FS_SIDE) if side else (LH_MAIN, FS_MAIN)
        numw = NUMW_SIDE if side else NUMW_MAIN
        edge = EDGE_SIDE if side else (BLUE if kind == "final" else EDGE_MAIN)
        ax.add_patch(FancyBboxPatch((x, y), w, h_of(rows, lh),
                                    boxstyle="round,pad=0,rounding_size=0.6",
                                    fc=TINT if kind == "final" else "white", ec=edge,
                                    lw=0.6 if side else (1.0 if kind == "final" else 0.8),
                                    zorder=3, clip_on=False, mutation_aspect=1))
        for i, row in enumerate(rows):
            number, text = row[0], row[1]
            yy = y + PADY + (i + 0.5) * lh
            if len(row) > 2:
                ax.text(x + PADX + INDENT, yy, text, fontsize=fs - 0.3, ha="left", va="center",
                        color=GREY, zorder=4, clip_on=False)
                continue
            if number is not None:
                ax.text(x + PADX + numw, yy, f"{number:,}",
                        fontsize=FS_FINAL if kind == "final" and i == 0 else fs,
                        fontweight="normal" if side else "bold", ha="right", va="center",
                        color=GREY if side else (BLUE if kind == "final" else INK),
                        zorder=4, clip_on=False)
            ax.text(x + PADX + numw + NUMGAP, yy, text, fontsize=fs, ha="left", va="center",
                    color=GREY if side else INK, zorder=4, clip_on=False)

    def trunk(y0, y1):
        ax.annotate("", xy=(CX, y1), xytext=(CX, y0), zorder=2,
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=EDGE_MAIN,
                                    mutation_scale=8, shrinkA=0, shrinkB=0))

    def branch(y):
        ax.annotate("", xy=(SX, y), xytext=(CX, y), zorder=2,
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color=EDGE_SIDE,
                                    mutation_scale=7, shrinkA=0, shrinkB=0))
        ax.plot([CX], [y], "o", ms=1.7, color=EDGE_MAIN, zorder=5, clip_on=False)

    for (y, main, side, h_main, h_side, mid, y_next) in zip(
            ys[:-1], mains[:-1], sides, hm[:-1], hs, mids, ys[1:]):
        box(MX, y, MW, main, "main")
        trunk(y + h_main, y_next)
        branch(mid)
        box(SX, mid - h_side / 2, SW, side, "side")
    box(MX, yD, MW, D, "final")

    ax.plot([MX, MX + MW], [y_rule, y_rule], lw=0.6, color=RULE, zorder=2, clip_on=False)
    ax.text(MX, y_note, "Studies contributed per patient", fontsize=FS_NOTE, va="center",
            ha="left", color=INK, clip_on=False)
    for i, k in enumerate(sorted(contrib, reverse=True)):
        yy = y_note + (i + 1) * LH_NOTE + 0.6
        ax.text(MX + NUMW_NOTE, yy, f"{contrib[k]:,}", fontsize=FS_NOTE, ha="right",
                va="center", color=INK, clip_on=False)
        ax.text(MX + NUMW_NOTE + NUMGAP, yy,
                f"patients contributed {k} stud{'y' if k == 1 else 'ies'}",
                fontsize=FS_NOTE, ha="left", va="center", color=GREY, clip_on=False)
    ax.text(MX, y_note + (len(contrib) + 1) * LH_NOTE + 1.8,
            f"Within-patient comparisons rest on the {n_rep} patients with repeated studies.",
            fontsize=FS_NOTE - 0.4, va="center", ha="left", color=GREY, style="italic",
            clip_on=False)
    save(fig, "Figure1a")
    flow = [("screened", n_screened), ("excluded_at_screening", n_excluded),
            ("enrolled", n_enrol), ("analysed", n_analysed), ("completed_all_visits", n_all),
            ("partial_follow_up", n_analysed - n_all), ("echocardiograms_obtained", len(echo)),
            ("failed_dimension_qc", len(echo) - len(qc)), ("echocardiograms_analysed", len(qc)),
            ("patients_analysed", n_pat), ("patients_with_repeated_studies", n_rep)] + \
           [(f"patients_contributing_{k}_studies", v) for k, v in sorted(contrib.items(), reverse=True)]
    pd.DataFrame(flow, columns=["stage", "n"]).to_csv(TABLES / "figure1a_flow.csv", index=False)


# ─────────────────────────────────────────────────────────────── panel b
def panel_b(baseline):
    """Histogram of eGFR at visit 1, with the KDIGO categories above it."""
    v = baseline["eGFR"].dropna().astype(float)
    bw = 5.0
    lo, hi = np.floor(v.min() / bw) * bw, np.ceil(v.max() / bw) * bw
    counts, edges = np.histogram(v, bins=np.arange(lo, hi + bw, bw))
    assert counts.sum() == len(v)
    fig = mm_figure(56.0, 40.0)
    ax = fig.add_axes([0.155, 0.185, 0.825, 0.715])
    ymax = counts.max()
    ax.bar(edges[:-1], counts, width=bw, align="edge", color=BLUE, edgecolor="white",
           lw=0.45, zorder=3)
    ax.set_ylim(0, ymax * 1.34)
    ax.set_yticks(np.arange(0, ymax + 0.1, 5 if ymax <= 30 else 10))
    y0, y1, y_n = ymax * 1.19, ymax * 1.31, ymax * 1.09
    stages = []
    for a, b, name, tint in STAGES:
        n = int(((v >= a) & (v < b)).sum())
        if n == 0:
            continue
        aa, bb = max(a, lo), min(b, hi)
        stages.append({"stage": name, "egfr_lo": a, "egfr_hi": b, "patients": n,
                       "pct": round(100 * n / len(v), 1)})
        ax.add_patch(Rectangle((aa + 0.3, y0), (bb - aa) - 0.6, y1 - y0, fc=tint, ec="none",
                               zorder=2, clip_on=False))
        ax.text((aa + bb) / 2, (y0 + y1) / 2, name, fontsize=FS_TINY, ha="center",
                va="center", color=SUBTLE, zorder=4)
        ax.text((aa + bb) / 2, y_n, f"{n}", fontsize=FS_TINY, ha="center", va="center",
                color=INK, zorder=4)
    ax.set_xlim(lo - 2, hi + 2)
    ax.set_xlabel("eGFR (ml/min/1.73 m\u00b2)")
    ax.set_ylabel("Patients")
    ax.set_title("Kidney function at visit 1", loc="left", fontsize=FS_TITLE, color=INK, pad=4)
    style_axis(ax)
    save(fig, "Figure1b")
    pd.DataFrame({"bin_lo": edges[:-1], "bin_hi": edges[1:], "patients": counts}) \
        .to_csv(TABLES / "figure1b_histogram.csv", index=False)
    pd.DataFrame(stages).to_csv(TABLES / "figure1b_stages.csv", index=False)


# ─────────────────────────────────────────────────────────────── panel c
def panel_c(baseline):
    """Histogram of baseline LV mass index minus the LVH threshold, by sex."""
    s = baseline.dropna(subset=["LVMI"]).copy()
    s["delta"] = s["LVMI"] - s["male"].map(LVH_CUT)
    bw = 5.0
    lo, hi = np.floor(s["delta"].min() / bw) * bw, np.ceil(s["delta"].max() / bw) * bw
    bins = np.arange(lo, hi + bw, bw)
    men, _ = np.histogram(s.loc[s["male"] == 1, "delta"], bins=bins)
    women, _ = np.histogram(s.loc[s["male"] == 0, "delta"], bins=bins)
    assert men.sum() + women.sum() == len(s)
    ymax = (men + women).max()
    fig = mm_figure(56.0, 50.0)
    ax = fig.add_axes([0.155, 0.330, 0.800, 0.600])
    ax.axvspan(-NEAR_THRESHOLD, NEAR_THRESHOLD, color=ZONE, lw=0, zorder=0)
    ax.bar(bins[:-1], men, width=bw, align="edge", color=BLUE, edgecolor="white", lw=0.45,
           zorder=3, label=f"Men (n = {int(men.sum())})")
    ax.bar(bins[:-1], women, width=bw, align="edge", bottom=men, color=GOLD, edgecolor="white",
           lw=0.45, zorder=3, label=f"Women (n = {int(women.sum())})")
    ax.axvline(0, color=INK, lw=0.9, zorder=5)
    ax.set_ylim(0, ymax * 1.22)
    ax.text(0, ymax * 1.19, " LVH threshold", fontsize=FS_TINY, ha="left", va="center",
            color=INK, zorder=6)
    ax.set_xlim(lo - 2, hi + 2)
    ax.set_xlabel("LV mass index minus threshold (g/m\u00b2)")
    ax.set_ylabel("Patients")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, -0.33), fontsize=FS_ANNOT,
              handlelength=1.0, handletextpad=0.45, labelspacing=0.3, frameon=False)
    style_axis(ax)
    save(fig, "Figure1c")
    pd.DataFrame({"bin_lo": bins[:-1], "bin_hi": bins[1:], "men": men, "women": women}) \
        .to_csv(TABLES / "figure1c_histogram.csv", index=False)


# ─────────────────────────────────────────────────────────────── panel e
def _histogram(ax, values, cut, bw, xlabel, n_total):
    """One histogram with its cut-off marked; returns the bin counts."""
    v = values.dropna().astype(float)
    lo = np.floor(v.min() / bw) * bw
    bins = np.arange(lo, np.ceil(v.max() / bw) * bw + bw, bw)
    counts, _ = np.histogram(v, bins=bins)
    assert counts.sum() == len(v)
    hi = max(bins[-1], cut) + 2 * bw
    ymax = counts.max()
    ax.axvspan(cut, hi, color=ZONE, lw=0, zorder=0)
    ax.bar(bins[:-1], counts, width=bw, align="edge", color=BLUE, edgecolor="white", lw=0.4,
           zorder=3)
    ax.axvline(cut, color=INK, lw=0.9, zorder=5)
    ax.set_ylim(0, ymax * 1.24)
    ax.set_xlim(lo - bw, hi)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Patients")
    style_axis(ax)
    ax.text(cut + bw * 0.6, ymax * 1.12, f"{int((v > cut).sum())} of {len(v)} above {cut:g}",
            fontsize=FS_TINY, ha="left", va="center", color=INK, zorder=6)
    if n_total - len(v):
        ax.text(0.995, 1.12, f"{n_total - len(v)} not measured", transform=ax.transAxes,
                fontsize=FS_TINY, ha="right", va="center", color=SUBTLE)
    return pd.DataFrame({"bin_lo": bins[:-1], "bin_hi": bins[1:], "patients": counts})


def panel_e(baseline):
    """Histograms of E/e' and left atrial volume index at baseline."""
    fig = mm_figure(62.0, 48.0)
    top = fig.add_axes([0.149, 0.615, 0.736, 0.288])
    bottom = fig.add_axes([0.149, 0.140, 0.736, 0.288])
    ee = _histogram(top, baseline["E_e_avg"], 14.0, 0.5, f"E/e{PRIME} (average)", len(baseline))
    lavi = _histogram(bottom, baseline["LAVI"], 34.0, 2.0, "LA volume index (ml/m\u00b2)",
                      len(baseline))
    top.set_title("Diastolic function", loc="left", fontsize=FS_TITLE, color=INK, pad=4)
    save(fig, "Figure1e")
    ee.to_csv(TABLES / "figure1e_ee_histogram.csv", index=False)
    lavi.to_csv(TABLES / "figure1e_lavi_histogram.csv", index=False)


# ─────────────────────────────────────────────────────────────── panel d
GEOMETRY = [("Normal geometry", "Normal", BLUE), ("Concentric remodelling", "Concentric\nremodel.", "#A9C5E0"),
            ("Eccentric hypertrophy", "Eccentric\nhypertr.", GOLD),
            ("Concentric hypertrophy", "Concentric\nhypertr.", "#B04A42")]


def panel_d(qc):
    """Percentage of echocardiograms in each LV geometry category."""
    share = 100 * qc["geometry"].value_counts(normalize=True)
    values = [float(share.get(k, 0.0)) for k, _, _ in GEOMETRY]
    fig = mm_figure(56.0, 40.0)
    ax = fig.add_axes([0.155, 0.235, 0.825, 0.665])
    xs = np.arange(len(GEOMETRY))
    ax.bar(xs, values, width=0.62, color=[c for _, _, c in GEOMETRY], zorder=3)
    for x, v in zip(xs, values):
        ax.text(x, v + 1.5, f"{v:.1f}", ha="center", va="bottom", fontsize=FS_ANNOT, color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels([label for _, label, _ in GEOMETRY], fontsize=FS_TINY)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Echocardiograms (%)")
    ax.set_title("LV geometry", loc="left", fontsize=FS_TITLE, color=INK, pad=4)
    style_axis(ax)
    save(fig, "Figure1d")
    pd.DataFrame({"geometry": [k for k, _, _ in GEOMETRY],
                  "echocardiograms": [int(qc["geometry"].eq(k).sum()) for k, _, _ in GEOMETRY],
                  "pct": [round(v, 1) for v in values]}).to_csv(TABLES / "figure1d_geometry.csv", index=False)


# ─────────────────────────────────────────────────────────────── panel f
def panel_f(qc):
    """Geometry at visits 3 and 5 in participants normal at visit 1."""
    wide = qc.pivot_table(index="Patient_ID", columns="visit", values="geometry", aggfunc="first")
    wide = wide.dropna(subset=[1, 3, 5])
    base = wide[wide[1] == "Normal geometry"]
    ab3, ab5 = base[3] != "Normal geometry", base[5] != "Normal geometry"
    groups = [("Remained normal", int((~ab3 & ~ab5).sum()), BLUE),
              ("Late onset (visit 5)", int((~ab3 & ab5).sum()), GOLD),
              ("Persistent conversion", int((ab3 & ab5).sum()), "#B04A42"),
              ("Transient", int((ab3 & ~ab5).sum()), "#CFCFCF")]
    assert sum(n for _, n, _ in groups) == len(base)
    fig = mm_figure(60.0, 40.0)
    ax = fig.add_axes([0.395, 0.235, 0.545, 0.665])
    ys = np.arange(len(groups))[::-1]
    ax.barh(ys, [n for _, n, _ in groups], height=0.55, color=[c for _, _, c in groups], zorder=3)
    for y, (_, n, _) in zip(ys, groups):
        ax.text(n + 0.015 * len(base), y, f"{n}", va="center", ha="left", fontsize=FS_ANNOT, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels([label for label, _, _ in groups], fontsize=FS_ANNOT)
    ax.set_xlim(0, len(base) * 0.95)
    ax.set_xlabel(f"Patients (of {len(base)} normal at baseline)")
    ax.set_title("Change over 24 months", loc="left", fontsize=FS_TITLE, color=INK, pad=4)
    style_axis(ax)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=HAIR, lw=0.35, alpha=0.9)
    ax.tick_params(axis="y", length=0)
    save(fig, "Figure1f")
    pd.DataFrame([(label, n) for label, n, _ in groups], columns=["trajectory", "patients"]) \
        .to_csv(TABLES / "figure1f_trajectories.csv", index=False)


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    baseline = df[df["visit"] == 1]
    panel_a(df)
    panel_b(baseline)
    panel_c(baseline)
    panel_d(qc)
    panel_e(baseline)
    panel_f(qc)


if __name__ == "__main__":
    main()
