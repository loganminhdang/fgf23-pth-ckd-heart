"""Figures 1 to 4 and Supplementary Figures S1 and S2, as submitted.

  Figure 1   participant flow (a), drawn here; panels b to f from figure1.py
  Figure 2   panels a to d, drawn with the panel functions of figure2.py
  Figure 3   from Tables S3a to S3c (eGFR not shown in panel a)
  Figure 4   spline curves (a, b), Tables S4a to S4c (c, d) and PTH relative
             to the upper reference limit (e)
  Figure S1  panels a and b from figureS1.py, stacked
  Figure S2  Table S4d (a) and the spline contrasts (b)

Values are read from outputs/, so the analysis scripts, figure1.py and
figureS1.py run first (run_all.sh does this). Written to
outputs/figures/submission/ as TIFF (500 dpi), PNG, PDF and SVG. In the SVG
copies the pasted panels are replaced by their own SVGs (vector_svg.py), so
every line and label can be edited.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from matplotlib.patches import Rectangle

from kipanel import (ACCENT, FS_ANNOT, FS_LAB, FS_TICK, HAIR, INK, LW_AXIS, LW_DATA,
                     LW_REF, NEUTRAL, OUT, PRIMARY, SUBTLE, TABLES, forest, frame,
                     mm_figure, plt)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure2  # noqa: E402
from vector_svg import replace_pasted_panels  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
import prep  # noqa: E402

DEST = OUT / "submission"
PASTED = []   # panels pasted into the current figure, for its vector SVG copy
DPI = 500
TITLE_COLOUR = "#2E3192"
FS_PANEL_TITLE = 8.2
FS_LETTER = 11.0

# Colours of the Figure 4c and 4d panels (figures/style.py).
GREY, LIGHT_GREY, GRID, DASH, TEXT = "#4A4A4A", "#8A8F94", "#D8D8D8", "#7A7A7A", "#111111"
BAND, STRIP = "#DCE7F2", "#C3CCD6"


def axes_mm(fig, size, x0, x1, top, bottom):
    """Axes spanning x0-x1 mm from the left and top-bottom mm from the top."""
    w, h = size
    return fig.add_axes([x0 / w, 1 - bottom / h, (x1 - x0) / w, (bottom - top) / h])


def restyle(ax, ytick=None, xtick=None, xlabel=None, ylabel=None):
    """Set font sizes after a panel function has drawn into ax."""
    if ytick:
        ax.tick_params(axis="y", labelsize=ytick)
    if xtick:
        ax.tick_params(axis="x", labelsize=xtick)
    if xlabel:
        ax.xaxis.label.set_size(xlabel)
    if ylabel:
        ax.yaxis.label.set_size(ylabel)


def title(ax, text, size=None):
    """Centred panel title in the submitted colour, replacing any left title."""
    ax.set_title("", loc="left")
    if text:
        ax.set_title(text, loc="center", fontsize=size or FS_PANEL_TITLE,
                     color=TITLE_COLOUR, pad=5)


def letter(fig, size, label, x, y):
    """Panel letter with its top-left corner at x, y mm from the top left."""
    fig.text(x / size[0], 1 - y / size[1], label, fontsize=FS_LETTER,
             fontweight="bold", va="top", ha="left", color="black")


def ytick_marks(ax):
    """Short tick marks on the y axis."""
    ax.tick_params(axis="y", length=2.0, width=LW_AXIS, color=INK)


def write(fig, name):
    """Save SVG, PDF, PNG and TIFF (500 dpi) to outputs/figures/submission/."""
    DEST.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf"):
        fig.savefig(DEST / f"{name}.{extension}")
    fig.savefig(DEST / f"{name}.png", dpi=DPI)
    plt.close(fig)
    if PASTED:   # the SVG copy gets the panels' own vector SVGs instead of their pictures
        replace_pasted_panels(DEST / f"{name}.svg", PASTED)
        PASTED.clear()
    with Image.open(DEST / f"{name}.png") as im:
        im.convert("RGB").save(DEST / f"{name}.tif", dpi=(DPI, DPI), compression="tiff_lzw")
    print(f"{name}: written to {DEST}")


# --------------------------------------------------------------------------- Figure 1
RED, BOX_EDGE, TINT = "#D7262B", "#3A3A3A", "#EEF2F6"


def flow_panel(ax, flow):
    """Participant flow (coordinates in mm, y downwards)."""
    n = dict(zip(flow["stage"], flow["n"]))
    ax.set_xlim(0, 78)
    ax.set_ylim(84, 0)
    ax.axis("off")
    main = [(2, f"{n['screened']}", "Patients assessed for eligibility", False),
            (18.5, f"{n['enrolled']}", "Patients with moderate-to-advanced\nnon-dialysis CKD enrolled", False),
            (35, f"{n['echocardiograms_obtained']}", "Echocardiograms obtained\nat visits 1, 3, and 5", False),
            (51.5, f"{n['echocardiograms_analysed']}",
             f"Echocardiograms analyzed\nfrom {n['patients_analysed']} patients", True)]
    x0, w, h = 2.0, 40.0, 7.0
    for y, num, text, final in main:
        ax.add_patch(Rectangle((x0, y), w, h, fc=TINT if final else "white",
                               ec="#2C5F8A" if final else BOX_EDGE, lw=0.8))
        ax.text(x0 + 2.0, y + h / 2, num, fontsize=8.5, va="center", ha="left",
                color="#2C5F8A" if final else INK)
        ax.text(x0 + 9.0, y + h / 2, text, fontsize=5.2, va="center", ha="left", color=INK,
                linespacing=1.25)
    cx = x0 + w / 2
    side = [(8.05, 11.4, [f"{n['excluded_at_screening']}   Were excluded", "Prevalent LVH",
                          "Established cardiac disease", "Hypertension predating CKD",
                          "eGFR outside 15\u201390"]),
            (27.45, 5.6, [f"{n['completed_all_visits']}   Completed all 5 visits",
                          f"  {n['partial_follow_up']}   Had partial follow-up"]),
            (43.05, 7.4, [f"{n['failed_dimension_qc']}   Excluded owing to",
                          "internally inconsistent", "wall measurements"])]
    for (y_box, h_box, lines), (y_top, _, _, _), (y_next, _, _, _) in zip(side, main[:-1], main[1:]):
        ax.annotate("", xy=(cx, y_next), xytext=(cx, y_top + h),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=BOX_EDGE, mutation_scale=6))
        y_mid = y_box + h_box / 2
        ax.plot([cx], [y_mid], "o", ms=1.6, color=BOX_EDGE)
        ax.annotate("", xy=(48.0, y_mid), xytext=(cx, y_mid),
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color=RED, mutation_scale=5))
        ax.add_patch(Rectangle((48.0, y_box), 29.5, h_box, fc="white", ec=RED, lw=0.7))
        ax.text(49.5, y_mid, "\n".join(lines), fontsize=4.6, va="center", ha="left", color=RED,
                linespacing=1.3)
    k = {c: n[f"patients_contributing_{c}_studies"] for c in (3, 2, 1)}
    ax.text(x0, 63.5, "Echocardiographic assessments per patient", fontsize=4.6, color=INK)
    for i, c in enumerate((3, 2, 1)):
        ax.text(x0 + 4.0, 66.0 + 2.1 * i, f"{k[c]}", fontsize=4.2, ha="right", color=INK)
        word = "assessment" if c == 1 else "assessments"
        ax.text(x0 + 5.0, 66.0 + 2.1 * i, f"patients underwent {c} {word}", fontsize=4.2,
                color=SUBTLE)
    ax.text(x0, 73.5, f"{n['patients_with_repeated_studies']} patients with repeated "
            "echocardiographic assessments.", fontsize=4.4, color=INK)


def paste(fig, size, path, x0, top, width):
    """Place a panel PNG at x0, top (mm from the top left), scaled to width mm."""
    with Image.open(path) as im:
        rgb = np.asarray(im.convert("RGB"))
    height = width * rgb.shape[0] / rgb.shape[1]
    ax = axes_mm(fig, size, x0, x0 + width, top, top + height)
    ax.imshow(rgb, interpolation="lanczos")
    ax.axis("off")
    PASTED.append((Path(path).with_suffix(".svg"), x0, top, width))
    return ax


def figure_1():
    """Figure 1: flow chart (a) and the panels drawn by figure1.py (b to f)."""
    flow = pd.read_csv(TABLES / "figure1a_flow.csv")
    panels = OUT
    size = (150.0, 136.0)
    fig = mm_figure(*size)
    ax_a = axes_mm(fig, size, 4.0, 82.0, 4.0, 88.0)
    flow_panel(ax_a, flow)
    paste(fig, size, panels / "Figure1b.png", 92.0, 4.0, 52.0)
    paste(fig, size, panels / "Figure1c.png", 92.0, 44.0, 52.0)
    paste(fig, size, panels / "Figure1d.png", 2.0, 94.0, 46.0)
    paste(fig, size, panels / "Figure1e.png", 49.0, 92.0, 51.0)
    paste(fig, size, panels / "Figure1f.png", 101.0, 94.0, 48.0)
    for label, (x, y) in zip("abcdef", [(1.0, 1.0), (88.0, 1.0), (88.0, 42.5), (1.0, 92.0),
                                         (49.0, 92.0), (100.0, 92.0)]):
        letter(fig, size, label, x, y)
    write(fig, "Figure 1")


# --------------------------------------------------------------------------- Figure 2
def figure_2(qc):
    """Figure 2: panels a to d, one column."""
    size = (90.1, 165.0)
    fig = mm_figure(*size)
    axes = [axes_mm(fig, size, 18.9, 77.1, 6.2, 30.6),
            axes_mm(fig, size, 25.8, 81.5, 47.5, 71.2),
            axes_mm(fig, size, 20.9, 79.2, 88.3, 112.2),
            axes_mm(fig, size, 25.4, 81.3, 129.2, 154.0)]
    figure2.panel_a(axes[0], qc)
    figure2.panel_b(axes[1], qc)
    figure2.panel_c(axes[2], qc)
    figure2.panel_d(axes[3])
    for ax, text in zip(axes, ["", "Across cardiac phenotypes",
                               "Across visits and patient comparison", "Sensitivity to assay"]):
        title(ax, text, size=7.1)
    restyle(axes[0], ytick=6.2, xtick=6.6, xlabel=6.2, ylabel=6.2)
    for ax in axes[1:]:
        ytick_marks(ax)
        restyle(ax, ytick=5.6, xtick=5.2, xlabel=5.9)
    axes[2].set_xticks(range(-4, 5))
    axes[2].set_xticklabels([f"{v}".replace("-", "\u2212") for v in range(-4, 5)])
    for label, (x, y) in zip("abcd", [(1.0, 2.0), (2.0, 43.6), (2.0, 84.4), (2.0, 125.3)]):
        letter(fig, size, label, x, y)
    write(fig, "Figure 2")


# --------------------------------------------------------------------------- Figure 3
MARKERS_C = ["PTH", "Phosphate", "Calcium", "Albumin"]
ROWS_C = {"Septum, diastole": 5.2, "Septum, systole": 4.4, "Posterior wall, diastole": 3.6,
          "LV diameter, diastole": 2.4, "LV diameter, systole": 1.6}


def figure_3():
    """Figure 3: markers and LV mass (a), PTH versus FGF23 (b), dimensions (c)."""
    size = (97.1, 142.7)
    fig = mm_figure(*size)

    a = pd.read_csv(TABLES / "04_table_s3a.csv")
    a = a[a["Marker"] != "eGFR"].sort_values("estimate", ascending=False)
    rows = [(r.Marker, r.estimate, r.ci_low, r.ci_high,
             ACCENT if r.Marker == "FGF23" else PRIMARY if r.Marker == "PTH" else NEUTRAL)
            for r in a.itertuples()][::-1]
    ax_a = axes_mm(fig, size, 31.8, 65.7, 7.4, 38.9)
    forest(ax_a, rows, "Change in LV mass index (g/m\u00b2) per SD", "")
    ax_a.set_xticks([-5.0, -2.5, 0.0, 2.5, 5.0])
    ax_a.set_xticklabels(["\u22125.0", "\u22122.5", "0.0", "2.5", "5.0"])

    b = pd.read_csv(TABLES / "04_figure3b.csv")
    rows = [(r.Observations, r.difference, r.ci_low, r.ci_high, PRIMARY) for r in b.itertuples()]
    ax_b = axes_mm(fig, size, 35.8, 69.6, 52.3, 89.6)
    forest(ax_b, rows, "Difference in effect, PTH minus FGF23 (g/m\u00b2 per SD)", "")
    for ax in (ax_a, ax_b):
        ytick_marks(ax)
        restyle(ax, ytick=6.0, xtick=6.0, xlabel=6.2)

    c = pd.read_csv(TABLES / "04_table_s3c.csv")
    c = c[c["Outcome"].isin(ROWS_C)]
    limit = max(0.8, 1.05 * float(c[["ci_low", "ci_high"]].abs().max().max()))
    axes_c = []
    for j, marker in enumerate(MARKERS_C):
        x0 = 26.9 + j * 16.4
        ax = axes_mm(fig, size, x0, x0 + 13.6, 105.3, 134.1)
        axes_c.append(ax)
        colour = PRIMARY if marker == "PTH" else NEUTRAL
        for x in (-0.5, 0.5):
            ax.axvline(x, color=HAIR, lw=0.35, zorder=0)
        ax.axvline(0, ls=(0, (2.5, 2)), lw=LW_REF, color=SUBTLE, zorder=1)
        for r in c[c["Marker"] == marker].itertuples():
            y = ROWS_C[r.Outcome]
            ax.plot([r.ci_low, r.ci_high], [y, y], color=colour, lw=LW_DATA,
                    solid_capstyle="butt", zorder=2)
            for x in (r.ci_low, r.ci_high):
                ax.plot([x, x], [y - 0.13, y + 0.13], color=colour, lw=LW_AXIS + 0.15, zorder=2)
            ax.plot([r.estimate], [y], "o", color=colour, ms=2.6, mec="white", mew=0.45, zorder=3)
        ax.set_xlim(-limit, limit)
        ax.set_xticks([-0.5, 0, 0.5])
        ax.set_xticklabels(["\u22120.5", "0", "0.5"])
        ax.set_ylim(1.0, 5.9)
        ax.set_yticks(list(ROWS_C.values()))
        ax.set_yticklabels(list(ROWS_C) if j == 0 else [])
        ax.set_title(marker, loc="center", fontsize=5.3, color=INK, pad=3)
        ax.spines["left"].set_visible(False)
        frame(ax, "", None)
        ytick_marks(ax)
        restyle(ax, ytick=4.6, xtick=5.4)
    fig.text((26.9 + 89.7) / 2 / size[0], 1 - 139.4 / size[1], "Change (mm) per SD of marker",
             ha="center", va="center", fontsize=6.2, color=INK)

    for label, (x, y) in zip("abc", [(2.0, 2.0), (2.0, 47.5), (2.0, 95.0)]):
        letter(fig, size, label, x, y)
    write(fig, "Figure 3")


# --------------------------------------------------------------------------- Figure 4
ROWS_4C = [("Septum, diastole", "Septum, diastole", 6.5),
           ("Septum, systole", "Septum, systole", 5.6),
           ("Posterior wall, diastole", "Posterior wall, diastole", 4.7),
           ("LV diameter, diastole", "Diameter, diastole", 3.1),
           ("LV diameter, systole", "Diameter, systole", 2.2)]
ROWS_4D = [("E/e' (average)", "E/e\u2032 (average)", 5.0), ("e' septal", "e\u2032 septal", 4.2),
           ("e' lateral", "e\u2032 lateral", 3.4), ("LA volume index", "LA volume index", 2.6),
           ("E/A (derived)", "E/A", 1.8)]


def forest_frame(ax):
    """Axis style of the Figure 4 and Supplementary Figure S2 forest plots."""
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.spines["bottom"].set_color(GREY)
    ax.xaxis.grid(True, color=GRID, lw=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=FS_TICK, colors=TEXT, width=0.6, length=2.2, color=GREY)
    ax.tick_params(axis="y", length=0, labelsize=FS_TICK, colors=TEXT)


def spline_panel(ax, curves, strips, panel, ylim, xlim):
    """Figure 4a or 4b: spline curve with its 95% band, and the number of
    echocardiograms across PTH as a strip along the bottom."""
    c, s = curves[curves["figure"] == panel], strips[strips["figure"] == panel]
    span = ylim[1] - ylim[0]
    bottom = ylim[0] - 0.16 * span
    ax.set_ylim(bottom, ylim[1])
    ax.set_xlim(*xlim)
    ax.set_xscale("log")
    heights = 0.105 * span * s["count"] / s["count"].max()
    ax.bar(10 ** s["log10_pth_left"], heights,
           width=10 ** s["log10_pth_right"] - 10 ** s["log10_pth_left"],
           bottom=bottom + 0.012 * span, align="edge", color=STRIP, lw=0, zorder=2)
    ax.fill_between(c["pth_pg_ml"], c["ci_low"], c["ci_high"], color=BAND, lw=0, zorder=3)
    ax.plot(c["pth_pg_ml"], c["adjusted_lvmi"], color=PRIMARY, lw=LW_DATA + 0.15, zorder=4)
    ax.set_yticks([t for t in (80, 100) if ylim[0] <= t <= ylim[1]])
    frame(ax, "PTH (pg/ml, log scale)", None,
          ylabel="Adjusted LV mass index (g/m\u00b2)", grid="y")


def figure_4():
    """Figure 4: spline curves (a, b), dimensions (c), diastolic function (d),
    PTH relative to the upper reference limit (e)."""
    curves = pd.read_csv(TABLES / "05_figure4ab_curves.csv")
    strips = pd.read_csv(TABLES / "05_figure4ab_strips.csv")
    low, high = curves["ci_low"].min(), curves["ci_high"].max()
    pad = 0.07 * (high - low)
    ylim = (low - pad, high + pad)
    xlim = (10 ** (strips["log10_pth_left"].min() - 0.04),
            10 ** (strips["log10_pth_right"].max() + 0.04))
    every = pd.read_csv(TABLES / "09_table_s4a.csv").set_index("Outcome")
    normal = pd.read_csv(TABLES / "09_table_s4b.csv").set_index("Outcome")
    diastolic = pd.read_csv(TABLES / "06_table_s4c.csv")
    diastolic = diastolic[(diastolic["Exposure"] == "PTH")
                          & (diastolic["Covariates"] == "primary")].set_index("Outcome")

    uln = pd.read_csv(TABLES / "05_pth_uln_categories.csv")

    size = (160.0, 158.0)
    fig = mm_figure(*size)
    ax_a = axes_mm(fig, size, 31.1, 68.0, 11.5, 40.1)
    ax_b = axes_mm(fig, size, 105.3, 140.1, 10.9, 39.8)
    ax_c = axes_mm(fig, size, 32.8, 75.5, 60.8, 98.8)
    ax_d = axes_mm(fig, size, 98.3, 146.5, 61.4, 94.1)

    for ax, panel in ((ax_a, "4a"), (ax_b, "4b")):
        spline_panel(ax, curves, strips, panel, ylim, xlim)
        restyle(ax, ytick=8.8, xtick=8.8, xlabel=7.0, ylabel=7.0)

    # c: PTH and each dimension, all echocardiograms and normal geometry.
    ax = ax_c
    ax.axvline(0, ls=(0, (3, 3)), lw=0.9, color=DASH, zorder=1)
    offset = 0.19
    for i, (outcome, _, y) in enumerate(ROWS_4C):
        f, g = every.loc[outcome], normal.loc[outcome]
        ax.errorbar(f.estimate, y + offset, xerr=[[f.estimate - f.ci_low], [f.ci_high - f.estimate]],
                    fmt="o", ms=3.6, mfc="white", mec=LIGHT_GREY, mew=1.0, elinewidth=1.25,
                    capsize=1.8, capthick=1.25, color=LIGHT_GREY, zorder=3)
        ax.errorbar(g.estimate, y - offset, xerr=[[g.estimate - g.ci_low], [g.ci_high - g.estimate]],
                    fmt="o", ms=3.6, mew=0, elinewidth=1.25, capsize=1.8, capthick=1.25,
                    color=PRIMARY, zorder=4)
        if i == 0:
            ax.text(f.ci_high + 0.12, y + offset, "All observations", fontsize=5.4,
                    color=LIGHT_GREY, va="center", ha="left", zorder=5)
            ax.text(g.ci_high + 0.12, y - offset, "Normal LV geometry", fontsize=5.4,
                    color=PRIMARY, va="center", ha="left", zorder=5)
    lows = pd.concat([every["ci_low"], normal["ci_low"]])
    highs = pd.concat([every["ci_high"], normal["ci_high"]])
    ax.set_xlim(min(-1.15, lows.min() - 0.1), max(3.05, highs.max() + 0.1))
    ax.set_xticks([-1, 0, 1, 2, 3])
    ax.set_xticklabels(["\u22121", "0", "1", "2", "3"])
    ax.set_yticks([y for _, _, y in ROWS_4C])
    ax.set_yticklabels([label for _, label, _ in ROWS_4C])
    ax.set_ylim(1.5, 7.2)
    ax.set_xlabel("Change (mm) per 10-fold PTH", color=TEXT)
    forest_frame(ax)
    restyle(ax, ytick=6.0, xtick=6.0, xlabel=6.4)

    # d: PTH and diastolic function (primary covariates).
    ax = ax_d
    ax.axvline(0, ls=(0, (3, 3)), lw=0.9, color=DASH, zorder=1)
    for outcome, _, y in ROWS_4D:
        r = diastolic.loc[outcome]
        ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]],
                    fmt="o", ms=3.8, mew=0, elinewidth=1.4, capsize=2.0, capthick=1.4,
                    color=GREY, zorder=3)
    ax.set_yticks([y for _, _, y in ROWS_4D])
    ax.set_yticklabels([label for _, label, _ in ROWS_4D])
    ax.set_ylim(1.3, 5.5)
    ax.set_xticks([-1.0, -0.5, 0.0, 0.5])
    ax.set_xticklabels(["\u22121.0", "\u22120.5", "0.0", "0.5"])
    ax.set_xlabel("Standardized change per 10-fold PTH", color=TEXT)
    forest_frame(ax)
    restyle(ax, ytick=6.2, xtick=5.8, xlabel=6.0)

    # e: LV mass index by PTH relative to the upper reference limit of the assay.
    ax_e = axes_mm(fig, size, 52.0, 118.0, 117.0, 147.0)
    ax = ax_e
    ax.axvline(0, ls=(0, (3, 3)), lw=0.9, color=DASH, zorder=1)
    names = {"Up to 1x": "\u2264 1\u00d7 ULN", "1-2x": "> 1\u20132\u00d7 ULN",
             "2-4x": "> 2\u20134\u00d7 ULN", "Above 4x": "> 4\u00d7 ULN"}
    ys = list(range(len(uln), 0, -1))
    for y, r in zip(ys, uln.itertuples()):
        if pd.isna(r.ci_low):
            ax.plot(0, y, "o", ms=3.8, mfc="white", mec=GREY, mew=1.0, zorder=3)
        else:
            ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]],
                        fmt="o", ms=3.8, mew=0, elinewidth=1.4, capsize=2.0, capthick=1.4,
                        color=PRIMARY, zorder=3)
        upper = "" if pd.isna(r.pth_to_pg_ml) or r.pth_to_pg_ml == "" else f"{float(r.pth_to_pg_ml):.0f}"
        span = (f"\u2264 {upper}" if r.Index == 0 else
                (f"> {r.pth_from_pg_ml:.0f}" if not upper else f"{r.pth_from_pg_ml:.0f}\u2013{upper}"))
        ax.text(1.03, y, f"{span} pg/ml", transform=ax.get_yaxis_transform(),
                fontsize=5.6, color=TEXT, va="center", ha="left")
        ax.text(1.27, y, f"{r.n} ({r.patients})", transform=ax.get_yaxis_transform(),
                fontsize=5.6, color=TEXT, va="center", ha="left")
    for x_text, header in ((1.03, "PTH"), (1.27, "Echocardiograms\n(patients)")):
        ax.text(x_text, len(uln) + 0.45, header, transform=ax.get_yaxis_transform(),
                fontsize=5.6, color=SUBTLE, va="bottom", ha="left", linespacing=1.2)
    ax.set_yticks(ys)
    ax.set_yticklabels([names[g] for g in uln["group"]])
    ax.set_ylim(0.4, len(uln) + 1.0)
    lo = min(-6.0, uln["ci_low"].min() - 1.0)
    hi = max(26.0, uln["ci_high"].max() + 1.0)
    ax.set_xlim(lo, hi)
    ax.set_xticks([-5, 0, 5, 10, 15, 20, 25])
    ax.set_xticklabels(["\u22125", "0", "5", "10", "15", "20", "25"])
    ax.set_xlabel("Difference in LV mass index (g/m\u00b2)", color=TEXT)
    forest_frame(ax)
    restyle(ax, ytick=6.2, xtick=5.8, xlabel=6.0)

    for ax, text in zip((ax_a, ax_b, ax_c, ax_d, ax_e),
                        ("All echocardiograms", "Normal LV geometry only",
                         "Anatomical localization", "Diastolic function",
                         "PTH relative to the upper reference limit")):
        title(ax, text, size=9.0)
    for label, (x, y) in zip("abcde", [(14.0, 3.0), (89.0, 3.0), (6.0, 52.0), (88.0, 52.0),
                                        (6.0, 108.0)]):
        letter(fig, size, label, x, y)
    write(fig, "Figure 4")


# --------------------------------------------------------------------------- Supplementary Figure S1
def figure_s1():
    """Panels a and b of figures/figureS1.py, stacked, as in the Supplementary Information."""
    size = (114.4, 118.8)
    fig = mm_figure(*size)
    top_a, width = 4.0, 110.0
    with Image.open(OUT / "FigureS1a.png") as im:
        height_a = width * im.height / im.width
    paste(fig, size, OUT / "FigureS1a.png", 2.0, top_a, width)
    paste(fig, size, OUT / "FigureS1b.png", 2.0, top_a + height_a + 4.0, width)
    letter(fig, size, "a", 1.0, 1.0)
    letter(fig, size, "b", 1.0, top_a + height_a + 2.0)
    write(fig, "Supplementary Figure S1")


# --------------------------------------------------------------------------- Supplementary Figure S2
S2A_ROWS = [("Primary model", "Primary model"),
            ("Prespecified covariates", "Prespecified covariates only"),
            ("Visit-adjusted", "Adjusted for study visit"),
            ("Visit 1", "Visit 1 only"), ("Visit 3", "Visit 3 only"), ("Visit 5", "Visit 5 only"),
            ("Patient means", "Patient averages"),
            ("Between patients", "Between-patient differences"),
            ("Within patients", "Within-patient change"),
            ("Volume proxies adjusted", "Adjusted for weight, albumin, calcium"),
            ("Adjusted for FGF23, phosphate and calcium", "Adjusted for FGF23, phosphate, calcium"),
            ("Non-hypertrophic only", "Observations without LVH"),
            ("Normal LV geometry", "Normal LV geometry")]
S2B_TICKS = [65, 100, 200, 300, 600]


def figure_s2():
    """Supplementary Figure S2: PTH across model specifications (a) and LV mass
    index across PTH relative to 65 pg/ml (b)."""
    size = (134.1, 65.5)
    fig = mm_figure(*size)
    spec = pd.read_csv(TABLES / "09_table_s4d.csv").set_index("Specification")
    ax = axes_mm(fig, size, 37.0, 80.0, 9.5, 57.0)
    ys = np.arange(len(S2A_ROWS))[::-1]
    ax.axvline(0, ls=(0, (2, 2)), lw=0.5, color=DASH, zorder=1)
    for y, (key, _) in zip(ys, S2A_ROWS):
        r = spec.loc[key]
        colour = PRIMARY if key in ("Primary model", "Normal LV geometry") else "#6E6E6E"
        ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]],
                    fmt="o", ms=2.6, mew=0.4, mec="white", elinewidth=0.9, capsize=1.6,
                    capthick=0.9, color=colour, zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels([label for _, label in S2A_ROWS])
    ax.set_ylim(-0.6, len(S2A_ROWS) - 0.4)
    ax.set_xlim(min(-1.0, spec["ci_low"].min() - 0.8), max(27.0, spec["ci_high"].max() + 0.8))
    ax.set_xticks([0, 5, 10, 15, 20, 25])
    ax.set_xlabel("Change in LV mass index (g/m\u00b2) per 10-fold PTH", color=TEXT)
    forest_frame(ax)
    restyle(ax, ytick=5.6, xtick=5.2, xlabel=5.6)
    ax.set_title("PTH and LVMI across analyses", fontsize=6.6, fontweight="bold",
                 color=TITLE_COLOUR, pad=5)

    curve = pd.read_csv(TABLES / "05_spline_curve.csv")
    contrasts = pd.read_csv(TABLES / "05_spline_contrasts.csv")
    histogram = pd.read_csv(TABLES / "05_pth_histogram.csv")
    x = curve["log10_pth"].to_numpy()
    ax = axes_mm(fig, size, 95.0, 127.0, 11.5, 44.5)
    ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.axhline(0, ls=(0, (3, 3)), lw=0.6, color=DASH, zorder=1)
    ax.fill_between(x, curve["ci_low"], curve["ci_high"], color=PRIMARY, alpha=0.14, lw=0, zorder=2)
    ax.plot(x, curve["difference"], color=PRIMARY, lw=1.2, zorder=3)
    for _, r in contrasts.iterrows():
        ax.errorbar(np.log10(r.pth_pg_ml), r.estimate,
                    yerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]],
                    fmt="o", ms=2.8, mew=0, elinewidth=0.9, color=PRIMARY, zorder=4)
    ax.plot(np.log10(65), 0, "o", ms=2.8, mfc="white", mec=PRIMARY, mew=0.8, zorder=4)
    xlim = (x[0], np.log10(S2B_TICKS[-1]) + 0.01)   # the axis runs to the last tick
    ax.set_xticks(np.log10(S2B_TICKS))
    ax.set_xlim(*xlim)
    ax.set_ylabel("Difference in LV mass index (g/m\u00b2)\nvs PTH 65 pg/ml", fontsize=5.6, color=TEXT)
    ax.tick_params(axis="both", labelsize=5.4, colors=TEXT, width=0.6, length=2.2, color=GREY)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
        ax.spines[side].set_color(GREY)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_title("LV mass index across PTH value ranges", fontsize=6.3, fontweight="bold",
                 color=TITLE_COLOUR, pad=7)
    below = axes_mm(fig, size, 95.0, 127.0, 46.0, 50.5)
    below.bar(histogram["log10_pth_left"], histogram["count"],
              width=histogram["log10_pth_right"] - histogram["log10_pth_left"],
              align="edge", color=LIGHT_GREY, alpha=0.55)
    below.set_yticks([])
    below.set_xticks(np.log10(S2B_TICKS))
    below.set_xlim(*xlim)
    below.set_xticklabels([str(v) for v in S2B_TICKS], fontsize=5.4, color=TEXT)
    below.set_xlabel("PTH (pg/ml, log scale)", fontsize=5.8, color=TEXT)
    below.tick_params(axis="x", width=0.6, length=2.2, color=GREY)
    for spine in below.spines.values():
        spine.set_visible(False)
    letter(fig, size, "a", 4.0, 5.0)
    letter(fig, size, "b", 82.0, 7.0)
    write(fig, "Supplementary Figure S2")


def main():
    df = prep.load_visits()
    _, qc = prep.echo_qc(df)
    figure_1()
    figure_2(qc)
    figure_3()
    figure_4()
    figure_s1()
    figure_s2()


if __name__ == "__main__":
    main()
