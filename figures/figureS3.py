"""Supplementary Figure S3: FGFR4 and FGF23 detection across groups of the
human heart atlas (all nuclei).

Drawn from the tables written by single_cell/03_analyse.py
(outputs/single_cell/tables, or a folder given as the first argument):
  a  by_study.csv       contributing study
  b  by_chemistry.csv   library chemistry
  c  by_region.csv      cardiac region
  d  by_disease.csv     disease state
Each panel shows the percentage of nuclei with an FGFR4 transcript (circles)
and, for comparison, an FGF23 transcript (squares; FGF23 is mainly made in
bone), groups in alphabetical order from the bottom. Written to outputs/figures/submission/ as TIFF (500 dpi), PNG, PDF
and SVG.
"""

import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator

from figure5 import TITLE_COLOUR, letter, mm_axes, write
from kipanel import TABLES, mm_figure

SIZE = (187.96, 103.988)          # mm; 3700 x 2047 pixels at 500 dpi
GENES = {"FGFR4": dict(marker="o", ms=2.9, color="#A8322D", zorder=4),
         "FGF23": dict(marker="s", ms=2.5, color="#8F9194", zorder=3)}
LABELS = {"FGFR4": "FGFR4", "FGF23": "FGF23 (negative reference)"}
TICK_TEXT, AXIS, GRID = "#3A3A3C", "#6B6C6E", "#DBDCDD"
FS_Y, FS_X, FS_XLAB, FS_TITLE, FS_LEGEND, FS_LETTER = 5.0, 5.8, 6.8, 8.0, 5.6, 12.0

# file, title, axes (left, top, right, bottom in mm), title centre and baseline (mm), letter (mm),
# x-axis tick spacing (% of nuclei)
PANELS = [
    ("by_study.csv", "By constituent study", (27.69, 14.43, 91.03, 40.30), (47.2, 11.75), ("a", 6.4, 10.25), 0.05),
    ("by_chemistry.csv", "By library chemistry", (117.45, 14.43, 180.75, 40.30), (144.5, 10.28), ("b", 102.15, 10.96), 0.05),
    ("by_region.csv", "By cardiac region", (27.69, 60.55, 91.03, 86.42), (45.25, 57.88), ("c", 6.65, 59.90), 0.1),
    ("by_disease.csv", "By disease state", (119.99, 62.59, 183.34, 88.55), (144.75, 58.49), ("d", 102.3, 59.67), 0.05),
]


def wrap(label):
    """Labels longer than 30 characters are set on two lines."""
    return textwrap.fill(label, 26) if len(label) > 30 else label


def panel(fig, table, title, box, title_at, letter_at, tick_step):
    """One panel: percentage of nuclei with FGFR4 and FGF23, by group."""
    group = table.columns[0]
    t = table[table["gene"].isin(GENES)]
    det = t.pivot(index=group, columns="gene", values="detection_pct")
    groups = sorted(det.index)
    ys = np.arange(len(groups))
    left, top, right, bottom = box
    ax = mm_axes(fig, SIZE, left, top, right - left, bottom - top)
    for gene, style in GENES.items():
        ax.plot(det.loc[groups, gene].to_numpy(float), ys, ls="none", mew=0, **style)
    xmax = float(np.nanmax(det.loc[groups, list(GENES)].to_numpy(float)))
    ax.set_xlim(right=1.1 * xmax)
    ax.xaxis.set_major_locator(MultipleLocator(tick_step))
    ax.set_yticks(ys)
    ax.set_yticklabels([wrap(g) for g in groups], fontsize=FS_Y, color=TICK_TEXT)
    ax.tick_params(axis="y", length=0, pad=3.2)
    ax.tick_params(axis="x", length=1.75, width=0.6, color=TICK_TEXT, pad=2.2,
                   labelsize=FS_X, labelcolor=TICK_TEXT)
    ax.set_xlabel("Nuclei detected (%)", fontsize=FS_XLAB, color=TICK_TEXT, labelpad=3.1)
    ax.xaxis.label.set_x(0.5 - 1.3 / (right - left))     # as submitted, 1.3 mm left of centre
    ax.grid(axis="x", color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.spines["bottom"].set_color(AXIS)
    fig.text(title_at[0] / SIZE[0], 1 - title_at[1] / SIZE[1], title, fontsize=FS_TITLE,
             fontweight="bold", color=TITLE_COLOUR, ha="center", va="baseline")
    letter(fig, SIZE, *letter_at, fontsize=FS_LETTER)


def main():
    tables = Path(sys.argv[1]) if len(sys.argv) > 1 else TABLES / "single_cell" / "tables"
    fig = mm_figure(*SIZE)
    for name, title, box, title_at, letter_at, tick_step in PANELS:
        panel(fig, pd.read_csv(tables / name), title, box, title_at, letter_at, tick_step)
    handles = [fig.axes[0].plot([], [], ls="none", mew=0, **style)[0] for style in GENES.values()]
    fig.legend(handles, [LABELS[g] for g in GENES], loc="center", ncol=2, frameon=False,
               bbox_to_anchor=(89.5 / SIZE[0], 1 - 98.45 / SIZE[1]), fontsize=FS_LEGEND,
               handlelength=1.6, handletextpad=0.25, columnspacing=1.85, borderpad=0, labelcolor=TICK_TEXT)
    write(fig, "Supplementary Figure S3")


if __name__ == "__main__":
    main()
