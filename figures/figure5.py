"""Figure 5: FGF receptor and Klotho transcripts in the human heart atlas.

Drawn from the tables written by single_cell/03_analyse.py
(outputs/single_cell/tables, or a folder given as the first argument). Dot
size is the percentage of nuclei in which the transcript was detected; colour
is mean expression, on one scale per panel.
  a  Each gene by cardiac cell type (by_celltype.csv). Nuclei without an
     assigned cell type (atlas label "cell") are not shown.
  b  FGFR1, FGFR4 and FGF23 in ventricular cardiomyocytes by library chemistry,
     all studies pooled, then each study that used two chemistries
     (vcm_by_study_chemistry.csv).
  c  Ventricular cardiomyocytes in one study (Chaffin et al.) by disease state
     (vcm_chaffin_by_disease.csv).

The atlas is public (Broad Institute Single Cell Portal, SCP3689), so the
values drawn can be obtained from it directly. Written to
outputs/figures/submission/ as TIFF (500 dpi), PNG, PDF and SVG.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image

from kipanel import HAIR, INK, LW_AXIS, OUT, SUBTLE, TABLES, mm_figure, plt

DEST = OUT / "submission"
DPI = 500
TITLE_COLOUR = "#2E3192"
SIZE = (90.07, 150.07)            # mm; 1773 x 2954 pixels at 500 dpi

GENES = ["FGFR1", "FGFR2", "FGFR3", "KLB", "KL", "FGFR4", "FGF23"]
GENES_B = ["FGFR1", "FGFR4", "FGF23"]
CHEMISTRY = "library_preparation_protocol__ontology_label"
DISEASE = "disease__ontology_label"
SHORT = {"regular ventricular cardiac myocyte": "ventricular cardiomyocyte",
         "regular atrial cardiac myocyte": "atrial cardiomyocyte",
         "vascular associated smooth muscle cell": "vascular smooth muscle",
         "endothelial cell of lymphatic vessel": "lymphatic endothelium",
         "mesothelial cell of epicardium": "epicardial mesothelium",
         "fibroblast of cardiac tissue": "cardiac fibroblast",
         "cardiac endothelial cell": "cardiac endothelium"}
DISEASES = {"normal": "Normal", "hypertrophic cardiomyopathy": "HCM",
            "dilated cardiomyopathy": "DCM"}
DOT_EC, SPINE = "#6A6A6A", "#4A4A4A"
FS_AB, FS_C, FS_TINY, FS_TITLE, FS_LETTER = 5.4, 6.0, 5.0, 7.2, 10.0
NAVY = LinearSegmentedColormap.from_list(
    "navy", ["#FBFCFE", "#DCE6F0", "#7FA5C4", "#2C5F8A", "#16324B"])


# --------------------------------------------------------------------------- layout helpers

def mm_axes(fig, size, left, top, width, height):
    """Axes at left, top (mm from the top left of the figure), width x height mm."""
    w, h = size
    return fig.add_axes([left / w, 1 - (top + height) / h, width / w, height / h])


def letter(fig, figsize, label, x, baseline, fontsize=FS_LETTER):
    """Bold panel letter, its baseline at x, baseline (mm from the top left)."""
    fig.text(x / figsize[0], 1 - baseline / figsize[1], label, fontsize=fontsize,
             fontweight="bold", va="baseline", ha="left", color="black")


def title(fig, text, centre, baseline):
    """Panel title, centred at centre mm, baseline at baseline mm from the top."""
    fig.text(centre / SIZE[0], 1 - baseline / SIZE[1], text, fontsize=FS_TITLE,
             color=TITLE_COLOUR, ha="center", va="baseline")


def write(fig, name):
    """Save SVG, PDF, PNG and TIFF (500 dpi) to outputs/figures/submission/."""
    DEST.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf"):
        fig.savefig(DEST / f"{name}.{extension}")
    fig.savefig(DEST / f"{name}.png", dpi=DPI)
    plt.close(fig)
    with Image.open(DEST / f"{name}.png") as im:
        im.convert("RGB").save(DEST / f"{name}.tif", dpi=(DPI, DPI), compression="tiff_lzw")
    print(f"{name}: written to {DEST}")


def frame(ax, genes, groups, fontsize):
    """Gene names along x, groups along y, light grid."""
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=fontsize, color=INK)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=fontsize, color=INK)
    ax.set_xlim(-0.7, len(genes) - 0.3)
    ax.set_ylim(-0.7, len(groups) - 0.3)
    ax.grid(color=HAIR, lw=0.3, alpha=0.85, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.0, width=LW_AXIS, pad=1.8, colors=INK)
    for spine in ax.spines.values():
        spine.set_linewidth(LW_AXIS)
        spine.set_color(SPINE)


def colourbar(fig, cax, norm, cmap, ticksize, ticklength=1.4, tickpad=1.8):
    """Colour scale for mean expression."""
    bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    bar.set_label("Mean expression", fontsize=FS_TINY, labelpad=2.5)
    bar.outline.set_linewidth(0.35)
    bar.outline.set_edgecolor("#8A8A8A")
    bar.ax.tick_params(width=0.35, length=ticklength, pad=tickpad, labelsize=ticksize)


def size_key(fig, ax, values, area, anchor, ncol, fontsize, **kw):
    """Legend of dot sizes (percentage of nuclei with the transcript)."""
    handles = [ax.scatter([], [], s=area(v), c="white", edgecolors=DOT_EC, linewidths=kw.pop("lw", 0.25))
               for v in values]
    legend = fig.legend(handles, [str(v) for v in values], title="Nuclei detected (%)",
                        loc="upper left", bbox_to_anchor=(anchor[0] / SIZE[0], 1 - anchor[1] / SIZE[1]),
                        ncol=ncol, frameon=False, fontsize=fontsize, **kw)
    legend.get_title().set_fontsize(FS_TINY)
    legend.get_title().set_color(SUBTLE)


# --------------------------------------------------------------------------- panels

def panel_a(fig, tables):
    """Dot plot: genes by cardiac cell type, ventricular then atrial cells first."""
    t = pd.read_csv(tables / "by_celltype.csv")
    t = t[t["gene"].isin(GENES) & (t["cell_type"] != "cell")]
    wide = {v: t.pivot(index="cell_type", columns="gene", values=v)
            for v in ("detection_pct", "mean_expr")}
    order = sorted(wide["detection_pct"].index,
                   key=lambda c: (0 if "ventricular" in c else 1 if "atrial" in c else 2, c))
    det, expr = wide["detection_pct"].loc[order, GENES], wide["mean_expr"].loc[order, GENES]
    ax = mm_axes(fig, SIZE, 27.72, 8.80, 45.80, 32.0)
    norm, cmap = Normalize(0, float(np.nanmax(expr.values))), plt.get_cmap("Blues")
    for j, g in enumerate(GENES):
        ax.scatter(np.full(len(order), j), np.arange(len(order)),
                   s=np.clip(det[g].to_numpy(float), 0, None) / 100.0 * 42.0,
                   c=[cmap(norm(v)) for v in expr[g].to_numpy(float)],
                   edgecolors=DOT_EC, linewidths=0.25, zorder=3)
    frame(ax, GENES, [SHORT.get(c, c) for c in order], FS_AB)
    colourbar(fig, mm_axes(fig, SIZE, 75.32, 16.80, 1.6, 16.0), norm, cmap, FS_TINY - 0.6)
    size_key(fig, ax, [1, 10, 25, 50], lambda v: v / 100.0 * 42.0, (63.9, 48.69), 4, FS_TINY,
             handlelength=1.0, handletextpad=0.4, columnspacing=1.2, labelspacing=0.3, borderpad=0.0)
    title(fig, "Receptor detection across cardiac cell types", 52.66, 7.42)
    letter(fig, SIZE, "a", 4.95, 7.78)


def panel_b(fig, tables):
    """Dot plot: ventricular cardiomyocytes by library chemistry."""
    t = pd.read_csv(tables / "vcm_by_study_chemistry.csv")
    t = t[t["gene"].isin(GENES_B)].copy()
    t["sum_expr"] = t["mean_expr"] * t["n_nuclei"]

    def pool(frame_, keys):
        """Combine groups: total nuclei, nuclei detected and mean expression."""
        g = frame_.groupby(keys, as_index=False)[["n_nuclei", "n_detected", "sum_expr"]].sum()
        g["detection_pct"] = 100.0 * g["n_detected"] / g["n_nuclei"]
        g["mean_expr"] = g["sum_expr"] / g["n_nuclei"]
        return g

    pooled = pool(t, [CHEMISTRY, "gene"]).rename(columns={CHEMISTRY: "row"})
    pooled["block"] = "chemistry"
    frames = [pooled]
    chemistries = t.groupby("study")[CHEMISTRY].nunique()
    spanning = sorted(chemistries[chemistries > 1].index)
    if spanning:
        w = pool(t[t["study"].isin(spanning)], ["study", CHEMISTRY, "gene"])
        w["row"] = w["study"] + ", " + w[CHEMISTRY]
        w["block"] = "within-study"
        frames.append(w)
    rows = pd.concat(frames, ignore_index=True)
    wide = {v: rows.pivot_table(index="row", columns="gene", values=v)
            for v in ("detection_pct", "mean_expr")}
    block = rows.groupby("row")["block"].first()
    order = sorted(block.index[block == "chemistry"]) + sorted(block.index[block == "within-study"])
    det, expr = wide["detection_pct"].loc[order, GENES_B], wide["mean_expr"].loc[order, GENES_B]
    ax = mm_axes(fig, SIZE, 36.96, 60.05, 28.05, 26.0)
    norm, cmap = Normalize(0, float(np.nanmax(expr.values))), plt.get_cmap("Blues")
    for j, g in enumerate(GENES_B):
        ax.scatter(np.full(len(order), j), np.arange(len(order)),
                   s=np.clip(det[g].to_numpy(float), 0, None) / 100.0 * 54.0,
                   c=[cmap(norm(v)) for v in expr[g].to_numpy(float)],
                   edgecolors=DOT_EC, linewidths=0.25, zorder=3)
    frame(ax, GENES_B, order, FS_AB)
    colourbar(fig, mm_axes(fig, SIZE, 66.81, 67.08, 1.6, 14.0), norm, cmap, FS_TINY - 0.6)
    size_key(fig, ax, [1, 10, 25, 50], lambda v: v / 100.0 * 54.0, (61.36, 87.57), 4, FS_TINY,
             handlelength=1.0, handletextpad=0.4, columnspacing=1.2, labelspacing=0.3, borderpad=0.0)
    title(fig, "Across library chemistries", 51.44, 58.71)
    letter(fig, SIZE, "b", 5.1, 59.25)


def panel_c(fig, tables):
    """Dot plot: ventricular cardiomyocytes in Chaffin et al. by disease state."""
    t = pd.read_csv(tables / "vcm_chaffin_by_disease.csv")
    t["disease"] = t[DISEASE].str.lower().str.strip().map(DISEASES)
    t = t.dropna(subset=["disease"])
    t = t[t["gene"].isin(GENES)]
    rows = ["Normal", "HCM", "DCM"]
    wide = {v: t.pivot_table(index="disease", columns="gene", values=v).loc[rows, GENES]
            for v in ("detection_pct", "mean_expr")}
    det, expr = wide["detection_pct"], wide["mean_expr"]
    norm = Normalize(0, float(np.nanmax(expr.values)))
    dmax = float(np.nanmax(det.values))
    area = lambda v: np.sqrt(np.asarray(v, float) / dmax) * 46 + 1.0  # noqa: E731
    ax = mm_axes(fig, SIZE, 23.09, 104.17, 48.62, 27.94)
    for j, g in enumerate(GENES):
        ax.scatter(np.full(len(rows), j), np.arange(len(rows)), s=area(det[g].to_numpy(float)),
                   c=[NAVY(norm(v)) for v in expr[g].to_numpy(float)],
                   edgecolors=DOT_EC, linewidths=0.28, zorder=3)
    frame(ax, GENES, rows, FS_C)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    colourbar(fig, mm_axes(fig, SIZE, 73.53, 109.57, 1.559, 17.15), norm, NAVY, 4.8, ticklength=1.2, tickpad=3.4)
    size_key(fig, ax, [1, 10, 30], area, (25.2, 140.88), 3, 5.6, lw=0.28,
             handletextpad=0.5, columnspacing=1.5, borderpad=0.0)
    title(fig, "Single study", 47.16, 102.5)
    letter(fig, SIZE, "c", 6.55, 103.67)


def main():
    tables = Path(sys.argv[1]) if len(sys.argv) > 1 else TABLES / "single_cell" / "tables"
    with plt.rc_context({"font.size": FS_AB, "text.color": INK, "axes.labelcolor": INK}):
        fig = mm_figure(*SIZE)
        panel_a(fig, tables)
        panel_b(fig, tables)
        panel_c(fig, tables)
        write(fig, "Figure 5")


if __name__ == "__main__":
    main()
