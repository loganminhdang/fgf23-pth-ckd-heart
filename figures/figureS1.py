"""Supplementary Figure S1: FGF23 assay batches and study visits (panels).

  a  Number of FGF23 measurements by assay batch and study visit
     (outputs/07_batch_by_visit.csv).
  b  FGF23 by assay batch, all 865 measurements, including values reported
     below the limit of quantification (outputs/07_batch_distribution.csv).

submission_figures.py stacks the two panels into the submitted figure.
"""

import pandas as pd

from style import GREY, GRID, TABLES, TEXT, forest_axis, plt, save

WIDTH = 6.35   # inches


def panel_a(counts, labels):
    """Heat map of measurements per batch and visit, with the counts written in."""
    M = counts.to_numpy(float)
    fig = plt.figure(figsize=(WIDTH, 3.6))
    ax = fig.add_axes([0.27, 0.13, 0.70, 0.76])
    ax.imshow(M, cmap="Blues", aspect="auto", vmin=0, vmax=M.max(), zorder=2)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j]:
                ax.text(j, i, str(int(M[i, j])), ha="center", va="center", fontsize=8.5,
                        color="white" if M[i, j] > 0.55 * M.max() else "#2F3B47", zorder=3)
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([f"Visit {v}" for v in counts.columns], fontsize=9.5, color=TEXT)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(labels, fontsize=8.5, color=TEXT)
    ax.set_xlabel("Study visit", fontsize=10, color=TEXT)
    ax.set_ylabel("Assay batch", fontsize=10, color=TEXT)
    ax.set_title("When samples were assayed", fontsize=10.5, color=TEXT, pad=6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    save(fig, "FigureS1a")


def panel_b(dist):
    """Box plots of log10 FGF23 by batch (whiskers at 1.5 x IQR, no outliers)."""
    boxes = [{"med": r.log10_med, "q1": r.log10_q1, "q3": r.log10_q3,
              "whislo": r.log10_whislo, "whishi": r.log10_whishi, "fliers": []}
             for r in dist.itertuples()]
    fig = plt.figure(figsize=(WIDTH, 3.0))
    ax = fig.add_axes([0.14, 0.36, 0.83, 0.54])
    ax.yaxis.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    bp = ax.bxp(boxes, positions=range(len(boxes)), widths=0.6, patch_artist=True,
                showfliers=False, zorder=3)
    for box in bp["boxes"]:
        box.set(facecolor="#AECBE3", edgecolor="#2F3B47", linewidth=0.8)
    for part in ("whiskers", "caps"):
        for line in bp[part]:
            line.set(color="#2F3B47", linewidth=0.8)
    for line in bp["medians"]:
        line.set(color="#16191C", linewidth=1.3)
    ax.set_xticks(range(len(boxes)))
    ax.set_xticklabels(dist["batch"], fontsize=7.5, color=TEXT, rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(["1", "10", "100", "1,000", "10,000", "100,000"], fontsize=9, color=TEXT)
    ax.set_ylim(-0.15, 5.1)
    ax.set_xlabel("Assay batch", fontsize=10, color=TEXT)
    ax.set_ylabel("FGF23 (pg/ml, log scale)", fontsize=10, color=TEXT)
    ax.set_title("FGF23 by assay batch", fontsize=10.5, color=TEXT, pad=6)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=3.0, width=0.8, color=GREY, labelsize=9, colors=TEXT)
    forest_axis(ax)
    save(fig, "FigureS1b")


def main():
    by_visit = pd.read_csv(TABLES / "07_batch_by_visit.csv", index_col=0)
    counts = by_visit[[str(v) for v in range(1, 6)]]
    counts.columns = range(1, 6)
    panel_a(counts, by_visit["batch"].tolist())
    panel_b(pd.read_csv(TABLES / "07_batch_distribution.csv"))


if __name__ == "__main__":
    main()
