"""Shared settings for the Supplementary Figure S1 panels: output folders,
colours and fonts."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs"
OUT = ROOT / "outputs" / "figures"

GREY = "#4A4A4A"        # markers, bars, spines
LIGHT_GREY = "#8A8F94"  # secondary series
BLUE = "#2C5F8A"        # PTH
GRID = "#D8D8D8"
DASH = "#7A7A7A"        # reference lines
TEXT = "#111111"
DPI = 390

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "svg.fonttype": "none",
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "font.size": 10, "axes.unicode_minus": True,
})


def save(fig, name):
    """Write the figure as SVG, PDF and PNG."""
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf"):
        fig.savefig(OUT / f"{name}.{extension}")
    fig.savefig(OUT / f"{name}.png", dpi=DPI)
    plt.close(fig)
    print(f"{name}: written to {OUT}")


def forest_axis(ax, spine=GREY):
    """Remove top, right and left spines and style the x axis."""
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.spines["bottom"].set_color(spine)
