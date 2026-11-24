"""Shared figure settings for Figures 1 to 5: sizes in mm, 6 to 7 pt type,
text kept editable in SVG and PDF, and helpers for forest plots."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from style import OUT, TABLES  # noqa: E402,F401

FS_TICK, FS_LAB, FS_TITLE, FS_ANNOT, FS_TINY = 6.0, 6.8, 7.2, 5.6, 5.0
LW_AXIS, LW_DATA, LW_REF = 0.45, 1.15, 0.55
INK, SUBTLE, HAIR = "#1A1A1A", "#7A7A7A", "#D4D4D4"
ACCENT, PRIMARY, NEUTRAL = "#A8322D", "#2C5F8A", "#5A5A5A"

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS_TICK, "axes.labelsize": FS_LAB, "axes.titlesize": FS_TITLE,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK, "legend.fontsize": FS_ANNOT,
    "axes.linewidth": LW_AXIS, "xtick.major.width": LW_AXIS, "ytick.major.width": LW_AXIS,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "xtick.major.pad": 1.8, "ytick.major.pad": 1.8, "axes.labelpad": 2.2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.bbox": None, "savefig.pad_inches": 0.0,
}
matplotlib.rcParams.update(RC)


def mm_figure(width, height):
    """A figure of width x height mm."""
    return plt.figure(figsize=(width / 25.4, height / 25.4))


def save(fig, name, png=True):
    """Write SVG and PDF, and a 600 dpi PNG."""
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf"):
        fig.savefig(OUT / f"{name}.{extension}")
    if png:
        fig.savefig(OUT / f"{name}.png", dpi=600)
    plt.close(fig)
    print(f"{name}: written to {OUT}")


def frame(ax, xlabel, title, ylabel=None, grid=None):
    """Axis labels, title, ticks and axis lines."""
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", fontsize=FS_TITLE, color=INK, pad=4)
    if grid:
        ax.grid(axis=grid, color=HAIR, lw=0.35, alpha=0.9)
    ax.tick_params(length=2.0, width=LW_AXIS, pad=1.8, colors=INK)
    for spine in ax.spines.values():
        spine.set_linewidth(LW_AXIS)
        spine.set_color("#4A4A4A")


def forest(ax, rows, xlabel, title, sep_after=(), highlight=()):
    """Forest plot. rows: (label, estimate, low, high, colour), top to bottom."""
    ys = np.arange(len(rows))[::-1]
    ax.axvline(0.0, color=SUBTLE, lw=LW_REF, ls=(0, (2.5, 2)), zorder=1)
    for y, (_, est, lo, hi, colour) in zip(ys, rows):
        ax.plot([lo, hi], [y, y], color=colour, lw=LW_DATA, solid_capstyle="butt", zorder=2)
        for x in (lo, hi):
            ax.plot([x, x], [y - 0.13, y + 0.13], color=colour, lw=LW_AXIS + 0.15, zorder=2)
        ax.plot([est], [y], "o", color=colour, ms=3.1, mec="white", mew=0.55, zorder=3)
    for i in sep_after:
        ax.axhline(ys[i] - 0.5, color=HAIR, lw=0.4, zorder=0)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    for label, tick in zip([r[0] for r in rows], ax.get_yticklabels()):
        if label in highlight:
            tick.set_color(ACCENT)
            tick.set_fontweight("bold")
    ax.grid(axis="x", color=HAIR, lw=0.35, alpha=0.9)
    ax.tick_params(axis="y", length=0, pad=2.5, labelsize=FS_ANNOT)
    ax.spines["left"].set_visible(False)
    frame(ax, xlabel, title)


def letters(fig, axes, labels, dx, dy):
    """Bold panel letters above the top-left corner of each panel."""
    for ax, label in zip(axes, labels):
        box = ax.get_position()
        fig.text(max(box.x0 - dx, 0.003), box.y1 + dy, label, fontsize=8.2,
                 fontweight="bold", va="bottom", ha="left", color=INK)
