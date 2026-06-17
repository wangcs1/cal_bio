"""Shared figure style for all C-part report figures.

Centralises the warm-paper palette used by ``plot_report_figures`` so that the
interpretability figures (ISM heatmaps, variant delta profiles) match the rest
of the report. Only depends on matplotlib.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

# Warm-paper palette (kept in sync with src/reports/plot_report_figures.py)
BG = "#fbf7ef"      # figure background
PANEL = "#fffdf8"   # axes background
DARK = "#202124"    # primary text / ticks
MUTED = "#6b6258"   # secondary text / notes
GRID = "#d9cfc0"    # grid lines
SPINE = "#817768"   # axis spines

# Line palette (shared with the model colours in plot_report_figures.py)
WT_COLOR = "#4C78A8"
MUT_COLOR = "#F58518"
ACCENT = "#E45756"   # markers / variant position

# Diverging colormap for in-silico-mutagenesis heatmaps
HEATMAP_CMAP = "RdBu_r"


def style_axes(ax: plt.Axes, *, frame: bool = False) -> None:
    """Apply the warm-paper look to an Axes.

    ``frame=False`` (default) hides the top/right spines for line/bar plots;
    ``frame=True`` keeps a thin four-sided frame, which reads better for
    heatmaps / image plots.
    """
    ax.set_facecolor(PANEL)
    ax.figure.set_facecolor(BG)
    if frame:
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color(SPINE)
            ax.spines[side].set_linewidth(0.8)
    else:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_color(SPINE)
        ax.spines["bottom"].set_color(SPINE)
    ax.tick_params(colors=DARK, labelsize=9)
    ax.title.set_color(DARK)
    ax.xaxis.label.set_color(DARK)
    ax.yaxis.label.set_color(DARK)


def style_colorbar(cbar) -> None:
    """Match a colorbar to the palette."""
    cbar.outline.set_edgecolor(SPINE)
    cbar.outline.set_linewidth(0.8)
    cbar.ax.tick_params(colors=DARK, labelsize=8)
    cbar.ax.yaxis.label.set_color(MUTED)


def save(fig: plt.Figure, path: Path, *, dpi: int = 220) -> None:
    """Save with the figure background preserved."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
