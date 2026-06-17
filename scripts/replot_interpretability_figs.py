"""Re-render the interpretability figures (ISM heatmaps, variant delta profile)
from the already-computed result CSVs, using the shared report style so they
match the rest of the figures. No model inference required.

Run with an interpreter that has matplotlib + pandas + numpy, e.g.::

    PYTHONPATH=. /home/hunter/anaconda3/bin/python scripts/replot_interpretability_figs.py
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.reports import fig_style

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "experiment_3" / "tables"
OUT_DIRS = [ROOT / "report_letax" / "images", ROOT / "results" / "experiment_3" / "figures"]
BASES = "ACGT"


def _save_all(fig: plt.Figure, filename: str) -> None:
    # save to every target dir (re-uses the same rendered figure)
    for d in OUT_DIRS[:-1]:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / filename, dpi=220, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig_style.save(fig, OUT_DIRS[-1] / filename)


def plot_ism(matrix_csv: Path, filename: str, case_name: str) -> None:
    matrix = pd.read_csv(matrix_csv)
    pivot = matrix.pivot(index="alt_base", columns="relative_position", values="importance").loc[list(BASES)]
    fig, ax = plt.subplots(figsize=(12, 3.4))
    fig_style.style_axes(ax, frame=True)
    vmax = max(0.05, float(np.nanmax(np.abs(pivot.to_numpy()))))
    image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap=fig_style.HEATMAP_CMAP, vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(BASES)), labels=list(BASES))
    tick_positions = np.linspace(0, pivot.shape[1] - 1, 11, dtype=int)
    ax.set_xticks(tick_positions, labels=[str(pivot.columns[i]) for i in tick_positions])
    ax.set_xlabel("Relative position")
    ax.set_ylabel("Alternative base")
    ax.set_title(f"In silico mutagenesis: {case_name}", fontsize=13, fontweight="bold", loc="left")
    cbar = fig.colorbar(image, ax=ax, label="Target probability drop", fraction=0.025, pad=0.012)
    fig_style.style_colorbar(cbar)
    fig.tight_layout()
    _save_all(fig, filename)


def _variant_position(profile: pd.DataFrame) -> int:
    """Look up the variant's relative position; fall back to the largest WT/Mut
    score divergence, then the window centre (0)."""
    vid = str(profile["variant_id"].iloc[0])
    clinvar = ROOT / "data" / "experiment_3" / "clinvar_splicing_variants.csv"
    if clinvar.exists():
        cv = pd.read_csv(clinvar)
        hit = cv[cv["variant_id"].astype(str) == vid]
        if not hit.empty and "relative_pos" in hit.columns:
            try:
                return int(hit["relative_pos"].iloc[0])
            except (TypeError, ValueError):
                pass
    wt = profile[profile["sequence_type"] == "WT"].set_index("relative_position")
    mut = profile[profile["sequence_type"] == "Mut"].set_index("relative_position")
    common = wt.index.intersection(mut.index)
    if len(common):
        delta = (mut.loc[common, ["donor_score", "acceptor_score"]]
                 - wt.loc[common, ["donor_score", "acceptor_score"]]).abs().max(axis=1)
        return int(delta.idxmax())
    return 0


def plot_delta_profile(profile_csv: Path, filename: str, title: str) -> None:
    profile = pd.read_csv(profile_csv)
    pos = _variant_position(profile)
    line_colors = {"WT": fig_style.WT_COLOR, "Mut": fig_style.MUT_COLOR}
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.2), sharex=True)
    for seq_type, group in profile.groupby("sequence_type"):
        color = line_colors.get(str(seq_type), fig_style.WT_COLOR)
        group = group.sort_values("relative_position")
        axes[0].plot(group["relative_position"], group["donor_score"], label=seq_type, linewidth=2, color=color)
        axes[1].plot(group["relative_position"], group["acceptor_score"], label=seq_type, linewidth=2, color=color)
    axes[0].set_ylabel("Donor score")
    axes[1].set_ylabel("Acceptor score")
    axes[1].set_xlabel("Relative position")
    for ax in axes:
        fig_style.style_axes(ax)
        ax.axvline(pos, color=fig_style.ACCENT, linestyle="--", linewidth=1.3, label="variant")
        ax.grid(alpha=0.30, color=fig_style.GRID)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=fig_style.DARK)
    fig.suptitle(title, color=fig_style.DARK, fontsize=13, fontweight="bold")
    fig.tight_layout()
    _save_all(fig, filename)


def main() -> None:
    plot_ism(TABLES / "donor_ism_matrix.csv", "ism_donor_heatmap.png", "donor")
    plot_ism(TABLES / "acceptor_ism_matrix.csv", "ism_acceptor_heatmap.png", "acceptor")
    plot_ism(TABLES / "hard_negative_ism_matrix.csv", "ism_hard_negative_heatmap.png", "hard_negative")
    plot_delta_profile(
        TABLES / "variant_delta_profile_clinvar_real_case.csv",
        "variant_delta_profile_clinvar_real_case.png",
        "WT vs Mut delta profile: real ClinVar case",
    )
    print("Re-plotted interpretability figures to:")
    for d in OUT_DIRS:
        print(f"  {d}")


if __name__ == "__main__":
    main()
