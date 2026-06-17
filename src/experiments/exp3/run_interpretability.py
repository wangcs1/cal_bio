from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.data.build_clinvar_variant_dataset import build_clinvar_smoke as build_clinvar_format_control, build_clinvar_variants
from src.data.build_splice_site_dataset import build_and_write_real
from src.experiments.exp3.run_variant_effect import saturation_matrix, train_models
from src.reports import fig_style
from src.utils import (
    BASES,
    EXP2_FIGURES_DIR,
    EXP2_TABLES_DIR,
    EXP3_FIGURES_DIR,
    EXP3_TABLES_DIR,
    PROJECT_ROOT,
    ensure_dirs,
    read_csv,
    exp3_data_file,
    scan_splice_scores,
    shared_split_file,
    write_dataframe,
)


def ensure_inputs() -> None:
    if not shared_split_file("test_pm200.csv").exists():
        build_and_write_real()
    if not exp3_data_file("clinvar_splicing_variants.csv").exists():
        build_clinvar_variants()


def plot_matrix(matrix: pd.DataFrame, out_path: Path, title: str) -> None:
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
    ax.set_title(title, fontsize=13, fontweight="bold", loc="left")
    cbar = fig.colorbar(image, ax=ax, label="Target probability drop", fraction=0.025, pad=0.012)
    fig_style.style_colorbar(cbar)
    fig.tight_layout()
    fig_style.save(fig, out_path)


def run_ism(out_tables: Path, out_figures: Path, random_state: int) -> None:
    test = read_csv(shared_split_file("test_pm200.csv"))
    models = train_models(random_state)
    cases = [
        ("donor", test[test["label"].astype(int) == 0].iloc[0], 0, "ism_donor_heatmap.png"),
        ("acceptor", test[test["label"].astype(int) == 1].iloc[0], 1, "ism_acceptor_heatmap.png"),
        (
            "hard_negative",
            test[(test["label"].astype(int) == 2) & test["negative_type"].astype(str).str.contains("hard")].iloc[0],
            2,
            "ism_hard_negative_heatmap.png",
        ),
    ]
    for model in models:
        model_slug = (
            model.name.lower()
            .replace(" ", "_")
            .replace("+", "plus")
            .replace("(", "")
            .replace(")", "")
            .replace("-", "_")
        )[:36]
        for case_name, row, target_class, filename in cases:
            matrix = saturation_matrix(model, str(row["sequence"]), target_class, flank=30)
            matrix.insert(0, "model", model.name)
            matrix.insert(1, "case", case_name)
            matrix.insert(2, "sample_id", row["sample_id"])
            write_dataframe(out_tables / f"{model_slug}_{case_name}_ism_matrix.csv", matrix)
            if model is models[0]:
                write_dataframe(out_tables / f"{case_name}_ism_matrix.csv", matrix)
                plot_matrix(matrix, out_figures / filename, f"In silico mutagenesis: {case_name}")


def plot_delta_profile(variant: pd.Series, out_path: Path, title: str) -> pd.DataFrame:
    wt = scan_splice_scores(str(variant["wt_sequence"]), flank=80)
    mut = scan_splice_scores(str(variant["mut_sequence"]), flank=80)
    wt["sequence_type"] = "WT"
    mut["sequence_type"] = "Mut"
    combined = pd.concat([wt, mut], ignore_index=True)
    line_colors = {"WT": fig_style.WT_COLOR, "Mut": fig_style.MUT_COLOR}
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.2), sharex=True)
    for seq_type, group in combined.groupby("sequence_type"):
        color = line_colors.get(str(seq_type), fig_style.WT_COLOR)
        axes[0].plot(group["relative_position"], group["donor_score"], label=seq_type, linewidth=2, color=color)
        axes[1].plot(group["relative_position"], group["acceptor_score"], label=seq_type, linewidth=2, color=color)
    axes[0].set_ylabel("Donor score")
    axes[1].set_ylabel("Acceptor score")
    axes[1].set_xlabel("Relative position")
    for ax in axes:
        fig_style.style_axes(ax)
        ax.axvline(int(variant["relative_pos"]), color=fig_style.ACCENT, linestyle="--", linewidth=1.3, label="variant")
        ax.grid(alpha=0.30, color=fig_style.GRID)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=fig_style.DARK)
    fig.suptitle(title, color=fig_style.DARK, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig_style.save(fig, out_path)
    combined.insert(0, "variant_id", variant["variant_id"])
    combined.insert(1, "variant_type", variant["variant_type"])
    return combined


def run_variant_profiles(out_tables: Path, out_figures: Path) -> None:
    variants = read_csv(exp3_data_file("clinvar_splicing_variants.csv"))
    positives = variants[variants["label"].astype(int) == 1]
    if positives.empty:
        raise ValueError("ClinVar interpretability requires at least one splice-altering positive variant.")
    case = positives.iloc[0]
    target_name = str(case.get("target_class_name", "splice"))
    profile = plot_delta_profile(
        case,
        out_figures / "variant_delta_profile_clinvar_real_case.png",
        f"WT vs Mut delta profile: real ClinVar {target_name} case",
    )
    write_dataframe(out_tables / "variant_delta_profile_clinvar_real_case.csv", profile)


def run_clinvar_delta_profile(out_tables: Path, out_figures: Path) -> None:
    clinvar = build_clinvar_format_control()
    case = clinvar[clinvar["label"].astype(int) == 1].iloc[0]
    profile = plot_delta_profile(
        case,
        out_figures / "variant_delta_profile_clinvar_format_control_case.png",
        "WT vs Mut delta profile: ClinVar format-control case",
    )
    profile["chrom"] = case["chrom"]
    profile["pos"] = case["pos"]
    profile["ref"] = case["ref"]
    profile["alt"] = case["alt"]
    profile["ref_sequence_window"] = case["wt_sequence"]
    profile["alt_sequence_window"] = case["mut_sequence"]
    wt = profile[profile["sequence_type"] == "WT"].set_index("relative_position")
    mut = profile[profile["sequence_type"] == "Mut"].set_index("relative_position")
    delta = (mut[["donor_score", "acceptor_score"]] - wt[["donor_score", "acceptor_score"]]).abs().max(axis=1)
    max_pos = int(delta.idxmax()) if not delta.empty else 0
    profile["max_abs_delta_position"] = max_pos
    write_dataframe(out_tables / "variant_delta_profile_clinvar_format_control_case.csv", profile)


def run(out_tables: Path, out_figures: Path, random_state: int = 42) -> None:
    ensure_dirs(out_tables, out_figures)
    ensure_inputs()
    run_ism(out_tables, out_figures, random_state)
    run_variant_profiles(out_tables, out_figures)
    run_clinvar_delta_profile(out_tables, out_figures)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run C-part interpretability analyses.")
    parser.add_argument("--tables", type=Path, default=EXP3_TABLES_DIR)
    parser.add_argument("--figures", type=Path, default=EXP3_FIGURES_DIR)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(args.tables, args.figures, args.seed)
    print("Interpretability outputs written.")


if __name__ == "__main__":
    main()
