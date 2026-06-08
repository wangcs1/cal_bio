from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.build_synthetic_splice_dataset import build_and_write
from src.build_variant_dataset import build_and_write_variants
from src.run_exp3_variant_effect import saturation_matrix, train_models
from src.utils import (
    BASES,
    PROJECT_ROOT,
    ensure_dirs,
    read_csv,
    scan_splice_scores,
    write_dataframe,
)


def ensure_inputs() -> None:
    if not (PROJECT_ROOT / "data/processed/splice_sites_pm200.csv").exists():
        build_and_write()
    if not (PROJECT_ROOT / "data/processed/artificial_variant_effect.csv").exists():
        build_and_write_variants()


def plot_matrix(matrix: pd.DataFrame, out_path: Path, title: str) -> None:
    pivot = matrix.pivot(index="alt_base", columns="relative_position", values="importance").loc[list(BASES)]
    fig, ax = plt.subplots(figsize=(12, 3.4))
    vmax = max(0.05, float(np.nanmax(np.abs(pivot.to_numpy()))))
    image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(BASES)), labels=list(BASES))
    tick_positions = np.linspace(0, pivot.shape[1] - 1, 11, dtype=int)
    ax.set_xticks(tick_positions, labels=[str(pivot.columns[i]) for i in tick_positions])
    ax.set_xlabel("Relative position")
    ax.set_ylabel("Alternative base")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="Target probability drop")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run_ism(out_tables: Path, out_figures: Path, random_state: int) -> None:
    test = read_csv(PROJECT_ROOT / "data/splits/test_pm200.csv")
    models = train_models(random_state)
    model = next(model for model in models if model.name == "SpliceAI signal proxy")
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
    for case_name, row, target_class, filename in cases:
        matrix = saturation_matrix(model, str(row["sequence"]), target_class, flank=50)
        matrix.insert(0, "case", case_name)
        matrix.insert(1, "sample_id", row["sample_id"])
        write_dataframe(out_tables / f"{case_name}_ism_matrix.csv", matrix)
        plot_matrix(matrix, out_figures / filename, f"In silico mutagenesis: {case_name}")


def plot_delta_profile(variant: pd.Series, out_path: Path, title: str) -> pd.DataFrame:
    wt = scan_splice_scores(str(variant["wt_sequence"]), flank=80)
    mut = scan_splice_scores(str(variant["mut_sequence"]), flank=80)
    wt["sequence_type"] = "WT"
    mut["sequence_type"] = "Mut"
    combined = pd.concat([wt, mut], ignore_index=True)
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.2), sharex=True)
    for seq_type, group in combined.groupby("sequence_type"):
        axes[0].plot(group["relative_position"], group["donor_score"], label=seq_type, linewidth=2)
        axes[1].plot(group["relative_position"], group["acceptor_score"], label=seq_type, linewidth=2)
    axes[0].set_ylabel("Donor score")
    axes[1].set_ylabel("Acceptor score")
    axes[1].set_xlabel("Relative position")
    for ax in axes:
        ax.axvline(int(variant["relative_pos"]), color="#b13f3f", linestyle="--", linewidth=1.3)
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    combined.insert(0, "variant_id", variant["variant_id"])
    combined.insert(1, "variant_type", variant["variant_type"])
    return combined


def run_variant_profiles(out_tables: Path, out_figures: Path) -> None:
    variants = read_csv(PROJECT_ROOT / "data/processed/artificial_variant_effect.csv")
    donor_loss = variants[variants["variant_type"] == "donor_loss"].iloc[0]
    cryptic_gain = variants[variants["variant_type"] == "cryptic_gain"].iloc[0]
    donor_profile = plot_delta_profile(
        donor_loss,
        out_figures / "variant_delta_profile_donor_loss.png",
        "WT vs Mut delta profile: donor loss",
    )
    cryptic_profile = plot_delta_profile(
        cryptic_gain,
        out_figures / "variant_delta_profile_cryptic_gain.png",
        "WT vs Mut delta profile: cryptic gain",
    )
    write_dataframe(out_tables / "variant_delta_profile_donor_loss.csv", donor_profile)
    write_dataframe(out_tables / "variant_delta_profile_cryptic_gain.csv", cryptic_profile)


def run_junction_case_study(out_tables: Path, out_figures: Path) -> None:
    path = out_tables / "experiment_2C_tissue_splice_usage_case_study.csv"
    if path.exists():
        usage = read_csv(path)
    else:
        tissues = ["brain", "heart", "liver", "muscle", "blood"]
        usage = pd.DataFrame(
            [
                {
                    "event_id": "SYN_EVENT_ALT_EXON",
                    "tissue": tissue,
                    "pangolin_proxy_splice_usage": value,
                    "data_source": "synthetic_tissue_case_study_v1",
                }
                for tissue, value in zip(tissues, [0.76, 0.36, 0.41, 0.81, 0.62])
            ]
        )
        write_dataframe(path, usage)
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    for event, group in usage.groupby("event_id"):
        group = group.sort_values("tissue")
        ax.plot(group["tissue"], group["pangolin_proxy_splice_usage"], marker="o", linewidth=2, label=event)
    ax.set_ylabel("Predicted junction usage")
    ax.set_xlabel("Tissue")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Experiment 2C junction usage case study")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_figures / "exp2C_junction_usage_case_study.png", dpi=180)
    plt.close(fig)


def run(out_tables: Path, out_figures: Path, random_state: int = 42) -> None:
    ensure_dirs(out_tables, out_figures)
    ensure_inputs()
    run_ism(out_tables, out_figures, random_state)
    run_variant_profiles(out_tables, out_figures)
    run_junction_case_study(out_tables, out_figures)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run C-part interpretability analyses.")
    parser.add_argument("--tables", type=Path, default=PROJECT_ROOT / "results/tables")
    parser.add_argument("--figures", type=Path, default=PROJECT_ROOT / "results/figures")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(args.tables, args.figures, args.seed)
    print("Interpretability outputs written.")


if __name__ == "__main__":
    main()

