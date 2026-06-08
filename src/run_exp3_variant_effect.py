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
from src.models.simple_splice_models import make_model_suite, zero_shot_embedding_distance
from src.utils import (
    BASES,
    PROJECT_ROOT,
    acceptor_consensus_score,
    binary_ranking_metrics,
    donor_consensus_score,
    ensure_dirs,
    mutate_base,
    read_csv,
    write_dataframe,
)


def ensure_inputs() -> None:
    if not (PROJECT_ROOT / "data/processed/splice_sites_pm200.csv").exists():
        build_and_write()
    if not (PROJECT_ROOT / "data/processed/artificial_variant_effect.csv").exists():
        build_and_write_variants()


def train_models(random_state: int) -> list[object]:
    train = read_csv(PROJECT_ROOT / "data/splits/train_pm200.csv")
    valid = read_csv(PROJECT_ROOT / "data/splits/valid_pm200.csv")
    train_full = pd.concat([train, valid], ignore_index=True)
    x_train = train_full["sequence"].astype(str).tolist()
    y_train = train_full["label"].astype(int).to_numpy()
    models = make_model_suite(random_state=random_state)
    for model in models:
        model.fit(x_train, y_train)
    return models


def delta_from_proba(row: pd.Series, wt_proba: np.ndarray, mut_proba: np.ndarray) -> float:
    target = int(row["target_class"])
    variant_type = str(row["variant_type"])
    if variant_type.endswith("loss"):
        return float(wt_proba[target] - mut_proba[target])
    if "gain" in variant_type:
        return float(mut_proba[target] - wt_proba[target])
    return float(np.max(np.abs(mut_proba[:2] - wt_proba[:2])))


def maxent_proxy_delta(row: pd.Series) -> float:
    wt = str(row["wt_sequence"])
    mut = str(row["mut_sequence"])
    target = int(row["target_class"])
    variant_type = str(row["variant_type"])
    if target == 0:
        wt_score = donor_consensus_score(wt)
        mut_score = donor_consensus_score(mut)
    elif target == 1:
        wt_score = acceptor_consensus_score(wt)
        mut_score = acceptor_consensus_score(mut)
    else:
        wt_score = max(donor_consensus_score(wt), acceptor_consensus_score(wt))
        mut_score = max(donor_consensus_score(mut), acceptor_consensus_score(mut))
    if variant_type.endswith("loss"):
        return float(wt_score - mut_score)
    if "gain" in variant_type:
        return float(mut_score - wt_score)
    return float(abs(mut_score - wt_score))


def score_variants(models: list[object], variants: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    wt_sequences = variants["wt_sequence"].astype(str).tolist()
    mut_sequences = variants["mut_sequence"].astype(str).tolist()

    for model in models:
        wt_proba = model.predict_proba(wt_sequences)
        mut_proba = model.predict_proba(mut_sequences)
        for idx, variant in variants.iterrows():
            rows.append(
                {
                    "variant_id": variant["variant_id"],
                    "variant_type": variant["variant_type"],
                    "label": int(variant["label"]),
                    "label_name": variant["label_name"],
                    "target_class": int(variant["target_class"]),
                    "model": model.name,
                    "impact_score": delta_from_proba(variant, wt_proba[idx], mut_proba[idx]),
                    "wt_donor": float(wt_proba[idx, 0]),
                    "wt_acceptor": float(wt_proba[idx, 1]),
                    "wt_non_splice": float(wt_proba[idx, 2]),
                    "mut_donor": float(mut_proba[idx, 0]),
                    "mut_acceptor": float(mut_proba[idx, 1]),
                    "mut_non_splice": float(mut_proba[idx, 2]),
                }
            )

    for mode, model_name in [
        ("rnafm", "RNA-FM zero-shot embedding distance"),
        ("rnabert", "RNABERT zero-shot token distance"),
    ]:
        for _, variant in variants.iterrows():
            rows.append(
                {
                    "variant_id": variant["variant_id"],
                    "variant_type": variant["variant_type"],
                    "label": int(variant["label"]),
                    "label_name": variant["label_name"],
                    "target_class": int(variant["target_class"]),
                    "model": model_name,
                    "impact_score": zero_shot_embedding_distance(
                        str(variant["wt_sequence"]), str(variant["mut_sequence"]), mode
                    ),
                    "wt_donor": np.nan,
                    "wt_acceptor": np.nan,
                    "wt_non_splice": np.nan,
                    "mut_donor": np.nan,
                    "mut_acceptor": np.nan,
                    "mut_non_splice": np.nan,
                }
            )

    for _, variant in variants.iterrows():
        rows.append(
            {
                "variant_id": variant["variant_id"],
                "variant_type": variant["variant_type"],
                "label": int(variant["label"]),
                "label_name": variant["label_name"],
                "target_class": int(variant["target_class"]),
                "model": "MaxEntScan consensus proxy",
                "impact_score": maxent_proxy_delta(variant),
                "wt_donor": np.nan,
                "wt_acceptor": np.nan,
                "wt_non_splice": np.nan,
                "mut_donor": np.nan,
                "mut_acceptor": np.nan,
                "mut_non_splice": np.nan,
            }
        )
    return pd.DataFrame(rows)


def summarize_metrics(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, group in scores.groupby("model"):
        metrics = binary_ranking_metrics(group["label"], group["impact_score"], k_fraction=0.1)
        rows.append({"model": model, **metrics, "variants": len(group)})
    return pd.DataFrame(rows).sort_values("auprc", ascending=False)


def plot_metric_bars(metrics: pd.DataFrame, out_dir: Path) -> None:
    for metric, filename, label in [
        ("auroc", "exp3_variant_auroc.png", "AUROC"),
        ("auprc", "exp3_variant_auprc.png", "AUPRC"),
    ]:
        fig, ax = plt.subplots(figsize=(8.2, 4.8))
        order = metrics.sort_values(metric)
        ax.barh(order["model"], order[metric], color="#4f7cac")
        ax.set_xlabel(label)
        ax.set_title(f"Experiment 3 artificial variant effect {label}")
        ax.set_xlim(0.0, 1.03)
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)


def plot_delta_boxplot(scores: pd.DataFrame, out_dir: Path) -> None:
    selected_models = [
        "SpliceAI signal proxy",
        "RNA-FM frozen k-mer + MLP",
        "RNABERT frozen token + MLP",
        "MaxEntScan consensus proxy",
    ]
    subset = scores[scores["model"].isin(selected_models)].copy()
    fig, axes = plt.subplots(1, len(selected_models), figsize=(14, 4.4), sharey=True)
    for ax, model in zip(axes, selected_models):
        group = subset[subset["model"] == model]
        values = [
            group[group["label"] == 0]["impact_score"].to_numpy(),
            group[group["label"] == 1]["impact_score"].to_numpy(),
        ]
        ax.boxplot(values, labels=["neutral", "splice altering"], showfliers=False)
        ax.set_title(model, fontsize=9)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Impact score")
    fig.suptitle("Experiment 3 delta score distribution")
    fig.tight_layout()
    fig.savefig(out_dir / "exp3_delta_score_boxplot.png", dpi=180)
    plt.close(fig)


def saturation_matrix(model: object, sequence: str, target_class: int, flank: int = 50) -> pd.DataFrame:
    c = len(sequence) // 2
    original = model.predict_proba([sequence])[0, target_class]
    rows = []
    for rel in range(-flank, flank + 1):
        idx = c + rel
        for base in BASES:
            if sequence[idx] == base:
                importance = 0.0
            else:
                mutated = mutate_base(sequence, idx, base)
                mut_score = model.predict_proba([mutated])[0, target_class]
                importance = float(original - mut_score)
            rows.append({"relative_position": rel, "alt_base": base, "importance": importance})
    return pd.DataFrame(rows)


def plot_saturation(matrix: pd.DataFrame, out_path: Path, title: str) -> None:
    pivot = matrix.pivot(index="alt_base", columns="relative_position", values="importance").loc[list(BASES)]
    fig, ax = plt.subplots(figsize=(12, 3.4))
    image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_yticks(range(len(BASES)), labels=list(BASES))
    tick_positions = np.linspace(0, pivot.shape[1] - 1, 11, dtype=int)
    ax.set_xticks(tick_positions, labels=[str(pivot.columns[i]) for i in tick_positions])
    ax.set_xlabel("Relative position")
    ax.set_ylabel("Alternative base")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="P(original class) drop")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run_saturation(models: list[object], out_tables: Path, out_figures: Path) -> None:
    test = read_csv(PROJECT_ROOT / "data/splits/test_pm200.csv")
    model = next(model for model in models if model.name == "SpliceAI signal proxy")
    donor = test[test["label"].astype(int) == 0].iloc[0]
    acceptor = test[test["label"].astype(int) == 1].iloc[0]
    donor_matrix = saturation_matrix(model, str(donor["sequence"]), 0)
    acceptor_matrix = saturation_matrix(model, str(acceptor["sequence"]), 1)
    write_dataframe(out_tables / "exp3_saturation_mutagenesis_donor_matrix.csv", donor_matrix)
    write_dataframe(out_tables / "exp3_saturation_mutagenesis_acceptor_matrix.csv", acceptor_matrix)
    plot_saturation(
        donor_matrix,
        out_figures / "exp3_saturation_mutagenesis_heatmap.png",
        "Saturation mutagenesis: donor probability sensitivity",
    )
    plot_saturation(
        acceptor_matrix,
        out_figures / "exp3_saturation_mutagenesis_acceptor_heatmap.png",
        "Saturation mutagenesis: acceptor probability sensitivity",
    )


def run(output_tables: Path, output_figures: Path, random_state: int = 42) -> dict[str, pd.DataFrame]:
    ensure_dirs(output_tables, output_figures)
    ensure_inputs()
    variants = read_csv(PROJECT_ROOT / "data/processed/artificial_variant_effect.csv")
    models = train_models(random_state=random_state)
    scores = score_variants(models, variants)
    metrics = summarize_metrics(scores)
    write_dataframe(output_tables / "experiment_3A_artificial_variant_scores.csv", scores)
    write_dataframe(output_tables / "experiment_3A_artificial_variant_metrics.csv", metrics)
    plot_metric_bars(metrics, output_figures)
    plot_delta_boxplot(scores, output_figures)
    run_saturation(models, output_tables, output_figures)
    return {"scores": scores, "metrics": metrics}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run C-part experiment 3: artificial splice variant effects.")
    parser.add_argument("--tables", type=Path, default=PROJECT_ROOT / "results/tables")
    parser.add_argument("--figures", type=Path, default=PROJECT_ROOT / "results/figures")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run(args.tables, args.figures, args.seed)
    print(outputs["metrics"].to_string(index=False))


if __name__ == "__main__":
    main()

