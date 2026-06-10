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

from src.data.build_clinvar_variant_dataset import build_clinvar_smoke
from src.data.build_synthetic_splice_dataset import build_and_write
from src.data.build_sqtl_variant_dataset import build_sqtl_smoke
from src.data.build_variant_dataset import build_and_write_variants
from src.models.maxentscan_wrapper import maxentscan_variant_delta
from src.models.mmsplice_wrapper import mmsplice_variant_delta
from src.models.pangolin_wrapper import pangolin_variant_delta
from src.models.simple_splice_models import make_model_suite, pseudo_likelihood_proxy_score, zero_shot_embedding_distance
from src.models.spliceai_wrapper import spliceai_variant_delta
from src.utils import (
    BASES,
    CONFIG_ROOT,
    EXP3_DATA_DIR,
    EXP3_FIGURES_DIR,
    EXP3_TABLES_DIR,
    PROJECT_ROOT,
    acceptor_consensus_score,
    binary_ranking_metrics,
    calibration_bins,
    donor_consensus_score,
    ensure_dirs,
    load_config,
    mutate_base,
    read_csv,
    exp3_data_file,
    shared_split_file,
    topk_enrichment_curve,
    write_dataframe,
)


def ensure_inputs() -> None:
    if not shared_split_file("train_pm200.csv").exists():
        build_and_write()
    if not exp3_data_file("artificial_variant_effect.csv").exists():
        build_and_write_variants()


def train_models(random_state: int) -> list[object]:
    train = read_csv(shared_split_file("train_pm200.csv"))
    valid = read_csv(shared_split_file("valid_pm200.csv"))
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
                    "ref_score": float(wt_proba[idx, int(variant["target_class"])]),
                    "alt_score": float(mut_proba[idx, int(variant["target_class"])]),
                    "delta_score": delta_from_proba(variant, wt_proba[idx], mut_proba[idx]),
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
            distance = zero_shot_embedding_distance(str(variant["wt_sequence"]), str(variant["mut_sequence"]), mode)
            rows.append(
                {
                    "variant_id": variant["variant_id"],
                    "variant_type": variant["variant_type"],
                    "label": int(variant["label"]),
                    "label_name": variant["label_name"],
                    "target_class": int(variant["target_class"]),
                    "model": model_name,
                    "ref_score": 0.0,
                    "alt_score": distance,
                    "delta_score": distance,
                    "impact_score": distance,
                    "wt_donor": np.nan,
                    "wt_acceptor": np.nan,
                    "wt_non_splice": np.nan,
                    "mut_donor": np.nan,
                    "mut_acceptor": np.nan,
                    "mut_non_splice": np.nan,
                }
            )
            ref_pll = pseudo_likelihood_proxy_score(str(variant["wt_sequence"]), mode)
            alt_pll = pseudo_likelihood_proxy_score(str(variant["mut_sequence"]), mode)
            delta = abs(alt_pll - ref_pll)
            rows.append(
                {
                    "variant_id": variant["variant_id"],
                    "variant_type": variant["variant_type"],
                    "label": int(variant["label"]),
                    "label_name": variant["label_name"],
                    "target_class": int(variant["target_class"]),
                    "model": model_name.replace("embedding distance", "pseudo-likelihood").replace("token distance", "pseudo-likelihood"),
                    "ref_score": ref_pll,
                    "alt_score": alt_pll,
                    "delta_score": delta,
                    "impact_score": delta,
                    "wt_donor": np.nan,
                    "wt_acceptor": np.nan,
                    "wt_non_splice": np.nan,
                    "mut_donor": np.nan,
                    "mut_acceptor": np.nan,
                    "mut_non_splice": np.nan,
                }
            )

    tool_fns = [
        ("MaxEntScan optional tool (proxy fallback)", maxentscan_variant_delta),
        ("MMSplice optional tool (proxy fallback)", mmsplice_variant_delta),
        ("SpliceAI optional real tool (proxy fallback)", spliceai_variant_delta),
        ("Pangolin optional tool (small case-study proxy)", pangolin_variant_delta),
    ]
    for _, variant in variants.iterrows():
        for model_name, fn in tool_fns:
            ref_score, alt_score, delta = fn(
                str(variant["wt_sequence"]),
                str(variant["mut_sequence"]),
                int(variant["target_class"]),
                str(variant["variant_type"]),
            )
            rows.append(
                {
                    "variant_id": variant["variant_id"],
                    "variant_type": variant["variant_type"],
                    "label": int(variant["label"]),
                    "label_name": variant["label_name"],
                    "target_class": int(variant["target_class"]),
                    "model": model_name,
                    "ref_score": ref_score,
                    "alt_score": alt_score,
                    "delta_score": delta,
                    "impact_score": delta,
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


def summarize_topk(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, group in scores.groupby("model"):
        curve = topk_enrichment_curve(group["label"], group["impact_score"], fractions=[0.02, 0.05, 0.10, 0.20])
        curve.insert(0, "model", model)
        rows.append(curve)
    return pd.concat(rows, ignore_index=True)


def summarize_variant_types(scores: pd.DataFrame) -> pd.DataFrame:
    top_models = [
        "RNABERT zero-shot token distance",
        "RNA-FM zero-shot embedding distance",
        "SpliceAI signal proxy",
        "MaxEntScan optional tool (proxy fallback)",
    ]
    rows = []
    for model, model_frame in scores[scores["model"].isin(top_models)].groupby("model"):
        for variant_type, group in model_frame.groupby("variant_type"):
            labels = (group["label_name"] == "splice_altering").astype(int).to_numpy()
            score = group["impact_score"].to_numpy()
            rows.append(
                {
                    "model": model,
                    "variant_type": variant_type,
                    "mean_score": float(np.mean(score)),
                    "median_score": float(np.median(score)),
                    "rows": len(group),
                    "positive_rate": float(np.mean(labels)),
                }
            )
    return pd.DataFrame(rows)


def plot_variant_type_summary(summary: pd.DataFrame, out_dir: Path) -> None:
    if summary.empty:
        return
    variant_order = summary["variant_type"].drop_duplicates().tolist()
    model_order = summary["model"].drop_duplicates().tolist()
    x = np.arange(len(variant_order))
    width = 0.8 / max(1, len(model_order))
    fig, ax = plt.subplots(figsize=(10, 4.8))
    for idx, model in enumerate(model_order):
        values = [
            float(summary[(summary["model"] == model) & (summary["variant_type"] == variant_type)]["mean_score"].iloc[0])
            if not summary[(summary["model"] == model) & (summary["variant_type"] == variant_type)].empty
            else 0.0
            for variant_type in variant_order
        ]
        ax.bar(x - 0.4 + width / 2 + idx * width, values, width, label=model)
    ax.set_xticks(x, variant_order, rotation=25, ha="right")
    ax.set_ylabel("Mean impact score")
    ax.set_title("Variant effect scores by perturbation type")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "variant_effect_stratified_by_type.png", dpi=180)
    plt.close(fig)


def plot_calibration(scores: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    selected = scores[scores["model"] == scores.groupby("model")["impact_score"].mean().idxmax()]
    if selected.empty:
        selected = scores
    model_name = str(selected["model"].iloc[0])
    bins = calibration_bins(selected["label"], selected["impact_score"], bins=8)
    bins.insert(0, "model", model_name)
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    ax.plot(bins["mean_score"], bins["positive_rate"], marker="o", linewidth=2)
    ax.set_xlabel("Mean impact score bin")
    ax.set_ylabel("Observed positive rate")
    ax.set_title(f"Calibration curve ({model_name}, n={len(selected)}, bins={len(bins)})")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "exp3_calibration_curve.png", dpi=180)
    plt.close(fig)
    return bins


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
        "MaxEntScan optional tool (proxy fallback)",
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
    test = read_csv(shared_split_file("test_pm200.csv"))
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
    variants = read_csv(exp3_data_file("artificial_variant_effect.csv"))
    models = train_models(random_state=random_state)
    scores = score_variants(models, variants)
    metrics = summarize_metrics(scores)
    topk = summarize_topk(scores)
    variant_summary = summarize_variant_types(scores)
    calibration = plot_calibration(scores, output_figures)
    write_dataframe(output_tables / "experiment_3A_artificial_variant_scores.csv", scores)
    write_dataframe(output_tables / "experiment_3A_artificial_variant_metrics.csv", metrics)
    write_dataframe(output_tables / "experiment_3A_topk_enrichment_curve.csv", topk)
    write_dataframe(output_tables / "variant_effect_stratified_by_type.csv", variant_summary)
    write_dataframe(output_tables / "experiment_3A_calibration_bins.csv", calibration)
    plot_metric_bars(metrics, output_figures)
    plot_delta_boxplot(scores, output_figures)
    plot_variant_type_summary(variant_summary, output_figures)
    run_saturation(models, output_tables, output_figures)
    clinvar = run_clinvar_smoke(output_tables, output_figures, models)
    sqtl = run_sqtl_case_study(output_tables, models)
    return {
        "scores": scores,
        "metrics": metrics,
        "topk": topk,
        "variant_summary": variant_summary,
        "clinvar": clinvar,
        "sqtl": sqtl,
    }


def run_clinvar_smoke(output_tables: Path, output_figures: Path, models: list[object]) -> pd.DataFrame:
    clinvar = build_clinvar_smoke()
    scores = score_variants(models, clinvar)
    metrics = summarize_metrics(scores)
    write_dataframe(output_tables / "experiment_3B_clinvar_smoke_scores.csv", scores)
    write_dataframe(output_tables / "experiment_3B_clinvar_smoke_metrics.csv", metrics)
    best = scores[scores["model"] == metrics.iloc[0]["model"]].copy()
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    best.boxplot(column="impact_score", by="label_name", ax=ax)
    ax.set_title("ClinVar smoke scores")
    ax.set_xlabel("Label")
    ax.set_ylabel("Impact score")
    fig.suptitle("")
    fig.tight_layout()
    fig.savefig(output_figures / "exp3_clinvar_smoke_scores.png", dpi=180)
    plt.close(fig)
    return metrics


def run_sqtl_case_study(output_tables: Path, models: list[object]) -> pd.DataFrame:
    sqtl = build_sqtl_smoke()
    scored = []
    model = next(model for model in models if model.name == "SpliceAI signal proxy")
    for _, row in sqtl.iterrows():
        wt = model.predict_proba([str(row["wt_sequence"])])[0]
        mut = model.predict_proba([str(row["mut_sequence"])])[0]
        delta = float(np.max(np.abs(mut[:2] - wt[:2])))
        scored.append({**row.to_dict(), "model_delta_score": delta})
    frame = pd.DataFrame(scored)
    write_dataframe(output_tables / "experiment_3C_sqtl_case_study.csv", frame)
    return frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run C-part experiment 3: artificial splice variant effects.")
    parser.add_argument("--config", type=Path, default=CONFIG_ROOT / "exp3_variant_effect.yaml")
    parser.add_argument("--tables", type=Path, default=EXP3_TABLES_DIR)
    parser.add_argument("--figures", type=Path, default=EXP3_FIGURES_DIR)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config) if args.config else {}
    outputs = run(
        Path(config.get("tables_dir", args.tables)),
        Path(config.get("figures_dir", args.figures)),
        int(config.get("seed", args.seed)),
    )
    print(outputs["metrics"].to_string(index=False))


if __name__ == "__main__":
    main()
