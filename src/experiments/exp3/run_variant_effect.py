from __future__ import annotations

import argparse
import os
import subprocess
import sys
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
from src.data.build_sqtl_variant_dataset import build_sqtl_smoke as build_sqtl_format_control
from src.experiments.exp1.common import make_model
from src.utils import (
    BASES,
    CONFIG_ROOT,
    EXP3_DATA_DIR,
    EXP3_FIGURES_DIR,
    EXP3_TABLES_DIR,
    PROJECT_ROOT,
    binary_ranking_metrics,
    calibration_bins,
    ensure_dirs,
    load_config,
    mutate_base,
    read_csv,
    shared_split_file,
    topk_enrichment_curve,
    write_dataframe,
)


REAL_MODEL_KEYS = ["cnn", "rnafm", "rnabert"]
EXTERNAL_TOOL_ENV_PYTHON = Path(os.environ.get("CALBIO_SPLICE_TOOLS_PYTHON", r"C:\Users\Wangcs\miniconda3\envs\calbio-splice-tools\python.exe"))
PANGOLIN_ENV_PYTHON = Path(os.environ.get("CALBIO_PANGOLIN_PYTHON", sys.executable))


def markdown_table(frame: pd.DataFrame, columns: list[str] | None = None, max_rows: int = 24) -> str:
    if frame.empty:
        return "_Not generated yet._"
    subset = frame.copy()
    if columns is not None:
        subset = subset[[column for column in columns if column in subset.columns]]
    subset = subset.head(max_rows)
    for column in subset.columns:
        if pd.api.types.is_float_dtype(subset[column]):
            subset[column] = subset[column].map(lambda value: "" if pd.isna(value) else f"{value:.4f}")
    return subset.to_markdown(index=False)


def write_report(
    path: Path,
    variant_path: Path,
    variant_counts: pd.DataFrame,
    metrics: pd.DataFrame,
    matched_metrics: pd.DataFrame,
    variant_summary: pd.DataFrame,
    clinvar: pd.DataFrame,
    sqtl: pd.DataFrame,
) -> None:
    ensure_dirs(path.parent)
    lines = [
        "# Experiment 3: Aberrant Splicing Variant Effects",
        "",
        f"Variant table: `{variant_path}`",
        "",
        "This run scores real ClinVar-labeled SNVs with real local models only. It includes the trained project models (CNN, RNA-FM frozen encoder, RNABERT frozen encoder) and external real splice tools (SpliceAI, Pangolin, MMSplice, MaxEntScan).",
        "External tools are executed through a Python 3.10 splice-tool environment and merged as real sequence-level variant scores.",
        "",
        "## Variant Set",
        "",
        markdown_table(variant_counts, ["variant_type", "label_name", "rows"]),
        "",
        "## 3A Real ClinVar Variant Ranking",
        "",
        markdown_table(metrics, ["model", "source", "auroc", "auprc", "top_k", "top_k_recall", "enrichment_at_k", "variants"]),
        "",
        "## 3B Distance-Matched ClinVar Ranking",
        "",
        "This diagnostic evaluates the same scores on an exact-distance-matched ClinVar subset to reduce the `closer to splice site = more likely pathogenic` shortcut.",
        "",
        markdown_table(matched_metrics, ["model", "source", "auroc", "auprc", "top_k", "top_k_recall", "enrichment_at_k", "variants"]),
        "",
        "## By Variant Type",
        "",
        markdown_table(variant_summary, ["model", "source", "variant_type", "mean_score", "median_score", "rows", "positive_rate"]),
        "",
        "## Format-Control Inputs",
        "",
        "The additional ClinVar subset and sQTL-style table are format-control checks, not main benchmark evidence.",
        "",
        "ClinVar subset metrics:",
        "",
        markdown_table(clinvar, ["model", "auroc", "auprc", "variants"]),
        "",
        f"sQTL-style rows: {len(sqtl)}",
        "",
        "Outputs:",
        "",
        "- `results/experiment_3/tables/experiment_3A_variant_scores.csv`",
        "- `results/experiment_3/tables/experiment_3A_variant_metrics.csv`",
        "- `results/experiment_3/tables/variant_effect_stratified_by_type.csv`",
        "- `results/experiment_3/tables/experiment_3B_format_control_metrics.csv`",
        "- `results/experiment_3/tables/experiment_3C_sqtl_case_study.csv`",
        "- `results/experiment_3/figures/exp3_variant_auroc.png`",
        "- `results/experiment_3/figures/exp3_variant_auprc.png`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def ensure_inputs(variant_path: Path) -> None:
    if not shared_split_file("train_pm200.csv").exists():
        build_and_write_real()
    if not variant_path.exists():
        build_clinvar_variants(output_path=variant_path)


def train_models(random_state: int) -> list[object]:
    train = read_csv(shared_split_file("train_pm200.csv"))
    valid = read_csv(shared_split_file("valid_pm200.csv"))
    train_full = pd.concat([train, valid], ignore_index=True)
    x_train = train_full["sequence"].astype(str).tolist()
    y_train = train_full["label"].astype(int).to_numpy()
    models = [make_model(model_key, random_state) for model_key in REAL_MODEL_KEYS]
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
                    "source": "trained_classifier",
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

    return pd.DataFrame(rows)


def run_external_tool_scores(variants_path: Path) -> pd.DataFrame:
    out_dir = EXP3_DATA_DIR / "external_tools"
    out_path = out_dir / "external_splice_tools_variant_scores.csv"
    python_exe = EXTERNAL_TOOL_ENV_PYTHON
    if not python_exe.exists():
        raise FileNotFoundError(
            f"External splice tool Python not found: {python_exe}. "
            "Set CALBIO_SPLICE_TOOLS_PYTHON to a Python environment with spliceai/mmsplice/maxentpy."
        )
    subprocess.run(
        [
            str(python_exe),
            str(PROJECT_ROOT / "scripts/run_real_external_tools.py"),
            "--variants",
            str(variants_path),
            "--out",
            str(out_dir),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    pangolin_out = out_dir / "pangolin_sequence_variant_scores.csv"
    subprocess.run(
        [
            str(PANGOLIN_ENV_PYTHON),
            str(PROJECT_ROOT / "scripts/run_real_pangolin_sequence.py"),
            "--variants",
            str(variants_path),
            "--out",
            str(out_dir),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    scores = pd.concat([read_csv(out_path), read_csv(pangolin_out)], ignore_index=True, sort=False)
    scores = scores.rename(
        columns={
            "ref_donor": "wt_donor",
            "ref_acceptor": "wt_acceptor",
            "ref_non_splice": "wt_non_splice",
            "alt_donor": "mut_donor",
            "alt_acceptor": "mut_acceptor",
            "alt_non_splice": "mut_non_splice",
        }
    )
    scores["delta_score"] = scores["impact_score"]
    target_cols = {0: "donor", 1: "acceptor", 2: "non_splice"}
    ref_scores = []
    alt_scores = []
    for _, row in scores.iterrows():
        suffix = target_cols[int(row["target_class"])]
        ref_scores.append(float(row[f"wt_{suffix}"]))
        alt_scores.append(float(row[f"mut_{suffix}"]))
    scores["ref_score"] = ref_scores
    scores["alt_score"] = alt_scores
    return scores


def summarize_metrics(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, group in scores.groupby("model"):
        metrics = binary_ranking_metrics(group["label"], group["impact_score"], k_fraction=0.1)
        source = str(group["source"].iloc[0]) if "source" in group.columns else "unknown"
        rows.append({"model": model, "source": source, **metrics, "variants": len(group)})
    return pd.DataFrame(rows).sort_values("auprc", ascending=False)


def summarize_distance_matched_metrics(scores: pd.DataFrame) -> pd.DataFrame:
    matched_path = EXP3_DATA_DIR / "clinvar_splicing_variants_distance_matched.csv"
    if not matched_path.exists():
        return pd.DataFrame()
    matched = read_csv(matched_path)
    keep_ids = set(matched["variant_id"].astype(str))
    subset = scores[scores["variant_id"].astype(str).isin(keep_ids)].copy()
    if subset.empty:
        return pd.DataFrame()
    return summarize_metrics(subset)


def summarize_topk(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, group in scores.groupby("model"):
        curve = topk_enrichment_curve(group["label"], group["impact_score"], fractions=[0.02, 0.05, 0.10, 0.20])
        curve.insert(0, "model", model)
        rows.append(curve)
    return pd.concat(rows, ignore_index=True)


def summarize_variant_types(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, model_frame in scores.groupby("model"):
        for variant_type, group in model_frame.groupby("variant_type"):
            labels = (group["label_name"] == "splice_altering").astype(int).to_numpy()
            score = group["impact_score"].to_numpy()
            rows.append(
                {
                    "model": model,
                    "source": str(model_frame["source"].iloc[0]) if "source" in model_frame.columns else "unknown",
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
    preferred_variants = ["donor_loss", "acceptor_loss", "donor_gain", "acceptor_gain", "neutral_far_snv"]
    available_variants = summary["variant_type"].drop_duplicates().tolist()
    variant_order = [variant for variant in preferred_variants if variant in available_variants]
    variant_order.extend([variant for variant in available_variants if variant not in variant_order])
    preferred_models = [
        "CNN baseline (PyTorch Conv1D)",
        "RNA-FM frozen encoder + MLP",
        "RNABERT frozen encoder + MLP",
        "SpliceAI real sequence model",
        "Pangolin real sequence model",
        "MMSplice real sequence model",
        "MaxEntScan real local score",
    ]
    available_models = summary["model"].drop_duplicates().tolist()
    model_order = [model for model in preferred_models if model in available_models]
    model_order.extend([model for model in available_models if model not in model_order])
    colors = {
        "CNN baseline (PyTorch Conv1D)": "#4C78A8",
        "RNA-FM frozen encoder + MLP": "#F58518",
        "RNABERT frozen encoder + MLP": "#54A24B",
        "SpliceAI real sequence model": "#B279A2",
        "Pangolin real sequence model": "#72B7B2",
        "MMSplice real sequence model": "#E45756",
        "MaxEntScan real local score": "#9D755D",
    }
    short_names = {
        "CNN baseline (PyTorch Conv1D)": "CNN baseline",
        "RNA-FM frozen encoder + MLP": "RNA-FM frozen",
        "RNABERT frozen encoder + MLP": "RNABERT frozen",
        "SpliceAI real sequence model": "SpliceAI real",
        "Pangolin real sequence model": "Pangolin real",
        "MMSplice real sequence model": "MMSplice real",
        "MaxEntScan real local score": "MaxEntScan real",
    }
    variant_labels = {
        "donor_loss": "Donor loss",
        "acceptor_loss": "Acceptor loss",
        "donor_gain": "Donor gain",
        "acceptor_gain": "Acceptor gain",
        "neutral_far_snv": "Neutral far SNV",
        "donor_clinvar_splice": "ClinVar donor splice",
        "acceptor_clinvar_splice": "ClinVar acceptor splice",
        "clinvar_benign_snv": "ClinVar benign SNV",
    }
    pivot = (
        summary.pivot_table(index="variant_type", columns="model", values="mean_score", aggfunc="mean")
        .reindex(variant_order)
        .fillna(0.0)
    )
    x = np.arange(len(variant_order))
    width = 0.8 / max(1, len(model_order))
    fig, (ax_bar, ax_line) = plt.subplots(
        2,
        1,
        figsize=(12.2, 8.1),
        gridspec_kw={"height_ratios": [2.25, 1.0], "hspace": 0.34},
        sharex=True,
    )
    fig.patch.set_facecolor("#fbf7ef")
    for ax in (ax_bar, ax_line):
        ax.set_facecolor("#fffdf8")
        ax.grid(axis="y", color="#d9cfc0", alpha=0.58, linewidth=0.8)
        ax.grid(axis="x", visible=False)
        for side in ["top", "right"]:
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_color("#817768")
        ax.spines["bottom"].set_color("#817768")
        ax.tick_params(colors="#202124", labelsize=9)

    max_value = float(np.nanmax(pivot[model_order].to_numpy())) if model_order else 0.0
    min_value = float(np.nanmin(pivot[model_order].to_numpy())) if model_order else 0.0
    label_threshold = max(0.02, max_value * 0.045)
    for idx, model in enumerate(model_order):
        values = pivot[model].to_numpy(dtype=float)
        positions = x - 0.4 + width / 2 + idx * width
        bars = ax_bar.bar(
            positions,
            values,
            width,
            label=short_names.get(model, model),
            color=colors.get(model, "#4C78A8"),
            edgecolor="#2f2a24",
            linewidth=0.45,
            alpha=0.94,
        )
        for bar, value in zip(bars, values):
            if abs(value) < label_threshold:
                continue
            va = "bottom" if value >= 0 else "top"
            y_offset = 0.012 if value >= 0 else -0.012
            ax_bar.text(
                bar.get_x() + bar.get_width() / 2,
                value + y_offset,
                f"{value:.2f}",
                ha="center",
                va=va,
                fontsize=7,
                color="#2f2a24",
                rotation=90,
            )

        positive_response = np.maximum(values, 0.0)
        scale = float(np.max(positive_response))
        normalized = positive_response / scale if scale > 0 else positive_response
        ax_line.plot(
            x,
            normalized,
            marker="o",
            markersize=4.6,
            linewidth=2.0,
            label=short_names.get(model, model),
            color=colors.get(model, "#4C78A8"),
            alpha=0.9 if scale > 0 else 0.45,
        )
    ax_bar.axhline(0, color="#2f2a24", linewidth=0.8, alpha=0.62)
    y_pad = max(0.035, (max_value - min_value) * 0.18)
    ax_bar.set_ylim(min(-0.03, min_value - y_pad * 0.35), max_value + y_pad)
    ax_bar.set_ylabel("Mean impact score", fontsize=10, fontweight="bold", color="#202124")
    ax_bar.set_title("Variant effect scores by perturbation type", fontsize=15, fontweight="bold", loc="left", color="#202124")
    ax_bar.text(
        0.01,
        0.93,
        "Bars compare mean delta scores across ClinVar target/label classes.",
        transform=ax_bar.transAxes,
        ha="left",
        va="top",
        color="#6b6258",
        fontsize=8.5,
        bbox={"facecolor": "#fffdf8", "edgecolor": "none", "alpha": 0.82, "pad": 2.5},
    )
    ax_line.set_ylabel("Within-model\nnormalized", fontsize=9, fontweight="bold", color="#202124")
    ax_line.set_ylim(-0.04, 1.08)
    ax_line.set_yticks([0.0, 0.5, 1.0])
    ax_line.text(
        0.01,
        0.94,
        "Lower panel normalizes each model to its own strongest positive response; use it for pattern, not magnitude.",
        transform=ax_line.transAxes,
        ha="left",
        va="top",
        color="#6b6258",
        fontsize=8.2,
        bbox={"facecolor": "#fffdf8", "edgecolor": "none", "alpha": 0.82, "pad": 2.5},
    )
    ax_line.set_xticks(x, [variant_labels.get(variant, variant) for variant in variant_order], rotation=18, ha="right")
    handles, labels = ax_bar.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=min(4, max(1, len(labels))),
        frameon=True,
        fontsize=8,
        bbox_to_anchor=(0.5, 0.005),
        facecolor="#fffdf8",
        edgecolor="#d9cfc0",
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.93, bottom=0.17)
    fig.savefig(out_dir / "variant_effect_stratified_by_type.png", dpi=220)
    plt.close(fig)


def plot_calibration(scores: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    metrics = summarize_metrics(scores)
    selected_model = str(metrics.sort_values("auprc", ascending=False).iloc[0]["model"])
    selected = scores[scores["model"] == selected_model]
    if selected.empty:
        selected = scores
    model_name = str(selected["model"].iloc[0])
    bins = calibration_bins(selected["label"], selected["impact_score"], bins=8)
    bins.insert(0, "model", model_name)
    base_rate = float(selected["label"].mean())
    short_names = {
        "CNN baseline (PyTorch Conv1D)": "CNN baseline",
        "RNA-FM frozen encoder + MLP": "RNA-FM frozen",
        "RNABERT frozen encoder + MLP": "RNABERT frozen",
        "SpliceAI real sequence model": "SpliceAI real",
        "Pangolin real sequence model": "Pangolin real",
        "MMSplice real sequence model": "MMSplice real",
        "MaxEntScan real local score": "MaxEntScan real",
    }
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    fig.patch.set_facecolor("#fbf7ef")
    ax.set_facecolor("#fffdf8")
    ax.plot(bins["mean_score"], bins["positive_rate"], marker="o", linewidth=2.4, color="#4C78A8")
    ax.axhline(base_rate, color="#E45756", linestyle="--", linewidth=1.4, label=f"Random baseline = {base_rate:.3f}")
    ax.set_xlabel("Mean impact score bin (delta score)")
    ax.set_ylabel("Observed splice-altering rate")
    ax.set_title(
        f"Calibration by delta-score bins ({short_names.get(model_name, model_name)})",
        fontsize=13,
        fontweight="bold",
        loc="left",
    )
    ax.text(
        0.02,
        0.94,
        f"n={len(selected)}, bins={len(bins)}; scores are current model/tool delta scores, not embedding-distance scores.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color="#6b6258",
        fontsize=8.4,
        bbox={"facecolor": "#fffdf8", "edgecolor": "none", "alpha": 0.82, "pad": 2.5},
    )
    ax.grid(axis="both", color="#d9cfc0", alpha=0.58, linewidth=0.8)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#817768")
    ax.spines["bottom"].set_color("#817768")
    ax.legend(frameon=True, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_dir / "exp3_calibration_curve.png", dpi=220)
    plt.close(fig)
    return bins


def plot_metric_bars(metrics: pd.DataFrame, out_dir: Path, auprc_baseline: float = 0.5) -> None:
    colors = {
        "CNN baseline (PyTorch Conv1D)": "#4C78A8",
        "RNA-FM frozen encoder + MLP": "#F58518",
        "RNABERT frozen encoder + MLP": "#54A24B",
        "SpliceAI real sequence model": "#B279A2",
        "Pangolin real sequence model": "#72B7B2",
        "MMSplice real sequence model": "#E45756",
        "MaxEntScan real local score": "#9D755D",
    }
    short_names = {
        "CNN baseline (PyTorch Conv1D)": "CNN baseline",
        "RNA-FM frozen encoder + MLP": "RNA-FM frozen",
        "RNABERT frozen encoder + MLP": "RNABERT frozen",
        "SpliceAI real sequence model": "SpliceAI real",
        "Pangolin real sequence model": "Pangolin real",
        "MMSplice real sequence model": "MMSplice real",
        "MaxEntScan real local score": "MaxEntScan real",
    }
    for metric, filename, label in [
        ("auroc", "exp3_variant_auroc.png", "AUROC"),
        ("auprc", "exp3_variant_auprc.png", "AUPRC"),
    ]:
        fig, ax = plt.subplots(figsize=(8.6, 5.1))
        fig.patch.set_facecolor("#fbf7ef")
        ax.set_facecolor("#fffdf8")
        order = metrics.sort_values(metric)
        labels = [short_names.get(model, model) for model in order["model"]]
        bars = ax.barh(labels, order[metric], color=[colors.get(model, "#4C78A8") for model in order["model"]], edgecolor="#2f2a24", linewidth=0.45)
        ax.set_xlabel(label)
        ax.set_title(f"Experiment 3 real ClinVar variant effect {label}", fontsize=13, fontweight="bold", loc="left")
        ax.set_xlim(0.0, 1.03)
        if metric == "auprc":
            ax.axvline(auprc_baseline, color="#E45756", linestyle="--", linewidth=1.4, label=f"Random baseline = {auprc_baseline:.3f}")
            ax.legend(frameon=True, fontsize=8, loc="lower right")
        elif metric == "auroc":
            ax.axvline(0.5, color="#E45756", linestyle="--", linewidth=1.4, label="Random baseline = 0.500")
            ax.legend(frameon=True, fontsize=8, loc="lower right")
        for bar, value in zip(bars, order[metric]):
            ax.text(float(value) + 0.012, bar.get_y() + bar.get_height() / 2, f"{float(value):.3f}", va="center", fontsize=8.5, color="#202124")
        ax.grid(axis="x", color="#d9cfc0", alpha=0.58, linewidth=0.8)
        for side in ["top", "right"]:
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_color("#817768")
        ax.spines["bottom"].set_color("#817768")
        fig.tight_layout()
        fig.savefig(out_dir / filename, dpi=220)
        plt.close(fig)


def plot_delta_boxplot(scores: pd.DataFrame, out_dir: Path) -> None:
    preferred = [
        "CNN baseline (PyTorch Conv1D)",
        "RNA-FM frozen encoder + MLP",
        "RNABERT frozen encoder + MLP",
        "SpliceAI real sequence model",
        "Pangolin real sequence model",
        "MMSplice real sequence model",
        "MaxEntScan real local score",
    ]
    available = scores["model"].drop_duplicates().tolist()
    selected_models = [model for model in preferred if model in available]
    selected_models.extend([model for model in available if model not in selected_models])
    subset = scores[scores["model"].isin(selected_models)].copy()
    short_names = {
        "CNN baseline (PyTorch Conv1D)": "CNN\nbaseline",
        "RNA-FM frozen encoder + MLP": "RNA-FM\nfrozen",
        "RNABERT frozen encoder + MLP": "RNABERT\nfrozen",
        "SpliceAI real sequence model": "SpliceAI\nreal",
        "Pangolin real sequence model": "Pangolin\nreal",
        "MMSplice real sequence model": "MMSplice\nreal",
        "MaxEntScan real local score": "MaxEntScan\nreal",
    }
    colors = ["#d8e8f5", "#fae1bd"]
    fig, axes = plt.subplots(1, len(selected_models), figsize=(2.35 * len(selected_models), 4.8), sharey=True)
    fig.patch.set_facecolor("#fbf7ef")
    if len(selected_models) == 1:
        axes = [axes]
    for ax, model in zip(axes, selected_models):
        ax.set_facecolor("#fffdf8")
        group = subset[subset["model"] == model]
        values = [
            group[group["label"] == 0]["impact_score"].to_numpy(),
            group[group["label"] == 1]["impact_score"].to_numpy(),
        ]
        boxes = ax.boxplot(values, labels=["neutral", "splice\naltering"], showfliers=False, patch_artist=True, widths=0.58)
        for patch, color in zip(boxes["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_edgecolor("#2f2a24")
            patch.set_alpha(0.95)
        for median in boxes["medians"]:
            median.set_color("#202124")
            median.set_linewidth(1.4)
        ax.set_title(short_names.get(model, model), fontsize=9.5, fontweight="bold")
        ax.grid(axis="y", color="#d9cfc0", alpha=0.58, linewidth=0.8)
        for side in ["top", "right"]:
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_color("#817768")
        ax.spines["bottom"].set_color("#817768")
        ax.tick_params(axis="x", labelsize=8)
    axes[0].set_ylabel("Impact score")
    fig.suptitle("Experiment 3 delta score distribution by current real model/tool", fontsize=14, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(out_dir / "exp3_delta_score_boxplot.png", dpi=220)
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
    model = models[0]
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


def run(
    output_tables: Path,
    output_figures: Path,
    random_state: int = 42,
    variant_path: Path = EXP3_DATA_DIR / "clinvar_splicing_variants.csv",
) -> dict[str, pd.DataFrame]:
    ensure_dirs(output_tables, output_figures)
    ensure_inputs(variant_path)
    variants = read_csv(variant_path)
    models = train_models(random_state=random_state)
    classifier_scores = score_variants(models, variants)
    external_scores = run_external_tool_scores(variant_path)
    scores = pd.concat([classifier_scores, external_scores], ignore_index=True, sort=False)
    metrics = summarize_metrics(scores)
    matched_metrics = summarize_distance_matched_metrics(scores)
    topk = summarize_topk(scores)
    variant_summary = summarize_variant_types(scores)
    calibration = plot_calibration(scores, output_figures)
    write_dataframe(output_tables / "experiment_3A_variant_scores.csv", scores)
    write_dataframe(output_tables / "experiment_3A_variant_metrics.csv", metrics)
    if not matched_metrics.empty:
        write_dataframe(output_tables / "experiment_3A_distance_matched_variant_metrics.csv", matched_metrics)
    write_dataframe(output_tables / "experiment_3A_topk_enrichment_curve.csv", topk)
    write_dataframe(output_tables / "variant_effect_stratified_by_type.csv", variant_summary)
    write_dataframe(output_tables / "experiment_3A_calibration_bins.csv", calibration)
    plot_metric_bars(metrics, output_figures, auprc_baseline=float(variants["label"].astype(int).mean()))
    plot_delta_boxplot(scores, output_figures)
    plot_variant_type_summary(variant_summary, output_figures)
    run_saturation(models, output_tables, output_figures)
    clinvar = run_clinvar_format_control(output_tables, output_figures, models)
    sqtl = run_sqtl_case_study(output_tables, models)
    variant_counts = (
        variants.groupby(["variant_type", "label_name"], as_index=False)
        .size()
        .rename(columns={"size": "rows"})
        .sort_values("variant_type")
    )
    write_report(PROJECT_ROOT / "reports/experiment_3.md", variant_path, variant_counts, metrics, matched_metrics, variant_summary, clinvar, sqtl)
    return {
        "scores": scores,
        "metrics": metrics,
        "topk": topk,
        "variant_summary": variant_summary,
        "clinvar": clinvar,
        "sqtl": sqtl,
    }


def run_clinvar_format_control(output_tables: Path, output_figures: Path, models: list[object]) -> pd.DataFrame:
    clinvar = build_clinvar_format_control()
    scores = score_variants(models, clinvar)
    metrics = summarize_metrics(scores)
    write_dataframe(output_tables / "experiment_3B_format_control_scores.csv", scores)
    write_dataframe(output_tables / "experiment_3B_format_control_metrics.csv", metrics)
    best = scores[scores["model"] == metrics.iloc[0]["model"]].copy()
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    best.boxplot(column="impact_score", by="label_name", ax=ax)
    ax.set_title("ClinVar format-control scores")
    ax.set_xlabel("Label")
    ax.set_ylabel("Impact score")
    fig.suptitle("")
    fig.tight_layout()
    fig.savefig(output_figures / "exp3_clinvar_format_control_scores.png", dpi=180)
    plt.close(fig)
    return metrics


def run_sqtl_case_study(output_tables: Path, models: list[object]) -> pd.DataFrame:
    sqtl = build_sqtl_format_control()
    scored = []
    model = models[0]
    for _, row in sqtl.iterrows():
        wt = model.predict_proba([str(row["wt_sequence"])])[0]
        mut = model.predict_proba([str(row["mut_sequence"])])[0]
        delta = float(np.max(np.abs(mut[:2] - wt[:2])))
        scored.append({**row.to_dict(), "model_delta_score": delta})
    frame = pd.DataFrame(scored)
    write_dataframe(output_tables / "experiment_3C_sqtl_case_study.csv", frame)
    return frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run C-part experiment 3: real ClinVar splice variant effects.")
    parser.add_argument("--config", type=Path, default=CONFIG_ROOT / "exp3_variant_effect.yaml")
    parser.add_argument("--tables", type=Path, default=EXP3_TABLES_DIR)
    parser.add_argument("--figures", type=Path, default=EXP3_FIGURES_DIR)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config) if args.config else {}
    variant_path = Path(config.get("variant_table", EXP3_DATA_DIR / "clinvar_splicing_variants.csv"))
    outputs = run(
        Path(config.get("tables_dir", args.tables)),
        Path(config.get("figures_dir", args.figures)),
        int(config.get("seed", args.seed)),
        variant_path=variant_path,
    )
    print(outputs["metrics"].to_string(index=False))


if __name__ == "__main__":
    main()
