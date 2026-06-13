from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

from src.utils import EXP1_TABLES_DIR, EXP2_TABLES_DIR, EXP3_TABLES_DIR, PROJECT_ROOT, ensure_dirs


REPORT_IMAGE_DIR = PROJECT_ROOT / "report_letax" / "images"

MODEL_COLORS = {
    "CNN motif baseline": "#4C78A8",
    "CNN baseline (PyTorch Conv1D)": "#4C78A8",
    "RNA-FM frozen k-mer + MLP": "#F58518",
    "RNA-FM frozen encoder + MLP": "#F58518",
    "RNA-FM zero-shot embedding distance": "#F58518",
    "RNABERT frozen token + MLP": "#54A24B",
    "RNABERT frozen encoder + MLP": "#54A24B",
    "RNABERT zero-shot token distance": "#54A24B",
    "SpliceAI signal proxy": "#E45756",
    "SpliceAI optional real tool (proxy fallback)": "#E45756",
    "MaxEntScan optional tool (proxy fallback)": "#8E6C8A",
    "Pangolin optional tool (small case-study proxy)": "#7F7F7F",
}
LABEL_COLORS = {"neutral": "#B9B9B9", "splice_altering": "#D95F5F"}
CLASS_ORDER = ["acceptor", "donor", "non_splice"]
CLASS_LABELS = {"acceptor": "Acceptor", "donor": "Donor", "non_splice": "Non-splice"}
MOTIF_ORDER = ["ESE_proxy", "ESS_proxy", "ISE_proxy", "ISS_proxy"]
VARIANT_ORDER = ["donor_loss", "acceptor_loss", "donor_gain", "acceptor_gain", "neutral_far_snv"]
VARIANT_LABELS = {
    "donor_loss": "Donor loss",
    "acceptor_loss": "Acceptor loss",
    "donor_gain": "Donor gain",
    "acceptor_gain": "Acceptor gain",
    "neutral_far_snv": "Neutral far SNV",
}
SELECTED_DELTA_MODELS = [
    "SpliceAI signal proxy",
    "RNA-FM frozen k-mer + MLP",
    "RNABERT frozen token + MLP",
    "MaxEntScan optional tool (proxy fallback)",
]
SELECTED_VARIANT_MODELS = [
    "RNABERT zero-shot token distance",
    "RNA-FM zero-shot embedding distance",
    "SpliceAI signal proxy",
    "MaxEntScan optional tool (proxy fallback)",
]
MAIN_EXP3_METRIC_MODELS = [
    "RNABERT zero-shot token distance",
    "RNA-FM zero-shot embedding distance",
    "RNABERT zero-shot pseudo-likelihood",
    "RNA-FM zero-shot pseudo-likelihood",
    "SpliceAI signal proxy",
    "CNN motif baseline",
    "RNABERT frozen token + MLP",
    "RNA-FM frozen k-mer + MLP",
]


def _set_style() -> None:
    sns.set_theme(
        context="paper",
        style="whitegrid",
        rc={
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#5A5A5A",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "grid.color": "#D8D8D8",
            "grid.linewidth": 0.7,
            "legend.frameon": False,
            "font.size": 9,
        },
    )


def _model_color(model: str) -> str:
    return MODEL_COLORS.get(model, "#4C78A8")


def _short_model_name(model: str) -> str:
    replacements = {
        "RNA-FM frozen k-mer + MLP": "RNA-FM fallback k-mer+signal",
        "RNA-FM frozen encoder + MLP": "RNA-FM fallback k-mer+signal",
        "RNABERT frozen token + MLP": "RNABERT fallback token+signal",
        "RNABERT frozen encoder + MLP": "RNABERT fallback token+signal",
        "SpliceAI signal proxy": "SpliceAI-style signal proxy",
        "SpliceAI optional real tool (proxy fallback)": "SpliceAI optional proxy",
        "MaxEntScan optional tool (proxy fallback)": "MaxEntScan proxy",
        "CNN motif baseline": "CNN motif baseline",
        "CNN baseline (PyTorch Conv1D)": "CNN motif baseline",
        "Pangolin optional tool (small case-study proxy)": "Pangolin small proxy",
    }
    return replacements.get(model, model)


def _wrapped_model_name(model: str, width: int = 26) -> str:
    return "\n".join(textwrap.wrap(_short_model_name(model), width=width, break_long_words=False))


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_exp1_macro_f1(metrics: pd.DataFrame, out_dir: Path) -> None:
    frame = metrics.copy().sort_values("macro_f1", ascending=True)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    labels = [_wrapped_model_name(model, width=28) for model in frame["model"]]
    colors = [_model_color(model) for model in frame["model"]]
    bars = ax.barh(labels, frame["macro_f1"], color=colors, alpha=0.9)
    ax.set_xlabel("Macro-F1")
    ax.set_title("Experiment 1 splice-site classification Macro-F1")
    ax.set_xlim(0.82, 1.01)
    ax.grid(axis="x", alpha=0.6)
    ax.grid(axis="y", visible=False)
    for bar, value in zip(bars, frame["macro_f1"]):
        ax.text(value + 0.004, bar.get_y() + bar.get_height() / 2, f"{value:.3f}", va="center", fontsize=8)
    _save(fig, out_dir / "experiment_1_macro_f1.png")


def plot_exp1_confusion_matrices(confusion: pd.DataFrame, out_dir: Path) -> None:
    models = confusion["model"].drop_duplicates().tolist()
    fig, axes = plt.subplots(1, len(models), figsize=(4.6 * len(models), 4.2), squeeze=False)
    order = ["donor", "acceptor", "non_splice"]
    for ax, model in zip(axes[0], models):
        sub = confusion[confusion["model"] == model]
        matrix = (
            sub.pivot_table(index="true_label_name", columns="pred_label_name", values="count", aggfunc="sum")
            .reindex(index=order, columns=order)
            .fillna(0)
            .astype(int)
        )
        sns.heatmap(
            matrix,
            annot=True,
            fmt="d",
            cmap="Blues",
            cbar=False,
            square=True,
            linewidths=0.5,
            linecolor="white",
            ax=ax,
        )
        ax.set_title(_wrapped_model_name(model, width=22), fontsize=9)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
    fig.suptitle("Experiment 1 confusion matrices", y=1.03, fontsize=12, fontweight="bold")
    _save(fig, out_dir / "experiment_1_confusion_matrices.png")


def plot_context_macro_f1(metrics: pd.DataFrame, out_dir: Path) -> None:
    frame = metrics.copy()
    frame["window_flank"] = frame["window_flank"].astype(int)
    fig, ax = plt.subplots(figsize=(8.4, 4.9))
    windows = sorted(frame["window_flank"].unique())
    markers = ["o", "s", "^", "D", "P"]
    linestyles = ["-", "-", "--", "-.", ":"]
    for idx, (model, group) in enumerate(frame.groupby("model", sort=False)):
        group = group.sort_values("window_flank")
        color = _model_color(model)
        label = _short_model_name(model)
        ax.plot(
            group["window_flank"],
            group["macro_f1"],
            marker=markers[idx % len(markers)],
            linestyle=linestyles[idx % len(linestyles)],
            linewidth=2.2,
            markersize=5,
            color=color,
            label=label,
        )
        end = group.iloc[-1]
        endpoint_offsets = {
            "CNN motif baseline": -2,
            "RNA-FM frozen k-mer + MLP": 0,
            "RNABERT frozen token + MLP": 0,
            "SpliceAI signal proxy": 0,
        }
        ax.annotate(
            label,
            xy=(end["window_flank"], end["macro_f1"]),
            xytext=(5, endpoint_offsets.get(model, 0)),
            textcoords="offset points",
            va="center",
            fontsize=8,
            color=color,
        )
    ax.set_title("Experiment 2A context ablation: Macro-F1")
    ax.set_ylabel("Macro-F1")
    ax.set_xlabel("Context flank (nt)")
    ax.set_xticks(windows, [f"{w} nt" for w in windows])
    ax.set_ylim(max(0.86, float(frame["macro_f1"].min()) - 0.008), 1.005)
    ax.grid(axis="y", alpha=0.7)
    ax.grid(axis="x", alpha=0.18)
    _save(fig, out_dir / "exp2A_context_macro_f1.png")


def plot_context_auprc(metrics: pd.DataFrame, out_dir: Path) -> None:
    frame = metrics.copy()
    frame["window_flank"] = frame["window_flank"].astype(int)
    fig, ax = plt.subplots(figsize=(8.4, 4.9))
    windows = sorted(frame["window_flank"].unique())
    markers = ["o", "s", "^", "D", "P"]
    linestyles = ["-", "-", "--", "-.", ":"]
    for idx, (model, group) in enumerate(frame.groupby("model", sort=False)):
        group = group.sort_values("window_flank")
        color = _model_color(model)
        label = _short_model_name(model)
        ax.plot(
            group["window_flank"],
            group["auprc"],
            marker=markers[idx % len(markers)],
            linestyle=linestyles[idx % len(linestyles)],
            linewidth=2.2,
            markersize=5,
            color=color,
            label=label,
        )
        end = group.iloc[-1]
        endpoint_offsets = {
            "RNA-FM frozen k-mer + MLP": 0,
            "RNABERT frozen token + MLP": 8,
            "SpliceAI signal proxy": -8,
        }
        ax.annotate(
            label,
            xy=(end["window_flank"], end["auprc"]),
            xytext=(5, endpoint_offsets.get(model, 0)),
            textcoords="offset points",
            va="center",
            fontsize=8,
            color=color,
        )
    ax.set_title("Experiment 2A context ablation: Macro AUPRC")
    ax.set_ylabel("Macro AUPRC")
    ax.set_xlabel("Context flank (nt)")
    ax.set_xticks(windows, [f"{w} nt" for w in windows])
    ax.set_ylim(max(0.885, float(frame["auprc"].min()) - 0.004), 1.0015)
    ax.grid(axis="y", alpha=0.7)
    ax.grid(axis="x", alpha=0.18)
    _save(fig, out_dir / "exp2A_context_auprc.png")


def plot_hard_negative_fpr(hard_metrics: pd.DataFrame, out_dir: Path) -> None:
    frame = hard_metrics.dropna(subset=["hard_negative_fpr"]).copy()
    frame = frame.sort_values("hard_negative_fpr", ascending=True)
    frame["model_short"] = frame["model"].map(lambda model: _wrapped_model_name(model, width=28))
    frame["false_positive_label"] = frame.apply(
        lambda row: f"{int(row['hard_negative_false_positives'])}/{int(row['hard_negative_rows'])}",
        axis=1,
    )
    fig, ax = plt.subplots(figsize=(7.7, 4.6))
    colors = [_model_color(model) for model in frame["model"]]
    bars = ax.barh(frame["model_short"], frame["hard_negative_fpr"], color=colors, alpha=0.9)
    ax.set_xlabel("Hard-negative false positive rate")
    ax.set_title("Experiment 2B GT/AG hard-negative FPR")
    ax.set_xlim(0.0, min(1.0, max(0.5, float(frame["hard_negative_fpr"].max()) * 1.25)))
    ax.grid(axis="x", alpha=0.65)
    ax.grid(axis="y", visible=False)
    for bar, (_, row) in zip(bars, frame.iterrows()):
        x = float(row["hard_negative_fpr"])
        ax.text(
            x + ax.get_xlim()[1] * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{x:.3f} ({row['false_positive_label']})",
            va="center",
            fontsize=8,
            color="#333333",
        )
    _save(fig, out_dir / "exp2B_hard_negative_fpr.png")


def plot_regulatory_motif_masking(summary: pd.DataFrame, out_dir: Path) -> None:
    frame = summary.copy()
    pivot = (
        frame.pivot_table(index="motif_group", columns="label_name", values="mean_delta", aggfunc="mean")
        .reindex(index=MOTIF_ORDER, columns=CLASS_ORDER)
        .fillna(0.0)
    )
    vmax = max(0.02, float(np.nanmax(np.abs(pivot.to_numpy()))))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    fig, ax = plt.subplots(figsize=(6.9, 4.7))
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    image = ax.imshow(pivot.to_numpy(), cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(pivot.columns)), [CLASS_LABELS[col] for col in pivot.columns])
    ax.set_yticks(range(len(pivot.index)), [idx.replace("_proxy", "") for idx in pivot.index])
    ax.set_xlabel("Predicted class")
    ax.set_ylabel("Masked motif group")
    ax.set_title("Regulatory motif masking effect on class probability")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = float(pivot.iloc[i, j])
            text_color = "white" if abs(value) > vmax * 0.55 else "#222222"
            ax.text(j, i, f"{value:+.3f}", ha="center", va="center", fontsize=8, color=text_color)
    cbar = fig.colorbar(image, ax=ax, shrink=0.86, pad=0.03)
    cbar.set_label("Mean delta: original - masked", weight="bold")
    ax.grid(False)
    _save(fig, out_dir / "regulatory_motif_masking_effect.png")


def _cohens_d(values_a: np.ndarray, values_b: np.ndarray) -> float:
    if len(values_a) < 2 or len(values_b) < 2:
        return float("nan")
    var_a = np.var(values_a, ddof=1)
    var_b = np.var(values_b, ddof=1)
    pooled = np.sqrt(((len(values_a) - 1) * var_a + (len(values_b) - 1) * var_b) / (len(values_a) + len(values_b) - 2))
    if pooled == 0:
        return float("nan")
    return float((np.mean(values_b) - np.mean(values_a)) / pooled)


def plot_delta_score_distribution(scores: pd.DataFrame, out_dir: Path) -> None:
    subset = scores[scores["model"].isin(SELECTED_DELTA_MODELS)].copy()
    subset["label_name"] = subset["label_name"].replace({"splice altering": "splice_altering"})
    fig, axes = plt.subplots(1, len(SELECTED_DELTA_MODELS), figsize=(12.5, 4.2), sharey=False)
    rng = np.random.default_rng(2026)
    for ax, model in zip(axes, SELECTED_DELTA_MODELS):
        group = subset[subset["model"] == model].copy()
        parts = []
        positions = [0, 1]
        for label in ["neutral", "splice_altering"]:
            values = group[group["label_name"] == label]["impact_score"].to_numpy(dtype=float)
            parts.append(values)
        violin = ax.violinplot(parts, positions=positions, widths=0.72, showmeans=False, showmedians=False, showextrema=False)
        for body, color in zip(violin["bodies"], [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]]):
            body.set_facecolor(color)
            body.set_edgecolor("#555555")
            body.set_alpha(0.34)
        box = ax.boxplot(
            parts,
            positions=positions,
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#111111", "linewidth": 1.2},
            boxprops={"linewidth": 0.8},
            whiskerprops={"linewidth": 0.8},
            capprops={"linewidth": 0.8},
        )
        for patch, color in zip(box["boxes"], [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]]):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        for pos, values, color in zip(positions, parts, [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]]):
            if len(values) > 120:
                values = rng.choice(values, size=120, replace=False)
            jitter = rng.normal(0, 0.045, size=len(values))
            ax.scatter(np.full(len(values), pos) + jitter, values, s=7, color=color, alpha=0.28, linewidths=0)
        neutral, altering = parts
        gap = float(np.median(altering) - np.median(neutral)) if len(neutral) and len(altering) else float("nan")
        effect = _cohens_d(neutral, altering)
        top = max(float(np.nanmax(group["impact_score"])), 0.01)
        bottom = min(0.0, float(np.nanmin(group["impact_score"])))
        pad = (top - bottom) * 0.16 if top > bottom else 0.05
        ax.set_ylim(bottom - pad * 0.25, top + pad)
        ax.text(
            0.02,
            0.98,
            f"Delta median={gap:.3f}\nCohen d={effect:.2f}\nn={len(group)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="#333333",
        )
        ax.axhline(0, color="#777777", linewidth=0.8)
        ax.set_xticks(positions, ["Neutral", "Splice\naltering"])
        ax.set_title(_short_model_name(model), color=_model_color(model), fontsize=10)
        ax.grid(axis="y", alpha=0.65)
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("Impact score")
    fig.suptitle("Experiment 3 delta score distribution", y=1.02, fontsize=12, fontweight="bold")
    _save(fig, out_dir / "exp3_delta_score_boxplot.png")


def plot_delta_score_distribution_neutral_zoom_preview(scores: pd.DataFrame, out_dir: Path) -> Path:
    subset = scores[scores["model"].isin(SELECTED_DELTA_MODELS)].copy()
    subset["label_name"] = subset["label_name"].replace({"splice altering": "splice_altering"})
    out_path = out_dir / "exp3_delta_score_boxplot_neutral_zoom_preview.png"
    fig, axes = plt.subplots(2, len(SELECTED_DELTA_MODELS), figsize=(12.8, 5.8), sharex="col")
    rng = np.random.default_rng(2026)
    for col, model in enumerate(SELECTED_DELTA_MODELS):
        group = subset[subset["model"] == model].copy()
        neutral = group[group["label_name"] == "neutral"]["impact_score"].to_numpy(dtype=float)
        altering = group[group["label_name"] == "splice_altering"]["impact_score"].to_numpy(dtype=float)
        full_ax = axes[0, col]
        zoom_ax = axes[1, col]

        parts = [neutral, altering]
        violin = full_ax.violinplot(parts, positions=[0, 1], widths=0.72, showmeans=False, showmedians=False, showextrema=False)
        for body, color in zip(violin["bodies"], [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]]):
            body.set_facecolor(color)
            body.set_edgecolor("#555555")
            body.set_alpha(0.34)
        box = full_ax.boxplot(
            parts,
            positions=[0, 1],
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#111111", "linewidth": 1.15},
            boxprops={"linewidth": 0.75},
            whiskerprops={"linewidth": 0.75},
            capprops={"linewidth": 0.75},
        )
        for patch, color in zip(box["boxes"], [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]]):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        for pos, values, color in zip([0, 1], parts, [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]]):
            sample = values
            if len(sample) > 100:
                sample = rng.choice(sample, size=100, replace=False)
            full_ax.scatter(np.full(len(sample), pos) + rng.normal(0, 0.045, len(sample)), sample, s=7, color=color, alpha=0.28, linewidths=0)

        neutral_high = float(np.nanpercentile(neutral, 98)) if len(neutral) else 0.01
        neutral_low = float(np.nanpercentile(neutral, 2)) if len(neutral) else 0.0
        zoom_top = max(0.01, neutral_high * 1.35)
        zoom_bottom = min(0.0, neutral_low - zoom_top * 0.08)
        full_ax.axhspan(zoom_bottom, zoom_top, color="#F2F2F2", zorder=0)
        full_ax.axhline(0, color="#777777", linewidth=0.8)
        top = max(float(np.nanmax(group["impact_score"])), 0.01)
        bottom = min(0.0, float(np.nanmin(group["impact_score"])))
        pad = (top - bottom) * 0.16 if top > bottom else 0.05
        full_ax.set_ylim(bottom - pad * 0.25, top + pad)
        full_ax.set_title(_short_model_name(model), color=_model_color(model), fontsize=10)
        full_ax.set_xticks([0, 1], [])
        full_ax.grid(axis="y", alpha=0.62)
        full_ax.grid(axis="x", visible=False)

        zoom_values = neutral
        jitter = rng.normal(0, 0.08, len(zoom_values))
        zoom_ax.scatter(
            np.full(len(zoom_values), 0) + jitter,
            zoom_values,
            s=13,
            color=LABEL_COLORS["neutral"],
            edgecolors="#777777",
            linewidths=0.25,
            alpha=0.62,
        )
        zoom_ax.boxplot(
            [zoom_values],
            positions=[0],
            widths=0.34,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#111111", "linewidth": 1.1},
            boxprops={"facecolor": LABEL_COLORS["neutral"], "alpha": 0.45, "linewidth": 0.8},
            whiskerprops={"linewidth": 0.8},
            capprops={"linewidth": 0.8},
        )
        zoom_ax.axhline(0, color="#777777", linewidth=0.8)
        zoom_ax.set_xlim(-0.55, 0.55)
        zoom_ax.set_ylim(zoom_bottom, zoom_top)
        zoom_ax.set_xticks([0], ["Neutral\nzoom"])
        zoom_ax.text(
            0.03,
            0.94,
            f"median={np.median(neutral):.4f}\n98th={neutral_high:.4f}",
            transform=zoom_ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.2,
            color="#333333",
        )
        zoom_ax.grid(axis="y", alpha=0.62)
        zoom_ax.grid(axis="x", visible=False)
    axes[0, 0].set_ylabel("Impact score")
    axes[1, 0].set_ylabel("Neutral-only zoom")
    fig.suptitle("Experiment 3 delta score distribution with neutral zoom preview", y=1.02, fontsize=12, fontweight="bold")
    _save(fig, out_path)
    return out_path


def plot_delta_score_distribution_neutral_emphasis_preview(scores: pd.DataFrame, out_dir: Path) -> Path:
    subset = scores[scores["model"].isin(SELECTED_DELTA_MODELS)].copy()
    subset["label_name"] = subset["label_name"].replace({"splice altering": "splice_altering"})
    out_path = out_dir / "exp3_delta_score_boxplot_neutral_emphasis_preview.png"
    fig, axes = plt.subplots(1, len(SELECTED_DELTA_MODELS), figsize=(12.5, 4.2), sharey=False)
    rng = np.random.default_rng(2026)
    for ax, model in zip(axes, SELECTED_DELTA_MODELS):
        group = subset[subset["model"] == model].copy()
        neutral = group[group["label_name"] == "neutral"]["impact_score"].to_numpy(dtype=float)
        altering = group[group["label_name"] == "splice_altering"]["impact_score"].to_numpy(dtype=float)
        parts = [neutral, altering]
        neutral_band_top = max(0.012, float(np.nanpercentile(neutral, 98)) * 1.45 if len(neutral) else 0.012)
        neutral_band_bottom = min(0.0, float(np.nanpercentile(neutral, 1)) if len(neutral) else 0.0)
        ax.axhspan(neutral_band_bottom, neutral_band_top, color="#F1F1F1", zorder=0)

        violin = ax.violinplot(parts, positions=[0, 1], widths=[0.9, 0.72], showmeans=False, showmedians=False, showextrema=False)
        for idx, (body, color) in enumerate(zip(violin["bodies"], [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]])):
            body.set_facecolor(color)
            body.set_edgecolor("#222222" if idx == 0 else "#555555")
            body.set_linewidth(1.2 if idx == 0 else 0.8)
            body.set_alpha(0.72 if idx == 0 else 0.34)

        box = ax.boxplot(
            parts,
            positions=[0, 1],
            widths=[0.34, 0.28],
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#111111", "linewidth": 1.35},
            boxprops={"linewidth": 1.0},
            whiskerprops={"linewidth": 1.0},
            capprops={"linewidth": 1.0},
        )
        for idx, (patch, color) in enumerate(zip(box["boxes"], [LABEL_COLORS["neutral"], LABEL_COLORS["splice_altering"]])):
            patch.set_facecolor(color)
            patch.set_edgecolor("#222222" if idx == 0 else "#555555")
            patch.set_alpha(0.9 if idx == 0 else 0.75)

        for pos, values, color, size, alpha in [
            (0, neutral, "#4F4F4F", 15, 0.55),
            (1, altering, LABEL_COLORS["splice_altering"], 7, 0.25),
        ]:
            sample = values
            if len(sample) > 180:
                sample = rng.choice(sample, size=180, replace=False)
            jitter = rng.normal(0, 0.07 if pos == 0 else 0.045, size=len(sample))
            ax.scatter(
                np.full(len(sample), pos) + jitter,
                sample,
                s=size,
                color=color,
                alpha=alpha,
                linewidths=0,
                zorder=3 if pos == 0 else 2,
            )

        neutral_98 = float(np.nanpercentile(neutral, 98)) if len(neutral) else 0.0
        ax.text(
            0.0,
            neutral_band_top,
            f"neutral 98th={neutral_98:.4f}",
            ha="center",
            va="bottom",
            fontsize=6.8,
            color="#444444",
        )
        gap = float(np.median(altering) - np.median(neutral)) if len(neutral) and len(altering) else float("nan")
        effect = _cohens_d(neutral, altering)
        top = max(float(np.nanmax(group["impact_score"])), 0.01)
        bottom = min(0.0, float(np.nanmin(group["impact_score"])))
        pad = (top - bottom) * 0.16 if top > bottom else 0.05
        ax.set_ylim(bottom - pad * 0.25, top + pad)
        ax.text(
            0.02,
            0.98,
            f"Delta median={gap:.3f}\nCohen d={effect:.2f}\nn={len(group)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="#333333",
        )
        ax.axhline(0, color="#777777", linewidth=0.8)
        ax.set_xticks([0, 1], ["Neutral", "Splice\naltering"])
        ax.set_title(_short_model_name(model), color=_model_color(model), fontsize=10)
        ax.grid(axis="y", alpha=0.65)
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("Impact score")
    fig.suptitle("Experiment 3 delta score distribution - neutral emphasis preview", y=1.02, fontsize=12, fontweight="bold")
    _save(fig, out_path)
    return out_path


def plot_variant_effect_by_type(summary: pd.DataFrame, out_dir: Path) -> None:
    frame = summary[summary["model"].isin(SELECTED_VARIANT_MODELS)].copy()
    frame["variant_type"] = pd.Categorical(frame["variant_type"], categories=VARIANT_ORDER, ordered=True)
    frame = frame.sort_values(["model", "variant_type"])
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 7.1), gridspec_kw={"height_ratios": [1.45, 1.0], "hspace": 0.34})

    panel = axes[0]
    model_order = [model for model in SELECTED_VARIANT_MODELS if model in set(frame["model"])]
    x = np.arange(len(VARIANT_ORDER))
    width = 0.72 / max(1, len(model_order))
    for idx, model in enumerate(model_order):
        sub = frame[frame["model"] == model].set_index("variant_type")
        values = [float(sub.loc[var, "mean_score"]) if var in sub.index else 0.0 for var in VARIANT_ORDER]
        offset = -0.36 + width / 2 + idx * width
        bars = panel.bar(x + offset, values, width=width, label=_short_model_name(model), color=_model_color(model), alpha=0.88)
        for bar, value in zip(bars, values):
            if abs(value) >= 0.2:
                panel.text(
                    bar.get_x() + bar.get_width() / 2,
                    value + 0.035 * np.sign(value if value else 1),
                    f"{value:.2f}",
                    ha="center",
                    va="bottom" if value >= 0 else "top",
                    fontsize=6.8,
                    rotation=90,
                )
    panel.axhline(0, color="#666666", linewidth=0.9)
    panel.set_xticks(x, [VARIANT_LABELS[v] for v in VARIANT_ORDER], rotation=20, ha="right")
    panel.set_ylabel("Mean impact score")
    panel.set_title("Variant effect scores by perturbation type")
    panel.legend(ncol=2, loc="upper right", fontsize=8)
    panel.grid(axis="y", alpha=0.65)
    panel.grid(axis="x", visible=False)

    normalized = []
    for model, group in frame.groupby("model", sort=False):
        low = float(group["mean_score"].min())
        high = float(group["mean_score"].max())
        denom = high - low
        for row in group.itertuples(index=False):
            normalized.append(
                {
                    "model": row.model,
                    "variant_type": row.variant_type,
                    "normalized": 0.0 if denom == 0 else (float(row.mean_score) - low) / denom,
                }
            )
    norm_frame = pd.DataFrame(normalized)
    dot_ax = axes[1]
    y_positions = {variant: idx for idx, variant in enumerate(VARIANT_ORDER)}
    for model in model_order:
        sub = norm_frame[norm_frame["model"] == model]
        dot_ax.plot(
            sub["normalized"],
            sub["variant_type"].map(y_positions),
            marker="o",
            linewidth=1.7,
            markersize=5,
            color=_model_color(model),
            label=_short_model_name(model),
        )
    dot_ax.set_yticks(list(y_positions.values()), [VARIANT_LABELS[v] for v in VARIANT_ORDER])
    dot_ax.set_xlim(-0.04, 1.04)
    dot_ax.set_xlabel("Within-model normalized effect")
    dot_ax.set_ylabel("Perturbation type")
    dot_ax.set_title("Relative response pattern within each method")
    dot_ax.grid(axis="x", alpha=0.65)
    dot_ax.grid(axis="y", alpha=0.18)
    dot_ax.invert_yaxis()
    _save(fig, out_dir / "variant_effect_stratified_by_type.png")


def plot_exp3_metric_bars(metrics: pd.DataFrame, out_dir: Path) -> None:
    for metric, filename, label in [
        ("auroc", "exp3_variant_auroc.png", "AUROC"),
        ("auprc", "exp3_variant_auprc.png", "AUPRC"),
    ]:
        frame = metrics[metrics["model"].isin(MAIN_EXP3_METRIC_MODELS)].copy()
        order_lookup = {model: idx for idx, model in enumerate(MAIN_EXP3_METRIC_MODELS)}
        frame["model_order"] = frame["model"].map(order_lookup)
        frame = frame.sort_values("model_order", ascending=False)
        fig, ax = plt.subplots(figsize=(9.4, 5.6))
        labels = [_wrapped_model_name(model, width=32) for model in frame["model"]]
        colors = [_model_color(model) for model in frame["model"]]
        bars = ax.barh(labels, frame[metric], color=colors, alpha=0.9)
        ax.set_xlabel(label)
        ax.set_title(f"Experiment 3 artificial variant effect {label}")
        ax.set_xlim(0.0, 1.03)
        ax.grid(axis="x", alpha=0.6)
        ax.grid(axis="y", visible=False)
        for bar, value in zip(bars, frame[metric]):
            ax.text(min(value + 0.01, 1.0), bar.get_y() + bar.get_height() / 2, f"{value:.3f}", va="center", fontsize=7.5)
        _save(fig, out_dir / filename)


def plot_score_bin_response_curve(bins: pd.DataFrame, out_dir: Path) -> None:
    if bins.empty:
        return
    model_name = str(bins["model"].iloc[0]) if "model" in bins.columns else "selected model"
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.plot(bins["mean_score"], bins["positive_rate"], marker="o", linewidth=2.2, color="#4C78A8")
    for _, row in bins.iterrows():
        ax.text(row["mean_score"], row["positive_rate"] + 0.035, f"n={int(row['rows'])}", ha="center", fontsize=7)
    ax.set_xlabel("Mean impact score bin")
    ax.set_ylabel("Observed splice-altering rate")
    ax.set_title(f"Score-bin response curve\n({_short_model_name(model_name)}, bins={len(bins)})")
    ax.set_ylim(-0.04, 1.08)
    ax.grid(alpha=0.35)
    _save(fig, out_dir / "exp3_calibration_curve.png")


def regenerate_report_figures(
    exp1_tables: Path = EXP1_TABLES_DIR,
    exp2_tables: Path = EXP2_TABLES_DIR,
    exp3_tables: Path = EXP3_TABLES_DIR,
    out_dir: Path = REPORT_IMAGE_DIR,
) -> None:
    _set_style()
    ensure_dirs(out_dir)
    plot_exp1_macro_f1(pd.read_csv(exp1_tables / "experiment_1_metrics.csv"), out_dir)
    plot_exp1_confusion_matrices(pd.read_csv(exp1_tables / "experiment_1_confusion_matrices.csv"), out_dir)
    context_metrics = pd.read_csv(exp2_tables / "experiment_2A_multiscale_context.csv")
    plot_context_macro_f1(context_metrics, out_dir)
    plot_context_auprc(context_metrics, out_dir)
    plot_hard_negative_fpr(pd.read_csv(exp2_tables / "experiment_2B_hard_negative.csv"), out_dir)
    plot_regulatory_motif_masking(pd.read_csv(exp2_tables / "regulatory_motif_masking.csv"), out_dir)
    plot_exp3_metric_bars(pd.read_csv(exp3_tables / "experiment_3A_artificial_variant_metrics.csv"), out_dir)
    plot_score_bin_response_curve(pd.read_csv(exp3_tables / "experiment_3A_calibration_bins.csv"), out_dir)
    plot_delta_score_distribution(pd.read_csv(exp3_tables / "experiment_3A_artificial_variant_scores.csv"), out_dir)
    plot_variant_effect_by_type(pd.read_csv(exp3_tables / "variant_effect_stratified_by_type.csv"), out_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Regenerate polished report figures from experiment result tables.")
    parser.add_argument("--exp1-tables", type=Path, default=EXP1_TABLES_DIR)
    parser.add_argument("--exp2-tables", type=Path, default=EXP2_TABLES_DIR)
    parser.add_argument("--exp3-tables", type=Path, default=EXP3_TABLES_DIR)
    parser.add_argument("--out-dir", type=Path, default=REPORT_IMAGE_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    regenerate_report_figures(args.exp1_tables, args.exp2_tables, args.exp3_tables, args.out_dir)
    print(f"Report figures written to {args.out_dir}")


if __name__ == "__main__":
    main()
