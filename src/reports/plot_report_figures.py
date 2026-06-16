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

from src.utils import EXP1_TABLES_DIR, EXP2_TABLES_DIR, EXP3_TABLES_DIR, PROJECT_ROOT, ensure_dirs


REPORT_IMAGE_DIR = PROJECT_ROOT / "report_letax" / "images"
MODEL_COLORS = {
    "CNN baseline (PyTorch Conv1D)": "#4C78A8",
    "RNA-FM frozen encoder + MLP": "#F58518",
    "RNABERT frozen encoder + MLP": "#54A24B",
    "SpliceAI real sequence model": "#B279A2",
    "Pangolin real sequence model": "#72B7B2",
    "MMSplice real sequence model": "#E45756",
    "MaxEntScan real local score": "#9D755D",
}
MODEL_SHORT = {
    "CNN baseline (PyTorch Conv1D)": "CNN\nbaseline",
    "RNA-FM frozen encoder + MLP": "RNA-FM\nfrozen encoder",
    "RNABERT frozen encoder + MLP": "RNABERT\nfrozen encoder",
    "SpliceAI real sequence model": "SpliceAI\nreal",
    "Pangolin real sequence model": "Pangolin\nreal",
    "MMSplice real sequence model": "MMSplice\nreal",
    "MaxEntScan real local score": "MaxEntScan\nreal",
}
BG = "#fbf7ef"
PANEL = "#fffdf8"
DARK = "#202124"
MUTED = "#6b6258"
GRID = "#d9cfc0"


def model_color(model: str) -> str:
    return MODEL_COLORS.get(model, "#4C78A8")


def short_model(model: str) -> str:
    return MODEL_SHORT.get(model, model)


def style_axes(ax: plt.Axes) -> None:
    ax.set_facecolor(PANEL)
    ax.figure.set_facecolor(BG)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#817768")
    ax.spines["bottom"].set_color("#817768")
    ax.tick_params(colors=DARK, labelsize=9)
    ax.title.set_color(DARK)
    ax.xaxis.label.set_color(DARK)
    ax.yaxis.label.set_color(DARK)


def add_note(ax: plt.Axes, text: str, x: float = 0.015, y: float = 0.965) -> None:
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=MUTED,
        fontsize=8.5,
        bbox={"facecolor": PANEL, "edgecolor": "none", "alpha": 0.82, "pad": 2.5},
    )


def save(fig: plt.Figure, path: Path) -> None:
    ensure_dirs(path.parent)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def add_value_labels(ax: plt.Axes, bars, values: pd.Series, fmt: str = "{:.3f}", pad: float = 0.004) -> None:
    x_max = ax.get_xlim()[1]
    for bar, value in zip(bars, values):
        ax.text(
            float(value) + pad,
            bar.get_y() + bar.get_height() / 2,
            fmt.format(float(value)),
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold",
            color=DARK,
            clip_on=False,
        )
    ax.set_xlim(ax.get_xlim()[0], max(x_max, float(values.max()) + pad * 8))


def plot_exp1_macro_f1(out_dir: Path) -> None:
    metrics = pd.read_csv(EXP1_TABLES_DIR / "experiment_1_metrics.csv").sort_values("macro_f1")
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    style_axes(ax)
    labels = [short_model(model) for model in metrics["model"]]
    bars = ax.barh(
        labels,
        metrics["macro_f1"],
        color=[model_color(model) for model in metrics["model"]],
        edgecolor="#2f2a24",
        linewidth=0.8,
        alpha=0.96,
    )
    ax.axvspan(0.98, 1.0, color="#d7ead3", alpha=0.55, zorder=0)
    ax.axvline(0.95, color="#8c8172", linestyle="--", linewidth=1.0, alpha=0.75)
    add_value_labels(ax, bars, metrics["macro_f1"], pad=0.0025)
    ax.set_xlabel("Macro-F1")
    ax.set_title("Experiment 1: splice-site classification", fontsize=15, fontweight="bold", loc="left")
    add_note(ax, "Zoomed axis highlights the real-model gap", x=0.64, y=0.965)
    ax.set_xlim(0.86, 1.006)
    ax.grid(axis="x", color=GRID, alpha=0.55, linewidth=0.8)
    ax.grid(axis="y", visible=False)
    save(fig, out_dir / "experiment_1_macro_f1.png")


def plot_exp2_context(out_dir: Path) -> None:
    metrics = pd.read_csv(EXP2_TABLES_DIR / "experiment_2A_multiscale_context.csv")
    for metric, filename, ylabel in [
        ("macro_f1", "exp2A_context_macro_f1.png", "Macro-F1"),
        ("auprc", "exp2A_context_auprc.png", "Macro AUPRC"),
    ]:
        fig, ax = plt.subplots(figsize=(8.4, 5.0))
        style_axes(ax)
        for model, group in metrics.groupby("model", sort=False):
            group = group.sort_values("window_flank")
            color = model_color(model)
            ax.plot(
                group["window_flank"],
                group[metric],
                marker="o",
                linewidth=2.6,
                markersize=6.5,
                color=color,
                label=short_model(model).replace("\n", " "),
            )
            ax.scatter(group["window_flank"], group[metric], s=70, color=color, edgecolor=PANEL, linewidth=1.0, zorder=3)
            end = group.iloc[-1]
            ax.annotate(
                short_model(model).replace("\n", " "),
                xy=(end["window_flank"], end[metric]),
                xytext=(8, 0),
                textcoords="offset points",
                va="center",
                fontsize=8.5,
                color=color,
                fontweight="bold",
            )
        if metric == "macro_f1":
            cnn = metrics[metrics["model"] == "CNN baseline (PyTorch Conv1D)"].sort_values("window_flank")
            rnabert = metrics[metrics["model"] == "RNABERT frozen encoder + MLP"].sort_values("window_flank")
            ax.fill_between(cnn["window_flank"], cnn[metric], rnabert[metric], color="#54A24B", alpha=0.10)
        ax.set_xlabel("Context flank (nt)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Experiment 2A: context ablation ({ylabel})", fontsize=15, fontweight="bold", loc="left")
        ax.set_xticks(sorted(metrics["window_flank"].unique()))
        lower = max(0.84, float(metrics[metric].min()) - 0.015)
        ax.set_ylim(lower, 1.01)
        ax.grid(axis="y", color=GRID, alpha=0.60, linewidth=0.8)
        ax.grid(axis="x", color=GRID, alpha=0.25, linewidth=0.8)
        add_note(ax, "Same split across windows; shaded band shows RNABERT-CNN gap")
        save(fig, out_dir / filename)


def plot_exp2_hard_negative(out_dir: Path) -> None:
    hard = pd.read_csv(EXP2_TABLES_DIR / "experiment_2B_hard_negative.csv").sort_values("hard_negative_fpr")
    fig, ax = plt.subplots(figsize=(8.4, 4.7))
    style_axes(ax)
    labels = [short_model(model) for model in hard["model"]]
    colors = ["#c84f3d" if model.startswith("CNN") else model_color(model) for model in hard["model"]]
    bars = ax.barh(labels, hard["hard_negative_fpr"], color=colors, edgecolor="#2f2a24", linewidth=0.8, alpha=0.96)
    for bar, (_, row) in zip(bars, hard.iterrows()):
        label = f"{row['hard_negative_fpr']:.3f}  ({int(row['hard_negative_false_positives'])}/{int(row['hard_negative_rows'])})"
        label_x = row["hard_negative_fpr"] + 0.012
        if row["model"] == "RNABERT frozen encoder + MLP":
            label_x = row["hard_negative_fpr"] + 0.035
        ax.text(
            label_x,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold",
            color=DARK,
        )
    cnn_fpr = float(hard.loc[hard["model"].str.startswith("CNN"), "hard_negative_fpr"].iloc[0])
    best = hard.iloc[0]
    reduction = 100.0 * (1.0 - float(best["hard_negative_fpr"]) / cnn_fpr)
    ax.annotate(
        f"{reduction:.0f}% lower FPR\nthan CNN",
        xy=(float(best["hard_negative_fpr"]) + 0.006, 0.0),
        xytext=(0.27, 0.42),
        arrowprops={"arrowstyle": "->", "color": "#2f6f3e", "linewidth": 1.4},
        va="center",
        fontsize=10,
        color="#2f6f3e",
        fontweight="bold",
    )
    ax.set_xlabel("Hard-negative false positive rate")
    ax.set_title("Experiment 2B: GT/AG hard-negative failure mode", fontsize=15, fontweight="bold", loc="left")
    add_note(ax, "Lower is better; labels show FPR and false positives / 67 hard negatives")
    ax.set_xlim(0.0, 0.6)
    ax.grid(axis="x", color=GRID, alpha=0.60, linewidth=0.8)
    ax.grid(axis="y", visible=False)
    save(fig, out_dir / "exp2B_hard_negative_fpr.png")


def plot_exp3_metrics(out_dir: Path) -> None:
    metrics = pd.read_csv(EXP3_TABLES_DIR / "experiment_3A_artificial_variant_metrics.csv")
    for metric, filename, label in [
        ("auroc", "exp3_variant_auroc.png", "AUROC"),
        ("auprc", "exp3_variant_auprc.png", "AUPRC"),
    ]:
        frame = metrics.sort_values(metric)
        fig, ax = plt.subplots(figsize=(8.8, 6.2))
        style_axes(ax)
        labels = [short_model(model) for model in frame["model"]]
        bars = ax.barh(
            labels,
            frame[metric],
            color=[model_color(model) for model in frame["model"]],
            edgecolor="#2f2a24",
            linewidth=0.8,
            alpha=0.96,
        )
        add_value_labels(ax, bars, frame[metric], pad=0.006)
        best = frame.iloc[-1]
        runner_up = frame.iloc[-2]
        margin = float(best[metric] - runner_up[metric])
        ax.annotate(
            f"lead +{margin:.3f}",
            xy=(float(best[metric]) + 0.028, len(labels) - 1),
            xytext=(float(best[metric]) + 0.055, len(labels) - 1 - 0.28),
            arrowprops={"arrowstyle": "->", "color": model_color(best["model"]), "linewidth": 1.3},
            va="center",
            fontsize=9,
            fontweight="bold",
            color=model_color(best["model"]),
        )
        ax.set_xlabel(label)
        ax.set_title(f"Experiment 3: artificial variant effect {label}", fontsize=15, fontweight="bold", loc="left")
        add_note(ax, "External splice tools are real sequence-model/tool runs, not fallback rows")
        lower = max(0.0, float(frame[metric].min()) - 0.06)
        upper = min(1.03, float(frame[metric].max()) + 0.12)
        ax.set_xlim(lower, upper)
        ax.grid(axis="x", color=GRID, alpha=0.60, linewidth=0.8)
        ax.grid(axis="y", visible=False)
        save(fig, out_dir / filename)


def plot_all(out_dir: Path = REPORT_IMAGE_DIR) -> None:
    ensure_dirs(out_dir)
    plot_exp1_macro_f1(out_dir)
    plot_exp2_context(out_dir)
    plot_exp2_hard_negative(out_dir)
    plot_exp3_metrics(out_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate report figures from real-model experiment outputs.")
    parser.add_argument("--out-dir", type=Path, default=REPORT_IMAGE_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plot_all(args.out_dir)
    print(f"Wrote report figures to {args.out_dir}")


if __name__ == "__main__":
    main()
