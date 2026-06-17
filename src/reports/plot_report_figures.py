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
    x_min = max(0.0, float(metrics["macro_f1"].min()) - 0.035)
    x_max = min(1.006, float(metrics["macro_f1"].max()) + 0.035)
    ax.axvspan(max(0.0, x_max - 0.02), x_max, color="#d7ead3", alpha=0.42, zorder=0)
    ax.axvline(float(metrics["macro_f1"].median()), color="#8c8172", linestyle="--", linewidth=1.0, alpha=0.75)
    add_value_labels(ax, bars, metrics["macro_f1"], pad=0.0025)
    ax.set_xlabel("Macro-F1")
    ax.set_title("Experiment 1: splice-site classification", fontsize=15, fontweight="bold", loc="left")
    add_note(ax, "Zoomed axis is computed from the current real-data scores", x=0.015, y=0.965)
    ax.set_xlim(x_min, x_max)
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
            note = "Same split across windows; shaded band shows RNABERT-CNN gap"
        else:
            note = "Same split across windows; labels mark each model at pm400"
        ax.set_xlabel("Context flank (nt)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Experiment 2A: context ablation ({ylabel})", fontsize=15, fontweight="bold", loc="left")
        ax.set_xticks(sorted(metrics["window_flank"].unique()))
        lower = max(0.0, float(metrics[metric].min()) - 0.025)
        upper = min(1.01, float(metrics[metric].max()) + 0.035)
        ax.set_ylim(lower, upper)
        ax.grid(axis="y", color=GRID, alpha=0.60, linewidth=0.8)
        ax.grid(axis="x", color=GRID, alpha=0.25, linewidth=0.8)
        add_note(ax, note)
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
    add_note(ax, "Lower is better; labels show FPR and false positives / 155 hard negatives")
    ax.set_xlim(0.0, 0.6)
    ax.grid(axis="x", color=GRID, alpha=0.60, linewidth=0.8)
    ax.grid(axis="y", visible=False)
    save(fig, out_dir / "exp2B_hard_negative_fpr.png")


def plot_exp3_metrics(out_dir: Path) -> None:
    metrics = pd.read_csv(EXP3_TABLES_DIR / "experiment_3A_variant_metrics.csv")
    for metric, filename, label in [
        ("auroc", "exp3_variant_auroc.png", "AUROC"),
        ("auprc", "exp3_variant_auprc.png", "AUPRC"),
    ]:
        baseline = 0.5
        frame = metrics.sort_values(metric).copy()
        frame["gain"] = frame[metric] - baseline
        fig, ax = plt.subplots(figsize=(8.8, 6.2))
        style_axes(ax)
        labels = [short_model(model) for model in frame["model"]]
        bars = ax.barh(
            labels,
            frame["gain"],
            left=baseline,
            color=[model_color(model) for model in frame["model"]],
            edgecolor="#2f2a24",
            linewidth=0.9,
            alpha=0.96,
        )
        for bar, (_, row) in zip(bars, frame.iterrows()):
            ax.text(
                float(row[metric]) + 0.006,
                bar.get_y() + bar.get_height() / 2,
                f"{row[metric]:.3f}  (+{row['gain']:.3f})",
                va="center",
                ha="left",
                fontsize=9,
                fontweight="bold",
                color=DARK,
                clip_on=False,
            )
        best = frame.iloc[-1]
        runner_up = frame.iloc[-2]
        margin = float(best[metric] - runner_up[metric])
        ax.annotate(
            f"lead +{margin:.3f}",
            xy=(float(best[metric]) + 0.012, len(labels) - 1),
            xytext=(float(best[metric]) + 0.055, len(labels) - 1 - 0.35),
            arrowprops={"arrowstyle": "->", "color": model_color(best["model"]), "linewidth": 1.3},
            va="center",
            fontsize=9,
            fontweight="bold",
            color=model_color(best["model"]),
        )
        ax.axvline(baseline, color="#E45756", linestyle="--", linewidth=1.7)
        ax.axvspan(0.0, baseline, color="#efe7dc", alpha=0.55, zorder=0)
        ax.text(
            baseline + 0.004,
            -0.58,
            "random baseline",
            color="#B6463C",
            fontsize=8.5,
            fontweight="bold",
            va="center",
        )
        ax.set_xlabel(f"{label} (bar length shows gain above random baseline)")
        ax.set_title(f"Experiment 3: real ClinVar variant effect {label}", fontsize=15, fontweight="bold", loc="left")
        add_note(ax, "Scores are sorted by full ClinVar performance; labels show absolute score and gain over 0.5")
        upper = min(1.03, float(frame[metric].max()) + 0.10)
        ax.set_xlim(0.48, upper)
        ax.grid(axis="x", color=GRID, alpha=0.60, linewidth=0.8)
        ax.grid(axis="y", visible=False)
        save(fig, out_dir / filename)


def plot_exp3_calibration(out_dir: Path) -> None:
    bins = pd.read_csv(EXP3_TABLES_DIR / "experiment_3A_calibration_bins.csv").sort_values("bin")
    fig, ax = plt.subplots(figsize=(8.8, 5.1))
    style_axes(ax)
    x = bins["bin"].astype(int).to_numpy()
    rates = bins["positive_rate"].to_numpy()
    rows = bins["rows"].to_numpy()
    ax2 = ax.twinx()
    ax2.set_facecolor("none")
    ax2.bar(x, rows, width=0.72, color="#d8cbb9", alpha=0.42, edgecolor="none", label="Rows per bin")
    ax2.set_ylim(0, max(rows) * 1.65)
    ax2.set_ylabel("Rows per bin", color=MUTED)
    ax2.tick_params(colors=MUTED, labelsize=8)
    for side in ["top", "right"]:
        ax2.spines[side].set_visible(False)
    ax.axhspan(0.5, 1.0, color="#d7ead3", alpha=0.20, zorder=0)
    ax.axhline(0.5, color="#E45756", linestyle="--", linewidth=1.5, label="Random baseline = 0.5")
    ax.plot(x, rates, color=MODEL_COLORS["RNA-FM frozen encoder + MLP"], linewidth=3.0, marker="o", markersize=7.5)
    ax.scatter(x, rates, s=95, color=MODEL_COLORS["RNA-FM frozen encoder + MLP"], edgecolor=PANEL, linewidth=1.2, zorder=4)
    for xi, rate, n in zip(x, rates, rows):
        ax.text(xi, rate + 0.028, f"{rate:.2f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=DARK)
        ax.text(xi, 0.225, f"n={int(n)}", ha="center", va="center", fontsize=7.8, color=MUTED)
    ax.set_xticks(x)
    ax.set_xlabel("Impact-score bin (low to high)")
    ax.set_ylabel("Observed splice-altering rate")
    ax.set_ylim(0.20, 1.02)
    ax.set_title("Experiment 3: ClinVar response by delta-score bins", fontsize=15, fontweight="bold", loc="left")
    add_note(ax, "RNA-FM score bins; higher bins should contain more splice-altering variants")
    ax.grid(axis="y", color=GRID, alpha=0.62, linewidth=0.8)
    ax.grid(axis="x", color=GRID, alpha=0.22, linewidth=0.8)
    save(fig, out_dir / "exp3_calibration_curve.png")


def plot_exp3_variant_type_summary(out_dir: Path) -> None:
    summary = pd.read_csv(EXP3_TABLES_DIR / "variant_effect_stratified_by_type.csv")
    label_map = {
        "acceptor_clinvar_splice": "Acceptor\nsplice",
        "donor_clinvar_splice": "Donor\nsplice",
        "clinvar_benign_snv": "Benign\nSNV",
    }
    order = [
        "RNA-FM frozen encoder + MLP",
        "Pangolin real sequence model",
        "CNN baseline (PyTorch Conv1D)",
        "RNABERT frozen encoder + MLP",
        "SpliceAI real sequence model",
        "MaxEntScan real local score",
        "MMSplice real sequence model",
    ]
    classes = ["acceptor_clinvar_splice", "donor_clinvar_splice", "clinvar_benign_snv"]
    pivot = (
        summary.pivot_table(index="model", columns="variant_type", values="mean_score", aggfunc="mean")
        .reindex(order)
        .reindex(columns=classes)
        .fillna(0.0)
    )
    separation = ((pivot["acceptor_clinvar_splice"] + pivot["donor_clinvar_splice"]) / 2.0 - pivot["clinvar_benign_snv"]).sort_values()
    fig, (ax_heat, ax_sep) = plt.subplots(
        1,
        2,
        figsize=(13.2, 5.9),
        gridspec_kw={"width_ratios": [1.32, 1.0], "wspace": 0.48},
    )
    style_axes(ax_heat)
    style_axes(ax_sep)
    matrix = pivot.to_numpy()
    image = ax_heat.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0.0, vmax=max(0.31, float(matrix.max())))
    ax_heat.set_yticks(np.arange(len(pivot.index)))
    ax_heat.set_yticklabels([short_model(model).replace("\n", " ") for model in pivot.index], fontsize=8.5)
    ax_heat.set_xticks(np.arange(len(classes)))
    ax_heat.set_xticklabels([label_map[name] for name in classes], fontsize=9, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            color = "white" if value > 0.18 else DARK
            ax_heat.text(j, i, f"{value:.2f}", ha="center", va="center", color=color, fontsize=8.5, fontweight="bold")
    ax_heat.set_title("Mean impact score by ClinVar class", fontsize=14, fontweight="bold", loc="left")
    add_note(ax_heat, "Darker cells indicate stronger average model response", x=0.02, y=1.08)
    cbar = fig.colorbar(image, ax=ax_heat, fraction=0.046, pad=0.035)
    cbar.set_label("Mean impact score", color=MUTED)
    cbar.ax.tick_params(colors=MUTED, labelsize=8)
    labels = [short_model(model) for model in separation.index]
    colors = [model_color(model) for model in separation.index]
    bars = ax_sep.barh(labels, separation.values, color=colors, edgecolor="#2f2a24", linewidth=0.8)
    ax_sep.tick_params(axis="y", labelsize=8.2)
    for bar, value in zip(bars, separation.values):
        ax_sep.text(value + 0.006, bar.get_y() + bar.get_height() / 2, f"+{value:.2f}", va="center", ha="left", fontsize=8.5, fontweight="bold")
    ax_sep.axvline(0.0, color="#817768", linewidth=1.0)
    ax_sep.set_xlabel("Mean positive score minus benign score")
    ax_sep.set_title("Class separation", fontsize=14, fontweight="bold", loc="left")
    add_note(ax_sep, "Higher means clearer splice-vs-benign separation", x=0.02, y=0.97)
    ax_sep.set_xlim(0.0, max(0.32, float(separation.max()) + 0.05))
    ax_sep.grid(axis="x", color=GRID, alpha=0.60, linewidth=0.8)
    ax_sep.grid(axis="y", visible=False)
    fig.suptitle("Experiment 3: ClinVar variant-type response profile", x=0.02, y=1.03, ha="left", fontsize=16, fontweight="bold", color=DARK)
    save(fig, out_dir / "variant_effect_stratified_by_type.png")


def plot_all(out_dir: Path = REPORT_IMAGE_DIR) -> None:
    ensure_dirs(out_dir)
    plot_exp1_macro_f1(out_dir)
    plot_exp2_context(out_dir)
    plot_exp2_hard_negative(out_dir)
    plot_exp3_metrics(out_dir)
    plot_exp3_calibration(out_dir)
    plot_exp3_variant_type_summary(out_dir)


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
