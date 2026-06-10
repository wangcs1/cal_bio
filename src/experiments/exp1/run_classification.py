from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.experiments.exp1.evaluate import evaluate_many
from src.experiments.exp1.common import parse_model_list, resolve_data_dir
from src.experiments.exp1.train import train_many
from src.utils import EXP1_FIGURES_DIR, EXP1_TABLES_DIR, LABELS, PROJECT_ROOT, ensure_dirs, write_dataframe


def plot_metrics(metrics: pd.DataFrame, out_dir: Path) -> None:
    ensure_dirs(out_dir)
    order = metrics.sort_values("macro_f1", ascending=False)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.barh(order["model"], order["macro_f1"], color="#3b6ea8")
    ax.set_xlabel("Macro-F1")
    ax.set_title("Experiment 1 splice-site classification")
    ax.set_xlim(0.0, 1.0)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "experiment_1_macro_f1.png", dpi=180)
    plt.close(fig)


def plot_confusion(confusion: pd.DataFrame, out_dir: Path) -> None:
    ensure_dirs(out_dir)
    models = confusion["model"].unique().tolist()
    fig, axes = plt.subplots(1, len(models), figsize=(4.2 * len(models), 3.8), squeeze=False)
    for ax, model in zip(axes[0], models):
        sub = confusion[confusion["model"] == model]
        matrix = sub.pivot(index="true_label", columns="pred_label", values="count").reindex(index=[0, 1, 2], columns=[0, 1, 2]).fillna(0)
        values = matrix.to_numpy(dtype=float)
        row_sums = values.sum(axis=1, keepdims=True)
        normalized = values / row_sums.clip(min=1.0)
        image = ax.imshow(normalized, cmap="Blues", vmin=0.0, vmax=1.0)
        ax.set_title(model)
        ax.set_xticks([0, 1, 2], [LABELS[idx] for idx in [0, 1, 2]], rotation=35, ha="right")
        ax.set_yticks([0, 1, 2], [LABELS[idx] for idx in [0, 1, 2]])
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        for i in range(3):
            for j in range(3):
                ax.text(j, i, str(int(values[i, j])), ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    fig.savefig(out_dir / "experiment_1_confusion_matrices.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(path: Path, data_dir: Path, train_summary: pd.DataFrame, test_metrics: pd.DataFrame) -> None:
    ensure_dirs(path.parent)
    lines = [
        "# Experiment 1: Splice-site three-class classification",
        "",
        f"Data directory: `{data_dir}`",
        "",
        "This run uses the available local Python environment. Because torch/RNA foundation model",
        "dependencies are not installed in this environment, RNA-FM and RNABERT are implemented as",
        "frozen-representation style k-mer/token proxies with lightweight classifier heads.",
        "",
        "## Validation summary",
        "",
        train_summary[["model", "train_rows", "valid_rows", "macro_f1", "accuracy", "auroc", "auprc"]].to_markdown(index=False),
        "",
        "## Test summary",
        "",
        test_metrics[["model", "rows", "accuracy", "macro_f1", "auroc", "auprc", "donor_f1", "acceptor_f1", "non_splice_f1"]].to_markdown(index=False),
        "",
        "Outputs:",
        "",
        "- `results/experiment_1/tables/experiment_1_metrics.csv`",
        "- `results/experiment_1/tables/experiment_1_confusion_matrices.csv`",
        "- `results/experiment_1/figures/experiment_1_macro_f1.png`",
        "- `results/experiment_1/figures/experiment_1_confusion_matrices.png`",
        "- `results/checkpoints/experiment_1/*.joblib`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(
    data_dir: Path | None,
    tables_dir: Path,
    figures_dir: Path,
    checkpoint_dir: Path,
    report_path: Path,
    model_keys: list[str],
    seed: int,
    max_train_rows: int | None,
    max_valid_rows: int | None,
    max_test_rows: int | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    resolved_data_dir = resolve_data_dir(data_dir)
    ensure_dirs(tables_dir, figures_dir, checkpoint_dir)

    train_summary = train_many(model_keys, resolved_data_dir, checkpoint_dir, seed, max_train_rows, max_valid_rows)
    test_metrics, confusion = evaluate_many(model_keys, checkpoint_dir, resolved_data_dir, "test", max_test_rows, seed)

    write_dataframe(tables_dir / "experiment_1_train_valid_metrics.csv", train_summary)
    write_dataframe(tables_dir / "experiment_1_metrics.csv", test_metrics)
    write_dataframe(tables_dir / "experiment_1_confusion_matrices.csv", confusion)
    plot_metrics(test_metrics, figures_dir)
    plot_confusion(confusion, figures_dir)
    write_report(report_path, resolved_data_dir, train_summary, test_metrics)
    return test_metrics, confusion


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run B-part experiment 1 classification end to end.")
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--tables-dir", type=Path, default=EXP1_TABLES_DIR)
    parser.add_argument("--figures-dir", type=Path, default=EXP1_FIGURES_DIR)
    parser.add_argument("--checkpoint-dir", type=Path, default=PROJECT_ROOT / "results/checkpoints/experiment_1")
    parser.add_argument("--report-path", type=Path, default=PROJECT_ROOT / "reports/experiment_1.md")
    parser.add_argument("--models", default="cnn,rnafm,rnabert")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-train-rows", type=int, default=60000)
    parser.add_argument("--max-valid-rows", type=int, default=15000)
    parser.add_argument("--max-test-rows", type=int, default=30000)
    parser.add_argument("--full-data", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model_keys = parse_model_list(args.models)
    max_train_rows = None if args.full_data else args.max_train_rows
    max_valid_rows = None if args.full_data else args.max_valid_rows
    max_test_rows = None if args.full_data else args.max_test_rows
    metrics, _ = run(
        data_dir=args.data_dir,
        tables_dir=args.tables_dir,
        figures_dir=args.figures_dir,
        checkpoint_dir=args.checkpoint_dir,
        report_path=args.report_path,
        model_keys=model_keys,
        seed=args.seed,
        max_train_rows=max_train_rows,
        max_valid_rows=max_valid_rows,
        max_test_rows=max_test_rows,
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
