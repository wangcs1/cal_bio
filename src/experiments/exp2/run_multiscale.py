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

from src.data.build_synthetic_splice_dataset import DEFAULT_WINDOWS
from src.data.build_splice_site_dataset import build_and_write_real
from src.experiments.exp1.common import make_model
from src.utils import (
    CONFIG_ROOT,
    EXP2_FIGURES_DIR,
    EXP2_TABLES_DIR,
    LABELS,
    PROJECT_ROOT,
    SHARED_SPLIT_DIR,
    confusion_rows,
    ensure_dirs,
    hard_negative_fpr,
    hard_negative_details,
    load_config,
    multiclass_metrics,
    read_csv,
    shared_processed_file,
    shared_split_file,
    write_dataframe,
)


REAL_MODEL_KEYS = ["cnn", "rnafm", "rnabert"]


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


def write_report(path: Path, metrics: pd.DataFrame, hard: pd.DataFrame, rare: pd.DataFrame) -> None:
    ensure_dirs(path.parent)
    lines = [
        "# Experiment 2: Context Scale And Hard Negatives",
        "",
        "This run uses only the implemented real-model set: CNN baseline, RNA-FM frozen encoder, and RNABERT frozen encoder.",
        "RNA-FM/RNABERT require local pretrained weights under `models/hf/`; the report contains only real-model rows.",
        "",
        "## 2A Multi-Scale Context",
        "",
        markdown_table(
            metrics,
            [
                "window_flank",
                "sequence_length",
                "model",
                "accuracy",
                "macro_f1",
                "auroc",
                "auprc",
                "hard_negative_fpr",
            ],
        ),
        "",
        "## 2B Hard-Negative Stress Test",
        "",
        markdown_table(
            hard,
            [
                "model",
                "test_easy_macro_f1",
                "test_hard_macro_f1",
                "cross_gene_macro_f1",
                "hard_negative_fpr",
                "hard_negative_false_positives",
                "hard_negative_rows",
            ],
        ),
        "",
        "## Rare Motif Stress Test",
        "",
        "The rare-motif table is a synthetic stress test, not a claim about population-scale rare splice-site recall.",
        "",
        markdown_table(
            rare,
            ["model", "motif_type", "rows", "mean_target_probability", "accuracy", "macro_f1"],
        ),
        "",
        "Outputs:",
        "",
        "- `results/experiment_2/tables/experiment_2A_multiscale_context.csv`",
        "- `results/experiment_2/tables/experiment_2B_hard_negative.csv`",
        "- `results/experiment_2/tables/experiment_2B_rare_motif.csv`",
        "- `results/experiment_2/figures/exp2A_context_macro_f1.png`",
        "- `results/experiment_2/figures/exp2A_context_auprc.png`",
        "- `results/experiment_2/figures/exp2B_hard_negative_fpr.png`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def ensure_inputs(windows: list[int]) -> None:
    missing = [shared_split_file(f"train_pm{window}.csv") for window in windows]
    if any(not path.exists() for path in missing):
        build_and_write_real(windows=windows)


def fit_and_eval(train: pd.DataFrame, test: pd.DataFrame, random_state: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    confusion: list[dict[str, object]] = []
    x_train = train["sequence"].astype(str).tolist()
    y_train = train["label"].astype(int).to_numpy()
    x_test = test["sequence"].astype(str).tolist()
    y_test = test["label"].astype(int).to_numpy()

    for model_key in REAL_MODEL_KEYS:
        model = make_model(model_key, random_state)
        model.fit(x_train, y_train)
        proba = model.predict_proba(x_test)
        metrics = multiclass_metrics(y_test, proba)
        metrics["hard_negative_fpr"] = hard_negative_fpr(test, proba)
        rows.append({"model": model.name, **metrics})
        confusion.extend(confusion_rows(y_test, proba, model.name, "test"))
    return pd.DataFrame(rows), pd.DataFrame(confusion)


def run_multiscale(windows: list[int], random_state: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    ensure_inputs(windows)
    metric_rows: list[pd.DataFrame] = []
    confusion_rows_all: list[pd.DataFrame] = []

    for window in windows:
        train = read_csv(shared_split_file(f"train_pm{window}.csv"))
        valid = read_csv(shared_split_file(f"valid_pm{window}.csv"))
        test = read_csv(shared_split_file(f"test_pm{window}.csv"))
        train_full = pd.concat([train, valid], ignore_index=True)
        metrics, confusion = fit_and_eval(train_full, test, random_state=random_state)
        metrics.insert(0, "window_flank", window)
        metrics.insert(1, "sequence_length", 2 * window + 1)
        confusion.insert(0, "window_flank", window)
        metric_rows.append(metrics)
        confusion_rows_all.append(confusion)

    return pd.concat(metric_rows, ignore_index=True), pd.concat(confusion_rows_all, ignore_index=True)


def evaluate_subset(model, frame: pd.DataFrame) -> dict[str, float]:
    proba = model.predict_proba(frame["sequence"].astype(str).tolist())
    metrics = multiclass_metrics(frame["label"].astype(int).to_numpy(), proba)
    metrics["hard_negative_fpr"] = hard_negative_fpr(frame, proba)
    return metrics


def run_hard_negative(random_state: int, window: int = 200) -> pd.DataFrame:
    ensure_inputs([window])
    train = read_csv(shared_split_file(f"train_pm{window}.csv"))
    valid = read_csv(shared_split_file(f"valid_pm{window}.csv"))
    test = read_csv(shared_split_file(f"test_pm{window}.csv"))
    train_full = pd.concat([train, valid], ignore_index=True)
    positives = test[test["label"].astype(int).isin([0, 1])]
    easy = pd.concat(
        [positives, test[(test["label"].astype(int) == 2) & (test["negative_type"].astype(str) == "easy_random")]],
        ignore_index=True,
    )
    hard = pd.concat(
        [positives, test[(test["label"].astype(int) == 2) & (test["negative_type"].astype(str).str.contains("hard"))]],
        ignore_index=True,
    )

    rows: list[dict[str, object]] = []
    for model_key in REAL_MODEL_KEYS:
        model = make_model(model_key, random_state)
        model.fit(train_full["sequence"].astype(str).tolist(), train_full["label"].astype(int).to_numpy())
        easy_metrics = evaluate_subset(model, easy)
        hard_metrics = evaluate_subset(model, hard)
        cross_gene_metrics = evaluate_subset(model, test)
        rows.append(
            {
                "model": model.name,
                "test_easy_macro_f1": easy_metrics["macro_f1"],
                "test_easy_auprc": easy_metrics["auprc"],
                "test_hard_macro_f1": hard_metrics["macro_f1"],
                "test_hard_auprc": hard_metrics["auprc"],
                "cross_gene_macro_f1": cross_gene_metrics["macro_f1"],
                "cross_gene_auprc": cross_gene_metrics["auprc"],
                **hard_negative_details(test, model.predict_proba(test["sequence"].astype(str).tolist())),
                "test_easy_rows": len(easy),
                "test_hard_rows": len(hard),
                "cross_gene_rows": len(test),
            }
        )
    return pd.DataFrame(rows)


def run_rare_motif_case(random_state: int, output_tables: Path) -> pd.DataFrame:
    rare_path = shared_processed_file("rare_motif_splice_sites.csv")
    if not rare_path.exists():
        raise FileNotFoundError(
            f"Rare-motif stress-test table is required by the current report but was not found: {rare_path}. "
            "Regenerate processed data before running experiment 2; the benchmark does not write placeholder not_run rows."
        )
    rare = read_csv(rare_path)
    train = read_csv(shared_split_file("train_pm200.csv"))
    valid = read_csv(shared_split_file("valid_pm200.csv"))
    train_full = pd.concat([train, valid], ignore_index=True)
    rows = []
    for model_key in REAL_MODEL_KEYS:
        model = make_model(model_key, random_state)
        model.fit(train_full["sequence"].astype(str).tolist(), train_full["label"].astype(int).to_numpy())
        for motif_type, group in rare.groupby("motif_type"):
            proba = model.predict_proba(group["sequence"].astype(str).tolist())
            metrics = multiclass_metrics(group["label"].astype(int).to_numpy(), proba)
            target_prob = [float(proba[i, int(label)]) for i, label in enumerate(group["label"].astype(int).tolist())]
            rows.append(
                {
                    "model": model.name,
                    "motif_type": motif_type,
                    "rows": len(group),
                    "mean_target_probability": float(np.mean(target_prob)),
                    "accuracy": metrics["accuracy"],
                    "macro_f1": metrics["macro_f1"],
                    "data_source": "synthetic rare-motif small case study",
                }
            )
    result = pd.DataFrame(rows)
    write_dataframe(output_tables / "experiment_2B_rare_motif.csv", result)
    return result


def plot_multiscale(metrics: pd.DataFrame, out_dir: Path) -> None:
    ensure_dirs(out_dir)
    for metric, filename, ylabel in [
        ("macro_f1", "exp2A_context_macro_f1.png", "Macro-F1"),
        ("auprc", "exp2A_context_auprc.png", "Macro AUPRC"),
    ]:
        fig, ax = plt.subplots(figsize=(8, 5))
        for model, group in metrics.groupby("model"):
            group = group.sort_values("window_flank")
            ax.plot(group["window_flank"], group[metric], marker="o", linewidth=2, label=model)
        ax.set_xlabel("Context flank (nt)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Experiment 2A context ablation: {ylabel}")
        ax.set_xticks(sorted(metrics["window_flank"].unique()))
        ax.set_ylim(0.0, 1.03)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)


def plot_hard_negative(hard_metrics: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    order = hard_metrics.sort_values("hard_negative_fpr")
    ax.barh(order["model"], order["hard_negative_fpr"], color="#c45a44")
    ax.set_xlabel("Hard-negative false positive rate")
    ax.set_title("Experiment 2B GT/AG hard-negative FPR")
    ax.set_xlim(0.0, max(1.0, float(order["hard_negative_fpr"].max()) * 1.15))
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "exp2B_hard_negative_fpr.png", dpi=180)
    plt.close(fig)


def run(output_tables: Path, output_figures: Path, random_state: int = 42, windows: list[int] | None = None) -> dict[str, pd.DataFrame]:
    ensure_dirs(output_tables, output_figures)
    windows = DEFAULT_WINDOWS if windows is None else windows
    metrics, confusion = run_multiscale(windows, random_state=random_state)
    hard = run_hard_negative(random_state=random_state, window=200)
    rare = run_rare_motif_case(random_state=random_state, output_tables=output_tables)
    write_dataframe(output_tables / "experiment_2A_multiscale_context.csv", metrics)
    pivot = metrics.pivot_table(index="model", columns="window_flank", values="macro_f1").reset_index()
    write_dataframe(output_tables / "experiment_2A_multiscale_context_pivot_macro_f1.csv", pivot)
    write_dataframe(output_tables / "experiment_2A_confusion_matrices.csv", confusion)
    write_dataframe(output_tables / "experiment_2B_hard_negative.csv", hard)
    plot_multiscale(metrics, output_figures)
    plot_hard_negative(hard, output_figures)
    write_report(PROJECT_ROOT / "reports/experiment_2.md", metrics, hard, rare)
    return {"experiment_2A": metrics, "experiment_2B": hard, "rare_motif": rare}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run C-part experiment 2: context ablation and hard negatives.")
    parser.add_argument("--config", type=Path, default=CONFIG_ROOT / "exp2_multiscale.yaml")
    parser.add_argument("--tables", type=Path, default=EXP2_TABLES_DIR)
    parser.add_argument("--figures", type=Path, default=EXP2_FIGURES_DIR)
    parser.add_argument("--windows", type=int, nargs="+", default=DEFAULT_WINDOWS)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config) if args.config else {}
    windows = [int(item) for item in config.get("windows", args.windows)]
    if windows != DEFAULT_WINDOWS:
        print(f"Using configured windows: {windows}")
    outputs = run(
        Path(config.get("tables_dir", args.tables)),
        Path(config.get("figures_dir", args.figures)),
        int(config.get("seed", args.seed)),
        windows=windows,
    )
    print(outputs["experiment_2A"][["window_flank", "model", "macro_f1", "auprc"]].to_string(index=False))
    print(outputs["experiment_2B"].to_string(index=False))


if __name__ == "__main__":
    main()
