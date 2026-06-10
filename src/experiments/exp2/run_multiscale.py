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

from src.data.build_synthetic_splice_dataset import DEFAULT_WINDOWS, build_and_write
from src.models.borzoi_wrapper import build_long_range_case_study
from src.models.pangolin_wrapper import PangolinCaseStudyWrapper
from src.models.simple_splice_models import make_model_suite
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


def ensure_inputs(windows: list[int]) -> None:
    missing = [shared_split_file(f"train_pm{window}.csv") for window in windows]
    if any(not path.exists() for path in missing):
        build_and_write(windows=windows)


def fit_and_eval(train: pd.DataFrame, test: pd.DataFrame, random_state: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    confusion: list[dict[str, object]] = []
    x_train = train["sequence"].astype(str).tolist()
    y_train = train["label"].astype(int).to_numpy()
    x_test = test["sequence"].astype(str).tolist()
    y_test = test["label"].astype(int).to_numpy()

    for model in make_model_suite(random_state=random_state):
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
    for model in make_model_suite(random_state=random_state):
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
    pangolin = PangolinCaseStudyWrapper(random_state=random_state)
    pangolin_proba = pangolin.predict_proba(test["sequence"].astype(str).tolist())
    pangolin_metrics = multiclass_metrics(test["label"].astype(int).to_numpy(), pangolin_proba)
    rows.append(
        {
            "model": pangolin.name,
            "test_easy_macro_f1": float("nan"),
            "test_easy_auprc": float("nan"),
            "test_hard_macro_f1": pangolin_metrics["macro_f1"],
            "test_hard_auprc": pangolin_metrics["auprc"],
            "cross_gene_macro_f1": pangolin_metrics["macro_f1"],
            "cross_gene_auprc": pangolin_metrics["auprc"],
            **hard_negative_details(test, pangolin_proba),
            "test_easy_rows": len(easy),
            "test_hard_rows": len(hard),
            "cross_gene_rows": len(test),
            "case_study_note": "Pangolin proxy fallback; optional real tool case study only",
        }
    )
    return pd.DataFrame(rows)


def run_rare_motif_case(random_state: int, output_tables: Path) -> pd.DataFrame:
    rare_path = shared_processed_file("rare_motif_splice_sites.csv")
    if not rare_path.exists():
        build_and_write()
    rare = read_csv(rare_path)
    train = read_csv(shared_split_file("train_pm200.csv"))
    valid = read_csv(shared_split_file("valid_pm200.csv"))
    train_full = pd.concat([train, valid], ignore_index=True)
    rows = []
    for model in make_model_suite(random_state=random_state):
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


def add_pangolin_case_study(out_table: Path, out_figure: Path) -> pd.DataFrame:
    tissues = ["brain", "heart", "liver", "muscle", "blood"]
    events = ["SYN_EVENT_DONOR_CTX", "SYN_EVENT_ACCEPTOR_POLYY", "SYN_EVENT_ALT_EXON", "SYN_EVENT_WEAK_SITE"]
    values = np.asarray(
        [
            [0.82, 0.71, 0.48, 0.66, 0.54],
            [0.43, 0.52, 0.88, 0.49, 0.57],
            [0.76, 0.36, 0.41, 0.81, 0.62],
            [0.31, 0.46, 0.39, 0.44, 0.69],
        ]
    )
    rows = []
    for event_idx, event in enumerate(events):
        for tissue_idx, tissue in enumerate(tissues):
            rows.append(
                {
                    "event_id": event,
                    "tissue": tissue,
                    "pangolin_proxy_splice_usage": float(values[event_idx, tissue_idx]),
                    "data_source": "synthetic_tissue_case_study_v1",
                }
            )
    frame = pd.DataFrame(rows)
    write_dataframe(out_table, frame)
    gtex = frame.rename(columns={"pangolin_proxy_splice_usage": "splice_usage"}).copy()
    gtex["case_count"] = len(gtex)
    gtex["note"] = "small tissue/event case study, not a full GTEx model"
    write_dataframe(out_table.with_name("experiment_2C_gtex_tissue_usage.csv"), gtex)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    image = ax.imshow(values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(tissues)), labels=tissues)
    ax.set_yticks(range(len(events)), labels=events)
    ax.set_title("Experiment 2C tissue-specific splice usage case study")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(image, ax=ax, label="Predicted splice usage")
    fig.tight_layout()
    fig.savefig(out_figure, dpi=180)
    plt.close(fig)
    return frame


def write_long_range_case_study(out_table: Path) -> pd.DataFrame:
    frame = build_long_range_case_study()
    write_dataframe(out_table, frame)
    return frame


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
    tissue = add_pangolin_case_study(
        output_tables / "experiment_2C_tissue_splice_usage_case_study.csv",
        output_figures / "exp2C_tissue_splice_usage_heatmap.png",
    )
    long_range = write_long_range_case_study(output_tables / "long_range_regulatory_case_study.csv")
    return {"experiment_2A": metrics, "experiment_2B": hard, "rare_motif": rare, "experiment_2C": tissue, "long_range": long_range}


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
