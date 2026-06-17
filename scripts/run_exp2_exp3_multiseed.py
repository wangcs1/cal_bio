from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.experiments.exp2.run_multiscale import run_hard_negative, run_multiscale
from src.experiments.exp1.common import make_model
from src.experiments.exp3.run_variant_effect import score_variants, summarize_distance_matched_metrics, summarize_metrics
from src.utils import EXP2_TABLES_DIR, EXP3_DATA_DIR, EXP3_TABLES_DIR, ensure_dirs, read_csv, shared_split_file, write_dataframe


METRIC_COLUMNS_EXP2 = ["accuracy", "macro_f1", "auroc", "auprc", "hard_negative_fpr"]
METRIC_COLUMNS_EXP3 = ["auroc", "auprc", "top_k_recall", "enrichment_at_k"]


def summarize(frame: pd.DataFrame, group_cols: list[str], metric_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for keys, group in frame.groupby(group_cols, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {column: value for column, value in zip(group_cols, keys)}
        row["seeds"] = ",".join(str(seed) for seed in sorted(group["seed"].astype(int).unique()))
        row["runs"] = len(group)
        for metric in metric_cols:
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_std"] = float(group[metric].std(ddof=1)) if len(group) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def run_exp2_multiseed(seeds: list[int], windows: list[int], tables_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    runs: list[pd.DataFrame] = []
    hard_runs: list[pd.DataFrame] = []
    for seed in seeds:
        metrics, _confusion = run_multiscale(windows, random_state=seed)
        hard = run_hard_negative(random_state=seed, window=200)
        metrics.insert(0, "seed", seed)
        hard.insert(0, "seed", seed)
        runs.append(metrics)
        hard_runs.append(hard)
    all_metrics = pd.concat(runs, ignore_index=True)
    all_hard = pd.concat(hard_runs, ignore_index=True)
    summary = summarize(all_metrics, ["window_flank", "model"], METRIC_COLUMNS_EXP2)
    hard_summary = summarize(all_hard, ["model"], ["test_hard_macro_f1", "cross_gene_macro_f1", "hard_negative_fpr"])
    write_dataframe(tables_dir / "experiment_2_multiseed_metrics.csv", all_metrics)
    write_dataframe(tables_dir / "experiment_2_multiseed_summary.csv", summary)
    write_dataframe(tables_dir / "experiment_2_multiseed_hard_negative_metrics.csv", all_hard)
    write_dataframe(tables_dir / "experiment_2_multiseed_hard_negative_summary.csv", hard_summary)
    return summary, hard_summary


def run_exp3_multiseed(seeds: list[int], tables_dir: Path, _figures_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    variants = read_csv(EXP3_DATA_DIR / "clinvar_splicing_variants.csv")
    train = read_csv(shared_split_file("train_pm200.csv"))
    valid = read_csv(shared_split_file("valid_pm200.csv"))
    train_full = pd.concat([train, valid], ignore_index=True)
    x_train = train_full["sequence"].astype(str).tolist()
    y_train = train_full["label"].astype(int).to_numpy()
    metric_runs: list[pd.DataFrame] = []
    matched_runs: list[pd.DataFrame] = []
    for seed in seeds:
        models = [make_model(model_key, seed) for model_key in ["cnn", "rnafm", "rnabert"]]
        for model in models:
            model.fit(x_train, y_train)
        scores = score_variants(models, variants)
        metrics = summarize_metrics(scores)
        metrics.insert(0, "seed", seed)
        metric_runs.append(metrics)
        matched = summarize_distance_matched_metrics(scores)
        matched.insert(0, "seed", seed)
        matched_runs.append(matched)
    all_metrics = pd.concat(metric_runs, ignore_index=True)
    summary = summarize(all_metrics, ["model"], METRIC_COLUMNS_EXP3)
    write_dataframe(tables_dir / "experiment_3_multiseed_metrics.csv", all_metrics)
    write_dataframe(tables_dir / "experiment_3_multiseed_summary.csv", summary)
    if matched_runs:
        all_matched = pd.concat(matched_runs, ignore_index=True)
        matched_summary = summarize(all_matched, ["model"], METRIC_COLUMNS_EXP3)
        write_dataframe(tables_dir / "experiment_3_multiseed_distance_matched_metrics.csv", all_matched)
        write_dataframe(tables_dir / "experiment_3_multiseed_distance_matched_summary.csv", matched_summary)
    else:
        matched_summary = pd.DataFrame()
    return summary, matched_summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run experiments 2 and 3 across multiple seeds.")
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--windows", default="50,100,200,400")
    parser.add_argument("--skip-exp2", action="store_true")
    parser.add_argument("--skip-exp3", action="store_true")
    parser.add_argument("--tables-exp2", type=Path, default=EXP2_TABLES_DIR)
    parser.add_argument("--tables-exp3", type=Path, default=EXP3_TABLES_DIR)
    parser.add_argument("--figures-exp3", type=Path, default=PROJECT_ROOT / "results/experiment_3/figures")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    seeds = [int(item.strip()) for item in args.seeds.split(",") if item.strip()]
    windows = [int(item.strip()) for item in args.windows.split(",") if item.strip()]
    ensure_dirs(args.tables_exp2, args.tables_exp3, args.figures_exp3)
    if not args.skip_exp2:
        exp2_summary, hard_summary = run_exp2_multiseed(seeds, windows, args.tables_exp2)
        print("Experiment 2 multi-seed summary:")
        print(exp2_summary.to_string(index=False))
        print("Experiment 2 hard-negative multi-seed summary:")
        print(hard_summary.to_string(index=False))
    if not args.skip_exp3:
        exp3_summary, matched_summary = run_exp3_multiseed(seeds, args.tables_exp3, args.figures_exp3)
        print("Experiment 3 multi-seed summary:")
        print(exp3_summary.to_string(index=False))
        if not matched_summary.empty:
            print("Experiment 3 distance-matched multi-seed summary:")
            print(matched_summary.to_string(index=False))


if __name__ == "__main__":
    main()
