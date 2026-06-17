from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.experiments.exp1.common import parse_model_list, resolve_data_dir
from src.experiments.exp1.evaluate import evaluate_many
from src.experiments.exp1.train import train_many
from src.utils import EXP1_TABLES_DIR, ensure_dirs, write_dataframe


METRIC_COLUMNS = ["accuracy", "macro_f1", "auroc", "auprc", "hard_negative_fpr"]


def summarize(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model_key, model), group in metrics.groupby(["model_key", "model"], sort=False):
        row: dict[str, object] = {
            "model_key": model_key,
            "model": model,
            "seeds": ",".join(str(seed) for seed in sorted(group["seed"].astype(int).unique())),
            "runs": len(group),
        }
        for column in METRIC_COLUMNS:
            row[f"{column}_mean"] = float(group[column].mean())
            row[f"{column}_std"] = float(group[column].std(ddof=1)) if len(group) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows).sort_values("macro_f1_mean", ascending=False)


def run(
    seeds: list[int],
    model_keys: list[str],
    data_dir: Path | None,
    tables_dir: Path,
    checkpoint_root: Path,
    max_train_rows: int | None,
    max_valid_rows: int | None,
    max_test_rows: int | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    resolved_data_dir = resolve_data_dir(data_dir)
    ensure_dirs(tables_dir, checkpoint_root)
    metrics_parts = []
    train_parts = []
    for seed in seeds:
        checkpoint_dir = checkpoint_root / f"seed_{seed}"
        train_summary = train_many(model_keys, resolved_data_dir, checkpoint_dir, seed, max_train_rows, max_valid_rows)
        test_metrics, _confusion = evaluate_many(model_keys, checkpoint_dir, resolved_data_dir, "test", max_test_rows, seed)
        train_summary.insert(0, "seed", seed)
        test_metrics.insert(0, "seed", seed)
        train_parts.append(train_summary)
        metrics_parts.append(test_metrics)

    all_train = pd.concat(train_parts, ignore_index=True)
    all_metrics = pd.concat(metrics_parts, ignore_index=True)
    summary = summarize(all_metrics)
    write_dataframe(tables_dir / "experiment_1_multiseed_train_valid_metrics.csv", all_train)
    write_dataframe(tables_dir / "experiment_1_multiseed_metrics.csv", all_metrics)
    write_dataframe(tables_dir / "experiment_1_multiseed_summary.csv", summary)
    return all_metrics, summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run experiment-1 classifiers across multiple random seeds.")
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--models", default="cnn,rnafm,rnabert")
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--tables-dir", type=Path, default=EXP1_TABLES_DIR)
    parser.add_argument("--checkpoint-root", type=Path, default=PROJECT_ROOT / "results/checkpoints/experiment_1_multiseed")
    parser.add_argument("--max-train-rows", type=int, default=855)
    parser.add_argument("--max-valid-rows", type=int, default=120)
    parser.add_argument("--max-test-rows", type=int, default=285)
    parser.add_argument("--full-data", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    seeds = [int(item.strip()) for item in args.seeds.split(",") if item.strip()]
    model_keys = parse_model_list(args.models)
    max_train_rows = None if args.full_data else args.max_train_rows
    max_valid_rows = None if args.full_data else args.max_valid_rows
    max_test_rows = None if args.full_data else args.max_test_rows
    _metrics, summary = run(
        seeds=seeds,
        model_keys=model_keys,
        data_dir=args.data_dir,
        tables_dir=args.tables_dir,
        checkpoint_root=args.checkpoint_root,
        max_train_rows=max_train_rows,
        max_valid_rows=max_valid_rows,
        max_test_rows=max_test_rows,
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
