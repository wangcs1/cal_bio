from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.experiments.exp1.common import load_model, load_split, parse_model_list, resolve_data_dir
from src.utils import EXP1_TABLES_DIR, PROJECT_ROOT, confusion_rows, hard_negative_details, multiclass_metrics, write_dataframe


def evaluate_one(model_key: str, checkpoint_dir: Path, data_dir: Path, split: str, max_rows: int | None, seed: int):
    frame = load_split(data_dir, split, max_rows, seed)
    model = load_model(checkpoint_dir / f"{model_key}.joblib")
    proba = model.predict_proba(frame["sequence"].astype(str).tolist())
    metrics = multiclass_metrics(frame["label"].to_numpy(), proba)
    metrics.update(hard_negative_details(frame, proba))
    row = {
        "model_key": model_key,
        "model": model.name,
        "split": split,
        "rows": len(frame),
        **metrics,
    }
    confusion = pd.DataFrame(confusion_rows(frame["label"].to_numpy(), proba, model.name, split))
    return row, confusion


def evaluate_many(
    model_keys: list[str],
    checkpoint_dir: Path,
    data_dir: Path,
    split: str,
    max_rows: int | None,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    confusion_frames = []
    for model_key in model_keys:
        row, confusion = evaluate_one(model_key, checkpoint_dir, data_dir, split, max_rows, seed)
        rows.append(row)
        confusion_frames.append(confusion)
    return pd.DataFrame(rows), pd.concat(confusion_frames, ignore_index=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate experiment-1 splice-site classifiers.")
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--checkpoint-dir", type=Path, default=PROJECT_ROOT / "results/checkpoints/experiment_1")
    parser.add_argument("--tables-dir", type=Path, default=EXP1_TABLES_DIR)
    parser.add_argument("--models", default="cnn,rnafm,rnabert")
    parser.add_argument("--split", default="test")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-rows", type=int, default=30000)
    parser.add_argument("--full-data", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = resolve_data_dir(args.data_dir)
    model_keys = parse_model_list(args.models)
    max_rows = None if args.full_data else args.max_rows
    metrics, confusion = evaluate_many(model_keys, args.checkpoint_dir, data_dir, args.split, max_rows, args.seed)
    write_dataframe(args.tables_dir / f"experiment_1_{args.split}_metrics.csv", metrics)
    write_dataframe(args.tables_dir / f"experiment_1_{args.split}_confusion_matrices.csv", confusion)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
