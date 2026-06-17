from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

from src.experiments.exp1.common import load_split, make_model, parse_model_list, resolve_data_dir, save_model
from src.utils import EXP1_TABLES_DIR, PROJECT_ROOT, multiclass_metrics, write_dataframe


def train_one(
    model_key: str,
    data_dir: Path,
    checkpoint_dir: Path,
    seed: int,
    max_train_rows: int | None,
    max_valid_rows: int | None,
) -> dict[str, object]:
    train = load_split(data_dir, "train", max_train_rows, seed)
    valid = load_split(data_dir, "valid", max_valid_rows, seed)
    model = make_model(model_key, seed)

    started = time.time()
    train_sequences = train["sequence"].astype(str).tolist()
    train_labels = train["label"].to_numpy()
    valid_sequences = valid["sequence"].astype(str).tolist()
    valid_labels = valid["label"].to_numpy()
    try:
        model.fit(
            train_sequences,
            train_labels,
            validation_sequences=valid_sequences,
            validation_labels=valid_labels,
        )
    except TypeError:
        model.fit(train_sequences, train_labels)
    seconds = time.time() - started

    proba = model.predict_proba(valid_sequences)
    metrics = multiclass_metrics(valid_labels, proba)
    checkpoint = checkpoint_dir / f"{model_key}.joblib"
    save_model(checkpoint, model)

    return {
        "model_key": model_key,
        "model": model.name,
        "split": "valid",
        "train_rows": len(train),
        "valid_rows": len(valid),
        "fit_seconds": seconds,
        "checkpoint": str(checkpoint),
        "best_epoch": getattr(model, "best_epoch_", None),
        "best_valid_macro_f1": getattr(model, "best_valid_macro_f1_", None),
        "epochs_trained": getattr(model, "epochs_trained_", None),
        **metrics,
    }


def train_many(
    model_keys: list[str],
    data_dir: Path,
    checkpoint_dir: Path,
    seed: int,
    max_train_rows: int | None,
    max_valid_rows: int | None,
) -> pd.DataFrame:
    rows = [
        train_one(model_key, data_dir, checkpoint_dir, seed, max_train_rows, max_valid_rows)
        for model_key in model_keys
    ]
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train experiment-1 splice-site classifiers.")
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--checkpoint-dir", type=Path, default=PROJECT_ROOT / "results/checkpoints/experiment_1")
    parser.add_argument("--tables-dir", type=Path, default=EXP1_TABLES_DIR)
    parser.add_argument("--models", default="cnn,rnafm,rnabert")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-train-rows", type=int, default=855)
    parser.add_argument("--max-valid-rows", type=int, default=120)
    parser.add_argument("--full-data", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = resolve_data_dir(args.data_dir)
    max_train_rows = None if args.full_data else args.max_train_rows
    max_valid_rows = None if args.full_data else args.max_valid_rows
    model_keys = parse_model_list(args.models)
    summary = train_many(model_keys, data_dir, args.checkpoint_dir, args.seed, max_train_rows, max_valid_rows)
    write_dataframe(args.tables_dir / "experiment_1_train_valid_metrics.csv", summary)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
