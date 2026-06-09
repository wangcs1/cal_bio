from __future__ import annotations

from pathlib import Path

import joblib
import pandas as pd

from src.models.cnn import CNNBaselineClassifier
from src.models.rnabert_mlp import RNABERTMLPClassifier
from src.models.rnafm_mlp import RNAFMMLPClassifier
from src.utils import PROJECT_ROOT, ensure_dirs


MODEL_FACTORIES = {
    "cnn": CNNBaselineClassifier,
    "rnafm": RNAFMMLPClassifier,
    "rnabert": RNABERTMLPClassifier,
}


def resolve_data_dir(data_dir: Path | None = None) -> Path:
    if data_dir is not None:
        return data_dir
    root_splits = PROJECT_ROOT / "splits"
    if (root_splits / "train.csv").exists():
        return root_splits
    return PROJECT_ROOT / "data/splits"


def split_path(data_dir: Path, split: str) -> Path:
    direct = data_dir / f"{split}.csv"
    if direct.exists():
        return direct
    pm200 = data_dir / f"{split}_pm200.csv"
    if pm200.exists():
        return pm200
    raise FileNotFoundError(f"Could not find {split}.csv or {split}_pm200.csv in {data_dir}")


def load_split(data_dir: Path, split: str, max_rows: int | None, seed: int) -> pd.DataFrame:
    path = split_path(data_dir, split)
    frame = pd.read_csv(path, usecols=lambda col: col in {"sample_id", "chrom", "label", "sequence", "gene_id"})
    frame = frame.dropna(subset=["label", "sequence"]).copy()
    frame["label"] = frame["label"].astype(int)
    if max_rows is not None and len(frame) > max_rows:
        frame = stratified_sample(frame, max_rows=max_rows, seed=seed)
    return frame.reset_index(drop=True)


def stratified_sample(frame: pd.DataFrame, max_rows: int, seed: int) -> pd.DataFrame:
    if max_rows <= 0 or len(frame) <= max_rows:
        return frame
    parts = []
    base = max_rows // frame["label"].nunique()
    remainder = max_rows - base * frame["label"].nunique()
    for idx, (_, group) in enumerate(frame.groupby("label", sort=True)):
        n = min(len(group), base + (1 if idx < remainder else 0))
        parts.append(group.sample(n=n, random_state=seed + idx))
    sampled = pd.concat(parts, ignore_index=True)
    return sampled.sample(frac=1.0, random_state=seed).reset_index(drop=True)


def make_model(model_key: str, seed: int):
    try:
        return MODEL_FACTORIES[model_key](random_state=seed)
    except KeyError as exc:
        valid = ", ".join(sorted(MODEL_FACTORIES))
        raise ValueError(f"Unknown model {model_key!r}. Valid models: {valid}") from exc


def save_model(path: Path, model: object) -> None:
    ensure_dirs(path.parent)
    joblib.dump(model, path)


def load_model(path: Path) -> object:
    return joblib.load(path)


def parse_model_list(value: str) -> list[str]:
    models = [item.strip() for item in value.split(",") if item.strip()]
    unknown = sorted(set(models) - set(MODEL_FACTORIES))
    if unknown:
        valid = ", ".join(sorted(MODEL_FACTORIES))
        raise ValueError(f"Unknown models {unknown}. Valid models: {valid}")
    return models

