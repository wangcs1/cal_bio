from __future__ import annotations

from importlib import import_module
from pathlib import Path

import joblib
import pandas as pd

from src.data.build_synthetic_splice_dataset import build_and_write
from src.utils import PROJECT_ROOT, REQUIRED_SPLIT_COLUMNS, SHARED_SPLIT_DIR, ensure_dirs, validate_split_frame


MODEL_FACTORIES = {
    "cnn": "src.models.cnn:CNNBaselineClassifier",
    "rnafm": "src.models.rnafm_mlp:RNAFMMLPClassifier",
    "rnabert": "src.models.rnabert_mlp:RNABERTMLPClassifier",
    "spliceai": "src.models.spliceai_wrapper:SpliceAIThreeClassWrapper",
}


def resolve_data_dir(data_dir: Path | None = None) -> Path:
    if data_dir is not None:
        return data_dir
    return SHARED_SPLIT_DIR


def split_path(data_dir: Path, split: str) -> Path:
    if not data_dir.exists() or not (data_dir / "train.csv").exists():
        build_and_write()
    direct = data_dir / f"{split}.csv"
    if direct.exists():
        return direct
    pm200 = data_dir / f"{split}_pm200.csv"
    if pm200.exists():
        return pm200
    raise FileNotFoundError(f"Could not find {split}.csv or {split}_pm200.csv in {data_dir}")


def load_split(data_dir: Path, split: str, max_rows: int | None, seed: int) -> pd.DataFrame:
    path = split_path(data_dir, split)
    frame = pd.read_csv(path)
    validate_split_frame(frame, path)
    keep = sorted(REQUIRED_SPLIT_COLUMNS | {"motif_type", "data_source", "split"})
    frame = frame.loc[:, [col for col in keep if col in frame.columns]]
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
        factory_path = MODEL_FACTORIES[model_key]
    except KeyError as exc:
        valid = ", ".join(sorted(MODEL_FACTORIES))
        raise ValueError(f"Unknown model {model_key!r}. Valid models: {valid}") from exc
    module_name, class_name = factory_path.split(":", 1)
    factory = getattr(import_module(module_name), class_name)
    return factory(random_state=seed)


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
