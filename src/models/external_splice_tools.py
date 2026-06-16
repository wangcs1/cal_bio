from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import numpy as np


@dataclass
class ToolScore:
    donor: float
    acceptor: float
    non_splice: float


def _normalize_three(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    arr = np.clip(arr, 1e-8, None)
    arr = arr / arr.sum()
    if arr.shape != (3,):
        raise ValueError(f"Expected 3 scores, got shape {arr.shape}")
    return arr


def maxentscan_three_class(sequence: str) -> ToolScore:
    from maxentpy import maxent

    seq = sequence.upper().replace("U", "T")
    c = len(seq) // 2
    donor_seq = seq[max(0, c - 2) : c + 7]
    acceptor_seq = seq[max(0, c - 20) : c + 3]
    donor = float(maxent.score5(donor_seq)) if len(donor_seq) == 9 else 0.0
    acceptor = float(maxent.score3(acceptor_seq)) if len(acceptor_seq) == 23 else 0.0
    non = float(max(0.0, 1.0 - 0.15 * max(abs(donor), abs(acceptor))))
    return ToolScore(*_normalize_three(np.array([donor, acceptor, non], dtype=float)))


@lru_cache(maxsize=1)
def _load_spliceai_models():
    from keras.models import load_model
    from pkg_resources import resource_filename

    paths = [resource_filename("spliceai", f"models/spliceai{i}.h5") for i in range(1, 6)]
    return [load_model(path, compile=False) for path in paths]


def spliceai_three_class(sequence: str) -> ToolScore:
    from spliceai.utils import one_hot_encode

    seq = sequence.upper().replace("U", "T")
    min_len = 10001
    if len(seq) < min_len:
        left = (min_len - len(seq)) // 2
        right = min_len - len(seq) - left
        seq = "N" * left + seq + "N" * right
    models = _load_spliceai_models()
    x = one_hot_encode(seq)[None, :]
    preds = np.mean([model.predict(x, verbose=0) for model in models], axis=0)[0]
    if preds.ndim == 2 and preds.shape[0] == 1:
        preds = preds[0]
    if preds.ndim == 1:
        acceptor = float(preds[1])
        donor = float(preds[2])
    else:
        donor = float(preds[:, 2].max())
        acceptor = float(preds[:, 1].max())
    non = float(max(0.0, 1.0 - max(donor, acceptor)))
    return ToolScore(*_normalize_three(np.array([donor, acceptor, non], dtype=float)))


@lru_cache(maxsize=1)
def _load_mmsplice():
    from mmsplice import MMSplice

    return MMSplice()


def mmsplice_three_class(sequence: str) -> ToolScore:
    seq = sequence.upper().replace("U", "T")
    model = _load_mmsplice()
    scores = np.asarray(model.predict_on_seq(seq, overhang=(80, 80)), dtype=float)
    donor = float(np.max(scores[[3, 4]]))
    acceptor = float(np.max(scores[[0, 1]]))
    non = float(max(0.0, 1.0 - 0.12 * max(donor, acceptor)))
    return ToolScore(*_normalize_three(np.array([donor, acceptor, non], dtype=float)))


@lru_cache(maxsize=1)
def _load_pangolin_models():
    import torch
    from pangolin.model import AR, L, W, Pangolin
    from pkg_resources import resource_filename

    models = []
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    for i in [0, 2, 4, 6]:
        for j in range(1, 4):
            model = Pangolin(L, W, AR).to(device)
            weights_path = resource_filename("pangolin", f"models/final.{j}.{i}.3.v2")
            weights = torch.load(weights_path, map_location=device)
            model.load_state_dict(weights)
            model.eval()
            models.append(model)
    return models


def pangolin_three_class(sequence: str) -> ToolScore:
    from pangolin.pangolin import compute_score

    seq = sequence.upper().replace("U", "T")
    min_len = 10001
    if len(seq) < min_len:
        left = (min_len - len(seq)) // 2
        right = min_len - len(seq) - left
        seq = "N" * left + seq + "N" * right
    loss, gain = compute_score(seq, seq, "+", 50, _load_pangolin_models())
    signal = np.maximum(np.abs(loss), np.abs(gain))
    donor = float(np.max(signal))
    acceptor = float(np.mean(np.sort(signal)[-5:])) if len(signal) >= 5 else donor
    non = float(max(0.0, 1.0 - max(donor, acceptor)))
    return ToolScore(*_normalize_three(np.array([donor, acceptor, non], dtype=float)))
