from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_recall_fscore_support,
    roc_auc_score,
)
from sklearn.preprocessing import label_binarize


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = PROJECT_ROOT / "data"
SHARED_PROCESSED_DIR = DATA_ROOT / "shared/processed"
SHARED_SPLIT_DIR = DATA_ROOT / "shared/splits"
EXP3_DATA_DIR = DATA_ROOT / "experiment_3"
RESULTS_ROOT = PROJECT_ROOT / "results"
EXP1_RESULTS_DIR = RESULTS_ROOT / "experiment_1"
EXP2_RESULTS_DIR = RESULTS_ROOT / "experiment_2"
EXP3_RESULTS_DIR = RESULTS_ROOT / "experiment_3"
EXP1_TABLES_DIR = EXP1_RESULTS_DIR / "tables"
EXP1_FIGURES_DIR = EXP1_RESULTS_DIR / "figures"
EXP2_TABLES_DIR = EXP2_RESULTS_DIR / "tables"
EXP2_FIGURES_DIR = EXP2_RESULTS_DIR / "figures"
EXP3_TABLES_DIR = EXP3_RESULTS_DIR / "tables"
EXP3_FIGURES_DIR = EXP3_RESULTS_DIR / "figures"
BASES = "ACGT"
LABELS = {0: "donor", 1: "acceptor", 2: "non_splice"}
CLASS_ORDER = [0, 1, 2]
TRAIN_CHROMS = {f"chr{i}" for i in range(1, 17)}
VALID_CHROMS = {"chr17", "chr18"}
TEST_CHROMS = {"chr19", "chr20", "chr21", "chr22", "chrX"}
ALL_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def ensure_dirs(*paths: Path) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def first_existing_path(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


def shared_processed_file(name: str) -> Path:
    return first_existing_path(SHARED_PROCESSED_DIR / name, PROJECT_ROOT / "data/processed" / name)


def shared_split_file(name: str) -> Path:
    return first_existing_path(SHARED_SPLIT_DIR / name, PROJECT_ROOT / "data/splits" / name)


def exp3_data_file(name: str) -> Path:
    return first_existing_path(EXP3_DATA_DIR / name, PROJECT_ROOT / "data/processed" / name)


def stable_id(*parts: object, prefix: str = "id") -> str:
    text = "|".join(map(str, parts))
    digest = hashlib.sha1(text.encode("utf-8")).hexdigest()[:16]
    return f"{prefix}_{digest}"


def random_dna(length: int, rng: random.Random, gc: float = 0.42) -> str:
    gc_each = gc / 2.0
    at_each = (1.0 - gc) / 2.0
    alphabet = ["A", "C", "G", "T"]
    weights = [at_each, gc_each, gc_each, at_each]
    return "".join(rng.choices(alphabet, weights=weights, k=length))


def insert_motif(sequence: str, start: int, motif: str) -> str:
    chars = list(sequence)
    if start < 0 or start + len(motif) > len(chars):
        raise ValueError(f"Motif {motif!r} exceeds sequence bounds at {start}")
    chars[start : start + len(motif)] = motif
    return "".join(chars)


def mutate_base(sequence: str, index: int, alt: str) -> str:
    if index < 0 or index >= len(sequence):
        raise ValueError(f"Mutation index outside sequence: {index}")
    if alt not in BASES:
        raise ValueError(f"Alt must be one of {BASES}: {alt}")
    chars = list(sequence)
    chars[index] = alt
    return "".join(chars)


def choose_alt(ref: str, rng: random.Random, preferred: str | None = None) -> str:
    if preferred is not None and preferred != ref:
        return preferred
    choices = [base for base in BASES if base != ref]
    return rng.choice(choices)


def crop_center(sequence: str, flank: int) -> str:
    center = len(sequence) // 2
    return sequence[center - flank : center + flank + 1]


def split_name(chrom: str) -> str:
    if chrom in TRAIN_CHROMS:
        return "train"
    if chrom in VALID_CHROMS:
        return "valid"
    if chrom in TEST_CHROMS:
        return "test"
    return "unassigned"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    ensure_dirs(path.parent)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def save_json(path: Path, payload: dict[str, object]) -> None:
    ensure_dirs(path.parent)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def load_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def safe_roc_auc(y_true: np.ndarray, score: np.ndarray) -> float:
    try:
        if len(np.unique(y_true)) < 2:
            return float("nan")
        return float(roc_auc_score(y_true, score))
    except ValueError:
        return float("nan")


def safe_average_precision(y_true: np.ndarray, score: np.ndarray) -> float:
    try:
        if len(np.unique(y_true)) < 2:
            return float("nan")
        return float(average_precision_score(y_true, score))
    except ValueError:
        return float("nan")


def multiclass_metrics(y_true: Iterable[int], proba: np.ndarray) -> dict[str, float]:
    y = np.asarray(list(y_true), dtype=int)
    pred = np.argmax(proba, axis=1)
    precision, recall, macro_f1, _ = precision_recall_fscore_support(
        y, pred, average="macro", zero_division=0
    )
    _, _, per_class_f1, _ = precision_recall_fscore_support(
        y, pred, labels=CLASS_ORDER, zero_division=0
    )
    y_bin = label_binarize(y, classes=CLASS_ORDER)
    try:
        auroc = roc_auc_score(y_bin, proba, average="macro", multi_class="ovr")
    except ValueError:
        auroc = float("nan")
    try:
        auprc = average_precision_score(y_bin, proba, average="macro")
    except ValueError:
        auprc = float("nan")

    return {
        "accuracy": float(accuracy_score(y, pred)),
        "macro_precision": float(precision),
        "macro_recall": float(recall),
        "macro_f1": float(macro_f1),
        "auroc": float(auroc),
        "auprc": float(auprc),
        "donor_f1": float(per_class_f1[0]),
        "acceptor_f1": float(per_class_f1[1]),
        "non_splice_f1": float(per_class_f1[2]),
    }


def binary_ranking_metrics(y_true: Iterable[int], score: Iterable[float], k_fraction: float = 0.1) -> dict[str, float]:
    y = np.asarray(list(y_true), dtype=int)
    s = np.asarray(list(score), dtype=float)
    order = np.argsort(-s)
    k = max(1, int(math.ceil(len(y) * k_fraction)))
    top = y[order[:k]]
    positives = int(y.sum())
    top_positive = int(top.sum())
    base_rate = positives / len(y) if len(y) else float("nan")
    top_rate = top_positive / k if k else float("nan")
    top_k_recall = top_positive / positives if positives else float("nan")
    enrichment = top_rate / base_rate if base_rate and not math.isnan(base_rate) else float("nan")
    return {
        "auroc": safe_roc_auc(y, s),
        "auprc": safe_average_precision(y, s),
        "top_k": float(k),
        "top_k_recall": float(top_k_recall),
        "enrichment_at_k": float(enrichment),
    }


def hard_negative_fpr(frame: pd.DataFrame, proba: np.ndarray) -> float:
    if "negative_type" not in frame.columns:
        return float("nan")
    pred = np.argmax(proba, axis=1)
    mask = (frame["label"].astype(int) == 2) & frame["negative_type"].astype(str).str.contains("hard", na=False)
    if mask.sum() == 0:
        return float("nan")
    return float(np.mean(pred[mask.to_numpy()] != 2))


def confusion_rows(y_true: Iterable[int], proba: np.ndarray, model_name: str, split: str) -> list[dict[str, object]]:
    y = np.asarray(list(y_true), dtype=int)
    pred = np.argmax(proba, axis=1)
    matrix = confusion_matrix(y, pred, labels=CLASS_ORDER)
    rows: list[dict[str, object]] = []
    for i, true_label in enumerate(CLASS_ORDER):
        for j, pred_label in enumerate(CLASS_ORDER):
            rows.append(
                {
                    "model": model_name,
                    "split": split,
                    "true_label": true_label,
                    "true_label_name": LABELS[true_label],
                    "pred_label": pred_label,
                    "pred_label_name": LABELS[pred_label],
                    "count": int(matrix[i, j]),
                }
            )
    return rows


def sigmoid(value: float) -> float:
    value = max(-40.0, min(40.0, value))
    return 1.0 / (1.0 + math.exp(-value))


def gc_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def pyrimidine_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("C") + sequence.count("T")) / len(sequence)


def _match(sequence: str, index: int, base: str) -> float:
    if index < 0 or index >= len(sequence):
        return 0.0
    return 1.0 if sequence[index] == base else 0.0


def donor_consensus_score(sequence: str, center: int | None = None) -> float:
    c = len(sequence) // 2 if center is None else center
    terms = [
        (2.4, _match(sequence, c + 1, "G")),
        (2.4, _match(sequence, c + 2, "T")),
        (0.8, _match(sequence, c - 2, "C")),
        (0.7, _match(sequence, c - 1, "A")),
        (0.6, _match(sequence, c, "G")),
        (0.5, _match(sequence, c + 3, "A")),
        (0.5, _match(sequence, c + 4, "A")),
        (0.5, _match(sequence, c + 5, "G")),
        (0.4, _match(sequence, c + 6, "T")),
    ]
    raw = sum(weight * val for weight, val in terms) / sum(weight for weight, _ in terms)
    downstream = sequence[max(0, c + 1) : min(len(sequence), c + 80)]
    motif_density = downstream.count("GTA") + downstream.count("GTG")
    return raw + min(0.25, motif_density / 120.0)


def acceptor_consensus_score(sequence: str, center: int | None = None) -> float:
    c = len(sequence) // 2 if center is None else center
    ag = 0.5 * _match(sequence, c - 2, "A") + 0.5 * _match(sequence, c - 1, "G")
    poly = pyrimidine_fraction(sequence[max(0, c - 35) : max(0, c - 5)])
    branch = 1.0 if "TACTAAC" in sequence[max(0, c - 55) : max(0, c - 18)] else 0.0
    upstream_ag_penalty = 0.08 * sequence[max(0, c - 35) : max(0, c - 5)].count("AG")
    return max(0.0, 0.58 * ag + 0.32 * poly + 0.10 * branch - upstream_ag_penalty)


def engineered_signal_features(sequence: str) -> np.ndarray:
    c = len(sequence) // 2
    windows = [20, 50, 100, min(200, c)]
    values: list[float] = [
        donor_consensus_score(sequence),
        acceptor_consensus_score(sequence),
        gc_fraction(sequence),
        pyrimidine_fraction(sequence[max(0, c - 35) : max(0, c - 5)]),
        float(sequence[c + 1 : c + 3] == "GT") if c + 3 <= len(sequence) else 0.0,
        float(sequence[c - 2 : c] == "AG") if c - 2 >= 0 else 0.0,
    ]
    for flank in windows:
        sub = sequence[max(0, c - flank) : min(len(sequence), c + flank + 1)]
        values.extend(
            [
                gc_fraction(sub),
                sub.count("GT") / max(1, len(sub) - 1),
                sub.count("AG") / max(1, len(sub) - 1),
                sub.count("CT") / max(1, len(sub) - 1),
                sub.count("GA") / max(1, len(sub) - 1),
            ]
        )
    return np.asarray(values, dtype=float)


def spliceai_proxy_proba(sequence: str) -> np.ndarray:
    donor = sigmoid(8.0 * (donor_consensus_score(sequence) - 0.62))
    acceptor = sigmoid(8.0 * (acceptor_consensus_score(sequence) - 0.58))
    non = max(0.02, 1.15 - max(donor, acceptor) - 0.18 * min(donor, acceptor))
    raw = np.asarray([donor, acceptor, non], dtype=float)
    return raw / raw.sum()


def scan_splice_scores(sequence: str, flank: int = 80) -> pd.DataFrame:
    c = len(sequence) // 2
    rows = []
    for pos in range(max(3, c - flank), min(len(sequence) - 7, c + flank + 1)):
        rows.append(
            {
                "relative_position": pos - c,
                "donor_score": donor_consensus_score(sequence, pos),
                "acceptor_score": acceptor_consensus_score(sequence, pos),
            }
        )
    return pd.DataFrame(rows)


def kmer_vocabulary(k_values: Iterable[int]) -> list[str]:
    vocab: list[str] = []
    for k in k_values:
        current = [""]
        for _ in range(k):
            current = [prefix + base for prefix in current for base in BASES]
        vocab.extend(current)
    return vocab


def kmer_counts(sequence: str, vocab: list[str]) -> np.ndarray:
    values = np.zeros(len(vocab), dtype=float)
    seq = sequence.upper()
    for idx, kmer in enumerate(vocab):
        k = len(kmer)
        denom = max(1, len(seq) - k + 1)
        values[idx] = seq.count(kmer) / denom
    return values


def positional_base_features(sequence: str, positions: Iterable[int]) -> np.ndarray:
    c = len(sequence) // 2
    values: list[float] = []
    for rel in positions:
        idx = c + rel
        for base in BASES:
            values.append(1.0 if 0 <= idx < len(sequence) and sequence[idx] == base else 0.0)
    return np.asarray(values, dtype=float)


def write_dataframe(path: Path, frame: pd.DataFrame) -> None:
    ensure_dirs(path.parent)
    frame.to_csv(path, index=False, encoding="utf-8")
