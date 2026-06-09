from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from scipy import sparse
from sklearn.base import BaseEstimator, TransformerMixin

from src.utils import BASES, engineered_signal_features


BASE_TO_INDEX = {base: idx for idx, base in enumerate(BASES)}


def normalize_dna(sequence: object) -> str:
    return str(sequence).upper().replace("U", "T")


def normalize_sequences(sequences: Iterable[object]) -> list[str]:
    return [normalize_dna(sequence) for sequence in sequences]


class CenterOneHotTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, flank: int = 24) -> None:
        self.flank = flank

    def fit(self, x: Iterable[object], y: object | None = None) -> "CenterOneHotTransformer":
        return self

    def transform(self, x: Iterable[object]) -> sparse.csr_matrix:
        sequences = normalize_sequences(x)
        width = 2 * self.flank + 1
        n_features = width * len(BASES)
        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []

        for row_idx, sequence in enumerate(sequences):
            center = len(sequence) // 2
            start = center - self.flank
            for offset in range(width):
                seq_idx = start + offset
                if 0 <= seq_idx < len(sequence):
                    base_idx = BASE_TO_INDEX.get(sequence[seq_idx])
                    if base_idx is not None:
                        rows.append(row_idx)
                        cols.append(offset * len(BASES) + base_idx)
                        data.append(1.0)

        return sparse.csr_matrix((data, (rows, cols)), shape=(len(sequences), n_features), dtype=np.float32)


class PositionalBaseTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, positions: tuple[int, ...] = tuple(range(-20, 21))) -> None:
        self.positions = positions

    def fit(self, x: Iterable[object], y: object | None = None) -> "PositionalBaseTransformer":
        return self

    def transform(self, x: Iterable[object]) -> sparse.csr_matrix:
        sequences = normalize_sequences(x)
        n_features = len(self.positions) * len(BASES)
        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []

        for row_idx, sequence in enumerate(sequences):
            center = len(sequence) // 2
            for pos_idx, rel_pos in enumerate(self.positions):
                seq_idx = center + rel_pos
                if 0 <= seq_idx < len(sequence):
                    base_idx = BASE_TO_INDEX.get(sequence[seq_idx])
                    if base_idx is not None:
                        rows.append(row_idx)
                        cols.append(pos_idx * len(BASES) + base_idx)
                        data.append(1.0)

        return sparse.csr_matrix((data, (rows, cols)), shape=(len(sequences), n_features), dtype=np.float32)


class EngineeredSignalTransformer(BaseEstimator, TransformerMixin):
    def fit(self, x: Iterable[object], y: object | None = None) -> "EngineeredSignalTransformer":
        return self

    def transform(self, x: Iterable[object]) -> sparse.csr_matrix:
        values = np.vstack([engineered_signal_features(sequence) for sequence in normalize_sequences(x)])
        return sparse.csr_matrix(values.astype(np.float32))

