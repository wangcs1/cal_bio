from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from src.utils import (
    BASES,
    engineered_signal_features,
    kmer_counts,
    kmer_vocabulary,
    positional_base_features,
    spliceai_proxy_proba,
)


class SequenceClassifier:
    def fit(self, sequences: list[str], labels: np.ndarray) -> "SequenceClassifier":
        raise NotImplementedError

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        raise NotImplementedError


@dataclass
class FeatureLogisticClassifier(SequenceClassifier):
    name: str
    mode: str
    random_state: int = 42
    max_iter: int = 1000

    def __post_init__(self) -> None:
        self.vocab_3_4 = kmer_vocabulary([3, 4])
        self.vocab_3_5 = kmer_vocabulary([3, 5])
        self.scaler = StandardScaler()
        self.estimator = LogisticRegression(
            max_iter=self.max_iter,
            C=2.0,
            class_weight="balanced",
            solver="lbfgs",
            random_state=self.random_state,
        )

    def fit(self, sequences: list[str], labels: np.ndarray) -> "FeatureLogisticClassifier":
        x = self._features(sequences)
        self.scaler.fit(x)
        self.estimator.fit(self.scaler.transform(x), labels)
        return self

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        x = self._features(sequences)
        return self.estimator.predict_proba(self.scaler.transform(x))

    def _features(self, sequences: list[str]) -> np.ndarray:
        sequences = [seq.upper().replace("U", "T") for seq in sequences]
        if self.mode == "cnn":
            return np.vstack([self._one_hot(seq) for seq in sequences])
        if self.mode == "rnafm":
            return np.vstack(
                [
                    np.concatenate(
                        [
                            kmer_counts(seq, self.vocab_3_4),
                            engineered_signal_features(seq),
                        ]
                    )
                    for seq in sequences
                ]
            )
        if self.mode == "rnabert":
            return np.vstack(
                [
                    np.concatenate(
                        [
                            kmer_counts(seq, self.vocab_3_5),
                            positional_base_features(seq, range(-12, 13)),
                            engineered_signal_features(seq),
                        ]
                    )
                    for seq in sequences
                ]
            )
        raise ValueError(f"Unknown feature mode: {self.mode}")

    @staticmethod
    def _one_hot(sequence: str) -> np.ndarray:
        center = len(sequence) // 2
        sequence = sequence[max(0, center - 12) : min(len(sequence), center + 13)]
        values = np.zeros((len(sequence), len(BASES)), dtype=float)
        base_to_idx = {base: idx for idx, base in enumerate(BASES)}
        for pos, base in enumerate(sequence):
            idx = base_to_idx.get(base)
            if idx is not None:
                values[pos, idx] = 1.0
        return values.ravel()


class SpliceAIProxyClassifier(SequenceClassifier):
    def __init__(self, random_state: int = 42) -> None:
        self.name = "SpliceAI signal proxy"
        self.random_state = random_state
        self.estimator = RandomForestClassifier(
            n_estimators=180,
            max_depth=9,
            min_samples_leaf=2,
            class_weight="balanced_subsample",
            random_state=random_state,
            n_jobs=1,
        )

    def fit(self, sequences: list[str], labels: np.ndarray) -> "SpliceAIProxyClassifier":
        x = np.vstack([engineered_signal_features(seq) for seq in sequences])
        self.estimator.fit(x, labels)
        return self

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        x = np.vstack([engineered_signal_features(seq) for seq in sequences])
        model_proba = self.estimator.predict_proba(x)
        classes = list(self.estimator.classes_)
        aligned = np.zeros((len(sequences), 3), dtype=float)
        for col, label in enumerate(classes):
            aligned[:, int(label)] = model_proba[:, col]
        heuristic = np.vstack([spliceai_proxy_proba(seq) for seq in sequences])
        mixed = 0.72 * aligned + 0.28 * heuristic
        return mixed / mixed.sum(axis=1, keepdims=True)


def make_model_suite(random_state: int = 42) -> list[SequenceClassifier]:
    return [
        FeatureLogisticClassifier("CNN motif baseline", "cnn", random_state=random_state),
        FeatureLogisticClassifier("RNA-FM frozen k-mer + MLP", "rnafm", random_state=random_state),
        FeatureLogisticClassifier("RNABERT frozen token + MLP", "rnabert", random_state=random_state),
        SpliceAIProxyClassifier(random_state=random_state),
    ]


def zero_shot_embedding_distance(wt_sequence: str, mut_sequence: str, mode: str) -> float:
    if mode == "rnafm":
        vocab = kmer_vocabulary([3, 4])
        wt = np.concatenate([kmer_counts(wt_sequence, vocab), engineered_signal_features(wt_sequence)])
        mut = np.concatenate([kmer_counts(mut_sequence, vocab), engineered_signal_features(mut_sequence)])
    elif mode == "rnabert":
        vocab = kmer_vocabulary([3, 5])
        wt = np.concatenate(
            [
                kmer_counts(wt_sequence, vocab),
                positional_base_features(wt_sequence, range(-12, 13)),
                engineered_signal_features(wt_sequence),
            ]
        )
        mut = np.concatenate(
            [
                kmer_counts(mut_sequence, vocab),
                positional_base_features(mut_sequence, range(-12, 13)),
                engineered_signal_features(mut_sequence),
            ]
        )
    else:
        raise ValueError(f"Unknown zero-shot mode: {mode}")
    return float(np.linalg.norm(wt - mut))


def pseudo_likelihood_proxy_score(sequence: str, mode: str) -> float:
    features = engineered_signal_features(sequence)
    if mode == "rnafm":
        vocab = kmer_vocabulary([3, 4])
        kmer_signal = np.dot(kmer_counts(sequence, vocab), np.linspace(0.1, 0.4, len(vocab)))
    elif mode == "rnabert":
        vocab = kmer_vocabulary([3, 5])
        kmer_signal = np.dot(kmer_counts(sequence, vocab), np.linspace(0.1, 0.5, len(vocab)))
    else:
        raise ValueError(f"Unknown zero-shot mode: {mode}")
    return float(0.65 * features[0] + 0.65 * features[1] + 0.15 * features[2] + 0.02 * kmer_signal)
