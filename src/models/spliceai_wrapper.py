from __future__ import annotations

import numpy as np

from src.utils import spliceai_proxy_proba


class SpliceAIThreeClassWrapper:
    """Small-sample SpliceAI-compatible three-class baseline.

    The wrapper keeps the experiment interface runnable without requiring the
    external SpliceAI package. When a real SpliceAI workflow is available, this
    class can be replaced by one that fills the same `predict_proba` contract.
    """

    name = "SpliceAI optional real tool (proxy fallback)"

    def __init__(self, random_state: int = 42, center_radius: int = 50) -> None:
        self.random_state = random_state
        self.center_radius = center_radius

    def fit(self, sequences: list[str], labels: np.ndarray | None = None) -> "SpliceAIThreeClassWrapper":
        return self

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        return np.vstack([spliceai_proxy_proba(sequence) for sequence in sequences])


def spliceai_variant_delta(wt_sequence: str, mut_sequence: str, target_class: int, variant_type: str) -> tuple[float, float, float]:
    wt = spliceai_proxy_proba(wt_sequence)
    mut = spliceai_proxy_proba(mut_sequence)
    if variant_type.endswith("loss"):
        delta = float(wt[target_class] - mut[target_class])
    elif "gain" in variant_type:
        delta = float(mut[target_class] - wt[target_class])
    else:
        delta = float(np.max(np.abs(mut[:2] - wt[:2])))
    return float(wt[target_class]), float(mut[target_class]), delta
