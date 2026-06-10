from __future__ import annotations

import numpy as np

from src.utils import spliceai_proxy_proba


class PangolinCaseStudyWrapper:
    name = "Pangolin optional tool (small case-study proxy)"

    def __init__(self, random_state: int = 42) -> None:
        self.random_state = random_state

    def fit(self, sequences: list[str], labels: np.ndarray | None = None) -> "PangolinCaseStudyWrapper":
        return self

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        rows = []
        for sequence in sequences:
            base = spliceai_proxy_proba(sequence)
            # Pangolin is used here as a long-context case-study proxy; make it
            # slightly more conservative on splice calls than the SpliceAI proxy.
            adjusted = np.asarray([0.92 * base[0], 0.92 * base[1], min(1.0, 1.08 * base[2])])
            rows.append(adjusted / adjusted.sum())
        return np.vstack(rows)


def pangolin_variant_delta(wt_sequence: str, mut_sequence: str, target_class: int, variant_type: str) -> tuple[float, float, float]:
    wrapper = PangolinCaseStudyWrapper()
    wt, mut = wrapper.predict_proba([wt_sequence, mut_sequence])
    if variant_type.endswith("loss"):
        delta = float(wt[target_class] - mut[target_class])
    elif "gain" in variant_type:
        delta = float(mut[target_class] - wt[target_class])
    else:
        delta = float(np.max(np.abs(mut[:2] - wt[:2])))
    return float(wt[target_class]), float(mut[target_class]), delta
