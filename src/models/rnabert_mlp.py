from __future__ import annotations

from pathlib import Path

from src.models.foundation_backbones import FrozenRnaFoundationClassifier
from src.models.simple_splice_models import FeatureLogisticClassifier


class RNABERTProxyClassifier(FeatureLogisticClassifier):
    def __init__(self, random_state: int = 42) -> None:
        super().__init__("RNABERT proxy (token/k-mer + signal MLP)", "rnabert", random_state=random_state)


class RNABERTMLPClassifier(FrozenRnaFoundationClassifier):
    def __init__(self, random_state: int = 42, cache_dir: Path | None = None, batch_size: int | None = None) -> None:
        super().__init__("rnabert", random_state=random_state, cache_dir=cache_dir, batch_size=batch_size)
