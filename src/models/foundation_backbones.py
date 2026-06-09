from __future__ import annotations

import hashlib
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import torch
from safetensors.torch import load_file
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from multimolecule import (
    RnaBertConfig,
    RnaBertModel,
    RnaFmConfig,
    RnaFmModel,
    RnaTokenizer,
)

from src.models.feature_extractors import EngineeredSignalTransformer, normalize_sequences
from src.utils import PROJECT_ROOT, ensure_dirs


MODEL_SPECS = {
    "rnafm": {
        "name": "RNA-FM frozen encoder + MLP",
        "path": PROJECT_ROOT / "models/hf/rnafm",
        "config_cls": RnaFmConfig,
        "model_cls": RnaFmModel,
        "batch_size": 8,
        "max_length": 403,
    },
    "rnabert": {
        "name": "RNABERT frozen encoder + MLP",
        "path": PROJECT_ROOT / "models/hf/rnabert",
        "config_cls": RnaBertConfig,
        "model_cls": RnaBertModel,
        "batch_size": 64,
        "max_length": 403,
    },
}


class FrozenRnaFoundationClassifier:
    def __init__(
        self,
        model_key: str,
        random_state: int = 42,
        cache_dir: Path | None = None,
        batch_size: int | None = None,
    ) -> None:
        self.model_key = model_key
        self.spec = MODEL_SPECS[model_key]
        self.name = str(self.spec["name"])
        self.random_state = random_state
        self.cache_dir = cache_dir or PROJECT_ROOT / "results/embeddings/experiment_1"
        self.batch_size = batch_size or int(self.spec["batch_size"])
        self.max_length = int(self.spec["max_length"])
        self.signals = EngineeredSignalTransformer()
        self.scaler = StandardScaler(with_mean=False)
        self.estimator = LogisticRegression(
            max_iter=1200,
            C=1.5,
            class_weight="balanced",
            solver="lbfgs",
            random_state=random_state,
        )
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._tokenizer = None
        self._backbone = None

    def fit(self, sequences: list[str], labels: np.ndarray) -> "FrozenRnaFoundationClassifier":
        x = self._features(sequences)
        self.scaler.fit(x)
        self.estimator.fit(self.scaler.transform(x), labels)
        self._drop_backbone()
        return self

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        x = self._features(sequences)
        proba = self.estimator.predict_proba(self.scaler.transform(x))
        aligned = np.zeros((len(sequences), 3), dtype=float)
        for col, label in enumerate(self.estimator.classes_):
            aligned[:, int(label)] = proba[:, col]
        row_sum = aligned.sum(axis=1, keepdims=True)
        return np.divide(aligned, row_sum, out=np.full_like(aligned, 1.0 / 3.0), where=row_sum > 0)

    def _features(self, sequences: list[str]) -> sparse.csr_matrix:
        embeddings = self._embeddings(sequences)
        signal_features = self.signals.transform(sequences)
        return sparse.hstack([sparse.csr_matrix(embeddings), signal_features], format="csr", dtype=np.float32)

    def _embeddings(self, sequences: list[str]) -> np.ndarray:
        ensure_dirs(self.cache_dir)
        normalized = normalize_sequences(sequences)
        cache_path = self.cache_dir / f"{self.model_key}_{self._fingerprint(normalized)}.npy"
        if cache_path.exists():
            return np.load(cache_path)

        tokenizer, backbone = self._load_backbone()
        rna_sequences = [sequence.replace("T", "U") for sequence in normalized]
        chunks: list[np.ndarray] = []
        backbone.eval()
        with torch.no_grad():
            for start in range(0, len(rna_sequences), self.batch_size):
                batch = rna_sequences[start : start + self.batch_size]
                inputs = tokenizer(
                    batch,
                    return_tensors="pt",
                    padding=True,
                    truncation=True,
                    max_length=self.max_length,
                )
                inputs = {key: value.to(self.device) for key, value in inputs.items()}
                output = backbone(**inputs)
                hidden = output.last_hidden_state
                mask = inputs["attention_mask"].unsqueeze(-1).to(hidden.dtype)
                pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1.0)
                chunks.append(pooled.cpu().numpy().astype(np.float32))
        embeddings = np.vstack(chunks)
        np.save(cache_path, embeddings)
        return embeddings

    def _load_backbone(self):
        if self._tokenizer is not None and self._backbone is not None:
            return self._tokenizer, self._backbone

        model_path = Path(self.spec["path"])
        tokenizer = RnaTokenizer.from_pretrained(model_path)
        config = self.spec["config_cls"].from_pretrained(model_path)
        backbone = self.spec["model_cls"](config)
        state = load_file(model_path / "model.safetensors")
        backbone_state = {key[len("model.") :]: value for key, value in state.items() if key.startswith("model.")}
        missing, unexpected = backbone.load_state_dict(backbone_state, strict=False)
        unexpected = list(unexpected)
        critical_missing = [key for key in missing if not key.startswith("pooler.")]
        if critical_missing or unexpected:
            raise RuntimeError(
                f"Could not load {self.model_key} backbone cleanly. "
                f"critical_missing={critical_missing[:8]}, unexpected={unexpected[:8]}"
            )
        for param in backbone.parameters():
            param.requires_grad_(False)
        backbone.to(self.device)
        self._tokenizer = tokenizer
        self._backbone = backbone
        return tokenizer, backbone

    def _drop_backbone(self) -> None:
        self._tokenizer = None
        self._backbone = None

    def __getstate__(self) -> dict[str, object]:
        state = self.__dict__.copy()
        state["_tokenizer"] = None
        state["_backbone"] = None
        state["device"] = torch.device("cpu")
        return state

    @staticmethod
    def _fingerprint(sequences: list[str]) -> str:
        digest = hashlib.sha1()
        for sequence in sequences:
            digest.update(sequence.encode("ascii", errors="ignore"))
            digest.update(b"\0")
        digest.update(str(len(sequences)).encode("ascii"))
        return digest.hexdigest()[:16]
