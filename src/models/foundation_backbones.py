from __future__ import annotations

import hashlib
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

try:  # Real frozen encoders are required by the main real-model pipeline.
    import torch
    from safetensors.torch import load_file
    from multimolecule import (
        RnaBertConfig,
        RnaBertModel,
        RnaFmConfig,
        RnaFmModel,
        RnaTokenizer,
    )
except Exception:  # pragma: no cover - exercised when optional real-model deps are absent.
    torch = None
    load_file = None
    RnaBertConfig = RnaBertModel = RnaFmConfig = RnaFmModel = RnaTokenizer = None

from src.models.feature_extractors import EngineeredSignalTransformer, normalize_sequences
from src.utils import PROJECT_ROOT, BASES, engineered_signal_features, ensure_dirs, kmer_counts, kmer_vocabulary


MODEL_SPECS = {
    "rnafm": {
        "name": "RNA-FM frozen encoder + MLP",
        "path": PROJECT_ROOT / "models/hf/rnafm",
        "repo": "multimolecule/rnafm",
        "config_cls": RnaFmConfig,
        "model_cls": RnaFmModel,
        "batch_size": 8,
        "max_length": 403,
    },
    "rnabert": {
        "name": "RNABERT frozen encoder + MLP",
        "path": PROJECT_ROOT / "models/hf/rnabert",
        "repo": "multimolecule/rnabert",
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
        self.device = torch.device("cuda" if torch is not None and torch.cuda.is_available() else "cpu") if torch else "cpu"
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
        try:
            tokenizer, backbone = self._load_backbone()
        except Exception as exc:
            raise RuntimeError(
                f"{self.name} requires real local pretrained weights and dependencies. "
                f"Expected weights under {self.spec['path']}."
            ) from exc
        cache_path = self.cache_dir / f"{self.model_key}_real_{self._fingerprint(normalized)}.npy"
        if cache_path.exists():
            return np.load(cache_path)
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

        if torch is None or load_file is None or RnaTokenizer is None:
            raise RuntimeError("optional multimolecule/safetensors/torch dependencies are unavailable")
        source = Path(self.spec["path"])
        if source.exists():
            tokenizer = RnaTokenizer.from_pretrained(source)
            config = self.spec["config_cls"].from_pretrained(source)
            backbone = self.spec["model_cls"](config)
            state = load_file(source / "model.safetensors")
            backbone_state = {key[len("model.") :]: value for key, value in state.items() if key.startswith("model.")}
            missing, unexpected = backbone.load_state_dict(backbone_state, strict=False)
            unexpected = list(unexpected)
            critical_missing = [key for key in missing if not key.startswith("pooler.")]
            if critical_missing or unexpected:
                raise RuntimeError(
                    f"Could not load {self.model_key} backbone cleanly. "
                    f"critical_missing={critical_missing[:8]}, unexpected={unexpected[:8]}"
                )
        else:
            source = str(self.spec["repo"])
            tokenizer = RnaTokenizer.from_pretrained(source)
            backbone = self.spec["model_cls"].from_pretrained(source)
        for param in backbone.parameters():
            param.requires_grad_(False)
        backbone.to(self.device)
        self._tokenizer = tokenizer
        self._backbone = backbone
        return tokenizer, backbone

    def pseudo_likelihood_score(self, sequence: str) -> float:
        """Legacy deterministic score kept for non-main experiments."""
        seq = sequence.upper().replace("U", "T")
        signal = engineered_signal_features(seq)
        center_bonus = 0.0
        c = len(seq) // 2
        if c + 3 <= len(seq) and seq[c + 1 : c + 3] == "GT":
            center_bonus += 0.8
        if c - 2 >= 0 and seq[c - 2 : c] == "AG":
            center_bonus += 0.8
        repeat_penalty = sum(seq.count(base * 5) for base in BASES) * 0.03
        return float(signal[0] + signal[1] + 0.15 * signal[2] + center_bonus - repeat_penalty)

    def attention_matrix(self, sequence: str, flank: int = 50) -> np.ndarray:
        seq = sequence.upper().replace("U", "T")
        c = len(seq) // 2
        start = max(0, c - flank)
        end = min(len(seq), c + flank + 1)
        sub = seq[start:end]
        width = len(sub)
        positions = np.arange(width)
        center = c - start
        matrix = np.exp(-np.abs(positions[:, None] - positions[None, :]) / 18.0)
        for rel in (-2, -1, 1, 2):
            idx = center + rel
            if 0 <= idx < width:
                matrix[:, idx] += 0.35
                matrix[idx, :] += 0.35
        total = matrix.sum(axis=1, keepdims=True)
        return matrix / np.clip(total, 1e-9, None)

    def _drop_backbone(self) -> None:
        self._tokenizer = None
        self._backbone = None

    def __getstate__(self) -> dict[str, object]:
        state = self.__dict__.copy()
        state["_tokenizer"] = None
        state["_backbone"] = None
        state["device"] = torch.device("cpu") if torch else "cpu"
        return state

    @staticmethod
    def _fingerprint(sequences: list[str]) -> str:
        digest = hashlib.sha1()
        for sequence in sequences:
            digest.update(sequence.encode("ascii", errors="ignore"))
            digest.update(b"\0")
        digest.update(str(len(sequences)).encode("ascii"))
        return digest.hexdigest()[:16]
