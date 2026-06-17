from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
from sklearn.metrics import f1_score

from src.models.feature_extractors import BASE_TO_INDEX, normalize_dna


class _SequenceDataset(Dataset):
    def __init__(self, sequences: list[str], labels: np.ndarray | None = None) -> None:
        self.sequences = [normalize_dna(sequence) for sequence in sequences]
        self.x = self._encode(self.sequences)
        self.labels = None if labels is None else labels.astype(np.int64)

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int):
        if self.labels is None:
            return self.x[index]
        return self.x[index], torch.tensor(int(self.labels[index]), dtype=torch.long)

    @staticmethod
    def _encode(sequences: list[str]) -> torch.Tensor:
        width = max((len(sequence) for sequence in sequences), default=0)
        encoded = np.zeros((len(sequences), 4, width), dtype=np.float32)
        for row_idx, sequence in enumerate(sequences):
            for pos, base in enumerate(sequence):
                base_idx = BASE_TO_INDEX.get(base)
                if base_idx is not None:
                    encoded[row_idx, base_idx, pos] = 1.0
        return torch.from_numpy(encoded)


class _SpliceCNN(nn.Module):
    def __init__(self, n_classes: int = 3) -> None:
        super().__init__()
        self.branches = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Conv1d(4, 48, kernel_size=kernel_size, padding=kernel_size // 2),
                    nn.BatchNorm1d(48),
                    nn.ReLU(),
                    nn.Conv1d(48, 64, kernel_size=kernel_size, padding=kernel_size // 2),
                    nn.BatchNorm1d(64),
                    nn.ReLU(),
                    nn.AdaptiveMaxPool1d(1),
                )
                for kernel_size in (3, 5, 9)
            ]
        )
        self.classifier = nn.Sequential(
            nn.Dropout(0.25),
            nn.Linear(64 * 3 * 2, 128),
            nn.ReLU(),
            nn.Dropout(0.20),
            nn.Linear(128, n_classes),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        features = []
        center = x.shape[-1] // 2
        for branch in self.branches:
            h = branch[0](x)
            h = branch[1](h)
            h = branch[2](h)
            h = branch[3](h)
            h = branch[4](h)
            h = branch[5](h)
            features.append(h[:, :, center])
            features.append(torch.amax(h, dim=-1))
        return self.classifier(torch.cat(features, dim=1))


@dataclass
class CNNBaselineClassifier:
    name: str = "CNN baseline (PyTorch Conv1D)"
    random_state: int = 42
    epochs: int = 5
    batch_size: int = 256
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    early_stopping_patience: int = 2

    def __post_init__(self) -> None:
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = _SpliceCNN(n_classes=3).to(self.device)
        self.best_epoch_: int | None = None
        self.best_valid_macro_f1_: float | None = None
        self.epochs_trained_: int = 0

    def fit(
        self,
        sequences: list[str],
        labels: np.ndarray,
        validation_sequences: list[str] | None = None,
        validation_labels: np.ndarray | None = None,
    ) -> "CNNBaselineClassifier":
        torch.manual_seed(self.random_state)
        np.random.seed(self.random_state)
        self.model = _SpliceCNN(n_classes=3).to(self.device)
        self.best_epoch_ = None
        self.best_valid_macro_f1_ = None
        self.epochs_trained_ = 0
        dataset = _SequenceDataset(sequences, labels)
        generator = torch.Generator()
        generator.manual_seed(self.random_state)
        loader = DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=0,
            generator=generator,
        )
        counts = np.bincount(labels.astype(int), minlength=3).astype(np.float32)
        weights = counts.sum() / np.maximum(counts, 1.0)
        weights = weights / weights.mean()
        criterion = nn.CrossEntropyLoss(weight=torch.tensor(weights, dtype=torch.float32, device=self.device))
        optimizer = torch.optim.AdamW(
            self.model.parameters(),
            lr=self.learning_rate,
            weight_decay=self.weight_decay,
        )
        valid_loader = None
        valid_y = None
        if validation_sequences is not None and validation_labels is not None and len(validation_sequences) > 0:
            valid_dataset = _SequenceDataset(validation_sequences, validation_labels)
            valid_loader = DataLoader(valid_dataset, batch_size=self.batch_size, shuffle=False, num_workers=0)
            valid_y = validation_labels.astype(int)

        best_score = -np.inf
        best_state: dict[str, torch.Tensor] | None = None
        stale_epochs = 0
        for epoch in range(1, self.epochs + 1):
            self.model.train()
            for x, y in loader:
                x = x.to(self.device)
                y = y.to(self.device)
                optimizer.zero_grad(set_to_none=True)
                loss = criterion(self.model(x), y)
                loss.backward()
                optimizer.step()
            self.epochs_trained_ = epoch

            if valid_loader is None or valid_y is None:
                continue

            proba = self._predict_proba_loader(valid_loader)
            score = float(f1_score(valid_y, np.argmax(proba, axis=1), average="macro", zero_division=0))
            if score > best_score + 1e-12:
                best_score = score
                self.best_epoch_ = epoch
                self.best_valid_macro_f1_ = score
                best_state = {key: value.detach().cpu().clone() for key, value in self.model.state_dict().items()}
                stale_epochs = 0
            else:
                stale_epochs += 1
                if stale_epochs >= self.early_stopping_patience:
                    break

        if best_state is not None:
            self.model.load_state_dict(best_state)
        return self

    def predict_proba(self, sequences: list[str]) -> np.ndarray:
        dataset = _SequenceDataset(sequences)
        loader = DataLoader(dataset, batch_size=self.batch_size, shuffle=False, num_workers=0)
        return self._predict_proba_loader(loader)

    def _predict_proba_loader(self, loader: DataLoader) -> np.ndarray:
        chunks: list[np.ndarray] = []
        self.model.eval()
        with torch.no_grad():
            for batch in loader:
                x = batch[0] if isinstance(batch, (list, tuple)) else batch
                x = x.to(self.device)
                proba = torch.softmax(self.model(x), dim=1).cpu().numpy()
                chunks.append(proba)
        return np.vstack(chunks)
