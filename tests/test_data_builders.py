from __future__ import annotations

import pandas as pd

from src.data.build_synthetic_splice_dataset import build_dataset, write_windowed_outputs
from src.data.build_variant_dataset import build_variant_dataset
from src.utils import REQUIRED_SPLIT_COLUMNS


def test_synthetic_builder_keeps_required_columns(tmp_path):
    frame = build_dataset(samples_per_class=20, max_flank=400, seed=7)
    out = tmp_path / "processed"
    splits = tmp_path / "splits"
    write_windowed_outputs(frame, [200], out, splits)
    train = pd.read_csv(splits / "train.csv")
    assert REQUIRED_SPLIT_COLUMNS.issubset(train.columns)
    assert set(train["label_name"]).issubset({"donor", "acceptor", "non_splice"})


def test_variant_builder_splits_gain_types():
    frame = build_dataset(samples_per_class=80, max_flank=200, seed=11)
    variants = build_variant_dataset(frame, per_type=10, seed=11)
    assert {"donor_loss", "acceptor_loss", "neutral_far_snv"}.issubset(set(variants["variant_type"]))
    assert {"donor_gain", "acceptor_gain"}.intersection(set(variants["variant_type"]))
    assert {"wt_sequence", "mut_sequence", "target_class"}.issubset(variants.columns)
