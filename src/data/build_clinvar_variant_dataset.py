from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.data.build_variant_dataset import build_and_write_variants
from src.utils import EXP3_DATA_DIR, PROJECT_ROOT, exp3_data_file, read_csv, write_dataframe


def build_clinvar_smoke(output: Path = EXP3_DATA_DIR / "clinvar_splicing_variants_smoke.csv", rows: int = 12) -> pd.DataFrame:
    if not exp3_data_file("artificial_variant_effect.csv").exists():
        build_and_write_variants()
    variants = read_csv(exp3_data_file("artificial_variant_effect.csv"))
    positives = variants[variants["label"].astype(int) == 1].head(rows // 2).copy()
    controls = variants[variants["label"].astype(int) == 0].head(rows - len(positives)).copy()
    frame = pd.concat([positives, controls], ignore_index=True)
    frame["clinvar_label_source"] = frame["label"].map({1: "splice-related synthetic positive", 0: "matched neutral control"})
    frame["data_source"] = "clinvar_format_smoke_from_artificial_variants"
    raw_hint = PROJECT_ROOT / "data/raw/clinvar_smoke.csv"
    frame["raw_resource_hint"] = str(raw_hint) if raw_hint.exists() else "optional data/raw/clinvar_smoke.csv not required"
    write_dataframe(output, frame)
    return frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a small ClinVar-style splice variant smoke dataset.")
    parser.add_argument("--output", type=Path, default=EXP3_DATA_DIR / "clinvar_splicing_variants_smoke.csv")
    parser.add_argument("--rows", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame = build_clinvar_smoke(args.output, args.rows)
    print(frame[["variant_id", "variant_type", "label_name", "clinvar_label_source"]].to_string(index=False))


if __name__ == "__main__":
    main()
