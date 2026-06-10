from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.data.build_variant_dataset import build_and_write_variants
from src.utils import EXP3_DATA_DIR, exp3_data_file, read_csv, write_dataframe


def build_sqtl_smoke(output: Path = EXP3_DATA_DIR / "gtex_sqtl_variants_smoke.csv", rows: int = 10) -> pd.DataFrame:
    if not exp3_data_file("artificial_variant_effect.csv").exists():
        build_and_write_variants()
    variants = read_csv(exp3_data_file("artificial_variant_effect.csv")).head(rows).copy()
    tissues = ["brain", "heart", "liver", "muscle", "blood"]
    payload = []
    for idx, variant in variants.iterrows():
        direction = "increased_usage" if "gain" in str(variant["variant_type"]) else "decreased_usage"
        payload.append(
            {
                "variant_id": variant["variant_id"],
                "tissue": tissues[idx % len(tissues)],
                "target_gene": variant["gene_id"],
                "target_junction": f"{variant['gene_id']}:J{idx % 3}",
                "observed_effect_direction": direction,
                "model_delta_score": 0.0,
                "variant_type": variant["variant_type"],
                "wt_sequence": variant["wt_sequence"],
                "mut_sequence": variant["mut_sequence"],
                "data_source": "gtex_sqtl_smoke_case_study_proxy",
            }
        )
    frame = pd.DataFrame(payload)
    write_dataframe(output, frame)
    return frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a small GTEx/sQTL-style case-study variant table.")
    parser.add_argument("--output", type=Path, default=EXP3_DATA_DIR / "gtex_sqtl_variants_smoke.csv")
    parser.add_argument("--rows", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame = build_sqtl_smoke(args.output, args.rows)
    print(frame[["variant_id", "tissue", "target_gene", "observed_effect_direction"]].to_string(index=False))


if __name__ == "__main__":
    main()
