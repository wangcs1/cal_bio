from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd

from src.build_synthetic_splice_dataset import build_and_write
from src.utils import (
    BASES,
    PROJECT_ROOT,
    choose_alt,
    mutate_base,
    read_csv,
    stable_id,
    write_dataframe,
)


def ensure_splice_dataset(path: Path) -> pd.DataFrame:
    if not path.exists():
        build_and_write()
    return read_csv(path)


def add_variant_row(
    rows: list[dict[str, object]],
    source: pd.Series,
    variant_type: str,
    label: int,
    pos_index: int,
    alt: str,
    target_class: int,
    rng: random.Random,
) -> None:
    wt = str(source["sequence"])
    ref = wt[pos_index]
    if alt == ref:
        alt = choose_alt(ref, rng)
    mut = mutate_base(wt, pos_index, alt)
    c = len(wt) // 2
    rel = pos_index - c
    rows.append(
        {
            "variant_id": stable_id(source["sample_id"], variant_type, rel, ref, alt, prefix="var"),
            "sample_id": source["sample_id"],
            "chrom": source["chrom"],
            "pos": int(source["center"]) + rel,
            "relative_pos": rel,
            "strand": source["strand"],
            "ref": ref,
            "alt": alt,
            "variant_type": variant_type,
            "label": label,
            "label_name": "splice_altering" if label else "neutral",
            "wt_sequence": wt,
            "mut_sequence": mut,
            "target_class": target_class,
            "target_class_name": {0: "donor", 1: "acceptor", 2: "non_splice"}[target_class],
            "gene_id": source["gene_id"],
            "transcript_id": source["transcript_id"],
            "data_source": "synthetic_splice_benchmark_v1",
        }
    )


def find_one_snv_to_motif(sequence: str, rng: random.Random) -> tuple[int, str, int] | None:
    c = len(sequence) // 2
    candidates: list[tuple[int, str, int]] = []
    for start in range(max(0, c - 70), min(len(sequence) - 1, c + 70)):
        dinuc = sequence[start : start + 2]
        for motif, target_class in (("GT", 0), ("AG", 1)):
            mismatches = [i for i, (a, b) in enumerate(zip(dinuc, motif)) if a != b]
            if len(mismatches) == 1:
                idx = start + mismatches[0]
                if abs(idx - c) <= 2:
                    continue
                candidates.append((idx, motif[mismatches[0]], target_class))
    if not candidates:
        return None
    return rng.choice(candidates)


def build_variant_dataset(source: pd.DataFrame, per_type: int = 90, seed: int = 2026) -> pd.DataFrame:
    rng = random.Random(seed)
    frame = source[source["split"] == "test"].copy() if "split" in source.columns else source.copy()
    if frame.empty:
        frame = source.copy()

    rows: list[dict[str, object]] = []
    donors = frame[frame["label"].astype(int) == 0].sample(
        n=min(per_type, int((frame["label"].astype(int) == 0).sum())), random_state=seed
    )
    acceptors = frame[frame["label"].astype(int) == 1].sample(
        n=min(per_type, int((frame["label"].astype(int) == 1).sum())), random_state=seed + 1
    )
    non_splice = frame[frame["label"].astype(int) == 2].sample(
        n=min(per_type * 2, int((frame["label"].astype(int) == 2).sum())), random_state=seed + 2
    )

    for _, row in donors.iterrows():
        c = len(str(row["sequence"])) // 2
        add_variant_row(rows, row, "donor_loss", 1, c + 1, "A", 0, rng)

    for _, row in acceptors.iterrows():
        c = len(str(row["sequence"])) // 2
        add_variant_row(rows, row, "acceptor_loss", 1, c - 2, "T", 1, rng)

    cryptic_count = 0
    for _, row in non_splice.iterrows():
        found = find_one_snv_to_motif(str(row["sequence"]), rng)
        if found is None:
            continue
        pos_index, alt, target_class = found
        add_variant_row(rows, row, "cryptic_gain", 1, pos_index, alt, target_class, rng)
        cryptic_count += 1
        if cryptic_count >= per_type:
            break

    neutral_pool = frame.sample(n=min(per_type * 2, len(frame)), random_state=seed + 3)
    for _, row in neutral_pool.iterrows():
        seq = str(row["sequence"])
        c = len(seq) // 2
        far_positions = list(range(max(0, c - 180), max(0, c - 90))) + list(
            range(min(len(seq), c + 90), min(len(seq), c + 180))
        )
        pos_index = rng.choice(far_positions)
        add_variant_row(rows, row, "neutral_far_snv", 0, pos_index, choose_alt(seq[pos_index], rng), 2, rng)

    variants = pd.DataFrame(rows).sample(frac=1.0, random_state=seed).reset_index(drop=True)
    return variants


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build artificial splice variant effect data.")
    parser.add_argument("--input", type=Path, default=PROJECT_ROOT / "data/processed/splice_sites_pm200.csv")
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "data/processed/artificial_variant_effect.csv")
    parser.add_argument("--per-type", type=int, default=90)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def build_and_write_variants(
    input_path: Path = PROJECT_ROOT / "data/processed/splice_sites_pm200.csv",
    output_path: Path = PROJECT_ROOT / "data/processed/artificial_variant_effect.csv",
    per_type: int = 90,
    seed: int = 2026,
) -> pd.DataFrame:
    source = ensure_splice_dataset(input_path)
    variants = build_variant_dataset(source, per_type=per_type, seed=seed)
    write_dataframe(output_path, variants)
    summary = variants.groupby(["variant_type", "label_name"]).size().reset_index(name="rows")
    write_dataframe(output_path.with_name("artificial_variant_effect_summary.csv"), summary)
    return variants


def main() -> None:
    args = parse_args()
    variants = build_and_write_variants(args.input, args.output, args.per_type, args.seed)
    print(variants.groupby(["variant_type", "label_name"]).size().reset_index(name="rows").to_string(index=False))


if __name__ == "__main__":
    main()
