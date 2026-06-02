#!/usr/bin/env python3
"""
Split splice-site samples by chromosome to avoid same-gene leakage.

Default input is the ±200 nt dataset:
  data/processed/splice_sites_pm200.csv

Default output:
  data/splits/train.csv
  data/splits/valid.csv
  data/splits/test.csv
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


DEFAULT_SPLITS = {
    "train": {f"chr{i}" for i in range(1, 17)},
    "valid": {"chr17", "chr18"},
    "test": {"chr19", "chr20", "chr21", "chr22", "chrX"},
}
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def parse_chroms(value: str) -> set[str]:
    chroms: set[str] = set()
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        chroms.add(token if token.startswith("chr") else f"chr{token}")
    return chroms


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT_ROOT / "data/processed/splice_sites_pm200.csv",
        help="Processed splice-site CSV to split.",
    )
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "data/splits")
    parser.add_argument("--train-chroms", default=",".join(sorted(DEFAULT_SPLITS["train"])))
    parser.add_argument("--valid-chroms", default="chr17,chr18")
    parser.add_argument("--test-chroms", default="chr19,chr20,chr21,chr22,chrX")
    parser.add_argument(
        "--drop-unassigned",
        dest="drop_unassigned",
        action="store_true",
        default=True,
        help="Drop rows whose chromosome is not in any split.",
    )
    parser.add_argument(
        "--fail-on-unassigned",
        dest="drop_unassigned",
        action="store_false",
        help="Fail if any row belongs to a chromosome outside the configured splits.",
    )
    return parser.parse_args()


def split_rows(args: argparse.Namespace) -> tuple[dict[str, list[dict[str, str]]], Counter[str]]:
    split_chroms = {
        "train": parse_chroms(args.train_chroms),
        "valid": parse_chroms(args.valid_chroms),
        "test": parse_chroms(args.test_chroms),
    }
    chrom_to_split: dict[str, str] = {}
    for split_name, chroms in split_chroms.items():
        overlap = set(chrom_to_split).intersection(chroms)
        if overlap:
            raise ValueError(f"Chromosomes assigned to multiple splits: {sorted(overlap)}")
        chrom_to_split.update({chrom: split_name for chrom in chroms})

    rows_by_split: dict[str, list[dict[str, str]]] = {"train": [], "valid": [], "test": []}
    unassigned = Counter()

    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Input CSV has no header: {args.input}")
        for row in reader:
            chrom = row["chrom"]
            split_name = chrom_to_split.get(chrom)
            if split_name is None:
                unassigned[chrom] += 1
                continue
            rows_by_split[split_name].append(row)

    if unassigned and not args.drop_unassigned:
        preview = ", ".join(f"{chrom}:{count}" for chrom, count in unassigned.most_common(10))
        raise ValueError(f"Unassigned chromosomes found ({preview}); pass --drop-unassigned to ignore them.")

    return rows_by_split, unassigned


def write_split(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize(split_name: str, rows: list[dict[str, str]]) -> str:
    labels = Counter(row["label"] for row in rows)
    chroms = Counter(row["chrom"] for row in rows)
    label_text = ", ".join(f"{label}:{labels[label]}" for label in sorted(labels))
    chrom_text = ", ".join(f"{chrom}:{chroms[chrom]}" for chrom in sorted(chroms))
    return f"{split_name}: rows={len(rows):,}; labels={{{label_text}}}; chroms={{{chrom_text}}}"


def main() -> None:
    args = parse_args()
    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Input CSV has no header: {args.input}")
        fieldnames = reader.fieldnames

    rows_by_split, unassigned = split_rows(args)
    for split_name in ("train", "valid", "test"):
        write_split(args.out_dir / f"{split_name}.csv", rows_by_split[split_name], fieldnames)
        print(summarize(split_name, rows_by_split[split_name]))
    if unassigned:
        preview = ", ".join(f"{chrom}:{count}" for chrom, count in unassigned.most_common(10))
        print(f"dropped_unassigned: rows={sum(unassigned.values()):,}; chroms={{{preview}}}")


if __name__ == "__main__":
    main()
