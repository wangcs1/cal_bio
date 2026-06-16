from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd

from src.utils import (
    ALL_CHROMS,
    LABELS,
    PROJECT_ROOT,
    SHARED_PROCESSED_DIR,
    SHARED_SPLIT_DIR,
    crop_center,
    ensure_dirs,
    insert_motif,
    random_dna,
    split_name,
    stable_id,
    validate_split_frame,
    write_dataframe,
)


DEFAULT_WINDOWS = [50, 100, 200, 400]


def add_donor_context(sequence: str) -> str:
    c = len(sequence) // 2
    sequence = insert_motif(sequence, c - 2, "CAGGTAAGT")
    sequence = insert_motif(sequence, c - 44, "GAAGAA")
    sequence = insert_motif(sequence, c + 24, "TGTGTA")
    return sequence


def add_acceptor_context(sequence: str) -> str:
    c = len(sequence) // 2
    sequence = insert_motif(sequence, c - 35, "TTCTTTCTCTTTTCTTTCTT")
    sequence = insert_motif(sequence, c - 2, "AG")
    sequence = insert_motif(sequence, c - 49, "TACTAAC")
    sequence = insert_motif(sequence, c + 16, "GAAGAA")
    return sequence


def add_non_splice_context(sequence: str, rng: random.Random, hard: bool) -> tuple[str, str]:
    c = len(sequence) // 2
    if hard:
        if rng.random() < 0.5:
            sequence = insert_motif(sequence, c + 1, "GT")
            return sequence, "hard_gt"
        sequence = insert_motif(sequence, c - 2, "AG")
        return sequence, "hard_ag"
    sequence = insert_motif(sequence, c - 20, "ACACAC")
    sequence = insert_motif(sequence, c + 20, "CACACA")
    return sequence, "easy_random"


def make_row(label: int, idx: int, max_flank: int, rng: random.Random) -> dict[str, object]:
    chrom = rng.choice(ALL_CHROMS)
    strand = rng.choice(["+", "-"])
    gc = rng.uniform(0.36, 0.54)
    sequence = random_dna(2 * max_flank + 1, rng, gc=gc)
    negative_type = "positive"
    motif_type = "canonical"
    if label == 0:
        sequence = add_donor_context(sequence)
    elif label == 1:
        sequence = add_acceptor_context(sequence)
    else:
        hard = idx % 3 != 0
        sequence, motif_type = add_non_splice_context(sequence, rng, hard=hard)
        negative_type = "hard_gtag" if hard else "easy_random"

    split = split_name(chrom)
    gene_number = idx // 3
    gene_id = f"SYN_{split.upper()}_{chrom.replace('chr', 'C')}_GENE_{gene_number:04d}"
    transcript_id = f"{gene_id}.T1"
    center = 1_000_000 + idx * 17 + max_flank
    sample_id = stable_id(label, idx, chrom, strand, prefix="syn")
    return {
        "sample_id": sample_id,
        "chrom": chrom,
        "start": center - max_flank,
        "end": center + max_flank,
        "strand": strand,
        "center": center,
        "label": label,
        "label_name": LABELS[label],
        "negative_type": negative_type,
        "motif_type": motif_type,
        "sequence": sequence,
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "data_source": "synthetic_splice_benchmark_v1",
        "split": split,
    }


def build_dataset(samples_per_class: int, max_flank: int, seed: int) -> pd.DataFrame:
    rng = random.Random(seed)
    rows = []
    for label in (0, 1, 2):
        for idx in range(samples_per_class):
            rows.append(make_row(label, idx + label * samples_per_class, max_flank, rng))
    frame = pd.DataFrame(rows)
    frame = frame[frame["split"].isin(["train", "valid", "test"])].copy()
    frame = frame.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    return frame


def write_windowed_outputs(frame: pd.DataFrame, windows: list[int], out_dir: Path, split_dir: Path) -> None:
    ensure_dirs(out_dir, split_dir)
    for window in windows:
        current = frame.copy()
        current["start"] = current["center"].astype(int) - window
        current["end"] = current["center"].astype(int) + window
        current["sequence"] = current["sequence"].astype(str).map(lambda seq: crop_center(seq, window))
        validate_split_frame(current, f"synthetic window pm{window}")
        path = out_dir / f"splice_sites_pm{window}.csv"
        write_dataframe(path, current.drop(columns=[]))
        for split in ("train", "valid", "test"):
            write_dataframe(split_dir / f"{split}_pm{window}.csv", current[current["split"] == split])
        if window == 200:
            for split in ("train", "valid", "test"):
                write_dataframe(split_dir / f"{split}.csv", current[current["split"] == split])
            write_dataframe(split_dir / "cross_gene_test.csv", current[current["split"] == "test"])


def summarize(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for split, split_frame in frame.groupby("split"):
        for label, label_frame in split_frame.groupby("label_name"):
            rows.append({"split": split, "label_name": label, "rows": len(label_frame)})
    return pd.DataFrame(rows).sort_values(["split", "label_name"])


def build_rare_motif_dataset(source: pd.DataFrame, max_rows_per_type: int = 40) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    donor_source = source[source["label"].astype(int) == 0].head(max_rows_per_type)
    acceptor_source = source[source["label"].astype(int) == 1].head(max_rows_per_type)
    for idx, (_, row) in enumerate(donor_source.iterrows()):
        seq = str(row["sequence"])
        c = len(seq) // 2
        rare = insert_motif(seq, c + 1, "GC")
        current = row.to_dict()
        current.update(
            {
                "sample_id": stable_id(row["sample_id"], "rare_gc_ag", idx, prefix="rare"),
                "sequence": rare,
                "motif_type": "rare_GC-AG_donor",
                "data_source": "synthetic_rare_motif_case_study_v1",
            }
        )
        rows.append(current)
    for idx, (_, row) in enumerate(acceptor_source.iterrows()):
        seq = str(row["sequence"])
        c = len(seq) // 2
        rare = insert_motif(seq, c - 2, "AC")
        current = row.to_dict()
        current.update(
            {
                "sample_id": stable_id(row["sample_id"], "rare_at_ac", idx, prefix="rare"),
                "sequence": rare,
                "motif_type": "rare_AT-AC_acceptor",
                "data_source": "synthetic_rare_motif_case_study_v1",
            }
        )
        rows.append(current)
    return pd.DataFrame(rows).reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a local synthetic splice-site benchmark.")
    parser.add_argument("--samples-per-class", type=int, default=420)
    parser.add_argument("--max-flank", type=int, default=400)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--out-dir", type=Path, default=SHARED_PROCESSED_DIR)
    parser.add_argument("--split-dir", type=Path, default=SHARED_SPLIT_DIR)
    parser.add_argument("--windows", type=int, nargs="+", default=DEFAULT_WINDOWS)
    return parser.parse_args()


def build_and_write(
    samples_per_class: int = 420,
    max_flank: int = 400,
    seed: int = 2026,
    out_dir: Path = SHARED_PROCESSED_DIR,
    split_dir: Path = SHARED_SPLIT_DIR,
    windows: list[int] | None = None,
) -> pd.DataFrame:
    windows = DEFAULT_WINDOWS if windows is None else windows
    frame = build_dataset(samples_per_class=samples_per_class, max_flank=max_flank, seed=seed)
    write_dataframe(out_dir / "synthetic_splice_sites_master_pm400.csv", frame)
    write_windowed_outputs(frame, windows, out_dir, split_dir)
    write_dataframe(out_dir / "synthetic_splice_sites_summary.csv", summarize(frame))
    rare = build_rare_motif_dataset(frame, max_rows_per_type=40)
    rare["sequence"] = rare["sequence"].astype(str).map(lambda seq: crop_center(seq, 200))
    write_dataframe(out_dir / "rare_motif_splice_sites.csv", rare)
    return frame


def main() -> None:
    args = parse_args()
    frame = build_and_write(
        samples_per_class=args.samples_per_class,
        max_flank=args.max_flank,
        seed=args.seed,
        out_dir=args.out_dir,
        split_dir=args.split_dir,
        windows=args.windows,
    )
    counts = frame.groupby(["split", "label_name"]).size().reset_index(name="rows")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
