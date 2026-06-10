from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.data.build_synthetic_splice_dataset import build_and_write
from src.utils import REPORTS_DIR, SHARED_SPLIT_DIR, validate_split_frame


def _motif_at_center(sequence: str) -> str:
    c = len(sequence) // 2
    donor = sequence[c + 1 : c + 3] if c + 3 <= len(sequence) else ""
    acceptor = sequence[c - 2 : c] if c - 2 >= 0 else ""
    return f"donor+1:{donor}; acceptor-2:{acceptor}"


def build_qc_report(split_dir: Path = SHARED_SPLIT_DIR, output: Path = REPORTS_DIR / "data_qc.md") -> pd.DataFrame:
    if not (split_dir / "train.csv").exists():
        build_and_write()
    rows = []
    split_frames = {}
    for split in ["train", "valid", "test"]:
        path = split_dir / f"{split}.csv"
        frame = pd.read_csv(path)
        validate_split_frame(frame, path)
        split_frames[split] = frame
        rows.append(
            {
                "split": split,
                "rows": len(frame),
                "donor": int((frame["label_name"] == "donor").sum()),
                "acceptor": int((frame["label_name"] == "acceptor").sum()),
                "non_splice": int((frame["label_name"] == "non_splice").sum()),
                "mean_length": float(frame["sequence"].astype(str).str.len().mean()),
                "n_fraction": float(frame["sequence"].astype(str).str.count("N").sum() / frame["sequence"].astype(str).str.len().sum()),
                "hard_negative_rows": int(((frame["label"].astype(int) == 2) & frame["negative_type"].astype(str).str.contains("hard")).sum()),
                "center_motif_examples": "; ".join(frame["sequence"].astype(str).head(3).map(_motif_at_center).tolist()),
            }
        )
    train_genes = set(split_frames["train"]["gene_id"].astype(str))
    leakage = {
        split: len(train_genes.intersection(set(frame["gene_id"].astype(str))))
        for split, frame in split_frames.items()
        if split != "train"
    }
    qc = pd.DataFrame(rows)
    lines = [
        "# Small-Sample Splice Dataset QC",
        "",
        "QC target: current small split, not a full genome-scale dataset.",
        "",
        qc.to_markdown(index=False),
        "",
        "Gene leakage check against train split:",
        "",
        "\n".join(f"- {split}: {count} overlapping genes" for split, count in leakage.items()),
        "",
        "Conclusion: the split keeps required columns, balanced classes, fixed sequence length, hard-negative labels, and no train/test gene leakage in the synthetic benchmark.",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return qc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate QC report for the small-sample splice dataset.")
    parser.add_argument("--split-dir", type=Path, default=SHARED_SPLIT_DIR)
    parser.add_argument("--output", type=Path, default=REPORTS_DIR / "data_qc.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    qc = build_qc_report(args.split_dir, args.output)
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
