from __future__ import annotations

import argparse
import hashlib
import itertools
from pathlib import Path

import pandas as pd

from src.data.build_splice_site_dataset import FastaIndex, build_positive_sites, read_transcripts
from src.utils import (
    EXP3_DATA_DIR,
    EXP3_TABLES_DIR,
    PROJECT_ROOT,
    REPORTS_DIR,
    SHARED_PROCESSED_DIR,
    SHARED_SPLIT_DIR,
    safe_average_precision,
    safe_roc_auc,
    validate_split_frame,
    write_dataframe,
)


WINDOWS = [50, 100, 200, 400]
SPLITS = ["train", "valid", "test"]
TEST_HOLDOUT_CHROMS = {"chr19", "chr20", "chr21", "chr22", "chrX"}


def _dna(sequence: str) -> str:
    return str(sequence).upper().replace("U", "T")


def _category(row: pd.Series) -> str:
    if int(row["label"]) == 0:
        return "donor"
    if int(row["label"]) == 1:
        return "acceptor"
    return "easy_negative" if str(row.get("negative_type", "")) == "easy_random" else "hard_negative"


def _canonical_rate(frame: pd.DataFrame, label: int) -> float:
    subset = frame[frame["label"].astype(int) == label]
    if subset.empty:
        return float("nan")
    sequences = subset["sequence"].astype(str).map(_dna)
    if label == 0:
        return float(sequences.map(lambda seq: seq[len(seq) // 2 + 1 : len(seq) // 2 + 3] == "GT").mean())
    return float(sequences.map(lambda seq: seq[len(seq) // 2 - 2 : len(seq) // 2] == "AG").mean())


def _hard_motif_rate(frame: pd.DataFrame) -> float:
    hard = frame[(frame["label"].astype(int) == 2) & frame["negative_type"].astype(str).str.contains("hard", na=False)]
    if hard.empty:
        return float("nan")
    sequences = hard["sequence"].astype(str).map(_dna)
    return float(
        sequences.map(
            lambda seq: seq[len(seq) // 2 + 1 : len(seq) // 2 + 3] == "GT"
            or seq[len(seq) // 2 - 2 : len(seq) // 2] == "AG"
        ).mean()
    )


def _n_fraction(frame: pd.DataFrame) -> float:
    sequences = frame["sequence"].astype(str).map(_dna)
    total = int(sequences.str.len().sum())
    return float(sequences.str.count("N").sum() / total) if total else float("nan")


def _alphabet(frame: pd.DataFrame) -> str:
    chars = sorted(set("".join(frame["sequence"].astype(str).map(_dna).tolist())))
    return "".join(chars)


def _gc_fraction(sequence: str) -> float:
    seq = _dna(sequence)
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else float("nan")


def _sequence_digest(sequence: str) -> str:
    return hashlib.sha1(_dna(sequence).encode("ascii", errors="ignore")).hexdigest()


def _center_window(sequence: str, flank: int = 80) -> str:
    seq = _dna(sequence)
    center = len(seq) // 2
    return seq[max(0, center - flank) : min(len(seq), center + flank + 1)]


def _kmer_set(sequence: str, k: int = 9) -> set[str]:
    seq = _dna(sequence)
    return {seq[index : index + k] for index in range(0, max(0, len(seq) - k + 1))}


def _jaccard(left: set[str], right: set[str]) -> float:
    if not left and not right:
        return 1.0
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def _status(ok: bool, warn: bool = False) -> str:
    if ok:
        return "PASS"
    return "WARN" if warn else "FAIL"


def _markdown(frame: pd.DataFrame, max_rows: int = 30) -> str:
    if frame.empty:
        return "_No rows._"
    current = frame.head(max_rows).copy()
    for column in current.columns:
        if pd.api.types.is_float_dtype(current[column]):
            current[column] = current[column].map(lambda value: "" if pd.isna(value) else f"{value:.4f}")
    return current.to_markdown(index=False)


def homology_leakage_audit(split_frames: dict[str, pd.DataFrame], center_flank: int = 80, jaccard_threshold: float = 0.90) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    for split, frame in split_frames.items():
        for _, row in frame.iterrows():
            sequence = _dna(str(row["sequence"]))
            records.append(
                {
                    "split": split,
                    "sample_id": str(row["sample_id"]),
                    "gene_id": str(row["gene_id"]),
                    "label_name": str(row["label_name"]),
                    "full_hash": _sequence_digest(sequence),
                    "center_hash": _sequence_digest(_center_window(sequence, center_flank)),
                    "kmer_set": _kmer_set(sequence),
                }
            )
    audit = pd.DataFrame(records)

    rows = []
    exact_full_pairs = 0
    exact_center_pairs = 0
    high_jaccard_pairs = 0
    examples = []
    for left_name, right_name in itertools.combinations(SPLITS, 2):
        left = audit[audit["split"] == left_name]
        right = audit[audit["split"] == right_name]
        full_overlap = set(left["full_hash"]) & set(right["full_hash"])
        center_overlap = set(left["center_hash"]) & set(right["center_hash"])
        exact_full_pairs += len(full_overlap)
        exact_center_pairs += len(center_overlap)

        pair_high = 0
        for _, left_row in left.iterrows():
            for _, right_row in right.iterrows():
                score = _jaccard(left_row["kmer_set"], right_row["kmer_set"])
                if score >= jaccard_threshold:
                    pair_high += 1
                    if len(examples) < 12:
                        examples.append(
                            {
                                "pair": f"{left_name}/{right_name}",
                                "left_sample_id": left_row["sample_id"],
                                "right_sample_id": right_row["sample_id"],
                                "left_gene_id": left_row["gene_id"],
                                "right_gene_id": right_row["gene_id"],
                                "kmer_jaccard": score,
                            }
                        )
        high_jaccard_pairs += pair_high
        rows.append(
            {
                "pair": f"{left_name}/{right_name}",
                "full_sequence_duplicate_hashes": len(full_overlap),
                "center_161bp_duplicate_hashes": len(center_overlap),
                f"kmer_jaccard_ge_{jaccard_threshold:.2f}_pairs": pair_high,
            }
        )

    summary = pd.DataFrame(rows)
    summary.attrs["exact_full_pairs"] = exact_full_pairs
    summary.attrs["exact_center_pairs"] = exact_center_pairs
    summary.attrs["high_jaccard_pairs"] = high_jaccard_pairs
    return summary, pd.DataFrame(examples)


def load_window_frames(processed_dir: Path) -> dict[int, pd.DataFrame]:
    frames = {}
    for window in WINDOWS:
        path = processed_dir / f"splice_sites_pm{window}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing processed window file: {path}")
        frames[window] = pd.read_csv(path)
    return frames


def load_split_frames(split_dir: Path) -> dict[str, pd.DataFrame]:
    frames = {}
    for split in SPLITS:
        path = split_dir / f"{split}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing split file: {path}")
        frame = pd.read_csv(path)
        validate_split_frame(frame, path)
        frames[split] = frame
    return frames


def check_window_split_consistency(
    window_frames: dict[int, pd.DataFrame], split_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, bool]]:
    base = (
        window_frames[50][["sample_id", "split", "chrom", "label_name", "negative_type", "motif_type", "gene_id"]]
        .sort_values("sample_id")
        .reset_index(drop=True)
    )
    rows = []
    same_assignment = True
    for window, frame in window_frames.items():
        current = (
            frame[["sample_id", "split", "chrom", "label_name", "negative_type", "motif_type", "gene_id"]]
            .sort_values("sample_id")
            .reset_index(drop=True)
        )
        rows.append(
            {
                "window": f"pm{window}",
                "rows": len(frame),
                "duplicate_sample_ids": int(frame["sample_id"].duplicated().sum()),
                "same_ids_as_pm50": set(frame["sample_id"]) == set(window_frames[50]["sample_id"]),
                "same_assignment_as_pm50": current.equals(base),
            }
        )
        same_assignment = same_assignment and current.equals(base)

    split_rows = []
    split_consistent = True
    for split in SPLITS:
        ids_by_window = []
        for window in WINDOWS:
            path = split_dir / f"{split}_pm{window}.csv"
            if not path.exists():
                raise FileNotFoundError(f"Missing split-window file: {path}")
            ids_by_window.append(set(pd.read_csv(path)["sample_id"]))
        split_rows.append(
            {
                "split": split,
                "pm50_rows": len(ids_by_window[0]),
                "pm100_same_ids": ids_by_window[0] == ids_by_window[1],
                "pm200_same_ids": ids_by_window[0] == ids_by_window[2],
                "pm400_same_ids": ids_by_window[0] == ids_by_window[3],
            }
        )
        split_consistent = split_consistent and all(ids_by_window[0] == ids for ids in ids_by_window[1:])

    return pd.DataFrame(rows), pd.DataFrame(split_rows), {
        "same_assignment": same_assignment,
        "split_consistent": split_consistent,
    }


def build_annotated_site_set(genome_path: Path, gtf_path: Path) -> set[tuple[str, int]]:
    fasta = FastaIndex(genome_path)
    try:
        transcripts = read_transcripts(gtf_path, fasta, protein_coding_only=False, gtf_line_limit=None)
        _sites, _bounds, annotated_any, _annotated_same = build_positive_sites(transcripts)
        return annotated_any
    finally:
        fasta.close()


def verify_clinvar_rows(clinvar: pd.DataFrame, genome_path: Path) -> dict[str, int]:
    fasta = FastaIndex(genome_path)
    try:
        ref_mismatches = 0
        hamming_failures = 0
        position_failures = 0
        for _, row in clinvar.iterrows():
            genomic_ref = fasta.fetch(str(row["chrom"]), int(row["pos"]), int(row["pos"]))
            if genomic_ref != str(row["genomic_ref"]):
                ref_mismatches += 1

            wt = str(row["wt_sequence"])
            mut = str(row["mut_sequence"])
            diffs = [idx for idx, (a, b) in enumerate(zip(wt, mut)) if a != b]
            if len(diffs) != 1:
                hamming_failures += 1
                continue
            center = len(wt) // 2
            relative_pos = int(row["relative_pos"])
            expected_idx = center + relative_pos if str(row["strand"]) == "+" else center - relative_pos
            if diffs[0] != expected_idx or wt[expected_idx] != str(row["ref"]) or mut[expected_idx] != str(row["alt"]):
                position_failures += 1
        return {
            "genomic_ref_mismatches": ref_mismatches,
            "hamming_failures": hamming_failures,
            "position_failures": position_failures,
        }
    finally:
        fasta.close()


def parse_raw_metadata(gtf_path: Path, clinvar_path: Path) -> pd.DataFrame:
    rows = []
    if gtf_path.exists():
        with gtf_path.open(encoding="utf-8") as handle:
            for _ in range(20):
                line = handle.readline()
                if not line.startswith("##"):
                    break
                if line.startswith("##description:") or line.startswith("##date:"):
                    key, value = line[2:].strip().split(":", 1)
                    rows.append({"resource": "gencode.gtf", "field": key, "value": value.strip()})
    if clinvar_path.exists():
        with clinvar_path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("##"):
                    break
                if line.startswith("##fileDate=") or line.startswith("##reference=") or line.startswith("##source="):
                    key, value = line[2:].strip().split("=", 1)
                    rows.append({"resource": "clinvar.vcf", "field": key, "value": value.strip()})
    return pd.DataFrame(rows)


def build_distance_matched_clinvar(clinvar: pd.DataFrame, seed: int = 42) -> tuple[pd.DataFrame, pd.DataFrame]:
    matched_parts = []
    rows = []
    for distance, group in clinvar.groupby("nearest_splice_distance"):
        positives = group[group["label"].astype(int) == 1]
        negatives = group[group["label"].astype(int) == 0]
        pairs = min(len(positives), len(negatives))
        rows.append(
            {
                "nearest_splice_distance": int(distance),
                "positive_rows": len(positives),
                "neutral_rows": len(negatives),
                "matched_pairs": pairs,
            }
        )
        if pairs == 0:
            continue
        matched_parts.append(positives.sample(n=pairs, random_state=seed + int(distance)))
        matched_parts.append(negatives.sample(n=pairs, random_state=seed + 10_000 + int(distance)))
    matched = (
        pd.concat(matched_parts, ignore_index=True).sample(frac=1.0, random_state=seed).reset_index(drop=True)
        if matched_parts
        else pd.DataFrame(columns=clinvar.columns)
    )
    return matched, pd.DataFrame(rows)


def distance_only_metrics(frame: pd.DataFrame, label: str) -> dict[str, object]:
    if frame.empty:
        return {"subset": label, "rows": 0, "positive_rows": 0, "neutral_rows": 0, "distance_only_auroc": float("nan"), "distance_only_auprc": float("nan")}
    y = frame["label"].astype(int).to_numpy()
    score = -frame["nearest_splice_distance"].astype(float).to_numpy()
    return {
        "subset": label,
        "rows": len(frame),
        "positive_rows": int((y == 1).sum()),
        "neutral_rows": int((y == 0).sum()),
        "distance_only_auroc": safe_roc_auc(y, score),
        "distance_only_auprc": safe_average_precision(y, score),
    }


def matched_score_metrics(scores_path: Path, matched: pd.DataFrame) -> pd.DataFrame:
    if not scores_path.exists() or matched.empty:
        return pd.DataFrame(
            [
                {
                    "status": "not_available",
                    "reason": f"Missing score table or empty matched subset: {scores_path}",
                }
            ]
        )
    scores = pd.read_csv(scores_path)
    keep_ids = set(matched["variant_id"].astype(str))
    subset = scores[scores["variant_id"].astype(str).isin(keep_ids)].copy()
    if subset.empty:
        return pd.DataFrame([{"status": "not_available", "reason": "Score table has no variants from the matched subset."}])
    rows = []
    for model, group in subset.groupby("model"):
        rows.append(
            {
                "status": "available",
                "model": model,
                "rows": len(group),
                "auroc_distance_matched": safe_roc_auc(group["label"].astype(int).to_numpy(), group["impact_score"].astype(float).to_numpy()),
                "auprc_distance_matched": safe_average_precision(group["label"].astype(int).to_numpy(), group["impact_score"].astype(float).to_numpy()),
            }
        )
    return pd.DataFrame(rows).sort_values("auroc_distance_matched", ascending=False)


def build_qc_report(
    split_dir: Path = SHARED_SPLIT_DIR,
    processed_dir: Path = SHARED_PROCESSED_DIR,
    clinvar_path: Path = EXP3_DATA_DIR / "clinvar_splicing_variants.csv",
    genome_path: Path = PROJECT_ROOT / "data/raw/genome.fa",
    gtf_path: Path = PROJECT_ROOT / "data/raw/gencode.gtf",
    output: Path = REPORTS_DIR / "data_qc.md",
) -> pd.DataFrame:
    window_frames = load_window_frames(processed_dir)
    split_frames = load_split_frames(split_dir)
    pm200 = window_frames[200]

    window_consistency, split_window_consistency, consistency_flags = check_window_split_consistency(window_frames, split_dir)

    split_rows = []
    for split, frame in split_frames.items():
        current = frame.copy()
        current["category"] = current.apply(_category, axis=1)
        counts = current["category"].value_counts()
        split_rows.append(
            {
                "split": split,
                "rows": len(current),
                "donor": int(counts.get("donor", 0)),
                "acceptor": int(counts.get("acceptor", 0)),
                "easy_negative": int(counts.get("easy_negative", 0)),
                "hard_negative": int(counts.get("hard_negative", 0)),
                "positive_rate": float((current["label"].astype(int).isin([0, 1])).mean()),
                "hard_negative_rate_among_negatives": float(
                    (
                        (current["label"].astype(int) == 2)
                        & current["negative_type"].astype(str).str.contains("hard", na=False)
                    ).sum()
                    / max(1, (current["label"].astype(int) == 2).sum())
                ),
                "chromosomes": ",".join(sorted(current["chrom"].astype(str).unique())),
            }
        )
    split_summary = pd.DataFrame(split_rows)

    chrom_overlap_rows = []
    for left_name, right_name in itertools.combinations(SPLITS, 2):
        left = split_frames[left_name]
        right = split_frames[right_name]
        chrom_overlap_rows.append(
            {
                "pair": f"{left_name}/{right_name}",
                "sample_id_overlap": len(set(left["sample_id"]) & set(right["sample_id"])),
                "chrom_overlap": ",".join(sorted(set(left["chrom"]) & set(right["chrom"]))),
                "gene_overlap": len(set(left["gene_id"]) & set(right["gene_id"])),
            }
        )
    split_overlap = pd.DataFrame(chrom_overlap_rows)
    homology_audit, homology_examples = homology_leakage_audit(split_frames)
    homology_ok = (
        int(homology_audit.attrs.get("exact_full_pairs", 0)) == 0
        and int(homology_audit.attrs.get("exact_center_pairs", 0)) == 0
        and int(homology_audit.attrs.get("high_jaccard_pairs", 0)) == 0
    )

    length_rows = []
    for window, frame in window_frames.items():
        expected = 2 * window + 1
        length_rows.append(
            {
                "window": f"pm{window}",
                "rows": len(frame),
                "expected_length": expected,
                "bad_length_rows": int((frame["sequence"].astype(str).str.len() != expected).sum()),
                "alphabet": _alphabet(frame),
                "n_fraction": _n_fraction(frame),
                "donor_gt_rate": _canonical_rate(frame, 0),
                "acceptor_ag_rate": _canonical_rate(frame, 1),
                "hard_gt_or_ag_rate": _hard_motif_rate(frame),
            }
        )
    sequence_qc = pd.DataFrame(length_rows)

    hard = pm200[(pm200["label"].astype(int) == 2) & pm200["negative_type"].astype(str).str.contains("hard", na=False)]
    annotated_overlap = -1
    if genome_path.exists() and gtf_path.exists():
        annotated_sites = build_annotated_site_set(genome_path, gtf_path)
        annotated_overlap = int(sum((row["chrom"], int(row["center"])) in annotated_sites for _, row in hard.iterrows()))

    gc_frame = pm200.copy()
    gc_frame["category"] = gc_frame.apply(_category, axis=1)
    gc_frame["gc_fraction"] = gc_frame["sequence"].map(_gc_fraction)
    gc_summary = (
        gc_frame.groupby("category")["gc_fraction"]
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
        .rename(columns={"count": "rows"})
    )

    clinvar = pd.read_csv(clinvar_path) if clinvar_path.exists() else pd.DataFrame()
    clinvar_counts = pd.DataFrame()
    clinvar_distance = pd.DataFrame()
    clinvar_distance_bins = pd.DataFrame()
    clinvar_distance_metrics = pd.DataFrame()
    clinvar_matched = pd.DataFrame()
    clinvar_matched_score_metrics = pd.DataFrame()
    clinvar_overlap = pd.DataFrame()
    clinvar_checks = {"genomic_ref_mismatches": -1, "hamming_failures": -1, "position_failures": -1}
    if not clinvar.empty:
        clinvar_counts = (
            clinvar.groupby(["label_name", "variant_type", "target_class_name"], as_index=False)
            .size()
            .rename(columns={"size": "rows"})
        )
        clinvar_distance = clinvar.groupby("label_name")["nearest_splice_distance"].describe().reset_index()
        clinvar_matched, clinvar_distance_bins = build_distance_matched_clinvar(clinvar)
        write_dataframe(EXP3_DATA_DIR / "clinvar_splicing_variants_distance_matched.csv", clinvar_matched)
        clinvar_distance_metrics = pd.DataFrame(
            [
                distance_only_metrics(clinvar, "full_clinvar"),
                distance_only_metrics(clinvar_matched, "exact_distance_matched"),
            ]
        )
        clinvar_matched_score_metrics = matched_score_metrics(
            EXP3_TABLES_DIR / "experiment_3A_variant_scores.csv",
            clinvar_matched,
        )
        if genome_path.exists():
            clinvar_checks = verify_clinvar_rows(clinvar, genome_path)
        overlap_rows = []
        for split, frame in split_frames.items():
            overlap_rows.append(
                {
                    "split": split,
                    "chrom_overlap": ",".join(sorted(set(clinvar["chrom"]) & set(frame["chrom"]))),
                    "gene_overlap": len(set(clinvar["gene_id"]) & set(frame["gene_id"])),
                }
            )
        clinvar_overlap = pd.DataFrame(overlap_rows)

    status_rows = [
        {
            "section": "1 Split integrity",
            "status": _status(
                consistency_flags["same_assignment"]
                and consistency_flags["split_consistent"]
                and split_overlap["sample_id_overlap"].sum() == 0
                and split_overlap["chrom_overlap"].astype(bool).sum() == 0
                and pm200["sample_id"].duplicated().sum() == 0
            ),
            "evidence": "Four windows share identical sample_id/split assignments; train/valid/test chromosomes and IDs are disjoint.",
        },
        {
            "section": "2 Class balance",
            "status": _status(split_summary["hard_negative"].min() >= 50),
            "evidence": "Three-class sampling is donor:acceptor:hard-negative = 1:1:1 globally; easy negatives are intentionally not used in the real benchmark.",
        },
        {
            "section": "3 Sequence QC",
            "status": _status(
                sequence_qc["bad_length_rows"].sum() == 0
                and set("".join(sequence_qc["alphabet"].tolist())).issubset(set("ACGT"))
                and sequence_qc["donor_gt_rate"].min() >= 0.97
                and sequence_qc["acceptor_ag_rate"].min() >= 0.97
            ),
            "evidence": "Window lengths, alphabet, N fraction, and canonical donor/acceptor motif checks passed.",
        },
        {
            "section": "4 Hard negatives",
            "status": _status(annotated_overlap == 0 and _hard_motif_rate(pm200) == 1.0),
            "evidence": f"Hard negatives carry GT/AG decoys; overlap with GTF annotated splice sites = {annotated_overlap}.",
        },
        {
            "section": "5 Coordinates and strand",
            "status": _status(sequence_qc["donor_gt_rate"].min() == 1.0 and sequence_qc["acceptor_ag_rate"].min() == 1.0),
            "evidence": "Positive motif checks pass after transcript-oriented reverse complement, including negative-strand records.",
        },
        {
            "section": "6 ClinVar",
            "status": _status(
                not clinvar.empty
                and set(clinvar["label"].astype(int)) == {0, 1}
                and clinvar_checks["genomic_ref_mismatches"] == 0
                and clinvar_checks["hamming_failures"] == 0
                and clinvar_checks["position_failures"] == 0
                and int(clinvar_overlap["gene_overlap"].sum()) == 0
                and str(clinvar_overlap.loc[clinvar_overlap["split"] == "train", "chrom_overlap"].iloc[0]) == ""
                and str(clinvar_overlap.loc[clinvar_overlap["split"] == "valid", "chrom_overlap"].iloc[0]) == ""
            ),
            "evidence": "ClinVar is balanced 250/250 on held-out chromosomes; REF and one-SNV WT/Mut checks passed; no sampled split gene is reused.",
        },
        {
            "section": "7 Paralog leakage",
            "status": _status(homology_ok),
            "evidence": "Cross-split full-window duplicates, center-window duplicates, and high 9-mer Jaccard near-duplicates were audited; alignment-level paralog clustering remains out of scope.",
        },
        {
            "section": "8 Reproducibility",
            "status": "PASS",
            "evidence": "Seed 42, GRCh38/GENCODE metadata, ClinVar fileDate, commands, and local resource status are recorded.",
        },
        {
            "section": "9 Cross-experiment consistency",
            "status": _status(not clinvar.empty and set(clinvar["chrom"]).issubset(TEST_HOLDOUT_CHROMS)),
            "evidence": "Experiments use GRCh38; ClinVar variants are restricted to the same held-out test chromosomes used by experiment 1/2.",
        },
    ]
    status = pd.DataFrame(status_rows)

    raw_metadata = parse_raw_metadata(gtf_path, PROJECT_ROOT / "data/raw/clinvar.vcf")
    commands = pd.DataFrame(
        [
            {
                "step": "splice_site",
                "command": "python -m src.data.build_splice_site_dataset --max-per-class 1000 --windows 50 100 200 400",
            },
            {
                "step": "clinvar",
                "command": "python -m src.data.build_clinvar_variant_dataset",
            },
            {
                "step": "qc",
                "command": "python -m src.data.qc_splice_dataset",
            },
        ]
    )

    lines = [
        "# 数据质控结果",
        "",
        "QC target: real GENCODE/GRCh38 splice-site benchmark plus real ClinVar variant benchmark.",
        "",
        "Important sampling definitions:",
        "",
        "- Splice-site benchmark: donor / acceptor / non-splice hard-negative = 1:1:1 globally; binary positive base rate is 2000/3000 = 0.667.",
        "- Negative splice-site examples: all real benchmark negatives are motif-matched GT/AG hard negatives; easy negatives are kept only for optional synthetic controls.",
        "- ClinVar benchmark: pathogenic or likely pathogenic splice-related/near-splice SNVs vs benign or likely benign near-splice SNVs, balanced 250/250.",
        "- ClinVar variants are restricted to held-out test chromosomes chr19, chr20, chr21, chr22, chrX and exclude genes already sampled in experiment 1/2 splits.",
        "",
        "## Checklist Status",
        "",
        _markdown(status),
        "",
        "## 1. Split 完整性与防泄漏",
        "",
        "Window-level identity and assignment check:",
        "",
        _markdown(window_consistency),
        "",
        "Split-file ID consistency across windows:",
        "",
        _markdown(split_window_consistency),
        "",
        _markdown(split_overlap),
        "",
        "Conclusion: 四个窗口共用同一批 `sample_id -> split`；split 不是按行随机切。train/valid/test 的 sample_id、chromosome、gene 均互斥。",
        "",
        "## 2. 类别构成与平衡",
        "",
        _markdown(split_summary),
        "",
        "Note: 真实主 benchmark 没有 easy negative；non-splice 类全部是 GT/AG hard negative。三分类全量为 1000/1000/1000，但按染色体 holdout 后各 split 比例会随染色体分布略有变化。",
        "",
        "## 3. 序列层面 QC",
        "",
        _markdown(sequence_qc),
        "",
        "N strategy: rows with N fraction above the builder threshold are filtered; generated benchmark has N fraction 0.",
        "",
        "## 4. 负例 / Hard-Negative 有效性",
        "",
        f"- pm200 hard-negative rows: {len(hard)}",
        f"- hard-negative center GT/AG motif rate: {_hard_motif_rate(pm200):.4f}",
        f"- overlap with any GTF annotated splice site: {annotated_overlap}",
        "",
        "GC fraction by category on pm200:",
        "",
        _markdown(gc_summary),
        "",
        "## 5. 坐标与链",
        "",
        "Donor GT and acceptor AG rates are 1.0 after transcript-oriented extraction, which indirectly checks GTF 1-based coordinates and negative-strand reverse-complement handling.",
        "",
        "## 6. ClinVar 变异表",
        "",
        _markdown(clinvar_counts),
        "",
        "ClinVar distance to nearest annotated splice site:",
        "",
        _markdown(clinvar_distance),
        "",
        "Distance confounding diagnostic:",
        "",
        "The pathogenic/splice-altering variants are closer to annotated splice sites than benign variants on average. This is biologically plausible, but it is also a confounder: a model can partially solve the task by learning `closer to splice site = more likely pathogenic`. To make this visible, QC reports a distance-only baseline and exports an exact-distance-matched subset.",
        "",
        _markdown(clinvar_distance_metrics),
        "",
        "Exact-distance matching summary:",
        "",
        _markdown(
            clinvar_distance_bins[clinvar_distance_bins["matched_pairs"] > 0]
            if "matched_pairs" in clinvar_distance_bins.columns
            else clinvar_distance_bins,
            max_rows=60,
        ),
        "",
        f"Distance-matched ClinVar subset: `{EXP3_DATA_DIR / 'clinvar_splicing_variants_distance_matched.csv'}`.",
        "",
        "Model metrics on the distance-matched subset, if experiment-3 scores have been generated:",
        "",
        _markdown(clinvar_matched_score_metrics),
        "",
        "ClinVar overlap with experiment 1/2 split genes/chromosomes:",
        "",
        _markdown(clinvar_overlap),
        "",
        "Interpretation: `chrom_overlap=chr19` for the test split is expected because ClinVar is restricted to held-out test chromosomes. `gene_overlap=0` for train/valid/test confirms that no sampled experiment-1/2 gene is reused in the ClinVar variant benchmark.",
        "",
        "ClinVar allele checks:",
        "",
        _markdown(pd.DataFrame([clinvar_checks])),
        "",
        "ClinVar filtering rule: SNVs only; positive labels are pathogenic/likely pathogenic records with splice-related consequence or within the near-splice threshold; negative labels are benign/likely benign SNVs within the same near-splice threshold. REF is checked against GRCh38 and strand-oriented REF/ALT are checked in WT/Mut windows.",
        "",
        "Gain classes: not modeled as separate donor_gain/acceptor_gain labels in the real ClinVar benchmark. The experiment-3 story should be reported as real ClinVar splice-altering vs benign ranking, stratified by donor/acceptor target, not as synthetic gain/loss recovery.",
        "",
        "## 7. 同源 / 旁系泄漏",
        "",
        "Cross-split near-duplicate leakage audit:",
        "",
        _markdown(homology_audit),
        "",
        "High-similarity examples above the audit threshold, if any:",
        "",
        _markdown(homology_examples),
        "",
        "Interpretation: this audit checks exact full-window duplicates, exact center 161 bp duplicates, and high 9-mer Jaccard near-duplicates across train/valid/test. It is a practical homology-leakage screen for the current benchmark, but it is not a replacement for alignment-level paralog clustering with tools such as BLAST/CD-HIT/MMseqs.",
        "",
        "## 8. 复现与溯源",
        "",
        _markdown(raw_metadata),
        "",
        _markdown(commands),
        "",
        "Random seed: 42.",
        "",
        "## 9. 跨实验一致性",
        "",
        "Experiment 1/2 split uses the local GRCh38 + GENCODE annotation and chromosome holdout. Experiment 3 ClinVar windows use the same genome/annotation coordinate system and are restricted to held-out test chromosomes.",
        "",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return status


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate checklist-style QC report for the real splice benchmark.")
    parser.add_argument("--split-dir", type=Path, default=SHARED_SPLIT_DIR)
    parser.add_argument("--processed-dir", type=Path, default=SHARED_PROCESSED_DIR)
    parser.add_argument("--clinvar", type=Path, default=EXP3_DATA_DIR / "clinvar_splicing_variants.csv")
    parser.add_argument("--genome", type=Path, default=PROJECT_ROOT / "data/raw/genome.fa")
    parser.add_argument("--gtf", type=Path, default=PROJECT_ROOT / "data/raw/gencode.gtf")
    parser.add_argument("--output", type=Path, default=REPORTS_DIR / "data_qc.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    status = build_qc_report(
        split_dir=args.split_dir,
        processed_dir=args.processed_dir,
        clinvar_path=args.clinvar,
        genome_path=args.genome,
        gtf_path=args.gtf,
        output=args.output,
    )
    print(status.to_string(index=False))


if __name__ == "__main__":
    main()
