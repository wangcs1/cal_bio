#!/usr/bin/env python3
"""
Build balanced donor / acceptor / motif-matched non-splice datasets.

Input:
  data/raw/genome.fa
  data/raw/gencode.gtf

Output:
  data/shared/processed/splice_sites_pm{window}.csv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from src.utils import REQUIRED_SPLIT_COLUMNS, SHARED_PROCESSED_DIR


LABELS = {
    0: "donor",
    1: "acceptor",
    2: "non_splice",
}

RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STANDARD_EXPERIMENT_CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX"])
PROJECT_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Site:
    chrom: str
    center: int
    strand: str
    label: int
    gene_id: str
    transcript_id: str
    gene_type: str = ""
    transcript_type: str = ""


@dataclass(frozen=True)
class Candidate:
    chrom: str
    center: int
    strand: str
    gene_id: str
    transcript_id: str


class FastaIndex:
    def __init__(self, fasta_path: Path) -> None:
        self.fasta_path = fasta_path
        self.index_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
        if not self.index_path.exists():
            self._build_index()
        self.records = self._read_index()
        self.handle = fasta_path.open("rb")

    def close(self) -> None:
        self.handle.close()

    def chrom_length(self, chrom: str) -> int:
        return self.records[chrom][0]

    def has_chrom(self, chrom: str) -> bool:
        return chrom in self.records

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if start < 1 or end < start:
            raise ValueError(f"Invalid FASTA interval: {chrom}:{start}-{end}")
        length, offset, line_bases, line_width = self.records[chrom]
        if end > length:
            raise ValueError(f"Interval exceeds chromosome length: {chrom}:{start}-{end}")

        pos0 = start - 1
        remaining = end - start + 1
        chunks: list[bytes] = []
        while remaining:
            line_idx = pos0 // line_bases
            line_pos = pos0 % line_bases
            to_read = min(remaining, line_bases - line_pos)
            byte_pos = offset + line_idx * line_width + line_pos
            self.handle.seek(byte_pos)
            chunks.append(self.handle.read(to_read))
            pos0 += to_read
            remaining -= to_read
        return b"".join(chunks).decode("ascii").upper()

    def _build_index(self) -> None:
        print(f"Building FASTA index: {self.index_path}")
        rows: list[tuple[str, int, int, int, int]] = []
        name: str | None = None
        length = 0
        offset = 0
        line_bases = 0
        line_width = 0

        with self.fasta_path.open("rb") as handle:
            while True:
                line_start = handle.tell()
                line = handle.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    if name is not None:
                        rows.append((name, length, offset, line_bases, line_width))
                    name = line[1:].decode("ascii").strip().split()[0]
                    length = 0
                    offset = handle.tell()
                    line_bases = 0
                    line_width = 0
                    continue
                stripped = line.rstrip(b"\r\n")
                if not stripped:
                    continue
                if line_bases == 0:
                    line_bases = len(stripped)
                    line_width = len(line)
                    offset = line_start
                length += len(stripped)

        if name is not None:
            rows.append((name, length, offset, line_bases, line_width))

        with self.index_path.open("w", newline="") as out:
            for row in rows:
                out.write("\t".join(map(str, row)) + "\n")

    def _read_index(self) -> dict[str, tuple[int, int, int, int]]:
        records: dict[str, tuple[int, int, int, int]] = {}
        with self.index_path.open() as handle:
            for line in handle:
                chrom, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")
                records[chrom] = (int(length), int(offset), int(line_bases), int(line_width))
        return records


def parse_attributes(attribute_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attribute_text.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, value = item.partition(" ")
        attrs[key] = value.strip().strip('"')
    return attrs


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1].upper()


def natural_chrom_key(chrom: str) -> tuple[int, str]:
    token = chrom.removeprefix("chr")
    if token.isdigit():
        return (int(token), chrom)
    if token == "X":
        return (23, chrom)
    if token == "Y":
        return (24, chrom)
    if token in {"M", "MT"}:
        return (25, chrom)
    return (10_000, chrom)


def parse_chroms(chrom_text: str) -> set[str]:
    chroms: set[str] = set()
    for token in chrom_text.split(","):
        token = token.strip()
        if not token:
            continue
        chroms.add(token if token.startswith("chr") else f"chr{token}")
    return chroms


def make_sample_id(chrom: str, center: int, strand: str, gene_id: str, transcript_id: str) -> str:
    key = f"{chrom}|{center}|{strand}|{gene_id}|{transcript_id}"
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]
    return f"ss_{digest}"


def read_transcripts(
    gtf_path: Path,
    fasta: FastaIndex,
    protein_coding_only: bool,
    gtf_line_limit: int | None,
) -> dict[str, dict[str, object]]:
    transcripts: dict[str, dict[str, object]] = {}
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if gtf_line_limit is not None and line_number > gtf_line_limit:
                break
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue

            chrom, _, _, start_text, end_text, _, strand, _, attrs_text = fields
            if strand not in {"+", "-"} or not fasta.has_chrom(chrom):
                continue

            attrs = parse_attributes(attrs_text)
            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id")
            if not transcript_id or not gene_id:
                continue

            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))
            transcript_type = attrs.get("transcript_type", attrs.get("transcript_biotype", ""))
            if protein_coding_only and "protein_coding" not in {gene_type, transcript_type}:
                continue

            transcript = transcripts.setdefault(
                transcript_id,
                {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_type": gene_type,
                    "transcript_type": transcript_type,
                    "exons": [],
                },
            )
            transcript["exons"].append((int(start_text), int(end_text)))  # type: ignore[index, union-attr]

            if line_number % 1_000_000 == 0:
                print(f"Parsed {line_number:,} GTF lines; transcripts={len(transcripts):,}")

    return transcripts


def build_positive_sites(
    transcripts: dict[str, dict[str, object]],
) -> tuple[list[Site], dict[str, tuple[int, int]], set[tuple[str, int]], set[tuple[str, str, int]]]:
    raw_sites: list[Site] = []
    transcript_bounds: dict[str, tuple[int, int]] = {}

    for transcript_id, transcript in transcripts.items():
        exons = sorted(transcript["exons"])  # type: ignore[arg-type]
        if len(exons) < 2:
            continue

        chrom = transcript["chrom"]  # type: ignore[assignment]
        strand = transcript["strand"]  # type: ignore[assignment]
        gene_id = transcript["gene_id"]  # type: ignore[assignment]
        gene_type = transcript["gene_type"]  # type: ignore[assignment]
        transcript_type = transcript["transcript_type"]  # type: ignore[assignment]
        transcript_bounds[transcript_id] = (exons[0][0], exons[-1][1])

        for left_exon, right_exon in zip(exons, exons[1:]):
            if left_exon[1] >= right_exon[0]:
                continue
            if strand == "+":
                donor_center = left_exon[1]
                acceptor_center = right_exon[0]
            else:
                donor_center = right_exon[0]
                acceptor_center = left_exon[1]

            raw_sites.append(
                Site(chrom, donor_center, strand, 0, gene_id, transcript_id, gene_type, transcript_type)
            )
            raw_sites.append(
                Site(chrom, acceptor_center, strand, 1, gene_id, transcript_id, gene_type, transcript_type)
            )

    labels_by_site: dict[tuple[str, str, int], set[int]] = {}
    first_by_key: dict[tuple[str, str, int, int], Site] = {}
    for site in raw_sites:
        labels_by_site.setdefault((site.chrom, site.strand, site.center), set()).add(site.label)
        first_by_key.setdefault((site.chrom, site.strand, site.center, site.label), site)

    sites: list[Site] = []
    for key, site in first_by_key.items():
        chrom, strand, center, _ = key
        if len(labels_by_site[(chrom, strand, center)]) == 1:
            sites.append(site)

    annotated_any_strand = {(site.chrom, site.center) for site in sites}
    annotated_same_strand = {(site.chrom, site.strand, site.center) for site in sites}
    sites.sort(key=lambda item: (natural_chrom_key(item.chrom), item.center, item.strand, item.label))
    return sites, transcript_bounds, annotated_any_strand, annotated_same_strand


def oriented_window(
    fasta: FastaIndex,
    chrom: str,
    center: int,
    strand: str,
    flank: int,
) -> str | None:
    start = center - flank
    end = center + flank
    if start < 1 or end > fasta.chrom_length(chrom):
        return None
    sequence = fasta.fetch(chrom, start, end)
    if strand == "-":
        sequence = reverse_complement(sequence)
    return sequence


def extract_sequence(
    fasta: FastaIndex,
    chrom: str,
    center: int,
    strand: str,
    window: int,
    max_n_ratio: float,
    rna: bool,
) -> tuple[int, int, str] | None:
    start = center - window
    end = center + window
    sequence = oriented_window(fasta, chrom, center, strand, window)
    if sequence is None:
        return None
    if len(sequence) != 2 * window + 1:
        return None
    if sequence.count("N") / len(sequence) > max_n_ratio:
        return None
    if rna:
        sequence = sequence.replace("T", "U")
    return start, end, sequence


def candidate_has_motif(
    fasta: FastaIndex,
    chrom: str,
    center: int,
    strand: str,
) -> bool:
    sequence = oriented_window(fasta, chrom, center, strand, 2)
    if sequence is None or len(sequence) != 5:
        return False
    return sequence[3:5] == "GT" or sequence[0:2] == "AG"


def passes_window_filters(
    fasta: FastaIndex,
    chrom: str,
    center: int,
    strand: str,
    windows: list[int],
    max_n_ratio: float,
) -> bool:
    max_window = max(windows)
    sequence = oriented_window(fasta, chrom, center, strand, max_window)
    if sequence is None or len(sequence) != 2 * max_window + 1:
        return False

    midpoint = max_window
    for window in windows:
        sub_sequence = sequence[midpoint - window : midpoint + window + 1]
        if sub_sequence.count("N") / len(sub_sequence) > max_n_ratio:
            return False
    return True


def positive_has_canonical_motif(fasta: FastaIndex, site: Site) -> bool:
    sequence = oriented_window(fasta, site.chrom, site.center, site.strand, 2)
    if sequence is None or len(sequence) != 5:
        return False
    if site.label == 0:
        return sequence[3:5] == "GT"
    if site.label == 1:
        return sequence[0:2] == "AG"
    return False


def filter_positive_sites(
    sites: list[Site],
    fasta: FastaIndex,
    windows: list[int],
    max_n_ratio: float,
    allowed_chroms: set[str] | None,
    canonical_only: bool,
) -> list[Site]:
    filtered: list[Site] = []
    for site in sites:
        if allowed_chroms is not None and site.chrom not in allowed_chroms:
            continue
        if canonical_only and not positive_has_canonical_motif(fasta, site):
            continue
        if not passes_window_filters(fasta, site.chrom, site.center, site.strand, windows, max_n_ratio):
            continue
        filtered.append(site)
    return filtered


def filter_negative_candidates(
    candidates: list[Candidate],
    fasta: FastaIndex,
    windows: list[int],
    max_n_ratio: float,
) -> list[Candidate]:
    return [
        candidate
        for candidate in candidates
        if passes_window_filters(fasta, candidate.chrom, candidate.center, candidate.strand, windows, max_n_ratio)
    ]


def sample_negative_candidates(
    positives: list[Site],
    transcript_bounds: dict[str, tuple[int, int]],
    fasta: FastaIndex,
    annotated_any_strand: set[tuple[str, int]],
    annotated_same_strand: set[tuple[str, str, int]],
    max_window: int,
    target_count: int,
    radius: int,
    attempts_per_site: int,
    rng: random.Random,
) -> list[Candidate]:
    print(f"Sampling motif-matched negatives; target={target_count:,}")
    candidates: list[Candidate] = []
    seen = set()
    positive_pool = positives[:]
    if not positive_pool:
        return candidates

    passes = max(3, math.ceil(target_count / len(positive_pool)) * 4)
    for pass_idx in range(passes):
        rng.shuffle(positive_pool)
        for site in positive_pool:
            if len(candidates) >= target_count:
                return candidates

            transcript_start, transcript_end = transcript_bounds.get(site.transcript_id, (1, fasta.chrom_length(site.chrom)))
            local_start = max(transcript_start + max_window, site.center - radius, 1 + max_window)
            local_end = min(transcript_end - max_window, site.center + radius, fasta.chrom_length(site.chrom) - max_window)
            if local_start > local_end:
                local_start = max(transcript_start + max_window, 1 + max_window)
                local_end = min(transcript_end - max_window, fasta.chrom_length(site.chrom) - max_window)
            if local_start > local_end:
                continue

            for _ in range(attempts_per_site):
                center = rng.randint(local_start, local_end)
                key = (site.chrom, site.strand, center)
                if key in seen:
                    continue
                if (site.chrom, center) in annotated_any_strand or key in annotated_same_strand:
                    continue
                if not candidate_has_motif(fasta, site.chrom, center, site.strand):
                    continue
                candidates.append(Candidate(site.chrom, center, site.strand, site.gene_id, site.transcript_id))
                seen.add(key)
                break

        print(f"Negative pass {pass_idx + 1}/{passes}: {len(candidates):,}/{target_count:,}")

    return candidates


def select_balanced_examples(
    sites: Iterable[Site],
    negatives: Iterable[Candidate],
    max_per_class: int | None,
    rng: random.Random,
) -> tuple[list[Site], list[Candidate], int]:
    sites_by_label: dict[int, list[Site]] = {0: [], 1: []}
    for site in sites:
        sites_by_label[site.label].append(site)

    negative_list = list(negatives)
    class_count = min(len(sites_by_label[0]), len(sites_by_label[1]), len(negative_list))
    if max_per_class is not None:
        class_count = min(class_count, max_per_class)

    for label in (0, 1):
        rng.shuffle(sites_by_label[label])
    rng.shuffle(negative_list)

    selected_sites = sites_by_label[0][:class_count] + sites_by_label[1][:class_count]
    selected_negatives = negative_list[:class_count]
    return selected_sites, selected_negatives, class_count


def build_rows_for_window(
    sites: Iterable[Site],
    negatives: Iterable[Candidate],
    fasta: FastaIndex,
    window: int,
    rna: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for site in sites:
        extracted = extract_sequence(fasta, site.chrom, site.center, site.strand, window, 1.0, rna)
        if extracted is None:
            raise ValueError(f"Selected positive site failed extraction: {site}")
        start, end, sequence = extracted
        rows.append(
            {
                "sample_id": make_sample_id(site.chrom, site.center, site.strand, site.gene_id, site.transcript_id),
                "chrom": site.chrom,
                "start": start,
                "end": end,
                "strand": site.strand,
                "center": site.center,
                "label": site.label,
                "label_name": LABELS[site.label],
                "negative_type": "positive",
                "sequence": sequence,
                "gene_id": site.gene_id,
                "transcript_id": site.transcript_id,
            }
        )

    for candidate in negatives:
        extracted = extract_sequence(fasta, candidate.chrom, candidate.center, candidate.strand, window, 1.0, rna)
        if extracted is None:
            raise ValueError(f"Selected negative site failed extraction: {candidate}")
        start, end, sequence = extracted
        rows.append(
            {
                "sample_id": make_sample_id(
                    candidate.chrom,
                    candidate.center,
                    candidate.strand,
                    candidate.gene_id,
                    candidate.transcript_id,
                ),
                "chrom": candidate.chrom,
                "start": start,
                "end": end,
                "strand": candidate.strand,
                "center": candidate.center,
                "label": 2,
                "label_name": LABELS[2],
                "negative_type": "hard_gtag",
                "sequence": sequence,
                "gene_id": candidate.gene_id,
                "transcript_id": candidate.transcript_id,
            }
        )

    rows.sort(key=lambda row: str(row["sample_id"]))
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], full_columns: bool) -> None:
    if full_columns:
        fieldnames = [
            "sample_id",
            "chrom",
            "start",
            "end",
            "strand",
            "center",
            "label",
            "label_name",
            "negative_type",
            "sequence",
            "gene_id",
            "transcript_id",
        ]
    else:
        fieldnames = [
            "sample_id",
            "chrom",
            "center",
            "strand",
            "label",
            "sequence",
            "gene_id",
        ]
    missing = sorted(REQUIRED_SPLIT_COLUMNS - set(fieldnames))
    if missing:
        raise ValueError(
            "Output field list is missing required split columns. "
            "Use --full-columns or update fieldnames. Missing: " + ", ".join(missing)
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", type=Path, default=PROJECT_ROOT / "data/raw/genome.fa")
    parser.add_argument("--gtf", type=Path, default=PROJECT_ROOT / "data/raw/gencode.gtf")
    parser.add_argument("--out-dir", type=Path, default=SHARED_PROCESSED_DIR)
    parser.add_argument("--windows", type=int, nargs="+", default=[50, 100, 200, 400])
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-n-ratio", type=float, default=0.05)
    parser.add_argument("--negative-radius", type=int, default=5000)
    parser.add_argument("--negative-pool-ratio", type=float, default=1.5)
    parser.add_argument("--negative-attempts", type=int, default=80)
    parser.add_argument("--max-per-class", type=int, default=None)
    parser.add_argument("--protein-coding-only", action="store_true")
    parser.add_argument("--dna", action="store_true", help="Keep T instead of converting DNA to RNA U.")
    parser.add_argument(
        "--chroms",
        default=",".join(STANDARD_EXPERIMENT_CHROMS),
        help="Comma-separated chromosomes to keep. Defaults to chr1-chr22,chrX.",
    )
    parser.add_argument("--all-chroms", action="store_true", help="Keep all chromosomes and contigs.")
    parser.add_argument(
        "--keep-noncanonical-positives",
        action="store_true",
        help="Keep positive donor/acceptor sites without canonical GT/AG motifs.",
    )
    parser.add_argument(
        "--full-columns",
        action="store_true",
        default=True,
        help="Write traceability columns: start,end,label_name,transcript_id.",
    )
    parser.add_argument(
        "--minimal-columns",
        dest="full_columns",
        action="store_false",
        help="Write the legacy minimal column set; not recommended for C Part experiments.",
    )
    parser.add_argument("--gtf-line-limit", type=int, default=None, help=argparse.SUPPRESS)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)
    max_window = max(args.windows)
    allowed_chroms = None if args.all_chroms else parse_chroms(args.chroms)

    fasta = FastaIndex(args.genome)
    try:
        transcripts = read_transcripts(args.gtf, fasta, args.protein_coding_only, args.gtf_line_limit)
        print(f"Loaded transcripts with exons: {len(transcripts):,}")

        sites, transcript_bounds, annotated_any, annotated_same = build_positive_sites(transcripts)
        donor_count = sum(1 for site in sites if site.label == 0)
        acceptor_count = sum(1 for site in sites if site.label == 1)
        print(f"Annotated sites after de-duplication: donor={donor_count:,}, acceptor={acceptor_count:,}")

        sites = filter_positive_sites(
            sites,
            fasta,
            args.windows,
            args.max_n_ratio,
            allowed_chroms,
            not args.keep_noncanonical_positives,
        )
        donor_count = sum(1 for site in sites if site.label == 0)
        acceptor_count = sum(1 for site in sites if site.label == 1)
        print(
            "Positive sites after experiment filtering: "
            f"donor={donor_count:,}, acceptor={acceptor_count:,}"
        )

        base_target = min(donor_count, acceptor_count)
        if args.max_per_class is not None:
            base_target = min(base_target, args.max_per_class)
        negative_target = max(base_target, int(base_target * args.negative_pool_ratio))

        negatives = sample_negative_candidates(
            sites,
            transcript_bounds,
            fasta,
            annotated_any,
            annotated_same,
            max_window,
            negative_target,
            args.negative_radius,
            args.negative_attempts,
            rng,
        )
        print(f"Negative candidates before shared window filtering: {len(negatives):,}")
        negatives = filter_negative_candidates(negatives, fasta, args.windows, args.max_n_ratio)
        print(f"Negative candidates after shared window filtering: {len(negatives):,}")

        selected_sites, selected_negatives, class_count = select_balanced_examples(
            sites,
            negatives,
            args.max_per_class,
            rng,
        )
        print(f"Selected shared centers per class: {class_count:,}")

        for window in args.windows:
            rows = build_rows_for_window(
                selected_sites,
                selected_negatives,
                fasta,
                window,
                not args.dna,
            )
            out_path = args.out_dir / f"splice_sites_pm{window}.csv"
            write_csv(out_path, rows, args.full_columns)
            per_label = {label: sum(1 for row in rows if row["label"] == label) for label in LABELS}
            print(f"Wrote {out_path}: {len(rows):,} rows; labels={per_label}")
    finally:
        fasta.close()


if __name__ == "__main__":
    main()
