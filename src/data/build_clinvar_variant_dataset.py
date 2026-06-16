from __future__ import annotations

import argparse
import bisect
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from src.data.build_splice_site_dataset import (
    FastaIndex,
    Site,
    build_positive_sites,
    extract_sequence,
    parse_chroms,
    read_transcripts,
    reverse_complement,
    validate_raw_inputs,
)
from src.utils import EXP3_DATA_DIR, PROJECT_ROOT, stable_id, write_dataframe
from src.utils import SHARED_SPLIT_DIR


DEFAULT_GENOME = PROJECT_ROOT / "data/raw/genome.fa"
DEFAULT_GTF = PROJECT_ROOT / "data/raw/gencode.gtf"
DEFAULT_CLINVAR = PROJECT_ROOT / "data/raw/clinvar.vcf"
DEFAULT_OUTPUT = EXP3_DATA_DIR / "clinvar_splicing_variants.csv"
CLINVAR_DATA_SOURCE = "clinvar_real_splice_variants_v1"
DEFAULT_CLINVAR_CHROMS = "chr19,chr20,chr21,chr22,chrX"


@dataclass(frozen=True)
class ClinvarRecord:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alt: str
    clinvar_significance: str
    consequence: str
    info: dict[str, str]


def parse_info(info_text: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in info_text.split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        info[key] = value if sep else "true"
    return info


def normalize_label_text(value: str) -> str:
    return value.lower().replace("%20", "_").replace(" ", "_").replace("-", "_")


def is_conflicting(label_text: str) -> bool:
    text = normalize_label_text(label_text)
    return "conflict" in text or "uncertain" in text or "not_provided" in text


def is_pathogenic_label(label_text: str) -> bool:
    text = normalize_label_text(label_text)
    return "pathogenic" in text and "benign" not in text and not is_conflicting(text)


def is_benign_label(label_text: str) -> bool:
    text = normalize_label_text(label_text)
    return "benign" in text and "pathogenic" not in text and not is_conflicting(text)


def is_splice_related(record: ClinvarRecord) -> bool:
    joined = ";".join([record.consequence, *[f"{key}={value}" for key, value in record.info.items()]])
    return "splice" in normalize_label_text(joined)


def normalize_chrom(chrom: str, fasta: FastaIndex) -> str | None:
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom.removeprefix("chr"))
    else:
        candidates.append(f"chr{chrom}")
    for candidate in candidates:
        if fasta.has_chrom(candidate):
            return candidate
    return None


def iter_clinvar_snv_records(vcf_path: Path, fasta: FastaIndex):
    if not vcf_path.exists():
        raise FileNotFoundError(
            f"Real ClinVar dataset requires local VCF: {vcf_path}. "
            "Place ClinVar VCF at data/raw/clinvar.vcf."
        )

    with vcf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            raw_chrom, pos_text, raw_id, ref, alt, *_rest, info_text = fields[:8]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue
            chrom = normalize_chrom(raw_chrom, fasta)
            if chrom is None:
                continue
            info = parse_info(info_text)
            yield ClinvarRecord(
                chrom=chrom,
                pos=int(pos_text),
                variant_id=raw_id if raw_id and raw_id != "." else stable_id(chrom, pos_text, ref, alt, prefix="clinvar"),
                ref=ref.upper(),
                alt=alt.upper(),
                clinvar_significance=info.get("CLNSIG", ""),
                consequence=info.get("MC", info.get("MolecularConsequence", "")),
                info=info,
            )


def build_site_index(sites: list[Site]) -> dict[str, tuple[list[int], list[Site]]]:
    by_chrom: dict[str, list[Site]] = {}
    for site in sites:
        by_chrom.setdefault(site.chrom, []).append(site)

    index: dict[str, tuple[list[int], list[Site]]] = {}
    for chrom, chrom_sites in by_chrom.items():
        chrom_sites.sort(key=lambda site: site.center)
        index[chrom] = ([site.center for site in chrom_sites], chrom_sites)
    return index


def nearest_site(record: ClinvarRecord, site_index: dict[str, tuple[list[int], list[Site]]]) -> tuple[Site, int] | None:
    if record.chrom not in site_index:
        return None
    centers, sites = site_index[record.chrom]
    insertion = bisect.bisect_left(centers, record.pos)
    candidates = []
    for idx in (insertion - 1, insertion):
        if 0 <= idx < len(sites):
            site = sites[idx]
            candidates.append((abs(record.pos - site.center), site))
    if not candidates:
        return None
    distance, site = min(candidates, key=lambda item: item[0])
    return site, distance


def oriented_variant_index(site: Site, variant_pos: int, window: int) -> int:
    rel = variant_pos - site.center
    return window + rel if site.strand == "+" else window - rel


def oriented_allele(base: str, strand: str) -> str:
    return reverse_complement(base) if strand == "-" else base.upper()


def make_variant_row(
    record: ClinvarRecord,
    site: Site,
    distance: int,
    fasta: FastaIndex,
    window: int,
    label: int,
) -> dict[str, object] | None:
    genomic_ref = fasta.fetch(record.chrom, record.pos, record.pos)
    if genomic_ref != record.ref:
        return None

    extracted = extract_sequence(fasta, site.chrom, site.center, site.strand, window, 0.05, rna=False)
    if extracted is None:
        return None
    _start, _end, wt_sequence = extracted
    variant_index = oriented_variant_index(site, record.pos, window)
    if variant_index < 0 or variant_index >= len(wt_sequence):
        return None

    ref = oriented_allele(record.ref, site.strand)
    alt = oriented_allele(record.alt, site.strand)
    if wt_sequence[variant_index] != ref:
        return None
    chars = list(wt_sequence)
    chars[variant_index] = alt
    mut_sequence = "".join(chars)

    target_name = "donor" if site.label == 0 else "acceptor"
    variant_type = f"{target_name}_clinvar_splice" if label else "clinvar_benign_snv"
    return {
        "variant_id": stable_id(record.variant_id, site.chrom, site.center, site.strand, prefix="clinvar"),
        "clinvar_id": record.variant_id,
        "sample_id": stable_id(site.chrom, site.center, site.strand, site.gene_id, site.transcript_id, prefix="ss"),
        "chrom": site.chrom,
        "pos": record.pos,
        "relative_pos": record.pos - site.center,
        "strand": site.strand,
        "ref": ref,
        "alt": alt,
        "genomic_ref": record.ref,
        "genomic_alt": record.alt,
        "variant_type": variant_type,
        "label": label,
        "label_name": "splice_altering" if label else "neutral",
        "wt_sequence": wt_sequence,
        "mut_sequence": mut_sequence,
        "target_class": site.label,
        "target_class_name": target_name,
        "gene_id": site.gene_id,
        "transcript_id": site.transcript_id,
        "nearest_splice_distance": distance,
        "clinvar_significance": record.clinvar_significance,
        "clinvar_consequence": record.consequence,
        "clinvar_splice_related": is_splice_related(record),
        "data_source": CLINVAR_DATA_SOURCE,
    }


def balance_rows(rows: list[dict[str, object]], max_per_label: int | None, seed: int) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    positives = frame[frame["label"].astype(int) == 1]
    negatives = frame[frame["label"].astype(int) == 0]
    target = min(len(positives), len(negatives))
    if max_per_label is not None:
        target = min(target, max_per_label)
    if target == 0:
        raise ValueError("ClinVar builder found variants, but not both positive and benign labels after filtering.")
    sampled = pd.concat(
        [
            positives.sample(n=target, random_state=seed),
            negatives.sample(n=target, random_state=seed + 1),
        ],
        ignore_index=True,
    )
    return sampled.sample(frac=1.0, random_state=seed).reset_index(drop=True)


def load_excluded_split_genes(split_dir: Path | None) -> set[str]:
    if split_dir is None or not split_dir.exists():
        return set()
    genes: set[str] = set()
    for split in ("train", "valid", "test"):
        path = split_dir / f"{split}.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path, usecols=["gene_id"])
        genes.update(frame["gene_id"].astype(str))
    return genes


def build_clinvar_variants(
    input_path: Path = DEFAULT_CLINVAR,
    output_path: Path = DEFAULT_OUTPUT,
    genome_path: Path = DEFAULT_GENOME,
    gtf_path: Path = DEFAULT_GTF,
    window: int = 200,
    max_distance: int = 50,
    max_per_label: int | None = 250,
    seed: int = 42,
    stop_when_balanced: bool = True,
    chroms: str = DEFAULT_CLINVAR_CHROMS,
    all_chroms: bool = False,
    protein_coding_only: bool = False,
    exclude_split_genes_dir: Path | None = SHARED_SPLIT_DIR,
) -> pd.DataFrame:
    validate_raw_inputs(genome_path, gtf_path)
    allowed_chroms = None if all_chroms else parse_chroms(chroms)
    excluded_genes = load_excluded_split_genes(exclude_split_genes_dir)
    if excluded_genes:
        print(f"Excluding ClinVar variants nearest to {len(excluded_genes):,} genes already used in splice-site splits.")
    fasta = FastaIndex(genome_path)
    try:
        transcripts = read_transcripts(gtf_path, fasta, protein_coding_only, gtf_line_limit=None)
        sites, _bounds, _annotated_any, _annotated_same = build_positive_sites(transcripts)
        if allowed_chroms is not None:
            sites = [site for site in sites if site.chrom in allowed_chroms]
        site_index = build_site_index(sites)

        rows_by_label: dict[int, list[dict[str, object]]] = {0: [], 1: []}
        scanned = 0
        for record in iter_clinvar_snv_records(input_path, fasta):
            scanned += 1
            nearest = nearest_site(record, site_index)
            if nearest is None:
                continue
            site, distance = nearest
            if distance > max_distance and not is_splice_related(record):
                continue
            if is_pathogenic_label(record.clinvar_significance) and (is_splice_related(record) or distance <= max_distance):
                label = 1
            elif is_benign_label(record.clinvar_significance) and distance <= max_distance:
                label = 0
            else:
                continue
            row = make_variant_row(record, site, distance, fasta, window, label)
            if row is not None:
                if str(row["gene_id"]) in excluded_genes:
                    continue
                if max_per_label is None or len(rows_by_label[label]) < max_per_label:
                    rows_by_label[label].append(row)
                if (
                    stop_when_balanced
                    and max_per_label is not None
                    and len(rows_by_label[0]) >= max_per_label
                    and len(rows_by_label[1]) >= max_per_label
                ):
                    break

            if scanned % 250_000 == 0:
                print(
                    "Scanned ClinVar SNVs: "
                    f"{scanned:,}; positives={len(rows_by_label[1]):,}; negatives={len(rows_by_label[0]):,}"
                )

        rows = rows_by_label[0] + rows_by_label[1]
        frame = balance_rows(rows, max_per_label=max_per_label, seed=seed)
        if frame.empty:
            raise ValueError(
                "No ClinVar SNVs passed filtering. Check CLNSIG labels, splice-related annotations, "
                "REF/FASTA compatibility, and --max-distance."
            )
        write_dataframe(output_path, frame)
        summary = frame.groupby(["variant_type", "label_name"], as_index=False).size().rename(columns={"size": "rows"})
        write_dataframe(output_path.with_name(output_path.stem + "_summary.csv"), summary)
        return frame
    finally:
        fasta.close()


def build_clinvar_smoke(output: Path = EXP3_DATA_DIR / "clinvar_splicing_variants_smoke.csv", rows: int = 12) -> pd.DataFrame:
    return build_clinvar_variants(output_path=output, max_per_label=max(1, rows // 2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a real ClinVar splice-variant benchmark.")
    parser.add_argument("--input", type=Path, default=DEFAULT_CLINVAR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--genome", type=Path, default=DEFAULT_GENOME)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--window", type=int, default=200)
    parser.add_argument("--max-distance", type=int, default=50)
    parser.add_argument("--max-per-label", type=int, default=250)
    parser.add_argument(
        "--scan-all",
        dest="stop_when_balanced",
        action="store_false",
        help="Scan the full ClinVar VCF instead of stopping when both labels reach --max-per-label.",
    )
    parser.set_defaults(stop_when_balanced=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chroms", default=DEFAULT_CLINVAR_CHROMS)
    parser.add_argument("--all-chroms", action="store_true")
    parser.add_argument("--protein-coding-only", action="store_true")
    parser.add_argument(
        "--exclude-split-genes-dir",
        type=Path,
        default=SHARED_SPLIT_DIR,
        help="Directory containing train/valid/test split CSVs whose genes should be excluded from ClinVar.",
    )
    parser.add_argument(
        "--allow-split-gene-overlap",
        action="store_true",
        help="Do not exclude genes already used by experiment 1/2 split CSVs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame = build_clinvar_variants(
        input_path=args.input,
        output_path=args.output,
        genome_path=args.genome,
        gtf_path=args.gtf,
        window=args.window,
        max_distance=args.max_distance,
        max_per_label=args.max_per_label,
        seed=args.seed,
        stop_when_balanced=args.stop_when_balanced,
        chroms=args.chroms,
        all_chroms=args.all_chroms,
        protein_coding_only=args.protein_coding_only,
        exclude_split_genes_dir=None if args.allow_split_gene_overlap else args.exclude_split_genes_dir,
    )
    print(frame.groupby(["variant_type", "label_name"]).size().reset_index(name="rows").to_string(index=False))


if __name__ == "__main__":
    main()
