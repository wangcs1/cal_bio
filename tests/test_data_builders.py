from __future__ import annotations

import pandas as pd

from src.data.build_clinvar_variant_dataset import build_clinvar_variants
from src.data.build_splice_site_dataset import (
    FastaIndex,
    build_and_write_real,
    build_positive_sites,
    filter_positive_sites,
    oriented_window,
    read_transcripts,
)
from src.data.build_synthetic_splice_dataset import build_dataset, write_windowed_outputs
from src.data.build_variant_dataset import build_variant_dataset
from src.utils import REQUIRED_SPLIT_COLUMNS


def write_real_fixtures(tmp_path):
    def make_chrom(strand: str) -> str:
        chars = list(("GTAG" * 80)[:240])
        if strand == "+":
            chars[70] = "G"
            chars[71] = "T"
            chars[88] = "A"
            chars[89] = "G"
        else:
            chars[88] = "A"
            chars[89] = "C"
            chars[70] = "C"
            chars[71] = "T"
        return "".join(chars)

    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(
        "\n".join(
            [
                ">chr1",
                make_chrom("+"),
                ">chr17",
                make_chrom("-"),
                ">chr19",
                make_chrom("+"),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    gtf_path = tmp_path / "gencode.gtf"
    rows = []
    for chrom, strand, gene, tx in [
        ("chr1", "+", "GENE1", "TX1"),
        ("chr17", "-", "GENE2", "TX2"),
        ("chr19", "+", "GENE3", "TX3"),
    ]:
        attrs = f'gene_id "{gene}"; transcript_id "{tx}"; gene_type "protein_coding"; transcript_type "protein_coding";'
        rows.append(f"{chrom}\ttest\texon\t50\t70\t.\t{strand}\t.\t{attrs}")
        rows.append(f"{chrom}\ttest\texon\t91\t110\t.\t{strand}\t.\t{attrs}")
    gtf_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return fasta_path, gtf_path


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


def test_real_builder_extracts_stranded_canonical_sites_and_splits(tmp_path):
    fasta_path, gtf_path = write_real_fixtures(tmp_path)
    out = tmp_path / "processed"
    splits = tmp_path / "splits"
    outputs = build_and_write_real(
        genome_path=fasta_path,
        gtf_path=gtf_path,
        out_dir=out,
        split_dir=splits,
        windows=[5],
        max_per_class=3,
        negative_radius=40,
        negative_attempts=200,
        seed=7,
        chroms="chr1,chr17,chr19",
    )
    frame = pd.DataFrame(outputs[5])
    assert REQUIRED_SPLIT_COLUMNS.issubset(frame.columns)
    assert {"data_source", "motif_type", "split"}.issubset(frame.columns)
    assert set(frame["label_name"]) == {"donor", "acceptor", "non_splice"}
    assert set(pd.read_csv(splits / "train_pm5.csv")["chrom"]) == {"chr1"}
    assert set(pd.read_csv(splits / "valid_pm5.csv")["chrom"]) == {"chr17"}
    assert set(pd.read_csv(splits / "test_pm5.csv")["chrom"]) == {"chr19"}
    assert frame[frame["label"].astype(int) == 2]["center"].isin({70, 91}).sum() == 0


def test_real_builder_positive_filter_and_reverse_complement(tmp_path):
    fasta_path, gtf_path = write_real_fixtures(tmp_path)
    fasta = FastaIndex(fasta_path)
    try:
        transcripts = read_transcripts(gtf_path, fasta, protein_coding_only=False, gtf_line_limit=None)
        sites, _bounds, _annotated_any, _annotated_same = build_positive_sites(transcripts)
        filtered = filter_positive_sites(
            sites,
            fasta,
            windows=[5],
            max_n_ratio=0.0,
            allowed_chroms={"chr1", "chr17", "chr19"},
            canonical_only=True,
        )
        assert len(filtered) == 6
        minus_donor = next(site for site in filtered if site.chrom == "chr17" and site.label == 0)
        window = oriented_window(fasta, minus_donor.chrom, minus_donor.center, minus_donor.strand, 2)
        assert window is not None
        assert window[3:5] == "GT"
    finally:
        fasta.close()


def test_clinvar_builder_parses_real_snv_labels_and_negative_strand_complement(tmp_path):
    fasta_path, gtf_path = write_real_fixtures(tmp_path)
    vcf_path = tmp_path / "clinvar.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t71\tpos_plus\tG\tA\t.\t.\tCLNSIG=Pathogenic;MC=SO:0001575|splice_donor_variant",
                "chr17\t90\tpos_minus\tC\tA\t.\t.\tCLNSIG=Likely_pathogenic;MC=SO:0001575|splice_donor_variant",
                "chr1\t89\tbenign_plus\tA\tC\t.\t.\tCLNSIG=Benign;MC=SO:0001627|intron_variant",
                "chr17\t71\tbenign_minus\tC\tT\t.\t.\tCLNSIG=Likely_benign;MC=SO:0001627|intron_variant",
                "chr1\t72\tindel_skip\tT\tTA\t.\t.\tCLNSIG=Pathogenic;MC=SO:0001575|splice_donor_variant",
                "chr1\t73\tmulti_skip\tA\tC,G\t.\t.\tCLNSIG=Pathogenic;MC=SO:0001575|splice_donor_variant",
                "chr1\t74\tref_mismatch\tC\tA\t.\t.\tCLNSIG=Pathogenic;MC=SO:0001575|splice_donor_variant",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "clinvar_splicing_variants.csv"
    frame = build_clinvar_variants(
        input_path=vcf_path,
        output_path=output,
        genome_path=fasta_path,
        gtf_path=gtf_path,
        window=5,
        max_distance=5,
        max_per_label=2,
        seed=3,
        chroms="chr1,chr17,chr19",
    )
    assert output.exists()
    assert set(frame["label_name"]) == {"splice_altering", "neutral"}
    assert not {"indel_skip", "multi_skip", "ref_mismatch"}.intersection(set(frame["clinvar_id"]))
    minus = frame[frame["clinvar_id"] == "pos_minus"].iloc[0]
    assert minus["strand"] == "-"
    assert minus["ref"] == "G"
    assert minus["alt"] == "T"
    assert minus["wt_sequence"][int(5 - (90 - 91))] == "G"
