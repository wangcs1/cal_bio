from __future__ import annotations

import csv
import argparse
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a tiny ClinVar SNV smoke-test set.")
    parser.add_argument("--limit", type=int, default=5, help="Number of single-nucleotide variants to keep.")
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    vcf_meta: list[str] = []
    vcf_columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    vcf_rows: list[str] = []
    vcf_contigs: set[str] = set()
    existing_contigs: set[str] = set()
    with (PROJECT_ROOT / "data/raw/clinvar.vcf").open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##contig=<ID="):
                existing_contigs.add(line.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                vcf_meta.append(line)
                continue
            if line.startswith("##"):
                vcf_meta.append(line)
                continue
            if line.startswith("#CHROM"):
                vcf_columns = line
                continue
            chrom, pos, _vid, ref, alt, *_rest = line.rstrip("\n").split("\t")
            if len(ref) == 1 and len(alt) == 1 and "," not in alt:
                rows.append(
                    {
                        "CHROM": chrom if chrom.startswith("chr") else f"chr{chrom}",
                        "POS": pos,
                        "REF": ref,
                        "ALT": alt,
                    }
                )
                vcf_contigs.add(chrom)
                vcf_rows.append(line)
            if len(rows) >= args.limit:
                break
    out_path = PROJECT_ROOT / "data/raw/clinvar_smoke.csv"
    with out_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["CHROM", "POS", "REF", "ALT"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    vcf_path = PROJECT_ROOT / "data/raw/clinvar_smoke.vcf"
    with vcf_path.open("w", encoding="utf-8", newline="") as out:
        out.writelines(vcf_meta)
        for chrom in sorted(vcf_contigs - existing_contigs):
            out.write(f"##contig=<ID={chrom}>\n")
        out.write(vcf_columns)
        out.writelines(vcf_rows)
    print(f"wrote {out_path}: {len(rows)} rows")
    print(f"wrote {vcf_path}: {len(vcf_rows)} records")


if __name__ == "__main__":
    main()
