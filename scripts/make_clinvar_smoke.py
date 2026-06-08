from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    rows: list[dict[str, str]] = []
    with (PROJECT_ROOT / "data/raw/clinvar.vcf").open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
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
            if len(rows) >= 5:
                break
    out_path = PROJECT_ROOT / "data/raw/clinvar_smoke.csv"
    with out_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["CHROM", "POS", "REF", "ALT"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_path}: {len(rows)} rows")


if __name__ == "__main__":
    main()
