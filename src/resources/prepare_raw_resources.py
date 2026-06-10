from __future__ import annotations

import argparse
from pathlib import Path

from src.utils import PROJECT_ROOT, REPORTS_DIR, ensure_dirs


RAW_FILES = ["genome.fa", "gencode.gtf", "clinvar.vcf", "gtex_sqtl.tsv", "known_splice_events.tsv", "gencode.db"]


def write_resource_report(raw_dir: Path = PROJECT_ROOT / "data/raw", output: Path = REPORTS_DIR / "resource_setup.md") -> Path:
    ensure_dirs(output.parent)
    lines = [
        "# Optional Resource Setup",
        "",
        "The main C Part experiments use the small synthetic/proxy split and do not require raw genome files.",
        "`data/raw/` is reserved for optional real-data smoke tests and case studies.",
        "",
        "| Resource | Status | Purpose |",
        "| --- | --- | --- |",
    ]
    for name in RAW_FILES:
        path = raw_dir / name
        status = "present" if path.exists() else "optional / missing"
        lines.append(f"| `{name}` | {status} | optional real-resource smoke or case study |")
    lines.extend(
        [
            "",
            "Large raw files and `models/hf/` pretrained weights are intentionally excluded from git.",
            "Default experiment commands remain runnable without these resources.",
        ]
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write optional raw-resource setup notes.")
    parser.add_argument("--raw-dir", type=Path, default=PROJECT_ROOT / "data/raw")
    parser.add_argument("--output", type=Path, default=REPORTS_DIR / "resource_setup.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    path = write_resource_report(args.raw_dir, args.output)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
