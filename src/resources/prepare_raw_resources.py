from __future__ import annotations

import argparse
from pathlib import Path

from src.utils import PROJECT_ROOT, REPORTS_DIR, ensure_dirs


RAW_FILES = ["genome.fa", "gencode.gtf", "clinvar.vcf", "gtex_sqtl.tsv", "known_splice_events.tsv", "gencode.db"]


def write_resource_report(raw_dir: Path = PROJECT_ROOT / "data/raw", output: Path = REPORTS_DIR / "resource_setup.md") -> Path:
    ensure_dirs(output.parent)
    lines = [
        "# Real Resource Setup",
        "",
        "The main C Part experiments now use real local resources.",
        "`data/raw/genome.fa`, `data/raw/gencode.gtf`, and `data/raw/clinvar.vcf` are required for the default real-data pipeline.",
        "",
        "| Resource | Status | Purpose |",
        "| --- | --- | --- |",
    ]
    for name in RAW_FILES:
        path = raw_dir / name
        required = name in {"genome.fa", "gencode.gtf", "clinvar.vcf"}
        if path.exists():
            status = f"present ({path.stat().st_size:,} bytes)"
        elif required:
            status = "required / missing"
        else:
            status = "optional / missing"
        purpose = "default real benchmark" if required else "optional real-resource case study"
        lines.append(f"| `{name}` | {status} | {purpose} |")
    lines.extend(
        [
            "",
            "Large raw files and `models/hf/` pretrained weights are intentionally excluded from git.",
            "Default experiment commands fail fast when required real resources are missing.",
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
