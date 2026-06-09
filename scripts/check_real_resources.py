from __future__ import annotations

import importlib.util
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_FILES = [
    "genome.fa",
    "gencode.gtf",
    "clinvar.vcf",
    "gtex_sqtl.tsv",
    "known_splice_events.tsv",
]
PACKAGES = [
    "torch",
    "transformers",
    "esm",
    "multimolecule",
    "Bio",
    "pyfaidx",
    "pangolin",
    "spliceai",
    "tensorflow",
    "keras",
    "cyvcf2",
    "maxentpy",
    "mmsplice",
]


def main() -> None:
    raw = PROJECT_ROOT / "data/raw"
    print("Raw data:")
    for name in RAW_FILES:
        path = raw / name
        if path.exists():
            print(f"  OK {name}: {path.stat().st_size:,} bytes")
        else:
            print(f"  MISSING {name}")

    print("\nPython packages:")
    for package in PACKAGES:
        print(f"  {'OK' if importlib.util.find_spec(package) else 'MISSING'} {package}")

    try:
        import torch

        print("\nGPU:")
        print(f"  cuda_available={torch.cuda.is_available()}")
        print(f"  device_count={torch.cuda.device_count()}")
        if torch.cuda.is_available():
            print(f"  device_0={torch.cuda.get_device_name(0)}")
    except Exception as exc:  # pragma: no cover
        print(f"\nGPU check failed: {exc}")


if __name__ == "__main__":
    main()
