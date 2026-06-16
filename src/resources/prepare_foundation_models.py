from __future__ import annotations

import argparse
from pathlib import Path

from src.utils import PROJECT_ROOT


MODEL_DIRS = {
    "rnafm": PROJECT_ROOT / "models/hf/rnafm",
    "rnabert": PROJECT_ROOT / "models/hf/rnabert",
}


def check_foundation_models() -> list[dict[str, str]]:
    rows = []
    for name, path in MODEL_DIRS.items():
        rows.append(
            {
                "model": name,
                "path": str(path),
                "status": "present" if path.exists() else "missing; real-model pipeline will fail until downloaded",
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check optional local RNA-FM/RNABERT model directories.")
    parser.add_argument("--check", action="store_true", help="Only check local model directories.")
    return parser.parse_args()


def main() -> None:
    parse_args()
    for row in check_foundation_models():
        print(f"{row['model']}: {row['status']} ({row['path']})")


if __name__ == "__main__":
    main()
