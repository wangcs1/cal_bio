from __future__ import annotations

import argparse

from src.data.build_synthetic_splice_dataset import build_and_write
from src.data.build_variant_dataset import build_and_write_variants
from src.experiments.exp2.run_multiscale import run as run_exp2
from src.experiments.exp3.run_variant_effect import run as run_exp3
from src.experiments.exp2.run_full_context import run as run_full_context
from src.experiments.exp3.run_interpretability import run as run_interpretability
from src.utils import EXP2_FIGURES_DIR, EXP2_TABLES_DIR, EXP3_FIGURES_DIR, EXP3_TABLES_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the full C-part synthetic/proxy experiment pipeline.")
    parser.parse_args()

    print("[1/6] Building synthetic splice-site benchmark")
    build_and_write()
    print("[2/6] Running experiment 2")
    run_exp2(EXP2_TABLES_DIR, EXP2_FIGURES_DIR, random_state=42)
    print("[3/6] Building artificial variant dataset")
    build_and_write_variants()
    print("[4/6] Running experiment 3")
    run_exp3(EXP3_TABLES_DIR, EXP3_FIGURES_DIR, random_state=42)
    print("[5/6] Running interpretability analyses")
    run_interpretability(EXP3_TABLES_DIR, EXP3_FIGURES_DIR, random_state=42)
    print("[6/6] Running full context question experiments")
    run_full_context(EXP2_TABLES_DIR, EXP2_FIGURES_DIR)
    print("Done")


if __name__ == "__main__":
    main()
