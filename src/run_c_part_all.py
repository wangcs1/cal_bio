from __future__ import annotations

from src.build_synthetic_splice_dataset import build_and_write
from src.build_variant_dataset import build_and_write_variants
from src.run_exp2_multiscale import run as run_exp2
from src.run_exp3_variant_effect import run as run_exp3
from src.run_full_context_question import run as run_full_context
from src.run_interpretability import run as run_interpretability
from src.utils import PROJECT_ROOT


def main() -> None:
    tables = PROJECT_ROOT / "results/tables"
    figures = PROJECT_ROOT / "results/figures"
    print("[1/6] Building synthetic splice-site benchmark")
    build_and_write()
    print("[2/6] Running experiment 2")
    run_exp2(tables, figures, random_state=42)
    print("[3/6] Building artificial variant dataset")
    build_and_write_variants()
    print("[4/6] Running experiment 3")
    run_exp3(tables, figures, random_state=42)
    print("[5/6] Running interpretability analyses")
    run_interpretability(tables, figures, random_state=42)
    print("[6/6] Running full context question experiments")
    run_full_context(tables, figures)
    print("Done")


if __name__ == "__main__":
    main()
