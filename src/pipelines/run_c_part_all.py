from __future__ import annotations

import argparse

from src.data.build_clinvar_variant_dataset import build_clinvar_variants
from src.data.build_splice_site_dataset import build_and_write_real
from src.data.qc_splice_dataset import build_qc_report
from src.experiments.exp1.run_classification import run as run_exp1
from src.experiments.exp2.run_multiscale import run as run_exp2
from src.experiments.exp3.run_variant_effect import run as run_exp3
from src.experiments.exp3.run_interpretability import run as run_interpretability
from src.reports.write_c_part_report import write_report
from src.resources.prepare_raw_resources import write_resource_report
from src.utils import EXP1_FIGURES_DIR, EXP1_TABLES_DIR, EXP2_FIGURES_DIR, EXP2_TABLES_DIR, EXP3_FIGURES_DIR, EXP3_TABLES_DIR, PROJECT_ROOT


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the full C-part real-model experiment pipeline.")
    parser.parse_args()

    print("[1/9] Building real GENCODE/GRCh38 splice-site benchmark")
    build_and_write_real()
    print("[2/9] Running experiment 1")
    run_exp1(
        data_dir=None,
        tables_dir=EXP1_TABLES_DIR,
        figures_dir=EXP1_FIGURES_DIR,
        checkpoint_dir=PROJECT_ROOT / "results/checkpoints/experiment_1",
        report_path=PROJECT_ROOT / "reports/experiment_1.md",
        model_keys=["cnn", "rnafm", "rnabert"],
        seed=42,
        max_train_rows=None,
        max_valid_rows=None,
        max_test_rows=None,
    )
    print("[3/9] Running experiment 2")
    run_exp2(EXP2_TABLES_DIR, EXP2_FIGURES_DIR, random_state=42)
    print("[4/9] Building real ClinVar variant dataset")
    build_clinvar_variants()
    print("[5/9] Running experiment 3")
    run_exp3(EXP3_TABLES_DIR, EXP3_FIGURES_DIR, random_state=42)
    print("[6/9] Running interpretability analyses")
    run_interpretability(EXP3_TABLES_DIR, EXP3_FIGURES_DIR, random_state=42)
    print("[7/9] Writing QC/resource reports")
    build_qc_report()
    write_resource_report()
    print("[8/9] Writing combined C Part report")
    write_report()
    print("Done")


if __name__ == "__main__":
    main()
