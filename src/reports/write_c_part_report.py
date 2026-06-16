from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd

from src.utils import (
    EXP1_TABLES_DIR,
    EXP2_TABLES_DIR,
    EXP3_DATA_DIR,
    EXP3_TABLES_DIR,
    PROJECT_ROOT,
    REPORTS_DIR,
    SHARED_SPLIT_DIR,
    load_or_empty,
)


def markdown_table(frame: pd.DataFrame, max_rows: int = 12) -> str:
    if frame.empty:
        return "_Not generated yet._"
    subset = frame.head(max_rows).copy()
    for column in subset.columns:
        if pd.api.types.is_float_dtype(subset[column]):
            subset[column] = subset[column].map(lambda value: "" if pd.isna(value) else f"{value:.4f}")
    return subset.to_markdown(index=False)


def write_report(root: Path = PROJECT_ROOT) -> Path:
    split_counts = []
    for split in ["train", "valid", "test"]:
        path = SHARED_SPLIT_DIR / f"{split}.csv"
        if path.exists():
            frame = pd.read_csv(path)
            split_counts.append({"split": split, "rows": len(frame)})
    split_frame = pd.DataFrame(split_counts)

    exp1 = load_or_empty(EXP1_TABLES_DIR / "experiment_1_metrics.csv")
    exp2 = load_or_empty(EXP2_TABLES_DIR / "experiment_2B_hard_negative.csv")
    exp3 = load_or_empty(EXP3_TABLES_DIR / "experiment_3A_variant_metrics.csv")
    variants = load_or_empty(EXP3_DATA_DIR / "clinvar_splicing_variants_summary.csv")

    text = f"""# C Part Combined Report

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Scope

The main C Part pipeline uses the current small split, not a million-scale full dataset.
Experiments 1 and 2 use a real GENCODE/GRCh38 chromosome-holdout split with balanced
donor/acceptor/non-splice classes. This run reports real-model rows only. Experiments 1
and 2 use CNN, RNA-FM frozen encoder, and RNABERT frozen encoder. Experiment 3 additionally
includes real external splice tools: SpliceAI, Pangolin, MMSplice, and MaxEntScan.

## Data

{markdown_table(split_frame)}

Real ClinVar variant summary:

{markdown_table(variants)}

## Experiment 1

Splice-site donor/acceptor/non-splice classification with the real-model rows only.

{markdown_table(exp1)}

## Experiment 2

Multi-scale context, hard-negative, and rare-motif stress tests using the same real-model
set.

{markdown_table(exp2)}

## Experiment 3

<<<<<<< HEAD
Real ClinVar variant effect prediction plus small format-control case studies. The
main variant table includes both trained project models and real external splice tools;
no proxy/fallback rows are reported.
=======
Artificial variant effect prediction plus small ClinVar/sQTL-format smoke inputs. The
main artificial-variant table includes both trained project models and real external
splice-tool rows.
>>>>>>> 2fa1a6f0c0e7b80a5e169c5f05ddc43f5c3767f8

{markdown_table(exp3)}

## Reports

- `reports/experiment_1.md`
- `reports/experiment_2.md`
- `reports/experiment_3.md`
- `reports/data_qc.md`
- `reports/resource_setup.md`
"""

    out = REPORTS_DIR / "c_part_combined_report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out


def main() -> None:
    path = write_report()
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
