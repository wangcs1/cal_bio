# Main Pipeline Audit

This audit separates the current reported benchmark from historical helper code.

## Main Reported Pipeline

- Data QC: `python -m src.data.qc_splice_dataset`
- Experiment 1: `python -m src.experiments.exp1.run_classification --full-data`
- Experiment 1 multi-seed: `python scripts/run_exp1_multiseed.py --full-data`
- Experiment 2: `python -m src.experiments.exp2.run_multiscale`
- Experiment 3: `python -m src.experiments.exp3.run_variant_effect`
- Interpretability: `python -m src.experiments.exp3.run_interpretability`
- LaTeX sanity check: `python scripts/check_latex_report.py`

## Fail-Fast Checks

- RNA-FM and RNABERT require local pretrained weights under `models/hf/rnafm/` and `models/hf/rnabert/`.
- If local RNA-FM/RNABERT weights or dependencies are unavailable, `src.models.foundation_backbones` raises an error instead of using k-mer fallback features or remote downloads.
- Experiment 3 external tools require configured real tool environments. Missing environments raise errors.
- Experiment 2 rare-motif stress-test data must exist; missing data raises an error instead of writing placeholder `not_run` rows.

## What Is Not In The Main Tables

- `src/models/simple_splice_models.py` contains historical proxy models and is not imported by the current main experiment scripts.
- `src/experiments/exp2/run_full_context.py` contains historical proxy/case-study analyses and is not used by the current paper tables.
- Borzoi/AlphaGenome notes remain historical long-range case-study notes and are not part of the reported benchmark.
- Format-control ClinVar/sQTL subsets in experiment 3 are used only to verify input/output wiring; the main evidence is the 500-row ClinVar benchmark and the 326-row exact-distance-matched subset.

## Current Boundary

The current benchmark uses real GRCh38/GENCODE/ClinVar data, real local RNA-FM/RNABERT frozen encoders, and real external SpliceAI/Pangolin/MMSplice/MaxEntScan tool calls. It does not claim clinical-grade variant interpretation, genome-wide generalization, tissue-specific splicing modeling, or alignment-level paralog clustering.
