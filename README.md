# RNA Splice-Site Recognition and ClinVar Variant-Effect Benchmark

This repository contains a small but fully real-data benchmark for computational biology coursework. The current main results use GRCh38, GENCODE v49, and ClinVar to evaluate splice-site recognition, GT/AG hard-negative stress tests, and ClinVar splice-altering-vs-benign SNV ranking.

## Research Questions

1. Can models distinguish donor, acceptor, and GT/AG motif-matched hard-negative sites on a real GENCODE/GRCh38 chromosome-holdout split?
2. Does context-window length affect rejection of motif-matched negatives?
3. Can splice-site recognition signals be transferred to real ClinVar pathogenic-vs-benign variant-effect ranking, and how much of the signal remains after exact distance matching?

## Data

The splice-site benchmark is built by `src/data/build_splice_site_dataset.py` from `data/raw/genome.fa` and `data/raw/gencode.gtf`. The main benchmark samples donor / acceptor / hard-negative = 1000 / 1000 / 1000. The non-splice class is made of real-genome GT/AG motif-matched decoys, not easy random negatives.

| split | rows | donor | acceptor | hard-negative |
| --- | ---: | ---: | ---: | ---: |
| train | 2339 | 792 | 786 | 761 |
| valid | 230 | 79 | 67 | 84 |
| test | 431 | 129 | 147 | 155 |

The ClinVar benchmark is built by `src/data/build_clinvar_variant_dataset.py`, with 250 splice-altering SNVs and 250 benign/likely benign near-splice SNVs restricted to held-out test chromosomes. Because splice-altering variants are closer to annotated splice sites, the project also exports `data/experiment_3/clinvar_splicing_variants_distance_matched.csv` for confounding diagnostics.

QC is recorded in `reports/data_qc.md`. It includes split integrity, ClinVar REF checks, distance-only baselines, and a cross-split near-duplicate audit covering full-window duplicates, center 161 bp duplicates, and high 9-mer Jaccard near-duplicates.

## Models

Main experiment 1/2 models:

- CNN baseline
- RNA-FM frozen encoder + MLP
- RNABERT frozen encoder + MLP

Experiment 3 additionally includes real external tool outputs:

- SpliceAI
- Pangolin
- MMSplice
- MaxEntScan

The main benchmark is fail-fast: RNA-FM/RNABERT require local real weights under `models/hf/`, and external tools require their configured real tool environments. No k-mer fallback, proxy row, or placeholder metric is used in the main reported tables.

## Reproduction

```bash
python -m src.data.build_splice_site_dataset --max-per-class 1000 --windows 50 100 200 400
python -m src.data.build_clinvar_variant_dataset
python -m src.data.qc_splice_dataset
python -m src.experiments.exp1.run_classification --full-data
python scripts/run_exp1_multiseed.py --full-data
python -m src.experiments.exp2.run_multiscale
python -m src.experiments.exp3.run_variant_effect
python -m src.experiments.exp3.run_interpretability
python -m src.reports.write_c_part_report
python scripts/check_latex_report.py
```

## Key Results

Experiment 1 test-set classification:

| model | Accuracy | Macro-F1 | AUROC | AUPRC | Hard FPR |
| --- | ---: | ---: | ---: | ---: | ---: |
| CNN baseline | 0.8237 | 0.8248 | 0.9489 | 0.9070 | 0.3226 |
| RNA-FM frozen encoder + MLP | 0.7633 | 0.7670 | 0.9157 | 0.8402 | 0.3290 |
| RNABERT frozen encoder + MLP | 0.8190 | 0.8199 | 0.9298 | 0.8718 | 0.3226 |

Experiment 1 multi-seed summary is stored in `results/experiment_1/tables/experiment_1_multiseed_summary.csv`. CNN Macro-F1 is 0.8273 ± 0.0124 over seeds 42/43/44; RNA-FM/RNABERT are deterministic in the current frozen-encoder setup.

Experiment 3 full ClinVar ranking:

| model | AUROC | AUPRC |
| --- | ---: | ---: |
| RNA-FM frozen encoder + MLP | 0.7532 | 0.7839 |
| Pangolin real sequence model | 0.7023 | 0.7781 |
| CNN baseline | 0.6959 | 0.7374 |
| RNABERT frozen encoder + MLP | 0.6673 | 0.7326 |
| SpliceAI real sequence model | 0.6425 | 0.7289 |
| MaxEntScan real local score | 0.6165 | 0.6716 |
| MMSplice real sequence model | 0.5000 | 0.5000 |

Distance-matched ClinVar results drop for most models, showing that near-splice distance is an important confounder in the full benchmark.

## Reports

- Experiment reports: `reports/experiment_1.md`, `reports/experiment_2.md`, `reports/experiment_3.md`
- Combined report: `reports/c_part_combined_report.md`
- Data QC: `reports/data_qc.md`
- Main pipeline audit: `reports/main_pipeline_audit.md`
- LaTeX paper: `report_letax/njuthesis-sample.tex`
- Paper figures: `report_letax/images/`

## Limitations

This is a small real-data benchmark, not a clinical-grade or genome-wide evaluation. Experiment 2/3 still use seed 42 as the main run; ClinVar is limited to near-splice SNVs; distance matching diagnoses but does not eliminate every confounder; and the near-duplicate audit is not a substitute for alignment-level paralog clustering with BLAST/CD-HIT/MMseqs.
