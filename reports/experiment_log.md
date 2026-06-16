# Experiment Log

## 2026-06-16

- Command: `python -m src.pipelines.run_c_part_all`
- Input data:
- `data/shared/splits/train.csv`
- `data/shared/splits/valid.csv`
- `data/shared/splits/test.csv`
- `data/experiment_3/clinvar_splicing_variants.csv`
- Required raw resources:
- `data/raw/genome.fa`
- `data/raw/gencode.gtf`
- `data/raw/clinvar.vcf`
- Local model weights:
- `models/hf/rnafm/`
- `models/hf/rnabert/`
- Output directories:
- `results/experiment_1/`
- `results/experiment_2/`
- `results/experiment_3/`
- `reports/`
- Notes:
- Main experiments use a small real chromosome-holdout split after raw resources are provided.
- Current main results include only CNN baseline, RNA-FM frozen encoder + MLP, and RNABERT frozen encoder + MLP.
- RNA-FM/RNABERT loading fails fast if local real weights are unavailable.
