# Experiment Log

## 2026-06-10

- Command: `python -m src.pipelines.run_c_part_all`
- Input data:
  - `data/shared/splits/train.csv`
  - `data/shared/splits/valid.csv`
  - `data/shared/splits/test.csv`
  - `data/experiment_3/artificial_variant_effect.csv`
- Output directories:
  - `results/experiment_1/`
  - `results/experiment_2/`
  - `results/experiment_3/`
  - `reports/`
- Notes:
  - Main experiments use the small split train/valid/test = 855/120/285.
  - Raw genome/annotation/ClinVar/GTEx resources and local RNA-FM/RNABERT weights are optional.
  - Proxy/fallback rows are explicitly labeled in model names and reports.
