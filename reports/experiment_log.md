# Experiment Log

## 2026-06-16

- Command: `python -m src.pipelines.run_c_part_all`
- Input data:
- `data/shared/splits/train.csv`
- `data/shared/splits/valid.csv`
- `data/shared/splits/test.csv`
- `data/experiment_3/artificial_variant_effect.csv`
- Local model weights:
- `models/hf/rnafm/`
- `models/hf/rnabert/`
- Output directories:
- `results/experiment_1/`
- `results/experiment_2/`
- `results/experiment_3/`
- `reports/`
- Notes:
- Main experiments use the small split train/valid/test = 855/120/285.
- Current main results include only CNN baseline, RNA-FM frozen encoder + MLP, and RNABERT frozen encoder + MLP.
- RNA-FM/RNABERT loading fails fast if local real weights are unavailable.
