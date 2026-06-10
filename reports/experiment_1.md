# Experiment 1: Splice-site three-class classification

Data directory: `/mnt/d/CAL_BIO/data/shared/splits`

The main result uses the small split train/valid/test = 855/120/285.
RNA-FM/RNABERT rows are frozen-encoder models when local weights are present and
deterministic proxy embeddings otherwise. SpliceAI is an optional real-tool baseline
with a proxy fallback, evaluated only on the small test split.

## Validation summary

| model                                        |   train_rows |   valid_rows |   macro_f1 |   accuracy |    auroc |    auprc |
|:---------------------------------------------|-------------:|-------------:|-----------:|-----------:|---------:|---------:|
| CNN baseline (PyTorch Conv1D)                |          855 |          120 |   0.85903  |   0.891667 | 0.980611 | 0.958264 |
| RNA-FM frozen encoder + MLP                  |          855 |          120 |   0.981123 |   0.983333 | 0.999303 | 0.998251 |
| RNABERT frozen encoder + MLP                 |          855 |          120 |   1        |   1        | 1        | 1        |
| SpliceAI optional real tool (proxy fallback) |          855 |          120 |   0.870884 |   0.9      | 0.9936   | 0.986976 |

## Test summary

| model                                        |   rows |   accuracy |   macro_f1 |    auroc |    auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |   hard_negative_fpr |   hard_negative_rows |
|:---------------------------------------------|-------:|-----------:|-----------:|---------:|---------:|-----------:|--------------:|----------------:|--------------------:|---------------------:|
| CNN baseline (PyTorch Conv1D)                |    285 |   0.887719 |   0.88579  | 0.990841 | 0.981421 |   0.994764 |      0.853081 |        0.809524 |           0.477612  |                   67 |
| RNA-FM frozen encoder + MLP                  |    285 |   0.961404 |   0.962358 | 0.997192 | 0.994774 |   0.941799 |      1        |        0.945274 |           0.0746269 |                   67 |
| RNABERT frozen encoder + MLP                 |    285 |   0.968421 |   0.969041 | 0.998288 | 0.996833 |   0.958333 |      0.994475 |        0.954315 |           0.0895522 |                   67 |
| SpliceAI optional real tool (proxy fallback) |    285 |   0.873684 |   0.868013 | 0.998053 | 0.996211 |   0.900474 |      0.923077 |        0.780488 |           0.522388  |                   67 |

Outputs:

- `results/experiment_1/tables/experiment_1_metrics.csv`
- `results/experiment_1/tables/experiment_1_confusion_matrices.csv`
- `results/experiment_1/figures/experiment_1_macro_f1.png`
- `results/experiment_1/figures/experiment_1_confusion_matrices.png`
- `results/checkpoints/experiment_1/*.joblib`
