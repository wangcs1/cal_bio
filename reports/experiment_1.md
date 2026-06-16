# Experiment 1: Splice-site three-class classification

Data directory: `D:\CAL_BIO\data\shared\splits`

The main result uses the small split train/valid/test = 855/120/285.
RNA-FM/RNABERT rows require real local frozen encoder weights under `models/hf/`.
Only real-model rows are reported in this run.

## Validation summary

| model                         |   train_rows |   valid_rows |   macro_f1 |   accuracy |   auroc |    auprc |
|:------------------------------|-------------:|-------------:|-----------:|-----------:|--------:|---------:|
| CNN baseline (PyTorch Conv1D) |          855 |          120 |   0.85903  |   0.891667 | 0.98218 | 0.961394 |
| RNA-FM frozen encoder + MLP   |          855 |          120 |   0.944841 |   0.95     | 0.99556 | 0.99023  |
| RNABERT frozen encoder + MLP  |          855 |          120 |   1        |   1        | 1       | 1        |

## Test summary

| model                         | model_type     | backend                   |   rows |   accuracy |   macro_f1 |    auroc |    auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |   hard_negative_fpr |   hard_negative_rows |
|:------------------------------|:---------------|:--------------------------|-------:|-----------:|-----------:|---------:|---------:|-----------:|--------------:|----------------:|--------------------:|---------------------:|
| CNN baseline (PyTorch Conv1D) | baseline       | pytorch_conv1d            |    285 |   0.887719 |   0.88579  | 0.991394 | 0.982513 |   0.994764 |      0.853081 |        0.809524 |           0.477612  |                   67 |
| RNA-FM frozen encoder + MLP   | frozen encoder | local_pretrained_required |    285 |   0.950877 |   0.951663 | 0.993585 | 0.988255 |   0.938144 |      0.989011 |        0.927835 |           0.149254  |                   67 |
| RNABERT frozen encoder + MLP  | frozen encoder | local_pretrained_required |    285 |   0.992982 |   0.993059 | 0.999909 | 0.999825 |   0.994764 |      0.994413 |        0.99     |           0.0149254 |                   67 |

Outputs:

- `results/experiment_1/tables/experiment_1_metrics.csv`
- `results/experiment_1/tables/experiment_1_confusion_matrices.csv`
- `results/experiment_1/figures/experiment_1_macro_f1.png`
- `results/experiment_1/figures/experiment_1_confusion_matrices.png`
- `results/checkpoints/experiment_1/*.joblib`
