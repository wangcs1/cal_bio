# Experiment 1: Splice-site three-class classification

Data directory: `data\shared\splits`

The main result uses the current real chromosome-holdout small split; no row caps are applied in the full-data run.
RNA-FM/RNABERT rows require real local frozen encoder weights under `models/hf/`.
Only real-model rows are reported in this run.
CNN rows report the multi-seed seed=42 values (`experiment_1_multiseed_metrics.csv`); RNA-FM/RNABERT are deterministic, so single-run and multi-seed values coincide.

## Validation summary

| model                         |   train_rows |   valid_rows |   macro_f1 |   accuracy |    auroc |    auprc |
|:------------------------------|-------------:|-------------:|-----------:|-----------:|---------:|---------:|
| CNN baseline (PyTorch Conv1D) |         2339 |          230 |   0.848824 |   0.852174 | 0.965901 | 0.934968 |
| RNA-FM frozen encoder + MLP   |         2339 |          230 |   0.809933 |   0.808696 | 0.93486  | 0.878588 |
| RNABERT frozen encoder + MLP  |         2339 |          230 |   0.819815 |   0.821739 | 0.955022 | 0.916894 |

## Test summary

| model                         | model_type     | backend                   |   rows |   accuracy |   macro_f1 |    auroc |    auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |   hard_negative_fpr |   hard_negative_rows |
|:------------------------------|:---------------|:--------------------------|-------:|-----------:|-----------:|---------:|---------:|-----------:|--------------:|----------------:|--------------------:|---------------------:|
| CNN baseline (PyTorch Conv1D) | baseline       | pytorch_conv1d            |    431 |   0.821346 |   0.822698 | 0.948953 | 0.906649 |   0.898876 |      0.834951 |        0.734266 |            0.322581 |                  155 |
| RNA-FM frozen encoder + MLP   | frozen encoder | local_pretrained_required |    431 |   0.763341 |   0.766977 | 0.915666 | 0.840231 |   0.816794 |      0.810997 |        0.673139 |            0.329032 |                  155 |
| RNABERT frozen encoder + MLP  | frozen encoder | local_pretrained_required |    431 |   0.819026 |   0.819869 | 0.929839 | 0.871822 |   0.888889 |      0.833876 |        0.736842 |            0.322581 |                  155 |

Outputs:

- `results/experiment_1/tables/experiment_1_multiseed_metrics.csv`
- `results/experiment_1/tables/experiment_1_multiseed_summary.csv`
- `results/experiment_1/tables/experiment_1_confusion_matrices.csv`
- `results/experiment_1/figures/experiment_1_macro_f1.png`
- `results/experiment_1/figures/experiment_1_confusion_matrices.png`
- `results/checkpoints/experiment_1/*.joblib`
