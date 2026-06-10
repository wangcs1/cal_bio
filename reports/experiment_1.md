# Experiment 1: Splice-site three-class classification

Data directory: `data/shared/splits`

The main result uses the small split train/valid/test = 855/120/285.
RNA-FM/RNABERT rows are frozen-encoder models when local weights are present and
deterministic proxy embeddings otherwise. SpliceAI is an optional real-tool baseline
with a proxy fallback, evaluated only on the small test split.

## Validation summary

| model                                        |   train_rows |   valid_rows |   macro_f1 |   accuracy |   auroc |    auprc |
|:---------------------------------------------|-------------:|-------------:|-----------:|-----------:|--------:|---------:|
| CNN baseline (PyTorch Conv1D)                |          855 |          120 |   0.906517 |   0.925    | 0.98231 | 0.961331 |
| RNA-FM frozen encoder + MLP                  |          855 |          120 |   0.936268 |   0.941667 | 0.99502 | 0.989131 |
| RNABERT frozen encoder + MLP                 |          855 |          120 |   1        |   1        | 1       | 1        |
| SpliceAI optional real tool (proxy fallback) |          855 |          120 |   0.870884 |   0.9      | 0.9936  | 0.986976 |

## Test summary

| model                                        | model_type     | backend                            |   rows |   accuracy |   macro_f1 |    auroc |    auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |   hard_negative_fpr |   hard_negative_rows |
|:---------------------------------------------|:---------------|:-----------------------------------|-------:|-----------:|-----------:|---------:|---------:|-----------:|--------------:|----------------:|--------------------:|---------------------:|
| CNN baseline (PyTorch Conv1D)                | baseline       | pytorch_conv1d                     |    285 |   0.901754 |   0.900513 | 0.989931 | 0.98016  |   0.994764 |      0.869565 |        0.837209 |           0.41791   |                   67 |
| RNA-FM frozen encoder + MLP                  | frozen encoder | local_pretrained_or_proxy_fallback |    285 |   0.950877 |   0.951881 | 0.993841 | 0.988757 |   0.933333 |      0.994475 |        0.927835 |           0.149254  |                   67 |
| RNABERT frozen encoder + MLP                 | frozen encoder | local_pretrained_or_proxy_fallback |    285 |   0.992982 |   0.993059 | 0.999909 | 0.999825 |   0.994764 |      0.994413 |        0.99     |           0.0149254 |                   67 |
| SpliceAI optional real tool (proxy fallback) | proxy          | spliceai_signal_proxy              |    285 |   0.873684 |   0.868013 | 0.998053 | 0.996211 |   0.900474 |      0.923077 |        0.780488 |           0.522388  |                   67 |

Outputs:

- `results/experiment_1/tables/experiment_1_metrics.csv`
- `results/experiment_1/tables/experiment_1_confusion_matrices.csv`
- `results/experiment_1/figures/experiment_1_macro_f1.png`
- `results/experiment_1/figures/experiment_1_confusion_matrices.png`
- `results/checkpoints/experiment_1/*.joblib`
