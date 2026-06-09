# Experiment 1: Splice-site three-class classification

Data directory: `/home/hkqxhy/cal_bio/splits`

This run uses the installed CPU PyTorch environment (`torch 2.8.0+cpu`). The CNN baseline is
a real Conv1D neural network trained with PyTorch. RNA-FM and RNABERT are still implemented as
frozen-representation style k-mer/token proxies because `esm`, `transformers`, and
`multimolecule` model dependencies/weights are not installed locally.

## Validation summary

| model                         |   train_rows |   valid_rows |   macro_f1 |   accuracy |    auroc |    auprc |
|:------------------------------|-------------:|-------------:|-----------:|-----------:|---------:|---------:|
| CNN baseline (PyTorch Conv1D) |        60000 |        15000 |   0.865967 |   0.869267 | 0.966852 | 0.935381 |

## Test summary

| model                         |   rows |   accuracy |   macro_f1 |    auroc |    auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |
|:------------------------------|-------:|-----------:|-----------:|---------:|---------:|-----------:|--------------:|----------------:|
| CNN baseline (PyTorch Conv1D) |  30000 |   0.874000 |   0.871174 | 0.968419 | 0.938427 |   0.944290 |      0.873078 |        0.796155 |
| RNA-FM + MLP                  |  30000 |   0.845367 |   0.841689 | 0.943919 | 0.895642 |   0.911695 |      0.856929 |        0.756442 |
| RNABERT + MLP                 |  30000 |   0.869667 |   0.869450 | 0.958412 | 0.920519 |   0.930564 |      0.863340 |        0.814445 |

Outputs:

- `results/tables/experiment_1_metrics.csv`
- `results/tables/experiment_1_confusion_matrices.csv`
- `results/figures/experiment_1_macro_f1.png`
- `results/figures/experiment_1_confusion_matrices.png`
- `results/checkpoints/experiment_1/*.joblib`
