# Experiment 1: Real RNA Foundation Model Run

Data directory: `/home/hkqxhy/cal_bio/splits`

This run uses real MultiMolecule pretrained backbones downloaded from Hugging Face:

- RNA-FM: `multimolecule/rnafm`, loaded from `models/hf/rnafm/model.safetensors`
- RNABERT: `multimolecule/rnabert`, loaded from `models/hf/rnabert/model.safetensors`

The backbone encoders are frozen. Sequence embeddings are mean-pooled over non-padding tokens, concatenated with engineered splice-site signal features, and passed to a balanced logistic/MLP-style classifier head. Embeddings are cached under `results/embeddings/experiment_1/`.

CPU-only constraint: CUDA is not available, so full 60k/15k/30k RNA-FM/RNABERT inference would take much longer. The real foundation model run uses stratified samples of train 3000, valid 900, and test 900. The CNN full PyTorch run remains in `reports/experiment_1.md` and uses train 60000, valid 15000, test 30000.

## Environment

- Python environment: `/home/hkqxhy/ENTER/bin/python`
- torch: `2.8.0+cpu`; CUDA available: `False`
- transformers: `4.40.2`
- multimolecule: installed/importable
- numpy/scipy/scikit-learn: `1.26.4` / `1.13.1` / `1.6.1`
- Pillow: `11.3.0`

## Validation Summary

| model                       |   train_rows |   valid_rows |   macro_f1 |   accuracy |    auroc |    auprc |
|:----------------------------|-------------:|-------------:|-----------:|-----------:|---------:|---------:|
| RNA-FM frozen encoder + MLP |         3000 |          900 |   0.794977 |   0.796667 | 0.929194 | 0.860421 |

## Test Summary

| model                        |   rows |   accuracy |   macro_f1 |    auroc |    auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |
|:-----------------------------|-------:|-----------:|-----------:|---------:|---------:|-----------:|--------------:|----------------:|
| RNA-FM frozen encoder + MLP  |    900 |   0.792222 |   0.792175 | 0.920944 | 0.852817 |   0.885572 |      0.793277 |        0.697674 |
| RNABERT frozen encoder + MLP |    900 |   0.841111 |   0.838071 | 0.953783 | 0.916495 |   0.922314 |      0.849145 |        0.742754 |

## Outputs

- `results/tables_real_foundation/experiment_1_test_metrics.csv`
- `results/tables_real_foundation/experiment_1_test_confusion_matrices.csv`
- `results/figures_real_foundation/experiment_1_macro_f1.png`
- `results/figures_real_foundation/experiment_1_confusion_matrices.png`
- `results/checkpoints/real_foundation/rnafm.joblib`
- `results/checkpoints/real_foundation/rnabert.joblib`
- `models/hf/rnafm/`
- `models/hf/rnabert/`
