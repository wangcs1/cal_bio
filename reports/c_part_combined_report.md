# C Part Combined Report

Generated: 2026-06-17 19:42:14

## Scope

The main C Part pipeline uses the current small split, not a million-scale full dataset.
Experiments 1 and 2 use a real GENCODE/GRCh38 chromosome-holdout split with balanced
donor/acceptor/non-splice classes. This run reports real-model rows only. Experiments 1
and 2 use CNN, RNA-FM frozen encoder, and RNABERT frozen encoder. Experiment 3 additionally
includes real external splice tools: SpliceAI, Pangolin, MMSplice, and MaxEntScan.

## Data

| split   |   rows |
|:--------|-------:|
| train   |   2339 |
| valid   |    230 |
| test    |    431 |

Real ClinVar variant summary:

| variant_type            | label_name      |   rows |
|:------------------------|:----------------|-------:|
| acceptor_clinvar_splice | splice_altering |    140 |
| clinvar_benign_snv      | neutral         |    250 |
| donor_clinvar_splice    | splice_altering |    110 |

## Experiment 1

Splice-site donor/acceptor/non-splice classification with the real-model rows only.

| model_key   | model                         | model_type     | backend                   | split   |   rows |   accuracy |   macro_precision |   macro_recall |   macro_f1 |   auroc |   auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |   hard_negative_fpr |   hard_negative_rows |   hard_negative_count |   hard_negative_false_positives | hard_negative_filter                         |
|:------------|:------------------------------|:---------------|:--------------------------|:--------|-------:|-----------:|------------------:|---------------:|-----------:|--------:|--------:|-----------:|--------------:|----------------:|--------------------:|---------------------:|----------------------:|--------------------------------:|:---------------------------------------------|
| cnn         | CNN baseline (PyTorch Conv1D) | baseline       | pytorch_conv1d            | test    |    431 |     0.8213 |            0.8225 |         0.8284 |     0.8227 |  0.949  |  0.9066 |     0.8989 |        0.835  |          0.7343 |              0.3226 |                  155 |                   155 |                              50 | label == 2 and negative_type contains 'hard' |
| rnafm       | RNA-FM frozen encoder + MLP   | frozen encoder | local_pretrained_required | test    |    431 |     0.7633 |            0.7664 |         0.7677 |     0.767  |  0.9157 |  0.8402 |     0.8168 |        0.811  |          0.6731 |              0.329  |                  155 |                   155 |                              51 | label == 2 and negative_type contains 'hard' |
| rnabert     | RNABERT frozen encoder + MLP  | frozen encoder | local_pretrained_required | test    |    431 |     0.819  |            0.8196 |         0.8261 |     0.8199 |  0.9298 |  0.8718 |     0.8889 |        0.8339 |          0.7368 |              0.3226 |                  155 |                   155 |                              50 | label == 2 and negative_type contains 'hard' |

## Experiment 2

Multi-scale context, hard-negative, and rare-motif stress tests using the same real-model
set.

| model                         |   test_easy_macro_f1 |   test_easy_auprc |   test_hard_macro_f1 |   test_hard_auprc |   cross_gene_macro_f1 |   cross_gene_auprc |   hard_negative_fpr |   hard_negative_rows |   hard_negative_count |   hard_negative_false_positives | hard_negative_filter                         |   test_easy_rows |   test_hard_rows |   cross_gene_rows |
|:------------------------------|---------------------:|------------------:|---------------------:|------------------:|----------------------:|-------------------:|--------------------:|---------------------:|----------------------:|--------------------------------:|:---------------------------------------------|-----------------:|-----------------:|------------------:|
| CNN baseline (PyTorch Conv1D) |               0.6334 |            0.6651 |               0.8193 |            0.9049 |                0.8193 |             0.9049 |              0.3548 |                  155 |                   155 |                              55 | label == 2 and negative_type contains 'hard' |              276 |              431 |               431 |
| RNA-FM frozen encoder + MLP   |               0.6064 |            0.6639 |               0.7869 |            0.8439 |                0.7869 |             0.8439 |              0.3097 |                  155 |                   155 |                              48 | label == 2 and negative_type contains 'hard' |              276 |              431 |               431 |
| RNABERT frozen encoder + MLP  |               0.6264 |            0.6635 |               0.8149 |            0.8733 |                0.8149 |             0.8733 |              0.329  |                  155 |                   155 |                              51 | label == 2 and negative_type contains 'hard' |              276 |              431 |               431 |

## Experiment 3

Real ClinVar variant effect prediction plus small format-control case studies. The
main variant table includes both trained project models and real external splice tools;
all reported rows use real local models or real external tool outputs.

| model                         | source             |   auroc |   auprc |   top_k |   top_k_recall |   enrichment_at_k |   variants |
|:------------------------------|:-------------------|--------:|--------:|--------:|---------------:|------------------:|-----------:|
| RNA-FM frozen encoder + MLP   | trained_classifier |  0.7532 |  0.7839 |      50 |          0.196 |              1.96 |        500 |
| Pangolin real sequence model  | real_external_tool |  0.7023 |  0.7781 |      50 |          0.196 |              1.96 |        500 |
| CNN baseline (PyTorch Conv1D) | trained_classifier |  0.7019 |  0.7404 |      50 |          0.196 |              1.96 |        500 |
| RNABERT frozen encoder + MLP  | trained_classifier |  0.6673 |  0.7326 |      50 |          0.196 |              1.96 |        500 |
| SpliceAI real sequence model  | real_external_tool |  0.6425 |  0.7289 |      50 |          0.2   |              2    |        500 |
| MaxEntScan real local score   | real_external_tool |  0.6165 |  0.6716 |      50 |          0.196 |              1.96 |        500 |
| MMSplice real sequence model  | real_external_tool |  0.5    |  0.5    |      50 |          0.076 |              0.76 |        500 |

## Reports

- `reports/experiment_1.md`
- `reports/experiment_2.md`
- `reports/experiment_3.md`
- `reports/data_qc.md`
- `reports/resource_setup.md`
