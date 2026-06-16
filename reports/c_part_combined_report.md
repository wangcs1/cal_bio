# C Part Combined Report

Generated: 2026-06-17 02:17:40

## Scope

The main C Part pipeline uses the current small split, not a million-scale full dataset.
Default split size is train/valid/test = 855/120/285 when the synthetic builder is run with
the repository defaults. This run reports real-model rows only. Experiments 1 and 2 use
CNN, RNA-FM frozen encoder, and RNABERT frozen encoder. Experiment 3 additionally
includes real external splice tools: SpliceAI, Pangolin, MMSplice, and MaxEntScan.

## Data

| split   |   rows |
|:--------|-------:|
| train   |    855 |
| valid   |    120 |
| test    |    285 |

Artificial variant summary:

| variant_type    | label_name      |   rows |
|:----------------|:----------------|-------:|
| acceptor_gain   | splice_altering |     55 |
| acceptor_loss   | splice_altering |     90 |
| donor_gain      | splice_altering |     45 |
| donor_loss      | splice_altering |     90 |
| neutral_far_snv | neutral         |    180 |

## Experiment 1

Splice-site donor/acceptor/non-splice classification with the real-model rows only.

| model_key   | model                         | model_type     | backend                   | split   |   rows |   accuracy |   macro_precision |   macro_recall |   macro_f1 |   auroc |   auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |   hard_negative_fpr |   hard_negative_rows |   hard_negative_count |   hard_negative_false_positives | hard_negative_filter                         |
|:------------|:------------------------------|:---------------|:--------------------------|:--------|-------:|-----------:|------------------:|---------------:|-----------:|--------:|--------:|-----------:|--------------:|----------------:|--------------------:|---------------------:|----------------------:|--------------------------------:|:---------------------------------------------|
| cnn         | CNN baseline (PyTorch Conv1D) | baseline       | pytorch_conv1d            | test    |    285 |     0.8877 |            0.9111 |         0.8933 |     0.8858 |  0.9914 |  0.9825 |     0.9948 |        0.8531 |          0.8095 |              0.4776 |                   67 |                    67 |                              32 | label == 2 and negative_type contains 'hard' |
| rnafm       | RNA-FM frozen encoder + MLP   | frozen encoder | local_pretrained_required | test    |    285 |     0.9509 |            0.9516 |         0.9526 |     0.9517 |  0.9936 |  0.9883 |     0.9381 |        0.989  |          0.9278 |              0.1493 |                   67 |                    67 |                              10 | label == 2 and negative_type contains 'hard' |
| rnabert     | RNABERT frozen encoder + MLP  | frozen encoder | local_pretrained_required | test    |    285 |     0.993  |            0.9932 |         0.993  |     0.9931 |  0.9999 |  0.9998 |     0.9948 |        0.9944 |          0.99   |              0.0149 |                   67 |                    67 |                               1 | label == 2 and negative_type contains 'hard' |

## Experiment 2

Multi-scale context, hard-negative, and rare-motif stress tests using the same real-model
set.

| model                         |   test_easy_macro_f1 |   test_easy_auprc |   test_hard_macro_f1 |   test_hard_auprc |   cross_gene_macro_f1 |   cross_gene_auprc |   hard_negative_fpr |   hard_negative_rows |   hard_negative_count |   hard_negative_false_positives | hard_negative_filter                         |   test_easy_rows |   test_hard_rows |   cross_gene_rows |
|:------------------------------|---------------------:|------------------:|---------------------:|------------------:|----------------------:|-------------------:|--------------------:|---------------------:|----------------------:|--------------------------------:|:---------------------------------------------|-----------------:|-----------------:|------------------:|
| CNN baseline (PyTorch Conv1D) |               1      |            1      |               0.839  |            0.9983 |                0.8821 |             0.9987 |              0.4925 |                   67 |                    67 |                              33 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 |
| RNA-FM frozen encoder + MLP   |               0.9802 |            0.9987 |               0.9659 |            0.9911 |                0.9726 |             0.994  |              0.0746 |                   67 |                    67 |                               5 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 |
| RNABERT frozen encoder + MLP  |               0.9932 |            0.9997 |               0.9871 |            0.9997 |                0.9896 |             0.9998 |              0.0299 |                   67 |                    67 |                               2 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 |

## Experiment 3

Artificial variant effect prediction plus small ClinVar/sQTL-format smoke inputs. The
main artificial-variant table includes both trained project models and real external
splice-tool rows.

| model                         | source             |   auroc |   auprc |   top_k |   top_k_recall |   enrichment_at_k |   variants |
|:------------------------------|:-------------------|--------:|--------:|--------:|---------------:|------------------:|-----------:|
| CNN baseline (PyTorch Conv1D) | trained_classifier |  0.7137 |  0.8693 |      46 |         0.1643 |            1.6429 |        460 |
| RNABERT frozen encoder + MLP  | trained_classifier |  0.6828 |  0.8503 |      46 |         0.1643 |            1.6429 |        460 |
| Pangolin real sequence model  | real_external_tool |  0.6853 |  0.8368 |      46 |         0.1643 |            1.6429 |        460 |
| SpliceAI real sequence model  | real_external_tool |  0.6392 |  0.8271 |      46 |         0.1643 |            1.6429 |        460 |
| MaxEntScan real local score   | real_external_tool |  0.7268 |  0.789  |      46 |         0.1643 |            1.6429 |        460 |
| RNA-FM frozen encoder + MLP   | trained_classifier |  0.5962 |  0.7695 |      46 |         0.1643 |            1.6429 |        460 |
| MMSplice real sequence model  | real_external_tool |  0.6667 |  0.7    |      46 |         0.0964 |            0.9643 |        460 |

## Reports

- `reports/experiment_1.md`
- `reports/experiment_2.md`
- `reports/experiment_3.md`
- `reports/data_qc.md`
- `reports/resource_setup.md`
