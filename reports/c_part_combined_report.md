# C Part Combined Report

Generated: 2026-06-10 21:55:32

## Scope

The main C Part pipeline uses the current small split, not a million-scale full dataset.
Default split size is train/valid/test = 855/120/285 when the synthetic builder is run with
the repository defaults. Large files under `data/raw/` and pretrained weights under
`models/hf/` are optional resources for smoke tests or case studies.

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

Splice-site donor/acceptor/non-splice classification. The table includes proxy/fallback rows
and optional-real-tool rows where available.

| model_key   | model                                        | split   |   rows |   accuracy |   macro_precision |   macro_recall |   macro_f1 |   auroc |   auprc |   donor_f1 |   acceptor_f1 |   non_splice_f1 |   hard_negative_fpr |   hard_negative_rows |   hard_negative_false_positives | hard_negative_filter                         |
|:------------|:---------------------------------------------|:--------|-------:|-----------:|------------------:|---------------:|-----------:|--------:|--------:|-----------:|--------------:|----------------:|--------------------:|---------------------:|--------------------------------:|:---------------------------------------------|
| cnn         | CNN baseline (PyTorch Conv1D)                | test    |    285 |     0.8877 |            0.9111 |         0.8933 |     0.8858 |  0.9908 |  0.9814 |     0.9948 |        0.8531 |          0.8095 |              0.4776 |                   67 |                              32 | label == 2 and negative_type contains 'hard' |
| rnafm       | RNA-FM frozen encoder + MLP                  | test    |    285 |     0.9614 |            0.9625 |         0.9623 |     0.9624 |  0.9972 |  0.9948 |     0.9418 |        1      |          0.9453 |              0.0746 |                   67 |                               5 | label == 2 and negative_type contains 'hard' |
| rnabert     | RNABERT frozen encoder + MLP                 | test    |    285 |     0.9684 |            0.9688 |         0.9695 |     0.969  |  0.9983 |  0.9968 |     0.9583 |        0.9945 |          0.9543 |              0.0896 |                   67 |                               6 | label == 2 and negative_type contains 'hard' |
| spliceai    | SpliceAI optional real tool (proxy fallback) | test    |    285 |     0.8737 |            0.892  |         0.88   |     0.868  |  0.9981 |  0.9962 |     0.9005 |        0.9231 |          0.7805 |              0.5224 |                   67 |                              35 | label == 2 and negative_type contains 'hard' |

## Experiment 2

Multi-scale context, hard-negative, rare motif, tissue usage, regulatory motif, and topology
analyses. Pangolin/GTEx/Borzoi/AlphaGenome outputs are small case studies, not full benchmarks.

| model                                           |   test_easy_macro_f1 |   test_easy_auprc |   test_hard_macro_f1 |   test_hard_auprc |   cross_gene_macro_f1 |   cross_gene_auprc |   hard_negative_fpr |   hard_negative_rows |   hard_negative_false_positives | hard_negative_filter                         |   test_easy_rows |   test_hard_rows |   cross_gene_rows | case_study_note                                             |
|:------------------------------------------------|---------------------:|------------------:|---------------------:|------------------:|----------------------:|-------------------:|--------------------:|---------------------:|--------------------------------:|:---------------------------------------------|-----------------:|-----------------:|------------------:|:------------------------------------------------------------|
| CNN motif baseline                              |               0.9485 |            1      |               0.8347 |            0.8616 |                0.8733 |             0.894  |              0.4179 |                   67 |                              28 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 | nan                                                         |
| RNA-FM frozen k-mer + MLP                       |               0.9485 |            0.9911 |               0.9362 |            0.9847 |                0.9453 |             0.9876 |              0.1194 |                   67 |                               8 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 | nan                                                         |
| RNABERT frozen token + MLP                      |               1      |            1      |               0.987  |            1      |                0.9895 |             1      |              0.0448 |                   67 |                               3 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 | nan                                                         |
| SpliceAI signal proxy                           |               1      |            1      |               0.9957 |            1      |                0.9965 |             1      |              0.0149 |                   67 |                               1 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 | nan                                                         |
| Pangolin optional tool (small case-study proxy) |                      |                   |               0.8874 |            0.9969 |                0.8874 |             0.9969 |              0.4478 |                   67 |                              30 | label == 2 and negative_type contains 'hard' |              218 |              252 |               285 | Pangolin proxy fallback; optional real tool case study only |

## Experiment 3

Artificial variant effect prediction plus ClinVar/sQTL smoke case studies.

| model                                           |   auroc |   auprc |   top_k |   top_k_recall |   enrichment_at_k |   variants |
|:------------------------------------------------|--------:|--------:|--------:|---------------:|------------------:|-----------:|
| RNABERT zero-shot token distance                |  0.999  |  0.9993 |      46 |         0.1643 |            1.6429 |        460 |
| RNA-FM zero-shot embedding distance             |  0.9986 |  0.9991 |      46 |         0.1643 |            1.6429 |        460 |
| RNABERT zero-shot pseudo-likelihood             |  0.8873 |  0.944  |      46 |         0.1643 |            1.6429 |        460 |
| RNA-FM zero-shot pseudo-likelihood              |  0.8791 |  0.9394 |      46 |         0.1643 |            1.6429 |        460 |
| SpliceAI signal proxy                           |  0.8178 |  0.9048 |      46 |         0.1643 |            1.6429 |        460 |
| Pangolin optional tool (small case-study proxy) |  0.8339 |  0.8865 |      46 |         0.1643 |            1.6429 |        460 |
| SpliceAI optional real tool (proxy fallback)    |  0.8321 |  0.8814 |      46 |         0.1643 |            1.6429 |        460 |
| MMSplice optional tool (proxy fallback)         |  0.8286 |  0.8773 |      46 |         0.1643 |            1.6429 |        460 |
| MaxEntScan optional tool (proxy fallback)       |  0.8286 |  0.8773 |      46 |         0.1643 |            1.6429 |        460 |
| CNN motif baseline                              |  0.8196 |  0.8678 |      46 |         0.1643 |            1.6429 |        460 |
| RNABERT frozen token + MLP                      |  0.6829 |  0.8183 |      46 |         0.1607 |            1.6071 |        460 |
| RNA-FM frozen k-mer + MLP                       |  0.5958 |  0.7584 |      46 |         0.1643 |            1.6429 |        460 |

## Reports

- `reports/experiment_1.md`
- `reports/experiment_2.md`
- `reports/experiment_3.md`
- `reports/data_qc.md`
- `reports/resource_setup.md`
- `reports/model_cards.md`
- `reports/experiment_log.md`
