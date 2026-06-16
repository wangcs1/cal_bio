# Experiment 3: Aberrant Splicing Variant Effects

This run scores artificial SNVs with real local models only. It includes the trained project models (CNN, RNA-FM frozen encoder, RNABERT frozen encoder) and external real splice tools (SpliceAI, Pangolin, MMSplice, MaxEntScan).
External tools are not fallback/proxy rows; they are executed through a Python 3.10 splice-tool environment and merged as sequence-level variant scores.

## Variant Set

| variant_type    | label_name      |   rows |
|:----------------|:----------------|-------:|
| acceptor_gain   | splice_altering |     55 |
| acceptor_loss   | splice_altering |     90 |
| donor_gain      | splice_altering |     45 |
| donor_loss      | splice_altering |     90 |
| neutral_far_snv | neutral         |    180 |

## 3A Artificial Variant Ranking

| model                         | source             |   auroc |   auprc |   top_k |   top_k_recall |   enrichment_at_k |   variants |
|:------------------------------|:-------------------|--------:|--------:|--------:|---------------:|------------------:|-----------:|
| CNN baseline (PyTorch Conv1D) | trained_classifier |  0.7137 |  0.8693 |      46 |         0.1643 |            1.6429 |        460 |
| RNABERT frozen encoder + MLP  | trained_classifier |  0.6828 |  0.8503 |      46 |         0.1643 |            1.6429 |        460 |
| SpliceAI real sequence model  | real_external_tool |  0.6392 |  0.8271 |      46 |         0.1643 |            1.6429 |        460 |
| MaxEntScan real local score   | real_external_tool |  0.7268 |  0.789  |      46 |         0.1643 |            1.6429 |        460 |
| RNA-FM frozen encoder + MLP   | trained_classifier |  0.5962 |  0.7695 |      46 |         0.1643 |            1.6429 |        460 |
| MMSplice real sequence model  | real_external_tool |  0.3333 |  0.6087 |      46 |         0      |            0      |        460 |
| Pangolin real sequence model  | real_external_tool |  0.5    |  0.6087 |      46 |         0.0857 |            0.8571 |        460 |

## By Variant Type

| model                         | source             | variant_type    |   mean_score |   median_score |   rows |   positive_rate |
|:------------------------------|:-------------------|:----------------|-------------:|---------------:|-------:|----------------:|
| CNN baseline (PyTorch Conv1D) | trained_classifier | acceptor_gain   |      -0.0004 |        -0.0002 |     55 |               1 |
| CNN baseline (PyTorch Conv1D) | trained_classifier | acceptor_loss   |       0.2436 |         0.2396 |     90 |               1 |
| CNN baseline (PyTorch Conv1D) | trained_classifier | donor_gain      |      -0.0004 |         0      |     45 |               1 |
| CNN baseline (PyTorch Conv1D) | trained_classifier | donor_loss      |       0.1168 |         0.1158 |     90 |               1 |
| CNN baseline (PyTorch Conv1D) | trained_classifier | neutral_far_snv |       0.0008 |         0.0004 |    180 |               0 |
| MMSplice real sequence model  | real_external_tool | acceptor_gain   |       0      |         0      |     55 |               1 |
| MMSplice real sequence model  | real_external_tool | acceptor_loss   |       0      |         0      |     90 |               1 |
| MMSplice real sequence model  | real_external_tool | donor_gain      |       0      |         0      |     45 |               1 |
| MMSplice real sequence model  | real_external_tool | donor_loss      |       0      |         0      |     90 |               1 |
| MMSplice real sequence model  | real_external_tool | neutral_far_snv |       0.0055 |         0      |    180 |               0 |
| MaxEntScan real local score   | real_external_tool | acceptor_gain   |      -0.0121 |         0      |     55 |               1 |
| MaxEntScan real local score   | real_external_tool | acceptor_loss   |       0.2969 |         0      |     90 |               1 |
| MaxEntScan real local score   | real_external_tool | donor_gain      |      -0      |         0      |     45 |               1 |
| MaxEntScan real local score   | real_external_tool | donor_loss      |       0.0026 |         0      |     90 |               1 |
| MaxEntScan real local score   | real_external_tool | neutral_far_snv |       0      |         0      |    180 |               0 |
| Pangolin real sequence model  | real_external_tool | acceptor_gain   |       0      |         0      |     55 |               1 |
| Pangolin real sequence model  | real_external_tool | acceptor_loss   |       0      |         0      |     90 |               1 |
| Pangolin real sequence model  | real_external_tool | donor_gain      |       0      |         0      |     45 |               1 |
| Pangolin real sequence model  | real_external_tool | donor_loss      |       0      |         0      |     90 |               1 |
| Pangolin real sequence model  | real_external_tool | neutral_far_snv |       0      |         0      |    180 |               0 |
| RNA-FM frozen encoder + MLP   | trained_classifier | acceptor_gain   |       0.0031 |         0      |     55 |               1 |
| RNA-FM frozen encoder + MLP   | trained_classifier | acceptor_loss   |       0.0302 |         0.0014 |     90 |               1 |
| RNA-FM frozen encoder + MLP   | trained_classifier | donor_gain      |       0.0066 |         0      |     45 |               1 |
| RNA-FM frozen encoder + MLP   | trained_classifier | donor_loss      |       0.3745 |         0.318  |     90 |               1 |

## Smoke Inputs

ClinVar-format and sQTL-format smoke tables are derived from the artificial variant set to validate input/output plumbing.
They are not full ClinVar or GTEx benchmarks.

ClinVar-format smoke metrics:

| model                         |   auroc |   auprc |   variants |
|:------------------------------|--------:|--------:|-----------:|
| CNN baseline (PyTorch Conv1D) |  0.8611 |  0.8968 |         12 |
| RNABERT frozen encoder + MLP  |  0.3333 |  0.6146 |         12 |
| RNA-FM frozen encoder + MLP   |  0.3889 |  0.5758 |         12 |

sQTL-format rows: 10

Outputs:

- `results/experiment_3/tables/experiment_3A_artificial_variant_scores.csv`
- `results/experiment_3/tables/experiment_3A_artificial_variant_metrics.csv`
- `results/experiment_3/tables/variant_effect_stratified_by_type.csv`
- `results/experiment_3/tables/experiment_3B_clinvar_smoke_metrics.csv`
- `results/experiment_3/tables/experiment_3C_sqtl_case_study.csv`
- `results/experiment_3/figures/exp3_variant_auroc.png`
- `results/experiment_3/figures/exp3_variant_auprc.png`
