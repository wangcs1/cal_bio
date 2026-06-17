# Experiment 3: Aberrant Splicing Variant Effects

Variant table: `data\experiment_3\clinvar_splicing_variants.csv`

This run scores real ClinVar-labeled SNVs with real local models only. It includes the trained project models (CNN, RNA-FM frozen encoder, RNABERT frozen encoder) and external real splice tools (SpliceAI, Pangolin, MMSplice, MaxEntScan).
External tools are executed through a Python 3.10 splice-tool environment and merged as real sequence-level variant scores.

## Variant Set

| variant_type            | label_name      |   rows |
|:------------------------|:----------------|-------:|
| acceptor_clinvar_splice | splice_altering |    141 |
| clinvar_benign_snv      | neutral         |    250 |
| donor_clinvar_splice    | splice_altering |    109 |

## 3A Real ClinVar Variant Ranking

| model                         | source             |   auroc |   auprc |   top_k |   top_k_recall |   enrichment_at_k |   variants |
|:------------------------------|:-------------------|--------:|--------:|--------:|---------------:|------------------:|-----------:|
| RNA-FM frozen encoder + MLP   | trained_classifier |  0.7532 |  0.7839 |      50 |          0.196 |              1.96 |        500 |
| Pangolin real sequence model  | real_external_tool |  0.7023 |  0.7781 |      50 |          0.196 |              1.96 |        500 |
| CNN baseline (PyTorch Conv1D) | trained_classifier |  0.6959 |  0.7374 |      50 |          0.196 |              1.96 |        500 |
| RNABERT frozen encoder + MLP  | trained_classifier |  0.6673 |  0.7326 |      50 |          0.196 |              1.96 |        500 |
| SpliceAI real sequence model  | real_external_tool |  0.6425 |  0.7289 |      50 |          0.2   |              2    |        500 |
| MaxEntScan real local score   | real_external_tool |  0.6165 |  0.6716 |      50 |          0.196 |              1.96 |        500 |
| MMSplice real sequence model  | real_external_tool |  0.5    |  0.5    |      50 |          0.076 |              0.76 |        500 |

## 3B Distance-Matched ClinVar Ranking

This diagnostic evaluates the same scores on an exact-distance-matched ClinVar subset to reduce the `closer to splice site = more likely pathogenic` shortcut.
| model                         | source             |   auroc |   auprc |   top_k |   top_k_recall |   enrichment_at_k |   variants |
|:------------------------------|:-------------------|--------:|--------:|--------:|---------------:|------------------:|-----------:|
| RNA-FM frozen encoder + MLP   | trained_classifier |  0.6684 |  0.6559 |      33 |         0.135  |            1.3333 |        326 |
| Pangolin real sequence model  | real_external_tool |  0.5902 |  0.6416 |      33 |         0.1595 |            1.5758 |        326 |
| SpliceAI real sequence model  | real_external_tool |  0.5135 |  0.5693 |      33 |         0.135  |            1.3333 |        326 |
| CNN baseline (PyTorch Conv1D) | trained_classifier |  0.558  |  0.565  |      33 |         0.1104 |            1.0909 |        326 |
| MaxEntScan real local score   | real_external_tool |  0.4668 |  0.5496 |      33 |         0.1288 |            1.2727 |        326 |
| RNABERT frozen encoder + MLP  | trained_classifier |  0.5155 |  0.5442 |      33 |         0.1104 |            1.0909 |        326 |
| MMSplice real sequence model  | real_external_tool |  0.5    |  0.5    |      33 |         0.092  |            0.9091 |        326 |

## By Variant Type

| model                         | source             | variant_type            |   mean_score |   median_score |   rows |   positive_rate |
|:------------------------------|:-------------------|:------------------------|-------------:|---------------:|-------:|----------------:|
| CNN baseline (PyTorch Conv1D) | trained_classifier | acceptor_clinvar_splice |       0.1354 |         0.0021 |    141 |               1 |
| CNN baseline (PyTorch Conv1D) | trained_classifier | clinvar_benign_snv      |       0.0106 |         0.0004 |    250 |               0 |
| CNN baseline (PyTorch Conv1D) | trained_classifier | donor_clinvar_splice    |       0.0778 |         0.0021 |    109 |               1 |
| MMSplice real sequence model  | real_external_tool | acceptor_clinvar_splice |       0      |         0      |    141 |               1 |
| MMSplice real sequence model  | real_external_tool | clinvar_benign_snv      |       0      |         0      |    250 |               0 |
| MMSplice real sequence model  | real_external_tool | donor_clinvar_splice    |       0      |         0      |    109 |               1 |
| MaxEntScan real local score   | real_external_tool | acceptor_clinvar_splice |       0.1202 |         0      |    141 |               1 |
| MaxEntScan real local score   | real_external_tool | clinvar_benign_snv      |       0.0029 |         0      |    250 |               0 |
| MaxEntScan real local score   | real_external_tool | donor_clinvar_splice    |       0.1315 |         0      |    109 |               1 |
| Pangolin real sequence model  | real_external_tool | acceptor_clinvar_splice |       0.2122 |         0.0615 |    141 |               1 |
| Pangolin real sequence model  | real_external_tool | clinvar_benign_snv      |       0.0305 |         0.0146 |    250 |               0 |
| Pangolin real sequence model  | real_external_tool | donor_clinvar_splice    |       0.2853 |         0.0794 |    109 |               1 |
| RNA-FM frozen encoder + MLP   | trained_classifier | acceptor_clinvar_splice |       0.1592 |         0.0273 |    141 |               1 |
| RNA-FM frozen encoder + MLP   | trained_classifier | clinvar_benign_snv      |       0.0179 |         0.0003 |    250 |               0 |
| RNA-FM frozen encoder + MLP   | trained_classifier | donor_clinvar_splice    |       0.3015 |         0.0084 |    109 |               1 |
| RNABERT frozen encoder + MLP  | trained_classifier | acceptor_clinvar_splice |       0.2295 |         0.0145 |    141 |               1 |
| RNABERT frozen encoder + MLP  | trained_classifier | clinvar_benign_snv      |       0.0178 |         0.0035 |    250 |               0 |
| RNABERT frozen encoder + MLP  | trained_classifier | donor_clinvar_splice    |       0.2598 |         0.0082 |    109 |               1 |
| SpliceAI real sequence model  | real_external_tool | acceptor_clinvar_splice |       0.1885 |         0.0166 |    141 |               1 |
| SpliceAI real sequence model  | real_external_tool | clinvar_benign_snv      |       0.0293 |         0.0081 |    250 |               0 |
| SpliceAI real sequence model  | real_external_tool | donor_clinvar_splice    |       0.2978 |         0.0406 |    109 |               1 |

## Format-Control Inputs

The additional ClinVar subset and sQTL-style table are format-control checks, not main benchmark evidence.

ClinVar subset metrics:

| model                         |   auroc |   auprc |   variants |
|:------------------------------|--------:|--------:|-----------:|
| RNA-FM frozen encoder + MLP   |  0.7222 |  0.7583 |         12 |
| RNABERT frozen encoder + MLP  |  0.4167 |  0.6366 |         12 |
| CNN baseline (PyTorch Conv1D) |  0.4444 |  0.5742 |         12 |

sQTL-style rows: 10

Outputs:

- `results/experiment_3/tables/experiment_3A_variant_scores.csv`
- `results/experiment_3/tables/experiment_3A_variant_metrics.csv`
- `results/experiment_3/tables/variant_effect_stratified_by_type.csv`
- `results/experiment_3/tables/experiment_3B_format_control_metrics.csv`
- `results/experiment_3/tables/experiment_3C_sqtl_case_study.csv`
- `results/experiment_3/figures/exp3_variant_auroc.png`
- `results/experiment_3/figures/exp3_variant_auprc.png`
