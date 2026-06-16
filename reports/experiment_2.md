# Experiment 2: Context Scale And Hard Negatives

This run uses only the implemented real-model set: CNN baseline, RNA-FM frozen encoder, and RNABERT frozen encoder.
RNA-FM/RNABERT require local pretrained weights under `models/hf/`; the report contains only real-model rows.

## 2A Multi-Scale Context

|   window_flank |   sequence_length | model                         |   accuracy |   macro_f1 |   auroc |   auprc |   hard_negative_fpr |
|---------------:|------------------:|:------------------------------|-----------:|-----------:|--------:|--------:|--------------------:|
|             50 |               101 | CNN baseline (PyTorch Conv1D) |     0.9719 |     0.9719 |  1      |  1      |              0.1194 |
|             50 |               101 | RNA-FM frozen encoder + MLP   |     0.986  |     0.9863 |  0.9985 |  0.9973 |              0.0299 |
|             50 |               101 | RNABERT frozen encoder + MLP  |     0.9965 |     0.9966 |  1      |  1      |              0      |
|            100 |               201 | CNN baseline (PyTorch Conv1D) |     0.9263 |     0.9258 |  1      |  1      |              0.3134 |
|            100 |               201 | RNA-FM frozen encoder + MLP   |     0.9649 |     0.9655 |  0.9974 |  0.9952 |              0.0746 |
|            100 |               201 | RNABERT frozen encoder + MLP  |     0.9965 |     0.9966 |  0.9999 |  0.9999 |              0.0149 |
|            200 |               401 | CNN baseline (PyTorch Conv1D) |     0.8842 |     0.8821 |  0.9994 |  0.9988 |              0.4925 |
|            200 |               401 | RNA-FM frozen encoder + MLP   |     0.9719 |     0.9726 |  0.9968 |  0.994  |              0.0746 |
|            200 |               401 | RNABERT frozen encoder + MLP  |     0.9895 |     0.9896 |  0.9999 |  0.9998 |              0.0299 |
|            400 |               801 | CNN baseline (PyTorch Conv1D) |     0.8842 |     0.8821 |  0.9568 |  0.9095 |              0.4925 |
|            400 |               801 | RNA-FM frozen encoder + MLP   |     0.9579 |     0.9586 |  0.9953 |  0.9913 |              0.0746 |
|            400 |               801 | RNABERT frozen encoder + MLP  |     1      |     1      |  1      |  1      |              0      |

## 2B Hard-Negative Stress Test

| model                         |   test_easy_macro_f1 |   test_hard_macro_f1 |   cross_gene_macro_f1 |   hard_negative_fpr |   hard_negative_false_positives |   hard_negative_rows |
|:------------------------------|---------------------:|---------------------:|----------------------:|--------------------:|--------------------------------:|---------------------:|
| CNN baseline (PyTorch Conv1D) |               1      |               0.839  |                0.8821 |              0.4925 |                              33 |                   67 |
| RNA-FM frozen encoder + MLP   |               0.9802 |               0.9659 |                0.9726 |              0.0746 |                               5 |                   67 |
| RNABERT frozen encoder + MLP  |               0.9932 |               0.9871 |                0.9896 |              0.0299 |                               2 |                   67 |

## Rare Motif Stress Test

The rare-motif table is a synthetic stress test, not a claim about population-scale rare splice-site recall.

| model                         | motif_type          |   rows |   mean_target_probability |   accuracy |   macro_f1 |
|:------------------------------|:--------------------|-------:|--------------------------:|-----------:|-----------:|
| CNN baseline (PyTorch Conv1D) | rare_AT-AC_acceptor |     40 |                    0.3587 |      0.05  |     0.0476 |
| CNN baseline (PyTorch Conv1D) | rare_GC-AG_donor    |     40 |                    0.9153 |      1     |     1      |
| RNA-FM frozen encoder + MLP   | rare_AT-AC_acceptor |     40 |                    0.9711 |      1     |     1      |
| RNA-FM frozen encoder + MLP   | rare_GC-AG_donor    |     40 |                    0.7211 |      0.775 |     0.4366 |
| RNABERT frozen encoder + MLP  | rare_AT-AC_acceptor |     40 |                    0.9359 |      0.975 |     0.4937 |
| RNABERT frozen encoder + MLP  | rare_GC-AG_donor    |     40 |                    0.5851 |      0.55  |     0.3548 |

Outputs:

- `results/experiment_2/tables/experiment_2A_multiscale_context.csv`
- `results/experiment_2/tables/experiment_2B_hard_negative.csv`
- `results/experiment_2/tables/experiment_2B_rare_motif.csv`
- `results/experiment_2/figures/exp2A_context_macro_f1.png`
- `results/experiment_2/figures/exp2A_context_auprc.png`
- `results/experiment_2/figures/exp2B_hard_negative_fpr.png`
