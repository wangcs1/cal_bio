# Experiment 2: Context Scale And Hard Negatives

This run uses only the implemented real-model set: CNN baseline, RNA-FM frozen encoder, and RNABERT frozen encoder.
RNA-FM/RNABERT require local pretrained weights under `models/hf/`; the report contains only real-model rows.

## 2A Multi-Scale Context

|   window_flank |   sequence_length | model                         |   accuracy |   macro_f1 |   auroc |   auprc |   hard_negative_fpr |
|---------------:|------------------:|:------------------------------|-----------:|-----------:|--------:|--------:|--------------------:|
|             50 |               101 | CNN baseline (PyTorch Conv1D) |     0.8329 |     0.8336 |  0.9513 |  0.9106 |              0.3419 |
|             50 |               101 | RNA-FM frozen encoder + MLP   |     0.8005 |     0.8024 |  0.9158 |  0.8448 |              0.2774 |
|             50 |               101 | RNABERT frozen encoder + MLP  |     0.8121 |     0.8136 |  0.9285 |  0.8644 |              0.3097 |
|            100 |               201 | CNN baseline (PyTorch Conv1D) |     0.8121 |     0.8109 |  0.9471 |  0.9047 |              0.3871 |
|            100 |               201 | RNA-FM frozen encoder + MLP   |     0.761  |     0.7634 |  0.906  |  0.8126 |              0.3613 |
|            100 |               201 | RNABERT frozen encoder + MLP  |     0.8306 |     0.8323 |  0.9392 |  0.8885 |              0.2903 |
|            200 |               401 | CNN baseline (PyTorch Conv1D) |     0.8167 |     0.8169 |  0.9467 |  0.9045 |              0.3613 |
|            200 |               401 | RNA-FM frozen encoder + MLP   |     0.7842 |     0.7869 |  0.9181 |  0.8439 |              0.3097 |
|            200 |               401 | RNABERT frozen encoder + MLP  |     0.8144 |     0.8149 |  0.9308 |  0.8733 |              0.329  |
|            400 |               801 | CNN baseline (PyTorch Conv1D) |     0.819  |     0.8186 |  0.9466 |  0.9052 |              0.3742 |
|            400 |               801 | RNA-FM frozen encoder + MLP   |     0.7842 |     0.788  |  0.9095 |  0.8333 |              0.2903 |
|            400 |               801 | RNABERT frozen encoder + MLP  |     0.8469 |     0.8477 |  0.9354 |  0.8815 |              0.2645 |

## 2B Hard-Negative Stress Test

| model                         |   test_easy_macro_f1 |   test_hard_macro_f1 |   cross_gene_macro_f1 |   hard_negative_fpr |   hard_negative_false_positives |   hard_negative_rows |
|:------------------------------|---------------------:|---------------------:|----------------------:|--------------------:|--------------------------------:|---------------------:|
| CNN baseline (PyTorch Conv1D) |               0.6334 |               0.8193 |                0.8193 |              0.3548 |                              55 |                  155 |
| RNA-FM frozen encoder + MLP   |               0.6064 |               0.7869 |                0.7869 |              0.3097 |                              48 |                  155 |
| RNABERT frozen encoder + MLP  |               0.6264 |               0.8149 |                0.8149 |              0.329  |                              51 |                  155 |

## Rare Motif Stress Test

The rare-motif table is a synthetic stress test, not a claim about population-scale rare splice-site recall.

| model                         | motif_type          |   rows |   mean_target_probability |   accuracy |   macro_f1 |
|:------------------------------|:--------------------|-------:|--------------------------:|-----------:|-----------:|
| CNN baseline (PyTorch Conv1D) | rare_AT-AC_acceptor |     40 |                    0.1839 |      0.1   |     0.0909 |
| CNN baseline (PyTorch Conv1D) | rare_GC-AG_donor    |     40 |                    0.991  |      1     |     1      |
| RNA-FM frozen encoder + MLP   | rare_AT-AC_acceptor |     40 |                    0.1288 |      0.1   |     0.0606 |
| RNA-FM frozen encoder + MLP   | rare_GC-AG_donor    |     40 |                    0.3003 |      0.25  |     0.2    |
| RNABERT frozen encoder + MLP  | rare_AT-AC_acceptor |     40 |                    0.114  |      0.025 |     0.0163 |
| RNABERT frozen encoder + MLP  | rare_GC-AG_donor    |     40 |                    0.6402 |      0.725 |     0.4203 |

Outputs:

- `results/experiment_2/tables/experiment_2A_multiscale_context.csv`
- `results/experiment_2/tables/experiment_2B_hard_negative.csv`
- `results/experiment_2/tables/experiment_2B_rare_motif.csv`
- `results/experiment_2/figures/exp2A_context_macro_f1.png`
- `results/experiment_2/figures/exp2A_context_auprc.png`
- `results/experiment_2/figures/exp2B_hard_negative_fpr.png`
