# Small-Sample Splice Dataset QC

QC target: current small split, not a full genome-scale dataset.

| split   |   rows |   donor |   acceptor |   non_splice |   mean_length |   n_fraction |   hard_negative_rows | center_motif_examples                                                           |
|:--------|-------:|--------:|-----------:|-------------:|--------------:|-------------:|---------------------:|:--------------------------------------------------------------------------------|
| train   |    855 |     276 |        288 |          291 |           401 |            0 |                  195 | donor+1:GT; acceptor-2:AG; donor+1:GT; acceptor-2:CA; donor+1:GT; acceptor-2:CA |
| valid   |    120 |      49 |         42 |           29 |           401 |            0 |                   18 | donor+1:AT; acceptor-2:AG; donor+1:GC; acceptor-2:AG; donor+1:GT; acceptor-2:CA |
| test    |    285 |      95 |         90 |          100 |           401 |            0 |                   67 | donor+1:GT; acceptor-2:CA; donor+1:GT; acceptor-2:CA; donor+1:AC; acceptor-2:AG |

Gene leakage check against train split:

- valid: 0 overlapping genes
- test: 0 overlapping genes

Conclusion: the split keeps required columns, balanced classes, fixed sequence length, hard-negative labels, and no train/test gene leakage in the synthetic benchmark.
