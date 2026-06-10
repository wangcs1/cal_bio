# 实验二报告：多尺度上下文与 hard negative

同步日期：2026-06-11

## 实验目标

实验二评估剪接识别是否只依赖局部 `GT/AG` motif，还是需要更长程序列上下文、hard-negative 区分、rare motif、组织使用 case study 和 junction topology。

主实验只使用当前已经准备好的四个小样本窗口：

- `data/shared/processed/splice_sites_pm50.csv`
- `data/shared/processed/splice_sites_pm100.csv`
- `data/shared/processed/splice_sites_pm200.csv`
- `data/shared/processed/splice_sites_pm400.csv`

`pm1000`、5kb、10kb 不作为主实验；Borzoi/AlphaGenome 只作为 optional long-range case study。

## 2A：多尺度上下文

| flank | 模型 | Macro-F1 | AUPRC | Hard-negative FPR |
|--:|:--|--:|--:|--:|
| 50 | CNN motif baseline | 0.8733 | 0.8940 | 0.4179 |
| 50 | RNA-FM frozen k-mer + MLP | 0.9897 | 0.9998 | 0.0448 |
| 50 | RNABERT frozen token + MLP | 1.0000 | 1.0000 | 0.0000 |
| 50 | SpliceAI signal proxy | 0.9965 | 1.0000 | 0.0149 |
| 100 | CNN motif baseline | 0.8733 | 0.8940 | 0.4179 |
| 100 | RNA-FM frozen k-mer + MLP | 0.9761 | 0.9993 | 0.0746 |
| 100 | RNABERT frozen token + MLP | 1.0000 | 1.0000 | 0.0000 |
| 100 | SpliceAI signal proxy | 0.9965 | 1.0000 | 0.0149 |
| 200 | CNN motif baseline | 0.8733 | 0.8940 | 0.4179 |
| 200 | RNA-FM frozen k-mer + MLP | 0.9453 | 0.9876 | 0.1194 |
| 200 | RNABERT frozen token + MLP | 0.9895 | 1.0000 | 0.0448 |
| 200 | SpliceAI signal proxy | 0.9965 | 1.0000 | 0.0149 |
| 400 | CNN motif baseline | 0.8733 | 0.8940 | 0.4179 |
| 400 | RNA-FM frozen k-mer + MLP | 0.9588 | 0.9923 | 0.0896 |
| 400 | RNABERT frozen token + MLP | 0.9825 | 0.9984 | 0.0597 |
| 400 | SpliceAI signal proxy | 0.9965 | 1.0000 | 0.0149 |

## 2B：hard-negative 与 rare motif

Hard-negative 过滤条件为 `label == 2 and negative_type contains 'hard'`，共 67 条。

| 模型 | Hard Macro-F1 | Cross-gene Macro-F1 | Hard-negative FPR | FP / hard rows | 备注 |
|:--|--:|--:|--:|:--|:--|
| CNN motif baseline | 0.8347 | 0.8733 | 0.4179 | 28 / 67 | proxy 主链路 |
| RNA-FM frozen k-mer + MLP | 0.9362 | 0.9453 | 0.1194 | 8 / 67 | proxy/fallback |
| RNABERT frozen token + MLP | 0.9870 | 0.9895 | 0.0448 | 3 / 67 | proxy/fallback |
| SpliceAI signal proxy | 0.9957 | 0.9965 | 0.0149 | 1 / 67 | signal proxy |
| Pangolin optional tool | 0.8874 | 0.8874 | 0.4478 | 30 / 67 | small case-study proxy |

Rare motif 小样本 case study 各 motif 40 条。GC-AG donor 对多数模型较容易；AT-AC acceptor proxy 对 CNN 与 SpliceAI signal proxy 较困难。

| 模型 | motif | rows | mean target prob | Macro-F1 |
|:--|:--|--:|--:|--:|
| CNN motif baseline | rare_AT-AC_acceptor_proxy | 40 | 0.0026 | 0.0000 |
| CNN motif baseline | rare_GC-AG_donor | 40 | 0.9942 | 1.0000 |
| RNA-FM frozen k-mer + MLP | rare_AT-AC_acceptor_proxy | 40 | 0.9926 | 1.0000 |
| RNA-FM frozen k-mer + MLP | rare_GC-AG_donor | 40 | 0.5842 | 0.3750 |
| RNABERT frozen token + MLP | rare_AT-AC_acceptor_proxy | 40 | 0.9963 | 1.0000 |
| RNABERT frozen token + MLP | rare_GC-AG_donor | 40 | 0.9867 | 1.0000 |
| SpliceAI signal proxy | rare_AT-AC_acceptor_proxy | 40 | 0.3938 | 0.0000 |
| SpliceAI signal proxy | rare_GC-AG_donor | 40 | 0.2239 | 0.0000 |

## 2C：组织使用与长程 case study

GTEx-style tissue usage 当前是 200 行 synthetic/proxy case study，不是 GTEx 全量建模。`splice_usage` 范围为 0.2191 到 0.8959，均值 0.4946。

Borzoi/AlphaGenome 只保留两条 documented proxy long-range case：

| case | 模型 | 输入范围 | delta signal | 状态 |
|:--|:--|:--|--:|:--|
| LR_CASE_001 | Borzoi optional long-range case study | synthetic 10kb-centered locus sketch | -0.31 | real model not required |
| LR_CASE_002 | AlphaGenome optional long-range case study | synthetic regulatory-context sketch | 0.22 | API/weights optional |

## 结论

实验二支持“局部 motif 很重要但不充分”的结论。CNN motif baseline 在 hard-negative 上误报较高；RNABERT/RNA-FM proxy 和 SpliceAI signal proxy 更能利用上下文或更严格的 splice signal。Pangolin、GTEx-style tissue usage、Borzoi 和 AlphaGenome 均只作为小样本 case study，不作为大规模真实 benchmark。

## 输出文件

- `results/experiment_2/tables/experiment_2A_multiscale_context.csv`
- `results/experiment_2/tables/experiment_2B_hard_negative.csv`
- `results/experiment_2/tables/experiment_2B_rare_motif.csv`
- `results/experiment_2/tables/experiment_2C_gtex_tissue_usage.csv`
- `results/experiment_2/tables/long_range_regulatory_case_study.csv`
- `results/experiment_2/tables/regulatory_motif_masking.csv`
- `results/experiment_2/tables/junction_topology_ablation.csv`
- `results/experiment_2/figures/exp2A_context_macro_f1.png`
- `results/experiment_2/figures/exp2B_hard_negative_fpr.png`
- `results/experiment_2/figures/exp2C_tissue_splice_usage_heatmap.png`
