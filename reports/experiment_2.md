# 实验二报告：多尺度上下文与 hard negative 剪接建模

## 实验目标

本实验回答一个核心问题：模型是否只是识别局部 `GT/AG` motif，还是能够利用更长程的 intron/exon context、调控 motif、组织程序和 splice junction topology。

## 数据与输出位置

- 共用剪接位点数据：`data/shared/processed/`
- 共用划分文件：`data/shared/splits/`
- 实验二结果表：`results/experiment_2/tables/`
- 实验二图件：`results/experiment_2/figures/`

## 已完成分析

- 实验 2A：`±50 / ±100 / ±200 / ±400` 上下文长度消融。
- 实验 2B：GT/AG hard negative benchmark。
- 局部 motif sufficiency：motif-only、MaxEntScan local score、短上下文 one-hot。
- regulatory motif masking 与 motif rescue proxy。
- tissue-specific splice usage proxy case study。
- junction topology ablation 与 case study。

## 主要结果

### 多尺度上下文

在 `±50` 窗口下，代表性结果如下：

| 模型 | Macro-F1 | AUPRC | Hard-negative FPR |
| --- | ---: | ---: | ---: |
| CNN motif baseline | 0.873 | 0.894 | 0.418 |
| RNA-FM frozen k-mer + MLP | 0.990 | 1.000 | 0.045 |
| RNABERT frozen token + MLP | 1.000 | 1.000 | 0.000 |
| SpliceAI signal proxy | 0.996 | 1.000 | 0.015 |

结果显示，单纯局部 CNN/motif 特征在 hard negative 上误报较高；加入 k-mer、位置 token、剪接信号或任务专用 proxy 后，hard-negative FPR 明显下降。

### GT/AG hard negative

| 模型 | Easy Macro-F1 | Hard Macro-F1 | Cross-gene Macro-F1 | Hard-negative FPR |
| --- | ---: | ---: | ---: | ---: |
| CNN motif baseline | 0.948 | 0.835 | 0.873 | 0.418 |
| RNA-FM frozen k-mer + MLP | 0.948 | 0.936 | 0.945 | 0.119 |
| RNABERT frozen token + MLP | 1.000 | 0.987 | 0.990 | 0.045 |
| SpliceAI signal proxy | 1.000 | 0.996 | 0.996 | 0.015 |

hard negative 的差异说明 `GT/AG` 不是充分条件，模型需要利用 motif 周围上下文才能抑制非功能性 GT/AG 位点。

### 局部 motif sufficiency

motif-only GT/AG rule 的 hard-negative FPR 为 `1.000`，说明它会把所有含局部 canonical motif 的 hard negative 误判为剪接位点。MaxEntScan local score 与短上下文 one-hot 能降低误报，但仍明显弱于上下文/信号更完整的模型。

### 组织与 junction topology

tissue-specific usage proxy 和 junction topology case study 显示，组织程序、competing junctions、donor/acceptor degree、exon skipping 等信息可以解释单纯位点强度无法解释的 splice usage 差异。

## 关键文件

- `results/experiment_2/tables/experiment_2A_multiscale_context.csv`
- `results/experiment_2/tables/experiment_2B_hard_negative.csv`
- `results/experiment_2/tables/local_motif_sufficiency.csv`
- `results/experiment_2/tables/regulatory_motif_masking.csv`
- `results/experiment_2/tables/tissue_specific_usage_ablation.csv`
- `results/experiment_2/tables/junction_topology_ablation.csv`
- `results/experiment_2/figures/exp2A_context_macro_f1.png`
- `results/experiment_2/figures/exp2B_hard_negative_fpr.png`

## 结论

实验二支持 README 中的核心判断：剪接识别不能只靠中心 `GT/AG` motif。更长程上下文、调控 motif、组织程序与 junction topology 都会影响模型对剪接位点和 splice usage 的判断。
