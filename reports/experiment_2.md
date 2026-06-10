# 实验二报告：多尺度上下文与 hard negative

## 目标

回答剪接识别是否只依赖局部 `GT/AG` motif，还是需要更长程 intron/exon context、调控 motif、组织 case study 和 junction topology。

## 数据

默认只读取当前四个小样本窗口文件：

- `data/shared/processed/splice_sites_pm50.csv`
- `data/shared/processed/splice_sites_pm100.csv`
- `data/shared/processed/splice_sites_pm200.csv`
- `data/shared/processed/splice_sites_pm400.csv`

不把 `pm1000`、5kb、10kb 作为主实验。长程模型只保留为 optional case study。

## 已完成分析

- 2A：`pm50/pm100/pm200/pm400` context ablation。
- 2B：GT/AG hard-negative benchmark。
- 2B 扩展：rare motif small case study，例如 GC-AG donor proxy。
- 2C：Pangolin/GTEx-style tissue usage small case study。
- Full-context：local motif sufficiency、regulatory motif masking、tissue program、junction topology。
- Long-range：Borzoi/AlphaGenome optional case-study table。

## 结果边界

Pangolin、GTEx、Borzoi、AlphaGenome 均标注为 small case study 或 proxy fallback。它们用于说明可扩展方向，不作为大规模真实生物结论。

## 关键输出

- `results/experiment_2/tables/experiment_2A_multiscale_context.csv`
- `results/experiment_2/tables/experiment_2B_hard_negative.csv`
- `results/experiment_2/tables/experiment_2B_rare_motif.csv`
- `results/experiment_2/tables/experiment_2C_gtex_tissue_usage.csv`
- `results/experiment_2/tables/long_range_regulatory_case_study.csv`
- `results/experiment_2/tables/regulatory_motif_masking.csv`
- `results/experiment_2/tables/junction_topology_ablation.csv`

## 结论

实验二的证据链显示：`GT/AG` 是重要局部线索，但不是充分条件。hard-negative、rare motif、regulatory motif masking、tissue program 和 junction topology 都说明模型需要多尺度上下文。
