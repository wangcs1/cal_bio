# 实验三报告：异常剪接变异效应预测

## 实验目标

本实验将任务从剪接位点识别推进到 SNV 变异效应预测，评估模型是否能够识别 donor loss、acceptor loss、cryptic gain 与 neutral far SNV，并用 delta score、zero-shot distance 和可解释性图件分析变异影响。

## 数据与输出位置

- 人工变异数据：`data/experiment_3/artificial_variant_effect.csv`
- 实验三结果表：`results/experiment_3/tables/`
- 实验三图件：`results/experiment_3/figures/`

当前人工变异数据共 450 条，其中 donor_loss 90、acceptor_loss 90、cryptic_gain 90、neutral_far_snv 180。

## 模型与分数

实验三包含：

- 分类器 delta score：CNN、RNA-FM proxy、RNABERT proxy、SpliceAI signal proxy。
- zero-shot proxy：RNA-FM embedding distance、RNABERT token distance。
- 传统局部强度 proxy：MaxEntScan consensus proxy。
- 可解释性：saturation mutagenesis、in silico mutagenesis、variant delta profile。

## 主要结果

| 模型 | AUROC | AUPRC | Top-k recall | Enrichment@K |
| --- | ---: | ---: | ---: | ---: |
| RNABERT zero-shot token distance | 0.998 | 0.998 | 0.167 | 1.667 |
| RNA-FM zero-shot embedding distance | 0.997 | 0.998 | 0.167 | 1.667 |
| SpliceAI signal proxy | 0.838 | 0.913 | 0.167 | 1.667 |
| MaxEntScan consensus proxy | 0.841 | 0.884 | 0.167 | 1.667 |
| CNN motif baseline | 0.830 | 0.873 | 0.167 | 1.667 |
| RNABERT frozen token + MLP | 0.683 | 0.809 | 0.156 | 1.556 |
| RNA-FM frozen k-mer + MLP | 0.613 | 0.757 | 0.163 | 1.630 |

zero-shot proxy distance 对 donor_loss 和 acceptor_loss 的区分最强；SpliceAI signal proxy 与 MaxEntScan consensus proxy 也能较好地区分核心剪接 motif 破坏与 neutral far SNV。

## 分类型结果

从 `variant_effect_stratified_by_type.csv` 看，donor_loss 和 acceptor_loss 的平均影响分数明显高于 neutral far SNV。例如：

- RNA-FM zero-shot embedding distance：donor_loss 约 1.040，acceptor_loss 约 1.042，neutral_far_snv 约 0.010。
- RNABERT zero-shot token distance：donor_loss 约 1.755，acceptor_loss 约 1.757，neutral_far_snv 约 0.011。
- MaxEntScan consensus proxy：donor_loss / acceptor_loss 有明显正分数，neutral_far_snv 接近 0。

cryptic_gain 相比 donor/acceptor loss 更难，说明构造新剪接位点不仅依赖局部二核苷酸，还需要更完整的上下文建模。

## 可解释性结果

实验三输出了：

- donor / acceptor saturation mutagenesis heatmap。
- donor / acceptor / hard negative in silico mutagenesis heatmap。
- donor_loss 与 cryptic_gain 的 WT vs Mut delta profile。

这些图件用于展示变异是否降低原 splice site score，或在其他位置形成 cryptic splice peak。

## 关键文件

- `data/experiment_3/artificial_variant_effect.csv`
- `results/experiment_3/tables/experiment_3A_artificial_variant_metrics.csv`
- `results/experiment_3/tables/experiment_3A_artificial_variant_scores.csv`
- `results/experiment_3/tables/variant_effect_stratified_by_type.csv`
- `results/experiment_3/figures/exp3_variant_auroc.png`
- `results/experiment_3/figures/exp3_delta_score_boxplot.png`
- `results/experiment_3/figures/variant_delta_profile_donor_loss.png`
- `results/experiment_3/figures/variant_delta_profile_cryptic_gain.png`

## 结论

实验三说明当前 synthetic/proxy 变异效应链路已经可以区分核心剪接位点破坏与远端中性 SNV，并能给出可解释的 delta profile。后续应补齐 ClinVar 真实 benchmark、Pangolin/MMSplice/SpliceAI 真实打分，以及 GTEx sQTL 组织特异性 case study。
