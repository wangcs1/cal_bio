# 实验三报告：异常剪接变异效应预测

同步日期：2026-06-11

## 实验目标

实验三评估 SNV 是否导致 donor loss、acceptor loss、donor gain、acceptor gain 或 neutral far SNV，并用 delta score、zero-shot score、Top-k/enrichment、校准曲线和可解释性图件解释预测来源。

主实验使用：

- `data/experiment_3/artificial_variant_effect.csv`
- `data/experiment_3/clinvar_splicing_variants_smoke.csv`
- `data/experiment_3/gtex_sqtl_variants_smoke.csv`

这些数据均来自当前小样本 split 或 smoke/case study，不是全量 ClinVar、GTEx 或百万级变异 benchmark。

## 人工变异集

旧版合并 gain 分类已拆分为 donor/acceptor 两类：

| variant_type | label_name | rows |
|:--|:--|--:|
| acceptor_gain | splice_altering | 55 |
| acceptor_loss | splice_altering | 90 |
| donor_gain | splice_altering | 45 |
| donor_loss | splice_altering | 90 |
| neutral_far_snv | neutral | 180 |

合计 460 条人工变异，其中 splice-altering 280 条、neutral 180 条。

## 3A：人工变异效应预测

| 模型 | AUROC | AUPRC | top-k | Top-k recall | Enrichment@k | variants |
|:--|--:|--:|--:|--:|--:|--:|
| RNABERT zero-shot token distance | 0.9990 | 0.9993 | 46 | 0.1643 | 1.6429 | 460 |
| RNA-FM zero-shot embedding distance | 0.9986 | 0.9991 | 46 | 0.1643 | 1.6429 | 460 |
| RNABERT zero-shot pseudo-likelihood | 0.8875 | 0.9441 | 46 | 0.1643 | 1.6429 | 460 |
| RNA-FM zero-shot pseudo-likelihood | 0.8791 | 0.9394 | 46 | 0.1643 | 1.6429 | 460 |
| SpliceAI signal proxy | 0.8178 | 0.9048 | 46 | 0.1643 | 1.6429 | 460 |
| Pangolin optional tool (small case-study proxy) | 0.8339 | 0.8865 | 46 | 0.1643 | 1.6429 | 460 |
| SpliceAI optional real tool (proxy fallback) | 0.8321 | 0.8814 | 46 | 0.1643 | 1.6429 | 460 |
| MMSplice optional tool (proxy fallback) | 0.8286 | 0.8773 | 46 | 0.1643 | 1.6429 | 460 |
| MaxEntScan optional tool (proxy fallback) | 0.8286 | 0.8773 | 46 | 0.1643 | 1.6429 | 460 |
| CNN motif baseline | 0.8196 | 0.8678 | 46 | 0.1643 | 1.6429 | 460 |
| RNABERT frozen token + MLP | 0.6829 | 0.8183 | 46 | 0.1607 | 1.6071 | 460 |
| RNA-FM frozen k-mer + MLP | 0.5958 | 0.7584 | 46 | 0.1643 | 1.6429 | 460 |

Top-k/enrichment 曲线已生成 48 行，覆盖 12 个模型和 0.02、0.05、0.10、0.20 四个阈值。校准曲线当前选择 RNABERT zero-shot token distance，8 个 bin，每个 bin 57 到 58 条样本。

## 分变异类型结果

| 模型 | variant_type | mean score | median score | rows |
|:--|:--|--:|--:|--:|
| RNABERT zero-shot token distance | acceptor_loss | 1.7565 | 1.7564 | 90 |
| RNABERT zero-shot token distance | donor_loss | 1.7554 | 1.7554 | 90 |
| RNABERT zero-shot token distance | acceptor_gain | 0.2054 | 0.0219 | 55 |
| RNABERT zero-shot token distance | donor_gain | 0.2712 | 0.0205 | 45 |
| RNA-FM zero-shot embedding distance | acceptor_loss | 1.0418 | 1.0416 | 90 |
| RNA-FM zero-shot embedding distance | donor_loss | 1.0399 | 1.0399 | 90 |
| RNA-FM zero-shot embedding distance | acceptor_gain | 0.0320 | 0.0216 | 55 |
| RNA-FM zero-shot embedding distance | donor_gain | 0.0288 | 0.0198 | 45 |
| SpliceAI signal proxy | acceptor_loss | 0.5257 | 0.5261 | 90 |
| SpliceAI signal proxy | donor_loss | 0.6981 | 0.7092 | 90 |
| SpliceAI signal proxy | acceptor_gain | -0.0042 | 0.0000 | 55 |
| SpliceAI signal proxy | donor_gain | 0.0014 | 0.0003 | 45 |

loss 类变异在 zero-shot distance 和 SpliceAI signal proxy 上更容易被拉开；gain 类变异得分更弱，需要结合 case profile 和模型边界解释。

## 3B/3C：ClinVar 与 sQTL smoke

ClinVar smoke 只有 12 条，用于验证流程可运行：

| 模型 | AUROC | AUPRC | variants |
|:--|--:|--:|--:|
| RNABERT zero-shot token distance | 1.0000 | 1.0000 | 12 |
| RNA-FM zero-shot embedding distance | 1.0000 | 1.0000 | 12 |
| SpliceAI signal proxy | 0.9167 | 0.9167 | 12 |
| RNA-FM zero-shot pseudo-likelihood | 0.7222 | 0.8409 | 12 |
| RNABERT zero-shot pseudo-likelihood | 0.6667 | 0.8258 | 12 |
| Pangolin optional tool (small case-study proxy) | 0.7500 | 0.7500 | 12 |

sQTL case study 当前是小样本 tissue-specific 示例，不给出总体 GTEx 性能结论。`experiment_3C_sqtl_case_study.csv` 保留 `variant_id,tissue,target_gene,target_junction,observed_effect_direction,model_delta_score` 字段。

## 可解释性

已生成以下解释性输出：

- donor、acceptor、hard-negative 的 ISM 矩阵和 heatmap。
- RNA-FM/RNABERT attention proxy heatmap。
- donor_loss、donor_gain、acceptor_gain 的 variant delta profile。
- ClinVar smoke case 的 variant delta profile。

当前报告只使用拆分后的 donor_gain 和 acceptor_gain，不再使用旧的合并 gain 分类作为结果口径。

## 结论

实验三已经形成小样本闭环：人工变异构造、zero-shot/MLP/proxy 工具评分、Top-k/enrichment、校准、ClinVar smoke、sQTL case study 和解释性图件均已生成。当前结论限定在小样本 synthetic/proxy 与 smoke case study 范围内，不声称完成全量 ClinVar、GTEx 或真实外部工具大规模 benchmark。

## 输出文件

- `results/experiment_3/tables/experiment_3A_artificial_variant_scores.csv`
- `results/experiment_3/tables/experiment_3A_artificial_variant_metrics.csv`
- `results/experiment_3/tables/experiment_3A_topk_enrichment_curve.csv`
- `results/experiment_3/tables/variant_effect_stratified_by_type.csv`
- `results/experiment_3/figures/variant_effect_stratified_by_type.png`
- `results/experiment_3/figures/exp3_calibration_curve.png`
- `results/experiment_3/tables/experiment_3B_clinvar_smoke_metrics.csv`
- `results/experiment_3/tables/experiment_3C_sqtl_case_study.csv`
- `results/experiment_3/tables/variant_delta_profile_clinvar_smoke_case.csv`
