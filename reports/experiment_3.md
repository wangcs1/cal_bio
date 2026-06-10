# 实验三报告：异常剪接变异效应预测

## 目标

预测 SNV 是否导致 donor loss、acceptor loss、donor gain、acceptor gain，或属于 neutral far SNV，并用 delta score、zero-shot score 与可解释性图件解释影响区域。

## 数据

主实验使用 `data/experiment_3/artificial_variant_effect.csv`。该表由当前小样本 split 构造，不扩展为大规模变异集。ClinVar 与 GTEx/sQTL 只做 smoke/case study：

- `data/experiment_3/clinvar_splicing_variants_smoke.csv`
- `data/experiment_3/gtex_sqtl_variants_smoke.csv`

## 模型与分数

- 分类器 delta score：CNN、RNA-FM proxy/fallback、RNABERT proxy/fallback、SpliceAI signal proxy。
- Zero-shot：RNA-FM/RNABERT embedding distance 与 pseudo-likelihood proxy。
- Optional real-tool wrapper：SpliceAI、Pangolin、MMSplice、MaxEntScan 的小样本 proxy fallback。
- 可解释性：ISM、attention proxy、variant delta profile。

所有分数表均包含 `ref_score`、`alt_score`、`delta_score` 与 `impact_score`，缺少真实外部模型时不会阻塞主流程。

## 关键输出

- `results/experiment_3/tables/experiment_3A_artificial_variant_scores.csv`
- `results/experiment_3/tables/experiment_3A_artificial_variant_metrics.csv`
- `results/experiment_3/tables/experiment_3A_topk_enrichment_curve.csv`
- `results/experiment_3/figures/exp3_calibration_curve.png`
- `results/experiment_3/tables/experiment_3B_clinvar_smoke_metrics.csv`
- `results/experiment_3/tables/experiment_3C_sqtl_case_study.csv`
- `results/experiment_3/tables/variant_delta_profile_clinvar_smoke_case.csv`

## 结论

实验三已经从人工变异、小样本 ClinVar smoke、sQTL case study 和解释性 delta profile 四个角度形成闭环。当前结论限定在小样本 synthetic/proxy 与 smoke case study 范围内，不声称完成全量 ClinVar 或 GTEx benchmark。
