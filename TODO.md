# TODO Completion Audit: Real-Data Benchmark Revision

本轮整改目标：把 TODO 中此前未完成的验证、训练策略和审计项全部落地，并同步实验结果、图表和论文文字。

## 已完成

- [x] 数据主线已切换为真实 GRCh38/GENCODE v49 splice-site benchmark。
- [x] donor / acceptor / GT/AG hard-negative 全量采样为 1000 / 1000 / 1000。
- [x] pm50/pm100/pm200/pm400 共用同一批 sample_id -> split。
- [x] train/valid/test 使用 chromosome holdout，sample id、chromosome、gene id 互斥。
- [x] ClinVar benchmark 已生成：250 splice-altering vs 250 benign/likely benign near-splice SNV。
- [x] ClinVar 限制在 held-out test chromosomes，并排除实验一/二已采样 genes。
- [x] ClinVar 距离混杂已写入 reports/data_qc.md。
- [x] exact-distance-matched ClinVar 子集已导出：data/experiment_3/clinvar_splicing_variants_distance_matched.csv。
- [x] distance-only baseline 已写入 QC：full ClinVar AUROC 0.6475 / AUPRC 0.6917，distance-matched AUROC/AUPRC 0.5。
- [x] Paralog/homology 从“未做”升级为可执行近重复泄漏审计：跨 split full-window duplicate、center 161 bp duplicate、high 9-mer Jaccard near-duplicate 均通过；alignment-level paralog clustering 仍作为限制说明。
- [x] 数据版本已写入论文和 README：GRCh38 / GENCODE v49 / Ensembl 115 / ClinVar fileDate 2026-05-30 / seed 42。

## 实验完成情况

- [x] 实验一已在真实 full-data split 上重跑：python -m src.experiments.exp1.run_classification --full-data。
- [x] CNN baseline 已实现 validation Macro-F1 best epoch / early stopping 兼容接口；seed=42 最佳 epoch 为第 5 轮。
- [x] 实验一结果表、混淆矩阵、Macro-F1 图已更新并同步到 report_letax/images/。
- [x] 实验一 hard-negative FPR 分母统一为 155。
- [x] 实验一多随机种子已完成：python scripts/run_exp1_multiseed.py --full-data。
- [x] 多种子汇总已写入 results/experiment_1/tables/experiment_1_multiseed_summary.csv；CNN Macro-F1 = 0.8273 ± 0.0124，RNA-FM/RNABERT 在当前确定性 frozen-encoder 设置下 std 为 0。
- [x] 实验二已在真实数据上重跑：python -m src.experiments.exp2.run_multiscale。
- [x] 实验二多尺度上下文图、AUPRC 图、hard-negative FPR 图已更新。
- [x] 实验二 rare motif 已作为补充 synthetic stress test，以小表写入论文。
- [x] 实验三已使用真实 ClinVar 表：data/experiment_3/clinvar_splicing_variants.csv。
- [x] 实验三 full ClinVar 和 exact-distance-matched 指标已生成。
- [x] 实验三已接入真实外部工具输出：SpliceAI、Pangolin、MMSplice、MaxEntScan。
- [x] 实验三全量脚本已用长超时重跑成功：python -m src.experiments.exp3.run_variant_effect。
- [x] 实验三图表已更新：AUROC、AUPRC、calibration、variant type summary、真实 ClinVar delta profile。
- [x] 解释性分析已生成真实 ClinVar case：variant_delta_profile_clinvar_real_case.png。

## 论文与报告完成情况

- [x] 摘要、引言、RQ 已改为真实 benchmark 主线。
- [x] 数据与方法已写入真实 splice-site 构造、ClinVar 构造、距离混杂、AUPRC baseline。
- [x] 模型边界已改为真实本地权重 frozen encoder 与真实外部工具输出。
- [x] CNN early stopping / valid best epoch 已写入方法部分。
- [x] 实验一正文、表 3-1、图 3-1/3-2 已更新。
- [x] 实验二正文、表 4-1、rare motif 小表、图 4-1/4-2/4-3 已更新。
- [x] 实验三正文已从人工 gain/loss 主线改为真实 ClinVar splice-altering vs benign 主线。
- [x] 实验三已加入 full ClinVar 表和 exact-distance-matched 子集表。
- [x] 图 5-4 calibration 已改为真实 ClinVar delta score 分箱响应，不再沿用 zero-shot 旧量纲。
- [x] 讨论与限制已改为真实小样本 benchmark、near-splice SNV、距离混杂、多 seed 补充、近重复泄漏审计和 alignment-level paralog clustering 边界。
- [x] 结论已改成“能支持/不能支持”的真实 benchmark 版本。
- [x] 附录路径、复现命令、ClinVar Top-k 指标已更新。
- [x] README 已重写为当前真实 benchmark 版。
- [x] 主流程审计已写入 reports/main_pipeline_audit.md：明确主实验 fail-fast，无 k-mer fallback、proxy 主结果或 placeholder 指标。
- [x] Pangolin 引用已存在于 report_letax/njuthesis-sample.bib。
- [x] repo_exp1/repo_exp2/repo_exp3/repo_qc/repo_model_cards 本地产物引用已补入 bib。
- [x] report_letax/njuthesis-sample.tex 已按用户要求恢复 \codepath{...} 路径样式。

## 验证结果

- [x] python -m src.data.qc_splice_dataset 已通过：Step 1-9 全部 PASS；Step 7 Paralog leakage 为可执行近重复审计 PASS。
- [x] python -m compileall src scripts 已通过。
- [x] python -m pytest -q 已通过：8 passed。
- [x] python -m src.reports.write_c_part_report 已通过。
- [x] python scripts/check_latex_report.py 已通过：检查 LaTeX environment、图片文件、citation key、残留 codepath usage。
- [x] 当前环境仍未发现 latexmk/xelatex/tectonic；已尝试安装 tectonic，但提权审批未及时返回。因此本轮完成的是 LaTeX sanity check，不是 TeX 引擎生成 PDF。
- [x] 关键关键词扫描已完成：论文主文中 synthetic 仅作为补充 control / 非主结果边界说明出现。
- [x] LaTeX 图片路径已同步到 report_letax/images/。

## 当前主结论

- 真实 GT/AG hard-negative 是有效压力测试，三类模型仍有约三成 hard-negative FPR。
- 多尺度上下文影响模型表现，但趋势非单调。
- Full ClinVar 中 RNA-FM frozen encoder + MLP 排名最高；distance-matched 子集上整体下降，说明近剪接距离是重要混杂因素。
- 当前结果支持真实小样本排序诊断，不支持全基因组泛化、临床级变异解释、组织特异性剪接建模或 alignment-level paralog/homology 去重后的强结论。
