# TODO Completion Audit: Real-Data Benchmark Revision

本轮整改目标：把项目从早期 synthetic/proxy 叙事收敛为真实 GRCh38/GENCODE + ClinVar 小样本 benchmark，并保证实验、图表、报告和论文主张一致。

## 已完成

- [x] 数据主线切换为真实 GRCh38/GENCODE splice-site benchmark。
- [x] donor / acceptor / GT/AG hard-negative 全量采样为 1000 / 1000 / 1000。
- [x] pm50/pm100/pm200/pm400 共用同一批 `sample_id -> split`。
- [x] train/valid/test 使用 chromosome holdout，sample id、chromosome、gene id 互斥。
- [x] ClinVar benchmark 已生成：250 splice-altering vs 250 benign/likely benign near-splice SNV。
- [x] ClinVar 限制在 held-out test chromosomes，并排除实验一/二已采样 genes。
- [x] ClinVar 距离混杂已写入 `reports/data_qc.md`。
- [x] exact-distance-matched ClinVar 子集已导出：`data/experiment_3/clinvar_splicing_variants_distance_matched.csv`。
- [x] distance-only baseline 已写入 QC：full ClinVar AUROC 0.6475 / AUPRC 0.6917，distance-matched AUROC/AUPRC 0.5。
- [x] Paralog/homology clustering 未做，已作为 QC WARN 和论文限制说明。
- [x] 数据版本已写入论文和 README：GRCh38 / GENCODE v49 / Ensembl 115 / ClinVar fileDate 2026-05-30 / seed 42。

## 实验完成情况

- [x] 实验一已在真实 full-data split 上重跑：`python -m src.experiments.exp1.run_classification --full-data`。
- [x] 实验一结果表、混淆矩阵、Macro-F1 图已更新。
- [x] 实验一 hard-negative FPR 分母已统一为 155。
- [x] 实验二已在真实数据上重跑：`python -m src.experiments.exp2.run_multiscale`。
- [x] 实验二多尺度上下文图、AUPRC 图、hard-negative FPR 图已更新。
- [x] 实验二 rare motif 已作为 optional synthetic stress test，以小表写入论文。
- [x] 实验三已使用真实 ClinVar 表：`data/experiment_3/clinvar_splicing_variants.csv`。
- [x] 实验三 full ClinVar 和 exact-distance-matched 指标已生成。
- [x] 实验三已接入真实外部工具输出：SpliceAI、Pangolin、MMSplice、MaxEntScan。
- [x] 实验三图表已更新：AUROC、AUPRC、calibration、variant type summary、真实 ClinVar delta profile。
- [x] 解释性分析已生成真实 ClinVar case：`variant_delta_profile_clinvar_real_case.png`。

## 论文与报告完成情况

- [x] 摘要、引言、RQ 已改为真实 benchmark 主线。
- [x] 数据与方法已写入真实 splice-site 构造、ClinVar 构造、距离混杂、AUPRC baseline。
- [x] 模型边界已改为真实本地权重 frozen encoder 与真实外部工具输出。
- [x] 实验一正文、表 3-1、图 3-1/3-2 已更新。
- [x] 实验二正文、表 4-1、rare motif 小表、图 4-1/4-2/4-3 已更新。
- [x] 实验三正文已从人工 gain/loss 主线改为真实 ClinVar splice-altering vs benign 主线。
- [x] 实验三已加入 full ClinVar 表和 exact-distance-matched 子集表。
- [x] 图 5-4 calibration 已改为真实 ClinVar delta score 分箱响应，不再沿用 zero-shot 旧量纲。
- [x] 讨论与限制已改为真实小样本 benchmark、paralog 未去重、near-splice SNV、距离混杂、单 seed 等边界。
- [x] 结论已改为“能支持/不能支持”的真实 benchmark 版本。
- [x] 附录路径、复现命令、ClinVar Top-k 指标已更新。
- [x] README 已重写为当前真实 benchmark 版。
- [x] Pangolin 引用已存在于 `report_letax/njuthesis-sample.bib`。

## 验证结果

- [x] `python -m src.data.qc_splice_dataset` 已通过：Step 1-6、8-9 PASS，Step 7 paralog/homology clustering 为预期 WARN。
- [x] `python -m compileall src scripts` 已通过。
- [x] `python -m src.reports.write_c_part_report` 已通过。
- [x] 关键词扫描已完成；论文主文中 `synthetic` 仅作为 optional control / 非主结果边界说明出现。
- [x] LaTeX 图片路径已同步到 `report_letax/images/`。

## 未完成但已明确写入限制

- [ ] 多随机种子 mean ± std 未运行。本轮保留单 seed=42，并在论文限制中说明。
- [ ] CNN early stopping / valid 选最优轮未实现。本轮保留最后一轮 checkpoint，并在方法与限制中说明。
- [ ] Paralog/homology clustering 未实现，已作为 QC WARN 和论文限制。
- [ ] `pytest` 未运行成功：当前 Python 环境缺少 `pytest` 模块，不是测试用例失败。
- [ ] LaTeX 编译未完成：当前环境未发现 `latexmk` / `xelatex` / `tectonic` 命令。
- [ ] 实验三全量脚本二次重跑在 120 秒内超时；本轮使用此前已成功生成的真实 ClinVar 结果表和图件，未改动分数。

## 当前主结论

- 真实 GT/AG hard-negative 是有效压力测试，三类模型仍有约三成 hard-negative FPR。
- 多尺度上下文影响模型表现，但趋势非单调。
- Full ClinVar 上 RNA-FM frozen encoder + MLP 排名最高；distance-matched 子集上整体下降，说明近剪接距离是重要混杂因素。
- 当前结果支持真实小样本排序诊断，不支持全基因组泛化、临床级变异解释、组织特异性剪接建模或 paralog 去泄漏后的强结论。
