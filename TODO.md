# TODO: Real-Data Benchmark Revision Plan

项目定位已从“合成/proxy 沙盒”切换为“一个小而真实、边界清楚的 GENCODE/GRCh38 + ClinVar benchmark”。后续所有代码、结果和正文都要服务这个新故事：真实剪接位点、真实 GT/AG 诱饵 hard negative、真实 ClinVar 二分类、真实外部工具对比；限制也相应改成小样本、未去 paralog、近位点 ClinVar、距离混杂。

## Critical Path

1. 数据冻结与 QC 收尾。
2. 在真实数据上重跑实验一/二/三，生成新表和新图。
3. 按新数字重写正文。
4. 全文一致性检查。

不要先改正文再重跑结果；所有数字和图都会变。

---

## Step 0. 数据收尾与冻结

- [x] 真实 splice-site 数据已生成：GENCODE/GRCh38，donor / acceptor / hard-negative = 1000 / 1000 / 1000。
- [x] 四个窗口 `pm50/100/200/400` 共用同一批 `sample_id -> split`。
- [x] `train/valid/test` 按染色体 holdout，sample/gene/chromosome 互斥。
- [x] ClinVar 真实变异表已生成：250 pathogenic/likely pathogenic splice-altering vs 250 benign/likely benign。
- [x] ClinVar 已限制在 held-out test chromosomes，并排除实验一/二已采样 genes。
- [x] ClinVar 距离混杂已写入 QC：致病变异更贴近剪接位点。
- [x] 已导出 exact-distance-matched ClinVar 子集：`data/experiment_3/clinvar_splicing_variants_distance_matched.csv`。
- [x] 距离-only baseline 已写入 `reports/data_qc.md`。
- [x] Paralog/homology clustering 未做，已作为 `WARN` 写入 QC 限制。
- [ ] 冻结数据版本记录到正文/附录：
  - GRCh38 / GENCODE v50 / Ensembl 116。
  - ClinVar fileDate: 2026-06-15。
  - seed: 42。
  - 构建命令：
    - `python -m src.data.build_splice_site_dataset --max-per-class 1000 --windows 50 100 200 400`
    - `python -m src.data.build_clinvar_variant_dataset`
    - `python -m src.data.qc_splice_dataset`

验收：`reports/data_qc.md` 中 Step 1-6、8-9 为 PASS，Step 7 为已解释的 WARN。

---

## Step 1. 在真实数据上重跑实验

### Step 1.1 实验一：剪接位点三分类

- [ ] 运行真实数据上的实验一：

```bash
python -m src.experiments.exp1.run_classification --full-data
```

- [ ] 生成并检查：
  - `results/experiment_1/tables/experiment_1_metrics.csv`
  - `results/experiment_1/tables/experiment_1_confusion_matrices.csv`
  - `results/experiment_1/figures/experiment_1_macro_f1.png`
  - `results/experiment_1/figures/experiment_1_confusion_matrices.png`
  - `reports/experiment_1.md`
- [ ] 更新论文：
  - 新表 3-1。
  - 新图 3-1：Macro-F1。
  - 新图 3-2：混淆矩阵，只保留真实模型面板，不再有 fallback/proxy 面板。
- [ ] 检查 hard-negative FPR 分母和实际 false positives。
- [ ] 若时间允许，跑 3-5 个随机种子并汇总 mean±std。

### Step 1.2 实验二：上下文尺度与 hard negative

- [ ] 运行真实数据上的实验二：

```bash
python -m src.experiments.exp2.run_multiscale
```

- [ ] 生成并检查：
  - `results/experiment_2/tables/experiment_2A_multiscale_context.csv`
  - `results/experiment_2/tables/experiment_2B_hard_negative.csv`
  - `results/experiment_2/tables/experiment_2B_rare_motif.csv`
  - `results/experiment_2/figures/exp2A_context_macro_f1.png`
  - `results/experiment_2/figures/exp2A_context_auprc.png`
  - `results/experiment_2/figures/exp2B_hard_negative_fpr.png`
  - `reports/experiment_2.md`
- [ ] 更新论文：
  - 新表 4-1。
  - 新图 4-1 / 4-2 / 4-3。
  - rare motif 补成小表：motif × model。
- [ ] 明确说明：真实主 benchmark 中 non-splice 全部为 GT/AG hard negative，easy negative 只属于可选 synthetic control。
- [ ] 写清楚上下文结果是否非单调，以及哪个窗口在真实数据上最稳。
- [ ] 若时间允许，跑 3-5 个随机种子，尤其给 hard-negative FPR 加置信区间或 mean±std。

### Step 1.3 实验三：真实 ClinVar 变异效应

- [ ] 确认实验三使用真实 ClinVar 表：
  - `configs/exp3_variant_effect.yaml`
  - `data/experiment_3/clinvar_splicing_variants.csv`
- [ ] 运行真实数据上的实验三：

```bash
python -m src.experiments.exp3.run_variant_effect
```

- [ ] 运行解释性分析，delta profile 应换成真实 ClinVar case：

```bash
python -m src.experiments.exp3.run_interpretability
```

- [ ] 生成并检查：
  - `results/experiment_3/tables/experiment_3A_variant_scores.csv`
  - `results/experiment_3/tables/experiment_3A_variant_metrics.csv`
  - `results/experiment_3/tables/experiment_3A_distance_matched_variant_metrics.csv`
  - `results/experiment_3/tables/variant_effect_stratified_by_type.csv`
  - `results/experiment_3/figures/exp3_variant_auroc.png`
  - `results/experiment_3/figures/exp3_variant_auprc.png`
  - `results/experiment_3/figures/exp3_calibration_curve.png`
  - `results/experiment_3/figures/variant_delta_profile_clinvar_real_case.png`
  - `reports/experiment_3.md`
- [ ] 更新论文：
  - 新表 5-2。
  - 新图 5-1 / 5-2。
  - 图 5-3：delta profile 换真实 ClinVar case。
  - 图 5-4：分箱曲线改成基于分类器 delta，标题和横轴量纲不能沿用旧 zero-shot 的 0-1.75。
- [ ] 报告 ClinVar exact-distance-matched 子集上的 AUROC/AUPRC。
- [ ] 明确 AUPRC 基准线：ClinVar 250/250，因此 binary baseline = 0.5；splice-site 三分类 one-vs-rest/positive-vs-negative处按对应定义写清。
- [ ] 不再讲 synthetic donor_gain / acceptor_gain 主结论；故事改成真实 ClinVar splice-altering vs benign，按 donor/acceptor target 分层。
- [ ] 若时间允许，跑 3-5 个随机种子，给 ClinVar 指标加 mean±std 或 bootstrap CI。

---

## Step 2. 按章改写正文

### 摘要

- [ ] 删除 tissue usage / zero-shot / synthetic sandbox 等旧表述。
- [ ] 改成：真实 GENCODE/GRCh38 splice-site benchmark + ClinVar 小样本 benchmark。
- [ ] 方法摘要重写：真实位点、GT/AG hard negative、chromosome holdout、ClinVar pathogenic vs benign、distance-matched diagnostic。
- [ ] 结论语气升级但保持边界：真实小样本、限制清楚。

### 第 1 章：问题与 RQ

- [ ] RQ 保持，但解释更新：
  - RQ1：真实位点三分类。
  - RQ2：真实上下文与真实 hard negative，不再是合成自证。
  - RQ3：真实 ClinVar 变异效应排序。
- [ ] 更新“哪个实验回答哪个 RQ”的表述。
- [ ] 避免承诺 tissue usage 或 gain/loss 四类主分析。

### 第 2 章：数据、模型与评价

#### 2.1 真实 splice-site 数据构造

- [ ] 写清楚：
  - FASTA: GRCh38 / GENCODE v50 primary assembly。
  - GTF: GENCODE v50 annotation。
  - donor = intron 5' boundary。
  - acceptor = intron 3' boundary。
  - 负链反向互补。
  - canonical GT/AG 过滤。
  - non-splice = 不落在注释剪接位点上的真实 GT/AG 诱饵。
  - 采样比例 donor / acceptor / hard-negative = 1:1:1。
  - positive base rate = 2000/3000 = 0.667。
  - chromosome holdout: train chr1-16, valid chr17-18, test chr19-22+X。
- [ ] 表 2-1 替换为真实 split 细分计数：
  - donor / acceptor / easy-negative / hard-negative。
  - hard-FPR 分母。
  - 各 split positive rate。
- [ ] 明确 easy negative 不属于真实主 benchmark，只属于 optional synthetic control。

#### 2.2 ClinVar 构造

- [ ] 写清楚：
  - SNV only。
  - pathogenic / likely pathogenic + splice-related 或 near-splice → positive。
  - benign / likely benign + near-splice → neutral。
  - balanced 250/250。
  - target class = 最近 donor/acceptor。
  - REF 与 GRCh38 校验。
  - 负链 REF/ALT 转 transcript orientation。
  - ClinVar 限制在 held-out test chromosomes，并排除实验一/二已采样 genes。
- [ ] 主动写距离混杂：
  - splice_altering median distance ≈ 9.5。
  - neutral median distance ≈ 17。
  - distance-only AUROC ≈ 0.637。
  - exact-distance-matched subset: 318 rows, 159/159, distance-only AUROC = 0.5。
- [ ] 写明 exact-distance-matched AUROC 是诊断模型是否只学距离捷径。

#### 2.3 模型边界

- [ ] 删除 fallback/proxy 措辞。
- [ ] 模型写成真实可运行对象：
  - CNN baseline。
  - RNA-FM frozen encoder + MLP。
  - RNABERT frozen encoder + MLP。
  - 实验三外部工具：SpliceAI / Pangolin / MMSplice / MaxEntScan。
- [ ] 说明外部工具现在输入真实基因组序列窗口，比较更公平。

#### 2.4 训练协议

- [ ] 加入 valid 选择最优轮 / early stopping。
- [ ] 写清楚 checkpoint、seed、row caps 或 full-data 设置。
- [ ] 若跑多 seed，写 seed 列表。

#### 2.5 评价指标

- [ ] hard-FPR 写清分母：label=non_splice 且 hard_negative。
- [ ] AUPRC 基准线：
  - splice-site positive-vs-negative 时按实际 positive rate 0.667。
  - ClinVar 二分类基准线 0.5。
  - 多分类 macro AUPRC 需说明是 macro one-vs-rest。
- [ ] 若有多 seed，写 mean±std / CI。

### 第 3 章：实验一

- [ ] 使用真实重跑数字重写。
- [ ] 围绕真实 donor/acceptor/non-splice 三分类展开。
- [ ] hard-FPR 仍作为核心诊断。
- [ ] 说明 CNN 与 frozen encoder 差距是否显著；若只有单 seed，不要夸统计显著。
- [ ] 图 3-2 只保留真实模型混淆矩阵。

### 第 4 章：实验二

- [ ] 使用真实重跑数字重写。
- [ ] 讲“上下文非单调”或实际观察到的窗口趋势。
- [ ] hard-negative stress test 作为主证据。
- [ ] rare motif 改为小表，并明确是 synthetic control / optional stress test。
- [ ] 不再把合成上下文当作真实发现来源。

### 第 5 章：实验三

- [ ] 故事改成：真实 ClinVar 致病 vs 良性，按 donor/acceptor target 分层。
- [ ] 不再沿用 donor_loss / donor_gain / acceptor_loss / acceptor_gain 四类合成结论。
- [ ] 解释外部工具比较变公平：真实序列回到 SpliceAI/Pangolin 训练分布附近。
- [ ] 把上一版“强工具垫底”解释清楚：旧结果来自合成序列分布偏移，不作为新结论。
- [ ] 补距离混杂与 exact-distance-matched 诊断。
- [ ] 表 5-2 同时给 full ClinVar 与 distance-matched 指标，或至少在正文报告匹配子集 AUROC/AUPRC。
- [ ] 图 5-4 的标题、横轴、图注改成分类器 delta/calibration，不沿用 zero-shot 旧量纲。
- [ ] AUPRC 基准线标为 0.5。

### 第 6 章：讨论与限制

- [ ] 6.2 改写：
  - 删除 gain/loss 主讨论。
  - 改为 donor/acceptor target、hard-negative 失败模式、真实诱饵上下文。
- [ ] 6.3 限制整体换代：
  - 小样本真实 benchmark。
  - 未做 paralog/homology clustering。
  - ClinVar 限于近位点 SNV。
  - 距离混杂仍需匹配子集诊断。
  - ClinVar 样本集中在 held-out test chromosomes，外推性有限。
- [ ] 不再把“合成不代表真实”作为主限制，只保留为 synthetic control 的边界。

### 第 7 章：结论

- [ ] 重写“能支持/不能支持”清单：
  - 能支持：真实 GT/AG hard-negative 是有效压力测试；真实上下文有价值；ClinVar 上可做真实小样本排序诊断。
  - 不能支持：全基因组泛化、paralog 去泄漏后的结论、组织特异性剪接、gain 类全面评测。
- [ ] 摘要、引言、结论三处声称保持一致。

### 参考文献

- [ ] 补 Pangolin: Zeng & Li, Genome Biology 2022。
- [ ] 确认 SpliceAI、MMSplice、MaxEntScan、GENCODE、ClinVar 引用齐全。

---

## Step 3. 全文一致性与收尾检查

### 关键词清理

- [ ] 全文搜索并处理：

```bash
rg -n "synthetic|fallback|proxy|zero-shot|gain|loss|tissue usage|smoke" README.md reports report_letax
```

- [ ] 每处要么删除，要么明确标为 optional synthetic control / smoke plumbing，不作为主结果。

### 数字与图表一致性

- [ ] 表 3-1 与图 3-1/3-2 一致。
- [ ] 表 4-1 与图 4-1/4-2/4-3 一致。
- [ ] 表 5-2 与图 5-1/5-2/5-4 一致。
- [ ] 5.2、5.4 中“谁最高”的文字和表 5-2 对齐。
- [ ] 所有 AUPRC 处标清 baseline：
  - ClinVar = 0.5。
  - splice-site positive-vs-negative = 0.667。
  - macro AUPRC = macro one-vs-rest。

### 图注复查

- [ ] 图 3-2：真实模型混淆矩阵，不再说 fallback/proxy。
- [ ] 图 5-3：真实 ClinVar delta profile，不再说 artificial donor loss/gain。
- [ ] 图 5-4：分类器 delta/calibration，横轴量纲正确。

### 附录与路径

- [ ] 附录 A 路径更新：
  - `src/data/build_splice_site_dataset.py`
  - `src/data/build_clinvar_variant_dataset.py`
  - `src/data/qc_splice_dataset.py`
  - `data/experiment_3/clinvar_splicing_variants_distance_matched.csv`
- [ ] 命令更新为真实数据主线。
- [ ] 说明 synthetic builders 仅用于 optional control。

### 最终验证命令

- [ ] 数据 QC：

```bash
python -m src.data.qc_splice_dataset
```

- [ ] 单元测试：

```bash
python -m pytest -q
```

- [ ] 报告一致性搜索：

```bash
rg -n "synthetic|fallback|proxy|zero-shot|gain|loss|tissue usage|smoke" README.md reports report_letax
```

- [ ] LaTeX 编译通过，图路径无缺失。

---

## Done Definition

完成标准：

- [ ] `reports/data_qc.md` 除 paralog/homology limitation 外无 FAIL。
- [ ] 三个实验全部在真实数据上重跑，并生成新表新图。
- [ ] 正文所有数字来自新结果。
- [ ] ClinVar 距离混杂和 distance-matched diagnostic 写入第 2 章和第 5 章。
- [ ] Paralog 未去重写入限制。
- [ ] 全文不再把 synthetic/proxy/fallback 当主结果。
- [ ] 摘要、引言、结论三个层级的主张一致。
