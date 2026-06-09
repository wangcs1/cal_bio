# C Part 工作总结：剪接位点识别的多尺度上下文实验

本仓库完成了 C Part 的实验代码、数据构建、结果表格和可视化输出。核心问题是：剪接位点识别究竟主要依赖局部 `GT/AG` motif，还是需要更长程的 intron/exon context、剪接增强子/沉默子、组织特异性调控和 splice junction topology？

当前结果支持一个明确结论：`GT/AG` 是剪接识别的重要局部信号，但不是充分条件。模型要稳定地区分真实 donor、acceptor 和带有 `GT/AG` 的 hard negative，需要引入局部强度、上下文长度、调控 motif、组织程序和 junction 拓扑信息。

## 数据

主要实验使用可复现的合成剪接 benchmark，并结合真实资源与真实模型 smoke test 做依赖验证。

- 三分类剪接位点数据：`data/processed/synthetic_splice_sites_master_pm400.csv`
  - 任务：donor / acceptor / non-splice 三分类。
  - 划分：`train_pm400.csv` 855 条，`valid_pm400.csv` 120 条，`test_pm400.csv` 285 条，`cross_gene_test.csv` 285 条。
  - 包含普通 non-splice 和带局部 `GT/AG` 的 hard negative，用于检验 motif-only 模型是否会误报。
- 变异效应数据：`data/processed/artificial_variant_effect.csv`
  - 共 450 个 SNV，其中 `donor_loss` 90、`acceptor_loss` 90、`cryptic_gain` 90、`neutral_far_snv` 180。
  - 任务：预测 splice-altering 与 neutral，并按变异类型分析模型响应。
- 组织调控 case-study：`results/tables/tissue_specific_usage_case_study.csv`
  - 共 200 行，覆盖 40 个事件、5 个组织和 5 类组织调控程序。
  - 用于检验单纯位点强度是否足够解释组织间 splice usage 差异。
- Junction topology case-study：`results/tables/junction_topology_case_study.csv`
  - 共 42 条 junction，包含 donor/acceptor strength、regulatory density、donor/acceptor degree、competing junctions 和 exon skipping 标记。
- 真实资源与真实模型 smoke test：
  - 本地准备了 `data/raw/genome.fa`、`data/raw/gencode.gtf`、ClinVar VCF、GTEx/sQTL 与 known splice event 小样本。
  - 输出在 `results/real_smoke/`，包括 RNA-FM、RNABERT embedding smoke test，SpliceAI、Pangolin、MaxEntScan、MMSplice 的小样本调用结果。
  - 大型 raw 文件不作为必要提交内容，仓库提交重点是可复现实验代码、派生数据与结果。

## 方法

实验覆盖从局部 motif 到多尺度上下文的多个层级。

- 局部 motif 与局部强度：
  - `Motif-only GT/AG rule`
  - `MaxEntScan local score + motif`
  - short-context one-hot logistic 模型
- 上下文长度消融：
  - 对 `+/-10, +/-50, +/-100, +/-200, +/-400` flank 长度分别训练/评估。
  - 使用 RNA-FM-style k-mer + signal、RNABERT-style token + position、SpliceAI signal proxy 三类轻量代理模型。
- 调控 motif 分析：
  - 构造 ESE/ISE/ESS/ISS proxy motif 组。
  - 通过 motif masking 和 rescue proxy 观察模型目标类别概率变化。
- 变异效应预测：
  - 比较 RNABERT/RNA-FM zero-shot distance proxy、SpliceAI signal proxy 和 MaxEntScan consensus proxy 对 donor loss、acceptor loss、cryptic gain、neutral SNV 的响应。
- 组织特异性调控：
  - 比较 `site_sequence_proxy`、`site_plus_tissue_label`、`site_plus_tissue_program`。
  - 重点检验组织程序和组织交互是否解释 splice usage。
- Junction topology：
  - 比较 `site_strength_only`、`topology_only`、`site_plus_topology`。
  - 检验 competing junctions、junction degree、exon skipping topology 是否提供超出位点强度的信息。

主流程入口为：

```bash
python -m src.run_c_part_all
```

真实资源与模型 smoke test 可单独运行：

```bash
python scripts/run_real_model_smoke.py
```

## 关键结果

### 1. GT/AG motif 不是充分条件

结果表：`results/tables/local_motif_sufficiency.csv`

| 模型 | Macro-F1 | Hard negative FPR |
|---|---:|---:|
| Motif-only GT/AG rule | 0.710 | 1.000 |
| MaxEntScan local score + motif | 0.972 | 0.104 |
| Short-context one-hot +/-10 | 0.853 | 0.418 |

Motif-only 规则会把所有 hard `GT/AG` negative 误判为剪接位点，说明 canonical dinucleotide 只能作为必要局部线索，不能单独完成可靠识别。

### 2. 更长上下文显著降低 hard negative 误报

结果表：`results/tables/context_length_ablation_full.csv`

代表性结果：

| 模型 | +/-10 Macro-F1 | +/-50 Macro-F1 | +/-10 Hard FPR | +/-50 Hard FPR |
|---|---:|---:|---:|---:|
| RNA-FM-style k-mer + signal | 0.787 | 0.990 | 0.373 | 0.045 |
| RNABERT-style token + position | 0.789 | 1.000 | 0.388 | 0.000 |
| SpliceAI signal proxy | 0.886 | 0.996 | 0.478 | 0.015 |

在当前 benchmark 中，从 `+/-10` 增加到 `+/-50` 已经带来主要收益，说明剪接识别依赖 motif 附近更宽的序列上下文。

### 3. 增强子/沉默子 proxy motif 会改变模型判断

结果表：`results/tables/regulatory_motif_masking.csv`

代表性 masking 结果：

- donor 样本中 masking `ISE_proxy` 后，目标类别概率平均下降 `0.804`。
- donor 样本中 masking `ESE_proxy` 后，目标类别概率平均下降 `0.692`。
- acceptor 样本中 masking `ISS_proxy` 后，目标类别概率平均下降 `0.689`。

这说明模型判断不只来自中心 `GT/AG`，还受到周边调控 motif 的影响。

### 4. 变异效应预测能区分核心位点破坏与中性远端 SNV

结果表：`results/tables/variant_effect_stratified_by_type.csv`

代表性结果：

- RNABERT zero-shot token distance：`donor_loss` 均值 `1.755`，`acceptor_loss` 均值 `1.757`，`neutral_far_snv` 均值 `0.011`。
- RNA-FM zero-shot embedding distance：`donor_loss` 均值 `1.040`，`acceptor_loss` 均值 `1.042`，`neutral_far_snv` 均值 `0.010`。
- SpliceAI signal proxy：`donor_loss` 均值 `0.698`，`acceptor_loss` 均值 `0.526`，`neutral_far_snv` 均值 `0.001`。

核心剪接位点破坏的模型响应明显高于 neutral SNV；cryptic gain 对部分 proxy 模型更难，需要更真实的全上下文模型进一步验证。

### 5. 组织特异性需要事件-组织调控程序

结果表：`results/tables/tissue_specific_usage_ablation.csv`

| 模型 | MAE | R2 |
|---|---:|---:|
| site_sequence_proxy | 0.084 | 0.427 |
| site_plus_tissue_label | 0.085 | 0.409 |
| site_plus_tissue_program | 0.020 | 0.969 |

单纯加入 tissue label 并不够，真正有效的是事件对应的 tissue regulatory program 与 tissue 的交互。这支持“组织特异性剪接不是只由位点局部序列决定”的结论。

### 6. Junction topology 提供位点强度之外的信息

结果表：`results/tables/junction_topology_ablation.csv`

| 模型 | MAE | R2 |
|---|---:|---:|
| site_strength_only | 0.172 | 0.166 |
| topology_only | 0.085 | 0.790 |
| site_plus_topology | 0.075 | 0.828 |

加入 junction degree、competing junctions、exon skipping 等拓扑信息后，splice usage 预测明显改善。

## 主要输出

- 代码：
  - `src/run_c_part_all.py`
  - `src/run_full_context_question.py`
  - `build_variant_dataset.py`
  - `run_ablation.py`
  - `scripts/run_real_model_smoke.py`
- 结果表：
  - `results/tables/local_motif_sufficiency.csv`
  - `results/tables/context_length_ablation_full.csv`
  - `results/tables/regulatory_motif_masking.csv`
  - `results/tables/motif_rescue_proxy.csv`
  - `results/tables/variant_effect_stratified_by_type.csv`
  - `results/tables/tissue_specific_usage_ablation.csv`
  - `results/tables/junction_topology_ablation.csv`
- 图：
  - `results/figures/local_motif_sufficiency_hard_fpr.png`
  - `results/figures/context_length_vs_macro_f1.png`
  - `results/figures/context_length_vs_hard_fpr.png`
  - `results/figures/regulatory_motif_masking_effect.png`
  - `results/figures/variant_effect_stratified_by_type.png`
  - `results/figures/tissue_specific_usage_ablation.png`
  - `results/figures/junction_topology_ablation.png`
  - `results/figures/splice_junction_graph_case_study.png`

## 结论

当前 C Part 已经从多个角度回答了研究问题：

1. `GT/AG` motif 是必要局部信号，但不是充分条件。
2. 更长的 intron/exon context 能显著降低 hard negative 误报。
3. 增强子/沉默子 proxy motif 会改变模型类别概率，说明周边调控序列有贡献。
4. 组织特异性需要事件-组织调控程序，而不是简单 tissue label。
5. Junction topology 能解释位点强度无法解释的 splice usage 差异。

需要注意的是，组织调控和 junction topology 目前是可复现的 proxy/case-study 实验；真实 GTEx 全量组织模型、真实 sQTL 标签和真实 SpliceAI/Pangolin 全量 benchmark 可以作为后续扩展。
