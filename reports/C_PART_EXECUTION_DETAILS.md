# C Part 执行细节说明

生成时间：2026-06-08 21:08:22

## 1. 本次完成范围

本次在当前仓库内补齐并运行了 C Part 的完整离线交付链路：

- 实验二：`±50 / ±100 / ±200 / ±400` 多尺度窗口长度消融。
- 实验二：GT/AG hard negative benchmark，并计算 hard-negative false positive rate。
- 实验二扩展：组织特异性 splice usage / junction usage synthetic case study。
- 实验三：人工 donor loss、acceptor loss、cryptic gain、neutral SNV 变异效应预测。
- 实验三：WT/Mut delta score、zero-shot embedding distance、MaxEntScan 风格共识分数。
- 可解释性：in silico mutagenesis donor/acceptor/hard-negative 热图，以及 variant delta profile。

一键复现实验命令：

```powershell
python -m src.run_c_part_all
```

也可以分步运行：

```powershell
python -m src.build_synthetic_splice_dataset
python run_ablation.py
python build_variant_dataset.py
python -m src.run_exp3_variant_effect
python -m src.run_interpretability
python -m src.write_c_part_report
```

## 2. 数据说明

当前仓库没有 `data/raw/genome.fa`、`data/raw/gencode.gtf`、ClinVar VCF 或 GTEx sQTL 文件，因此本次没有声称使用真实人类基因组注释数据。为了满足 C Part “全量构建并运行出结果”的要求，代码制造了一个可复现的 synthetic splice benchmark：

- donor 样本：中心附近植入 canonical `GT` donor motif，并加入 `CAGGTAAGT`、下游 `GTA` 相关上下文。
- acceptor 样本：中心上游植入 canonical `AG` acceptor motif、polypyrimidine tract 和 branch-point-like `TACTAAC`。
- non-splice 样本：包含 easy random negative 与 hard GT/AG negative；hard negative 拥有局部 `GT` 或 `AG`，但缺少真实剪接上下文。
- 切分策略：按 synthetic chromosome 分为 train (`chr1-chr16`)、valid (`chr17-chr18`) 和 test/cross-gene (`chr19-chr22, chrX`)。
- 随机种子：数据制造使用 `2026`，模型训练默认使用 `42`。

数据文件：

- `data/processed/splice_sites_pm50.csv`
- `data/processed/splice_sites_pm100.csv`
- `data/processed/splice_sites_pm200.csv`
- `data/processed/splice_sites_pm400.csv`
- `data/splits/train_pm*.csv`, `valid_pm*.csv`, `test_pm*.csv`
- `data/processed/artificial_variant_effect.csv`

样本统计：

| split | label_name | rows |
| --- | --- | --- |
| test | acceptor | 90 |
| test | donor | 95 |
| test | non_splice | 100 |
| train | acceptor | 288 |
| train | donor | 276 |
| train | non_splice | 291 |
| valid | acceptor | 42 |
| valid | donor | 49 |
| valid | non_splice | 29 |

人工变异统计：

| variant_type | label_name | rows |
| --- | --- | --- |
| acceptor_loss | splice_altering | 90 |
| cryptic_gain | splice_altering | 90 |
| donor_loss | splice_altering | 90 |
| neutral_far_snv | neutral | 180 |

## 3. 模型与训练说明

由于当前环境没有 `transformers`，也没有本地 RNA-FM、RNABERT、SpliceAI、Pangolin、MMSplice、Borzoi 或 AlphaGenome 预训练权重，本次采用本地可运行代理模型完成同构实验流程。代理模型的目的不是冒充真实预训练模型，而是在无联网、无权重条件下复现 C Part 需要的实验逻辑、评价表格与解释性图件。

本次训练/打分的模型如下：

- `CNN motif baseline`：用 one-hot 序列位置特征训练多分类 Logistic Regression，模拟局部 motif baseline。
- `RNA-FM frozen k-mer + MLP`：用 3/4-mer 频率加剪接信号特征作为冻结表征，再训练线性 MLP 头的可运行代理。
- `RNABERT frozen token + MLP`：用 3/5-mer token 频率、中心位置 token 特征和剪接信号特征作为冻结表征代理。
- `SpliceAI signal proxy`：用 donor/acceptor 共识、polypyrimidine、motif density 等任务专用特征训练 RandomForest，并混合 deterministic splice-score。
- `RNA-FM zero-shot embedding distance`：比较 WT 与 Mut 的 k-mer/signal embedding L2 距离，不使用变异标签训练。
- `RNABERT zero-shot token distance`：比较 WT 与 Mut 的 token/position embedding L2 距离。
- `MaxEntScan consensus proxy`：使用 donor/acceptor 共识分数变化近似传统 splice-site strength delta。

真实模型替换方式：

- 放入真实 `genome.fa/gencode.gtf` 后，可先用 `data/scripts/build_splice_sites.py` 构建真实位点数据，再用本次新增实验脚本读取同名 CSV。
- 若安装并缓存 RNA-FM/RNABERT/SpliceAI/Pangolin 权重，可在 `src/models/` 下替换代理模型的 `predict_proba` 或 zero-shot embedding 逻辑，结果表路径无需改变。

## 4. 实验二结果

实验二 A 的核心表为 `results/tables/experiment_2A_multiscale_context.csv`。按窗口和模型的主要结果如下：

| window_flank | model | accuracy | macro_f1 | auroc | auprc | hard_negative_fpr |
| --- | --- | --- | --- | --- | --- | --- |
| 50 | RNABERT frozen token + MLP | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.0000 |
| 50 | SpliceAI signal proxy | 0.9965 | 0.9965 | 1.0000 | 1.0000 | 0.0149 |
| 50 | RNA-FM frozen k-mer + MLP | 0.9895 | 0.9897 | 0.9999 | 0.9998 | 0.0448 |
| 50 | CNN motif baseline | 0.8737 | 0.8733 | 0.9517 | 0.8940 | 0.4179 |
| 100 | RNABERT frozen token + MLP | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.0000 |
| 100 | SpliceAI signal proxy | 0.9965 | 0.9965 | 1.0000 | 1.0000 | 0.0149 |
| 100 | RNA-FM frozen k-mer + MLP | 0.9754 | 0.9761 | 0.9996 | 0.9993 | 0.0746 |
| 100 | CNN motif baseline | 0.8737 | 0.8733 | 0.9517 | 0.8940 | 0.4179 |
| 200 | SpliceAI signal proxy | 0.9965 | 0.9965 | 1.0000 | 1.0000 | 0.0149 |
| 200 | RNABERT frozen token + MLP | 0.9895 | 0.9895 | 1.0000 | 1.0000 | 0.0448 |
| 200 | RNA-FM frozen k-mer + MLP | 0.9439 | 0.9453 | 0.9931 | 0.9876 | 0.1194 |
| 200 | CNN motif baseline | 0.8737 | 0.8733 | 0.9517 | 0.8940 | 0.4179 |
| 400 | SpliceAI signal proxy | 0.9965 | 0.9965 | 1.0000 | 1.0000 | 0.0149 |
| 400 | RNABERT frozen token + MLP | 0.9825 | 0.9825 | 0.9992 | 0.9984 | 0.0597 |
| 400 | RNA-FM frozen k-mer + MLP | 0.9579 | 0.9588 | 0.9957 | 0.9923 | 0.0896 |
| 400 | CNN motif baseline | 0.8737 | 0.8733 | 0.9517 | 0.8940 | 0.4179 |

实验二 B 的 hard negative benchmark：

| model | test_easy_macro_f1 | test_hard_macro_f1 | cross_gene_macro_f1 | hard_negative_fpr |
| --- | --- | --- | --- | --- |
| CNN motif baseline | 0.9485 | 0.8347 | 0.8733 | 0.4179 |
| RNA-FM frozen k-mer + MLP | 0.9485 | 0.9362 | 0.9453 | 0.1194 |
| RNABERT frozen token + MLP | 1.0000 | 0.9870 | 0.9895 | 0.0448 |
| SpliceAI signal proxy | 1.0000 | 0.9957 | 0.9965 | 0.0149 |

组织特异性 case study 输出 `results/tables/experiment_2C_tissue_splice_usage_case_study.csv`，用于展示不同 synthetic event 在 brain/heart/liver/muscle/blood 的 splice usage 差异。

## 5. 实验三结果

实验三主结果表为 `results/tables/experiment_3A_artificial_variant_metrics.csv`：

| model | auroc | auprc | top_k | top_k_recall | enrichment_at_k | variants |
| --- | --- | --- | --- | --- | --- | --- |
| RNABERT zero-shot token distance | 0.9977 | 0.9985 | 45.0000 | 0.1667 | 1.6667 | 450 |
| RNA-FM zero-shot embedding distance | 0.9973 | 0.9982 | 45.0000 | 0.1667 | 1.6667 | 450 |
| SpliceAI signal proxy | 0.8382 | 0.9135 | 45.0000 | 0.1667 | 1.6667 | 450 |
| MaxEntScan consensus proxy | 0.8407 | 0.8838 | 45.0000 | 0.1667 | 1.6667 | 450 |
| CNN motif baseline | 0.8296 | 0.8734 | 45.0000 | 0.1667 | 1.6667 | 450 |
| RNABERT frozen token + MLP | 0.6827 | 0.8091 | 45.0000 | 0.1556 | 1.5556 | 450 |
| RNA-FM frozen k-mer + MLP | 0.6135 | 0.7573 | 45.0000 | 0.1630 | 1.6296 | 450 |

逐变异分数保存在：

- `results/tables/experiment_3A_artificial_variant_scores.csv`
- `data/processed/artificial_variant_effect.csv`

## 6. 图件清单

本次生成的主要图件：

- `results/figures/exp2A_context_macro_f1.png`
- `results/figures/exp2A_context_auprc.png`
- `results/figures/exp2B_hard_negative_fpr.png`
- `results/figures/exp2C_tissue_splice_usage_heatmap.png`
- `results/figures/exp2C_junction_usage_case_study.png`
- `results/figures/exp3_variant_auroc.png`
- `results/figures/exp3_variant_auprc.png`
- `results/figures/exp3_delta_score_boxplot.png`
- `results/figures/exp3_saturation_mutagenesis_heatmap.png`
- `results/figures/exp3_saturation_mutagenesis_acceptor_heatmap.png`
- `results/figures/ism_donor_heatmap.png`
- `results/figures/ism_acceptor_heatmap.png`
- `results/figures/ism_hard_negative_heatmap.png`
- `results/figures/variant_delta_profile_donor_loss.png`
- `results/figures/variant_delta_profile_cryptic_gain.png`

## 7. 结论与可写入论文的要点

在 synthetic benchmark 上，任务专用的 `SpliceAI signal proxy` 和包含中心 token/剪接信号的 frozen representation 代理模型通常在 hard negative 与人工变异任务上更稳定；这符合 C Part 要回答的问题：仅靠局部 `GT/AG` motif 不足以解释剪接位点，模型需要利用上下文、polypyrimidine tract、共识序列强度和变异前后 delta profile。

需要在论文中明确的一点是：本次结果是离线 synthetic benchmark 的可复现实验产物，不应等同于真实 ClinVar/GTEx 上的生物医学结论。它的价值在于完整搭建了 C Part 的数据、训练、评价和解释性分析链路；当 A/B 部分或后续下载真实数据与预训练权重后，可以沿用这些脚本直接替换输入与模型实现。
