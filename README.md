# RNA 剪接位点识别与 ClinVar 变异效应预测

本项目已经从早期演示性方案收敛为一个小规模但可追溯的真实数据 benchmark：基于 GRCh38、GENCODE v49 和 ClinVar，评估剪接位点识别、GT/AG hard-negative 压力测试，以及 ClinVar splice-altering vs benign SNV 变异效应排序。

## 研究问题

1. 在真实 GENCODE/GRCh38 小样本条件下，模型能否区分 donor、acceptor 与 GT/AG hard-negative？
2. 不同上下文窗口长度是否影响模型拒绝 motif-matched negative 的能力？
3. 剪接识别信号能否转化为真实 ClinVar pathogenic-vs-benign 变异效应排序，并在 distance-matched 子集上保持区分度？

## 数据

主 splice-site benchmark 由 `src/data/build_splice_site_dataset.py` 从 `data/raw/genome.fa` 和 `data/raw/gencode.gtf` 构建。全量采样为 donor / acceptor / hard-negative = 1000 / 1000 / 1000；non-splice 类全部为真实基因组 GT/AG motif-matched decoy，不使用 easy negative 作为主结果。

主 split 使用 chromosome holdout：

| split | rows | donor | acceptor | hard-negative |
| --- | ---: | ---: | ---: | ---: |
| train | 2339 | 792 | 786 | 761 |
| valid | 230 | 79 | 67 | 84 |
| test | 431 | 129 | 147 | 155 |

ClinVar benchmark 由 `src/data/build_clinvar_variant_dataset.py` 构建，包含 250 条 splice-altering SNV 与 250 条 benign/likely benign near-splice SNV，并限制在 held-out test chromosomes。由于 splice-altering 变异更靠近注释剪接位点，项目同时导出 exact-distance-matched 子集 `data/experiment_3/clinvar_splicing_variants_distance_matched.csv` 用于距离混杂诊断。

当前 QC 记录：

- GENCODE: v49 / Ensembl 115
- ClinVar fileDate: 2026-05-30
- seed: 42
- QC 报告: `reports/data_qc.md`

## 模型

实验一和实验二使用：

- CNN baseline
- RNA-FM frozen encoder + MLP
- RNABERT frozen encoder + MLP

实验三额外接入真实外部工具输出：

- SpliceAI
- Pangolin
- MMSplice
- MaxEntScan

所有主表只报告真实本地模型或真实外部工具输出；RNA-FM/RNABERT 使用本地真实权重 frozen encoder。

## 复现命令

```bash
python -m src.data.build_splice_site_dataset --max-per-class 1000 --windows 50 100 200 400
python -m src.data.build_clinvar_variant_dataset
python -m src.data.qc_splice_dataset
python -m src.experiments.exp1.run_classification --full-data
python -m src.experiments.exp2.run_multiscale
python -m src.experiments.exp3.run_variant_effect
python -m src.experiments.exp3.run_interpretability
python -m src.reports.write_c_part_report
```

## 关键结果

实验一测试集三分类结果：

| model | Accuracy | Macro-F1 | AUROC | AUPRC | Hard FPR |
| --- | ---: | ---: | ---: | ---: | ---: |
| CNN baseline | 0.8190 | 0.8206 | 0.9488 | 0.9072 | 0.3226 |
| RNA-FM frozen encoder + MLP | 0.7633 | 0.7670 | 0.9157 | 0.8402 | 0.3290 |
| RNABERT frozen encoder + MLP | 0.8190 | 0.8199 | 0.9298 | 0.8718 | 0.3226 |

实验三 full ClinVar 排序结果：

| model | AUROC | AUPRC |
| --- | ---: | ---: |
| RNA-FM frozen encoder + MLP | 0.7532 | 0.7839 |
| Pangolin real sequence model | 0.7023 | 0.7781 |
| CNN baseline | 0.7018 | 0.7399 |
| RNABERT frozen encoder + MLP | 0.6673 | 0.7326 |
| SpliceAI real sequence model | 0.6425 | 0.7289 |
| MaxEntScan real local score | 0.6165 | 0.6716 |
| MMSplice real sequence model | 0.5000 | 0.5000 |

实验三 exact-distance-matched 子集显示多数模型分数下降，说明 full ClinVar 的一部分可分性来自近剪接距离混杂。

## 报告与论文

- 实验报告：`reports/experiment_1.md`、`reports/experiment_2.md`、`reports/experiment_3.md`
- 组合报告：`reports/c_part_combined_report.md`
- 数据 QC：`reports/data_qc.md`
- 论文 LaTeX：`report_letax/njuthesis-sample.tex`
- 论文图片：`report_letax/images/`

## 限制

当前项目是小样本真实 benchmark，而不是临床级全量评测。主要限制包括：未做 paralog/homology clustering，只使用随机种子 42，CNN 未使用 early stopping，ClinVar 仅覆盖 near-splice SNV，distance-matched 子集只能诊断并部分缓解距离混杂，不能替代更大规模严格匹配 benchmark。
