# 实验一报告：剪接位点三分类

## 实验目标

本实验评估模型能否根据中心窗口序列区分三类样本：donor、acceptor 与 non-splice。该实验对应 README 中的基础任务，是后续多尺度上下文分析和变异效应预测的分类器基础。

## 数据与输出位置

- 共用剪接位点数据：`data/shared/processed/`
- 共用划分文件：`data/shared/splits/`
- 实验一结果表：`results/experiment_1/tables/`
- 实验一图件：`results/experiment_1/figures/`
- 小样本真实 foundation model 结果：`results/experiment_1/tables/real_foundation/` 与 `results/experiment_1/figures/real_foundation/`

## 模型

当前实验一主结果包含：

- `CNN baseline (PyTorch Conv1D)`
- `RNA-FM + MLP`
- `RNABERT + MLP`

其中 `results/experiment_1/tables/real_foundation/` 保存真实 RNA-FM / RNABERT frozen encoder 小样本结果；主表中的 RNA-FM/RNABERT 结果需结合报告说明区分 proxy 与真实 encoder 运行。

## 主要结果

测试集主指标如下：

| 模型 | Accuracy | Macro-F1 | AUROC | AUPRC |
| --- | ---: | ---: | ---: | ---: |
| CNN baseline | 0.874 | 0.871 | 0.968 | 0.938 |
| RNA-FM + MLP | 0.845 | 0.842 | 0.944 | 0.896 |
| RNABERT + MLP | 0.870 | 0.869 | 0.958 | 0.921 |

CNN 与 RNABERT + MLP 在当前测试集上接近，CNN 的 Macro-F1 略高；RNA-FM + MLP 略低。三类 F1 中 non-splice 类整体低于 donor/acceptor，说明负样本尤其是 hard negative 仍是主要难点。

## 关键文件

- `results/experiment_1/tables/experiment_1_metrics.csv`
- `results/experiment_1/tables/experiment_1_confusion_matrices.csv`
- `results/experiment_1/figures/experiment_1_macro_f1.png`
- `results/experiment_1/figures/experiment_1_confusion_matrices.png`

## 结论

实验一已经完成剪接位点三分类的基础训练与评估链路。当前结果说明模型可以稳定识别 donor/acceptor，但 non-splice 与剪接 motif 相似时仍存在混淆。后续最重要的补齐项是将 SpliceAI 真实 baseline 纳入实验一主表，并进一步明确 proxy 结果与真实 foundation encoder 结果的边界。
