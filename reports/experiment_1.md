# 实验一报告：剪接位点三分类

同步日期：2026-06-11

## 实验目标

实验一评估 donor、acceptor、non-splice 三分类任务。主实验只使用当前小样本 split：

- `data/shared/splits/train.csv`：855 条
- `data/shared/splits/valid.csv`：120 条
- `data/shared/splits/test.csv`：285 条

`data/raw/` 中的完整 genome/GTF 和外部 foundation model 权重均为可选资源，不是本实验默认运行前提。

## 模型边界

| 模型 | 当前报告中的含义 |
|:--|:--|
| CNN baseline | PyTorch Conv1D 基线 |
| RNA-FM frozen encoder + MLP | 有本地权重时使用 frozen encoder；否则使用确定性 proxy/fallback embedding |
| RNABERT frozen encoder + MLP | 有本地权重时使用 frozen encoder；否则使用确定性 proxy/fallback embedding |
| SpliceAI optional real tool | 当前环境没有真实 `spliceai` 包时使用 `spliceai_signal_proxy`，不写成完整真实 SpliceAI benchmark |

hard-negative FPR 只在 `label == 2` 且 `negative_type` 包含 hard 的 67 条 test 样本上统计。

## 验证集结果

| 模型 | train | valid | Macro-F1 | Accuracy | AUROC | AUPRC |
|:--|--:|--:|--:|--:|--:|--:|
| CNN baseline (PyTorch Conv1D) | 855 | 120 | 0.9065 | 0.9250 | 0.9823 | 0.9613 |
| RNA-FM frozen encoder + MLP | 855 | 120 | 0.9363 | 0.9417 | 0.9950 | 0.9891 |
| RNABERT frozen encoder + MLP | 855 | 120 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| SpliceAI optional real tool (proxy fallback) | 855 | 120 | 0.8709 | 0.9000 | 0.9936 | 0.9870 |

## 测试集结果

| 模型 | 类型 | backend | rows | Accuracy | Macro-F1 | AUROC | AUPRC | Hard-negative FPR |
|:--|:--|:--|--:|--:|--:|--:|--:|--:|
| CNN baseline (PyTorch Conv1D) | baseline | `pytorch_conv1d` | 285 | 0.9018 | 0.9005 | 0.9899 | 0.9802 | 0.4179 |
| RNA-FM frozen encoder + MLP | frozen encoder | `local_pretrained_or_proxy_fallback` | 285 | 0.9509 | 0.9519 | 0.9938 | 0.9888 | 0.1493 |
| RNABERT frozen encoder + MLP | frozen encoder | `local_pretrained_or_proxy_fallback` | 285 | 0.9930 | 0.9931 | 0.9999 | 0.9998 | 0.0149 |
| SpliceAI optional real tool (proxy fallback) | proxy | `spliceai_signal_proxy` | 285 | 0.8737 | 0.8680 | 0.9981 | 0.9962 | 0.5224 |

## 结论

RNABERT fallback/frozen 路线在当前小样本 test 上取得最高 Macro-F1 与最低 hard-negative FPR。RNA-FM 也明显优于 CNN baseline。SpliceAI 行当前是 signal proxy fallback，可作为规则型信号参照，但不能表述为真实 SpliceAI 完整推理结论。

## 输出文件

- `results/experiment_1/tables/experiment_1_metrics.csv`
- `results/experiment_1/tables/experiment_1_train_valid_metrics.csv`
- `results/experiment_1/tables/experiment_1_confusion_matrices.csv`
- `results/experiment_1/figures/experiment_1_macro_f1.png`
- `results/experiment_1/figures/experiment_1_confusion_matrices.png`
- `results/experiment_1/tables/real_foundation/experiment_1_metrics.csv`
