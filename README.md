# 多尺度 RNA 剪接建模与异常剪接变异效应预测

本项目围绕 RNA 剪接位点识别、多尺度上下文建模和异常剪接变异效应预测，比较 CNN、RNA foundation model proxy/frozen encoder、SpliceAI/Pangolin/MMSplice/MaxEntScan 风格工具以及长程调控 case study 的信息贡献。

## 数据规模说明

当前主实验默认使用小样本 synthetic/proxy split，不使用百万级完整数据集：

- `data/shared/splits/train.csv`：855 条样本
- `data/shared/splits/valid.csv`：120 条样本
- `data/shared/splits/test.csv`：285 条样本

`data/raw/` 中的 `genome.fa`、`gencode.gtf`、ClinVar、GTEx/sQTL 等大文件只作为可选真实资源，用于 smoke test 或 case study，不是主实验运行前提。`models/hf/rnafm/` 和 `models/hf/rnabert/` 也是可选本地权重目录；缺少权重时主流程自动使用 proxy/fallback 表征。

## 主要实验

1. 实验一：donor / acceptor / non-splice 三分类。
2. 实验二：多窗口上下文、GT/AG hard negative、rare motif、调控 motif、tissue usage 和 junction topology。
3. 实验三：donor loss、acceptor loss、donor gain、acceptor gain、neutral SNV 的变异效应预测。
4. 可解释性：attention proxy、in silico mutagenesis、variant delta profile。

## 运行

轻量环境：

```bash
pip install -r requirements.txt
```

一键运行当前小样本主流程：

```bash
python -m src.pipelines.run_c_part_all
```

可选真实模型/真实资源依赖：

```bash
pip install -r requirements-real.txt
python -m src.resources.prepare_raw_resources
python -m src.resources.prepare_foundation_models --check
```

## 输出

- 实验报告：`reports/experiment_1.md`、`reports/experiment_2.md`、`reports/experiment_3.md`
- 数据 QC：`reports/data_qc.md`
- 资源说明：`reports/resource_setup.md`
- 模型卡：`reports/model_cards.md`
- 实验日志：`reports/experiment_log.md`
- 结果表和图：`results/experiment_1/`、`results/experiment_2/`、`results/experiment_3/`

## 结论边界

当前结果用于证明“小样本可复现的多尺度剪接建模链路”已经跑通：局部 `GT/AG` motif 不是充分条件，模型还需要上下文、调控 motif、组织程序和 junction topology。真实 ClinVar/GTEx 全量 benchmark 与大型长程模型推理保留为可选扩展，不在当前主结果中声称已完成。
