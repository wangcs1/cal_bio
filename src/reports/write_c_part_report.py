from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd

from src.utils import (
    EXP2_TABLES_DIR,
    EXP3_DATA_DIR,
    EXP3_TABLES_DIR,
    PROJECT_ROOT,
    SHARED_PROCESSED_DIR,
    SHARED_SPLIT_DIR,
    load_or_empty,
)


def markdown_table(frame: pd.DataFrame, columns: list[str], max_rows: int = 12) -> str:
    if frame.empty:
        return "_未生成_"
    subset = frame.loc[:, [col for col in columns if col in frame.columns]].head(max_rows).copy()
    for col in subset.columns:
        if pd.api.types.is_float_dtype(subset[col]):
            subset[col] = subset[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
        else:
            subset[col] = subset[col].map(lambda x: "" if pd.isna(x) else str(x).replace("|", r"\|"))
    header = "| " + " | ".join(map(str, subset.columns)) + " |"
    sep = "| " + " | ".join(["---"] * len(subset.columns)) + " |"
    rows = []
    for _, row in subset.iterrows():
        rows.append("| " + " | ".join(str(row[col]) for col in subset.columns) + " |")
    return "\n".join([header, sep, *rows])


def write_report(root: Path = PROJECT_ROOT) -> Path:
    exp2_tables = EXP2_TABLES_DIR
    exp3_tables = EXP3_TABLES_DIR
    processed = SHARED_PROCESSED_DIR
    splits = SHARED_SPLIT_DIR

    summary = load_or_empty(processed / "synthetic_splice_sites_summary.csv")
    variants_summary = load_or_empty(EXP3_DATA_DIR / "artificial_variant_effect_summary.csv")
    exp2a = load_or_empty(exp2_tables / "experiment_2A_multiscale_context.csv")
    exp2b = load_or_empty(exp2_tables / "experiment_2B_hard_negative.csv")
    exp3 = load_or_empty(exp3_tables / "experiment_3A_artificial_variant_metrics.csv")
    tissue = load_or_empty(exp2_tables / "experiment_2C_tissue_splice_usage_case_study.csv")
    real_smoke = root / "results/real_smoke"
    foundation_smoke = load_or_empty(real_smoke / "foundation_model_smoke.csv")
    maxent_mmsplice_smoke = load_or_empty(real_smoke / "maxentscan_mmsplice_smoke.csv")
    spliceai_smoke = load_or_empty(real_smoke / "spliceai_clinvar_smoke_summary.csv")

    best_exp2 = exp2a.sort_values(["window_flank", "macro_f1"], ascending=[True, False])
    best_exp3 = exp3.sort_values("auprc", ascending=False)

    real_resource_section = """## 0. 真实资源补充状态（2026-06-09 更新）

当前版本已经在 WSL2 / Ubuntu-22.04 环境下补齐并验证 C Part 所需真实资源和真实模型依赖：

- `data/raw/genome.fa`：已下载并解压，GENCODE v49 / GRCh38 primary assembly genome FASTA。
- `data/raw/gencode.gtf`：已下载并解压，GENCODE v49 primary assembly annotation GTF。
- `data/raw/clinvar.vcf`：已下载并解压，使用 NCBI ClinVar GRCh38 archive_2.0 `clinvar_20260530.vcf.gz`，共 4,434,969 条变异记录。
- `data/raw/gtex_sqtl.tsv` 与 `data/raw/known_splice_events.tsv`：已通过 GTEx Portal API 拉取小型真实 case study。
- `data/raw/gencode.db`：已由 `gencode.gtf` 构建，用于真实 Pangolin 运行。
- RNA-FM 与 RNABERT：已安装 `multimolecule`，并成功下载/加载 `multimolecule/rnafm` 和 `multimolecule/rnabert` 预训练权重；已在 GPU 上完成前向验证。
- SpliceAI：已在 WSL Python 3.10 环境安装 `spliceai + tensorflow + pysam`，并对真实 ClinVar smoke VCF 输出结果。
- Pangolin：已安装 GitHub `tkzeng/Pangolin`，命令行可用，并已用真实 `genome.fa + gencode.db + ClinVar smoke CSV` 跑通 GPU smoke test。
- MMSplice：已安装 `mmsplice + kipoi + cyvcf2`，并修正为 `numpy==1.26.4` 以解决 ABI 兼容；已加载 H5 权重并完成 CPU smoke 预测。
- MaxEntScan：已从 `kepbod/maxentpy` 源码安装 `maxentpy`，并完成 donor/acceptor 最小评分。
- GPU：本机可用 `NVIDIA GeForce RTX 5070 Ti`。

真实模型 smoke 输出：

- `results/real_smoke/spliceai_clinvar_smoke.vcf`
- `results/real_smoke/spliceai_clinvar_smoke_summary.csv`
- `results/real_smoke/pangolin_clinvar_smoke.csv`
- `results/real_smoke/foundation_model_smoke.csv`
- `results/real_smoke/maxentscan_mmsplice_smoke.csv`

详细状态见 `REAL_RESOURCE_STATUS.md`，可重复执行脚本见 `python -m src.resources.run_real_model_smoke`。
"""

    text = f"""# C Part 执行细节说明

生成时间：{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

{real_resource_section}

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
python -m src.pipelines.run_c_part_all
```

也可以分步运行：

```powershell
python -m src.data.build_synthetic_splice_dataset
python -m src.experiments.exp2.run_multiscale
python -m src.data.build_variant_dataset
python -m src.experiments.exp3.run_variant_effect
python -m src.experiments.exp3.run_interpretability
python -m src.reports.write_c_part_report
```

## 2. 数据说明

本仓库当前包含两类数据资产。

第一类是真实资源，位于 `data/raw/`。其中 `genome.fa`、`gencode.gtf`、`clinvar.vcf`、`gtex_sqtl.tsv`、`known_splice_events.tsv` 和 `gencode.db` 已经用于真实模型 smoke test 与 case-study 支撑；由于原始 genome/annotation/ClinVar 文件体积过大，仍通过 `.gitignore` 排除，不随代码仓库提交。

第二类是可提交、可复现的实验数据，位于 `data/shared/processed/`、`data/shared/splits/` 与 `data/experiment_3/`。为了保证 C Part 主实验能够在没有外部大文件的机器上完整复现，代码同时制造了一个 synthetic splice benchmark：

- donor 样本：中心附近植入 canonical `GT` donor motif，并加入 `CAGGTAAGT`、下游 `GTA` 相关上下文。
- acceptor 样本：中心上游植入 canonical `AG` acceptor motif、polypyrimidine tract 和 branch-point-like `TACTAAC`。
- non-splice 样本：包含 easy random negative 与 hard GT/AG negative；hard negative 拥有局部 `GT` 或 `AG`，但缺少真实剪接上下文。
- 切分策略：按 synthetic chromosome 分为 train (`chr1-chr16`)、valid (`chr17-chr18`) 和 test/cross-gene (`chr19-chr22, chrX`)。
- 随机种子：数据制造使用 `2026`，模型训练默认使用 `42`。

数据文件：

- `data/shared/processed/splice_sites_pm50.csv`
- `data/shared/processed/splice_sites_pm100.csv`
- `data/shared/processed/splice_sites_pm200.csv`
- `data/shared/processed/splice_sites_pm400.csv`
- `data/shared/splits/train_pm*.csv`, `valid_pm*.csv`, `test_pm*.csv`
- `data/experiment_3/artificial_variant_effect.csv`

样本统计：

{markdown_table(summary, ["split", "label_name", "rows"], max_rows=20)}

人工变异统计：

{markdown_table(variants_summary, ["variant_type", "label_name", "rows"], max_rows=20)}

## 3. 模型与训练说明

当前 WSL 环境已经具备真实模型依赖与权重验证：RNA-FM、RNABERT、SpliceAI、Pangolin、MMSplice 和 MaxEntScan 均完成最小可执行 smoke test。考虑到完整 ClinVar/GTEx 规模训练和 foundation model fine-tuning 的运行成本较高，主实验表格仍保留可快速复现的同构 benchmark；真实模型 smoke 输出作为“环境与模型可运行性证明”和论文 case-study 支撑。

本次训练/打分的模型如下：

- `CNN motif baseline`：用 one-hot 序列位置特征训练多分类 Logistic Regression，模拟局部 motif baseline。
- `RNA-FM frozen k-mer + MLP`：用 3/4-mer 频率加剪接信号特征作为冻结表征，再训练线性 MLP 头的可运行代理。
- `RNABERT frozen token + MLP`：用 3/5-mer token 频率、中心位置 token 特征和剪接信号特征作为冻结表征代理。
- `SpliceAI signal proxy`：用 donor/acceptor 共识、polypyrimidine、motif density 等任务专用特征训练 RandomForest，并混合 deterministic splice-score。
- `RNA-FM zero-shot embedding distance`：比较 WT 与 Mut 的 k-mer/signal embedding L2 距离，不使用变异标签训练。
- `RNABERT zero-shot token distance`：比较 WT 与 Mut 的 token/position embedding L2 距离。
- `MaxEntScan consensus proxy`：使用 donor/acceptor 共识分数变化近似传统 splice-site strength delta。

真实模型替换方式：

- 放入真实 `genome.fa/gencode.gtf` 后，可先用 `python -m src.data.build_splice_site_dataset` 构建真实位点数据，再用本次新增实验脚本读取同名 CSV。
- 若将主实验从 synthetic benchmark 切换到真实全量数据，可在 `src/models/` 下替换代理模型的 `predict_proba` 或 zero-shot embedding 逻辑，结果表路径无需改变。

真实模型 smoke 结果：

RNA foundation model 前向验证：

{markdown_table(foundation_smoke, ["model", "device", "input_tokens", "hidden_size", "embedding_mean", "embedding_std"], max_rows=10)}

MaxEntScan 与 MMSplice 最小评分：

{markdown_table(maxent_mmsplice_smoke, ["model", "input", "score_name", "score"], max_rows=10)}

SpliceAI ClinVar smoke 摘要：

{markdown_table(spliceai_smoke, ["chrom", "pos", "ref", "alt", "spliceai_info"], max_rows=10)}

## 4. 实验二结果

实验二 A 的核心表为 `results/experiment_2/tables/experiment_2A_multiscale_context.csv`。按窗口和模型的主要结果如下：

{markdown_table(best_exp2, ["window_flank", "model", "accuracy", "macro_f1", "auroc", "auprc", "hard_negative_fpr"], max_rows=20)}

实验二 B 的 hard negative benchmark：

{markdown_table(exp2b, ["model", "test_easy_macro_f1", "test_hard_macro_f1", "cross_gene_macro_f1", "hard_negative_fpr"], max_rows=20)}

组织特异性 case study 输出 `results/experiment_2/tables/experiment_2C_tissue_splice_usage_case_study.csv`，用于展示不同 synthetic event 在 brain/heart/liver/muscle/blood 的 splice usage 差异。

## 5. 实验三结果

实验三主结果表为 `results/experiment_3/tables/experiment_3A_artificial_variant_metrics.csv`：

{markdown_table(best_exp3, ["model", "auroc", "auprc", "top_k", "top_k_recall", "enrichment_at_k", "variants"], max_rows=20)}

逐变异分数保存在：

- `results/experiment_3/tables/experiment_3A_artificial_variant_scores.csv`
- `data/experiment_3/artificial_variant_effect.csv`

## 6. 图件清单

本次生成的主要图件：

- `results/experiment_2/figures/exp2A_context_macro_f1.png`
- `results/experiment_2/figures/exp2A_context_auprc.png`
- `results/experiment_2/figures/exp2B_hard_negative_fpr.png`
- `results/experiment_2/figures/exp2C_tissue_splice_usage_heatmap.png`
- `results/experiment_2/figures/exp2C_junction_usage_case_study.png`
- `results/experiment_3/figures/exp3_variant_auroc.png`
- `results/experiment_3/figures/exp3_variant_auprc.png`
- `results/experiment_3/figures/exp3_delta_score_boxplot.png`
- `results/experiment_3/figures/exp3_saturation_mutagenesis_heatmap.png`
- `results/experiment_3/figures/exp3_saturation_mutagenesis_acceptor_heatmap.png`
- `results/experiment_3/figures/ism_donor_heatmap.png`
- `results/experiment_3/figures/ism_acceptor_heatmap.png`
- `results/experiment_3/figures/ism_hard_negative_heatmap.png`
- `results/experiment_3/figures/variant_delta_profile_donor_loss.png`
- `results/experiment_3/figures/variant_delta_profile_cryptic_gain.png`

## 7. 结论与可写入论文的要点

在 synthetic benchmark 上，任务专用的 `SpliceAI signal proxy` 和包含中心 token/剪接信号的 frozen representation 代理模型通常在 hard negative 与人工变异任务上更稳定；这符合 C Part 要回答的问题：仅靠局部 `GT/AG` motif 不足以解释剪接位点，模型需要利用上下文、polypyrimidine tract、共识序列强度和变异前后 delta profile。

需要在论文中明确的一点是：本次结果是离线 synthetic benchmark 的可复现实验产物，不应等同于真实 ClinVar/GTEx 上的生物医学结论。它的价值在于完整搭建了 C Part 的数据、训练、评价和解释性分析链路；当 A/B 部分或后续下载真实数据与预训练权重后，可以沿用这些脚本直接替换输入与模型实现。
"""

    out_path = root / "C_PART_EXECUTION_DETAILS.md"
    out_path.write_text(text, encoding="utf-8")
    reports_path = root / "reports/C_PART_EXECUTION_DETAILS.md"
    reports_path.write_text(text, encoding="utf-8")
    return out_path


def main() -> None:
    path = write_report()
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
