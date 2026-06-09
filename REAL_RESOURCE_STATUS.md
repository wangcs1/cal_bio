# Real Data and Model Resource Status

更新时间：2026-06-09

## 已补充到本地的真实数据

这些文件位于 `data/raw/`，目录仍被 `.gitignore` 忽略，不会推送到远端仓库。

| 文件 | 状态 | 说明 |
| --- | --- | --- |
| `data/raw/genome.fa` | 已下载并解压 | GENCODE v49 / GRCh38 primary assembly genome FASTA，3,151,417,447 bytes |
| `data/raw/gencode.gtf` | 已下载并解压 | GENCODE v49 primary assembly annotation GTF，3,326,907,803 bytes |
| `data/raw/clinvar.vcf` | 已下载并解压 | NCBI ClinVar GRCh38 archive_2.0 `clinvar_20260530.vcf.gz`，4,434,969 条变异记录 |
| `data/raw/gtex_sqtl.tsv` | 已生成 | GTEx v8 Portal API 小型 sQTL case study，177 条 variant 记录 |
| `data/raw/known_splice_events.tsv` | 已生成 | GTEx v8 Portal API 小型 splice event case study，25 条 event 记录 |
| `data/raw/gencode.db` | 已生成 | Pangolin/gffutils annotation database，由 `gencode.gtf` 生成 |
| `data/raw/clinvar_smoke.csv` | 已生成 | Pangolin smoke-test 输入，5 条 ClinVar SNV |
| `data/raw/clinvar_smoke.vcf` | 已生成 | SpliceAI smoke-test 输入，5 条 ClinVar SNV |

说明：NCBI 当前周 `clinvar_20260606.vcf.gz` 在本机解压后只有 header，因此改用官方归档 `clinvar_20260530.vcf.gz`。

## WSL 已安装并验证的真实模型环境

环境：WSL2 / Ubuntu-22.04 / Python 3.10 venv `/home/wcs/.venvs/cal_bio_cpart`。

| 资源 | 状态 | 说明 |
| --- | --- | --- |
| GPU | 已可用 | `NVIDIA GeForce RTX 5070 Ti`，`torch.cuda.is_available() == True` |
| PyTorch | 已可用 | WSL venv 中安装 `torch 2.12.0`，CUDA 可用 |
| TensorFlow | 已可用 | `tensorflow 2.21.0`；SpliceAI 使用 CPU smoke，避免 5070 Ti PTX JIT 兼容问题 |
| RNA-FM | 已可用 | `multimolecule/rnafm` 权重已下载并在 GPU 前向成功 |
| RNABERT | 已可用 | `multimolecule/rnabert` 权重已下载并在 GPU 前向成功 |
| SpliceAI | 已可用 | `spliceai 1.3.1 + pysam + tensorflow` 已安装，ClinVar smoke VCF 已输出 |
| Pangolin | 已可用 | GitHub `tkzeng/Pangolin` 已安装，CLI 能运行并使用 GPU |
| MMSplice | 已可用 | `mmsplice + kipoi + cyvcf2` 已安装；固定 `numpy==1.26.4` 解决 cyvcf2 ABI |
| MaxEntScan | 已可用 | `maxentpy` 已从 `kepbod/maxentpy` 源码安装并完成 donor/acceptor 评分 |
| Biopython / pyfaidx / pyfastx | 已可用 | 可处理 FASTA/GTF |

## 已新增脚本

| 脚本 | 用途 |
| --- | --- |
| `scripts/fetch_gtex_sqtl_cases.py` | 从 GTEx Portal API 拉取小型真实 sQTL / splice event case study |
| `scripts/check_real_resources.py` | 检查 raw 文件、Python 包、GPU 状态 |
| `scripts/make_clinvar_smoke.py` | 从真实 ClinVar VCF 生成 Pangolin CSV 和 SpliceAI VCF smoke-test 输入 |
| `scripts/run_real_model_smoke.py` | 一键运行 SpliceAI、Pangolin、RNA-FM、RNABERT、MaxEntScan、MMSplice smoke test |

## 验证命令

Windows PowerShell 调 WSL：

```powershell
wsl -d Ubuntu-22.04 -- bash -lc "source /home/wcs/.venvs/cal_bio_cpart/bin/activate && cd /mnt/d/CAL_BIO && python scripts/check_real_resources.py"
wsl -d Ubuntu-22.04 -- bash -lc "source /home/wcs/.venvs/cal_bio_cpart/bin/activate && cd /mnt/d/CAL_BIO && export HF_ENDPOINT=https://hf-mirror.com && python scripts/run_real_model_smoke.py"
wsl -d Ubuntu-22.04 -- bash -lc "source /home/wcs/.venvs/cal_bio_cpart/bin/activate && cd /mnt/d/CAL_BIO && python -m src.run_c_part_all"
```

在 WSL 终端中：

```bash
cd /mnt/d/CAL_BIO
source /home/wcs/.venvs/cal_bio_cpart/bin/activate
export HF_ENDPOINT=https://hf-mirror.com
python scripts/check_real_resources.py
python scripts/run_real_model_smoke.py
python -m src.run_c_part_all
```

## 已生成的真实模型结果

| 文件 | 说明 |
| --- | --- |
| `results/real_smoke/spliceai_clinvar_smoke.vcf` | SpliceAI 对 5 条真实 ClinVar SNV 的 VCF 输出 |
| `results/real_smoke/spliceai_clinvar_smoke_summary.csv` | SpliceAI INFO 字段摘要 |
| `results/real_smoke/pangolin_clinvar_smoke.csv` | Pangolin 对 5 条真实 ClinVar SNV 的输出 |
| `results/real_smoke/foundation_model_smoke.csv` | RNA-FM/RNABERT 真实权重前向验证，含 embedding 形状和统计 |
| `results/real_smoke/maxentscan_mmsplice_smoke.csv` | MaxEntScan donor/acceptor 分数与 MMSplice H5 权重预测分数 |

## 仍建议后续扩展

1. 将当前 synthetic benchmark 的主训练数据切到真实 `genome.fa + gencode.gtf` 生成的数据集上重跑。
2. 用 `multimolecule/rnafm` 与 `multimolecule/rnabert` 替换当前代理 embedding，实现真正的 frozen encoder + MLP 全量训练。
3. 若课程时间允许，追加更大规模 ClinVar/GTEx benchmark；当前 smoke test 证明依赖与接口已跑通，但不是全量真实生物医学结论。
