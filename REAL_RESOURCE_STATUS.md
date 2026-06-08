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

说明：NCBI 当前周 `clinvar_20260606.vcf.gz` 在本机解压后只有 header，因此改用官方归档 `clinvar_20260530.vcf.gz`。

## 已安装或验证的真实模型环境

| 资源 | 状态 | 说明 |
| --- | --- | --- |
| GPU | 已可用 | `NVIDIA GeForce RTX 5070 Ti`，`torch.cuda.is_available() == True` |
| PyTorch | 已可用 | `torch 2.9.0+cu130` |
| RNA-FM | 已可用 | `multimolecule/rnafm` 权重已从 Hugging Face 下载并可加载 |
| RNABERT | 已可用 | `multimolecule/rnabert` 权重已从 Hugging Face 下载并可加载 |
| Pangolin | 已可用 | GitHub `tkzeng/Pangolin` 已安装，CLI 能运行并使用 GPU |
| Biopython / pyfaidx | 已可用 | 可处理 FASTA/GTF |
| SpliceAI | 未完成 | Windows + Python 3.12 下 `pysam` 构建失败，建议用 WSL/Linux/conda Linux 环境 |
| MMSplice | 未完成 | Windows + Python 3.12 下 `cyvcf2` 构建失败，建议用 WSL/Linux/conda Linux 环境 |
| MaxEntPy | 未完成 | 当前 pip 源没有可用 `maxentpy` 包；项目内仍保留 MaxEntScan consensus proxy |

## 已新增脚本

| 脚本 | 用途 |
| --- | --- |
| `scripts/fetch_gtex_sqtl_cases.py` | 从 GTEx Portal API 拉取小型真实 sQTL / splice event case study |
| `scripts/check_real_resources.py` | 检查 raw 文件、Python 包、GPU 状态 |
| `scripts/make_clinvar_smoke.py` | 从真实 ClinVar VCF 生成 Pangolin smoke-test CSV |

## 验证命令

```powershell
python scripts\fetch_gtex_sqtl_cases.py
python scripts\check_real_resources.py
python scripts\make_clinvar_smoke.py
pangolin data\raw\clinvar_smoke.csv data\raw\genome.fa data\raw\gencode.db results\real_smoke\pangolin_clinvar_smoke -c CHROM,POS,REF,ALT
```

Pangolin smoke test 已运行成功，并输出 `results/real_smoke/pangolin_clinvar_smoke.csv`。

## 仍建议后续补充

1. 用 WSL/Linux 环境安装并验证 SpliceAI 与 MMSplice。
2. 将当前 synthetic benchmark 的实验脚本切到真实 `data/raw/genome.fa + data/raw/gencode.gtf` 生成的新数据集上重跑。
3. 用 `multimolecule/rnafm` 与 `multimolecule/rnabert` 替换当前代理 embedding，实现真正的 frozen encoder + MLP。

