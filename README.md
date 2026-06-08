# 多尺度 RNA 剪接建模与异常剪接变异效应预测

英文题目：**Multi-scale RNA Splicing Modeling and Aberrant Splicing Variant Effect Prediction Based on RNA Foundation Models and Genomic Regulatory Models**

本项目围绕 RNA 剪接建模与异常剪接变异效应预测，评估不同类型序列模型是否真正学到了可迁移、可解释的剪接调控规则。项目不再只做普通的“RNA 模型调研”，而是将 RNA foundation model、剪接任务专用模型和长程基因组调控模型放在统一实验框架下比较，重点回答三个问题：

1. **剪接位点识别**：模型能否区分 donor site、acceptor site 和 non-splice site？
2. **多尺度上下文建模**：模型是否只依赖局部 GT/AG motif，还是能够利用更长程的 intron/exon context、hard negative motif、组织特异性 splice usage 等信息？
3. **异常剪接变异效应预测**：模型能否判断一个 SNV 是否会破坏正常剪接、激活 cryptic splice site，并进一步解释其可能影响的 splice junction 或组织特异性剪接模式？

项目最终希望形成一个从**剪接位点分类**到**多尺度上下文分析**再到**variant effect prediction**的完整实验链条，使课设不仅是模型跑分，而是更接近当前 AI for Science 中“从序列预测到分子机制解释”的研究范式。

---

## 0. 项目动机

高通量测序技术已经能够快速获得大规模基因组和转录组数据，但如何从原始核酸序列中解释调控规则、功能后果和疾病机制，仍然是计算生物学中的核心问题。RNA 剪接是连接基因组序列和成熟转录本功能的关键步骤，它不仅依赖 donor site、acceptor site、branch point、polypyrimidine tract 等局部信号，也受到外显子 / 内含子剪接增强子、剪接沉默子、RNA-binding proteins、组织特异性调控和长程上下文共同影响。

大量疾病相关变异并不直接改变蛋白编码序列，却可能通过破坏 canonical splice site、激活 cryptic splice site 或改变 alternative splicing pattern 影响基因功能。传统方法通常依赖局部 motif 或人工特征，难以系统建模长程上下文和组织特异性剪接调控。近年来，RNA-FM、RNABERT、RNAErnie 等 RNA foundation model 通过自监督学习从大量 RNA 序列中获得通用表示；SpliceAI、Pangolin 等剪接任务专用模型在 splice site 和 splice variant prediction 上表现突出；Borzoi、AlphaGenome 等长程基因组调控模型进一步尝试从 DNA 序列预测 RNA-seq coverage、splice junction、splice site usage、gene expression 等多种分子表型。

因此，本项目拟以 RNA 剪接为切入点，比较三类模型在剪接相关任务中的能力边界：通用 RNA foundation model 是否真正学到了可迁移的剪接语法？剪接任务专用模型相比 foundation model 的优势来自哪里？长程基因组调控模型能否把“是否影响剪接”的二分类问题推进到“影响哪个 splice junction、在哪个组织中影响更明显、产生什么分子后果”的机制解释层面？

---

## 1. 实验总览

| 实验 | 名称 | 目标 | 模型 | 必做程度 |
| --- | --- | --- | --- | --- |
| 实验一 | 剪接位点三分类 | 判断中心位置是 donor、acceptor 还是 non-splice | CNN、RNA-FM、RNABERT、SpliceAI | 必做 |
| 实验二 | 多尺度上下文与 hard negative 剪接建模 | 检验模型是否真正利用上下文，而非只识别 GT/AG motif | RNA-FM、RNABERT、SpliceAI、Pangolin、Borzoi/AlphaGenome case study | 必做 + 扩展 |
| 实验三 | 异常剪接变异效应预测 | 预测 SNV 是否造成 donor loss、acceptor loss、donor gain、acceptor gain | RNA-FM zero-shot、RNABERT zero-shot、SpliceAI、Pangolin、MMSplice、MaxEntScan、Borzoi/AlphaGenome case study | 加分重点 |
| 可解释性 | Attention / in silico mutagenesis / delta profile | 解释模型关注的剪接信号和突变敏感区域 | RNA-FM、RNABERT、CNN、SpliceAI | 建议做 |

---

## 2. 推荐项目结构

```text
cal_bio/
├── README.md
├── guidance.md
├── data/
│   ├── raw/
│   │   ├── genome.fa
│   │   ├── gencode.gtf
│   │   ├── clinvar.vcf
│   │   ├── gtex_sqtl.tsv                 # 可选：GTEx sQTL / splicing QTL
│   │   └── known_splice_events.tsv        # 可选：已知 alternative splicing events
│   ├── processed/
│   │   ├── splice_sites_pm50.csv
│   │   ├── splice_sites_pm100.csv
│   │   ├── splice_sites_pm200.csv
│   │   ├── splice_sites_pm400.csv
│   │   ├── splice_sites_pm1000.csv
│   │   ├── hard_negative_gtag.csv
│   │   ├── artificial_variant_effect.csv
│   │   ├── clinvar_splicing_variants.csv
│   │   └── tissue_splice_usage_cases.csv
│   └── splits/
│       ├── train.csv
│       ├── valid.csv
│       ├── test.csv
│       └── cross_gene_test.csv
├── src/
│   ├── build_splice_site_dataset.py
│   ├── build_hard_negative_dataset.py
│   ├── build_variant_dataset.py
│   ├── build_sqtl_dataset.py
│   ├── models/
│   │   ├── cnn.py
│   │   ├── rnafm_mlp.py
│   │   ├── rnabert_mlp.py
│   │   ├── spliceai_wrapper.py
│   │   ├── pangolin_wrapper.py
│   │   ├── mmsplice_wrapper.py
│   │   ├── maxentscan_wrapper.py
│   │   ├── borzoi_wrapper.py
│   │   └── alphagenome_case_study.md
│   ├── train.py
│   ├── evaluate.py
│   ├── run_exp1_classification.py
│   ├── run_exp2_multiscale.py
│   ├── run_exp3_variant_effect.py
│   ├── run_zero_shot_scoring.py
│   ├── run_interpretability.py
│   └── utils.py
├── configs/
│   ├── cnn_pm200.yaml
│   ├── rnafm_pm200.yaml
│   ├── rnabert_pm200.yaml
│   ├── exp2_multiscale.yaml
│   ├── exp2_hard_negative.yaml
│   ├── exp3_variant_effect.yaml
│   └── interpretability.yaml
├── results/
│   ├── tables/
│   ├── figures/
│   └── checkpoints/
└── reports/
    ├── experiment_1.md
    ├── experiment_2.md
    ├── experiment_3.md
    ├── model_cards.md
    └── experiment_log.md
```

最小可完成版本只需要实现：

```text
data/processed/splice_sites_pm200.csv
src/models/cnn.py
src/models/rnafm_mlp.py
src/models/rnabert_mlp.py
src/train.py
src/evaluate.py
src/run_exp1_classification.py
src/run_exp2_multiscale.py
src/run_exp3_variant_effect.py
results/
```

高级扩展可以逐步补充 Pangolin、MMSplice、MaxEntScan、Borzoi / AlphaGenome case study。

---

## 3. 环境准备

### 3.1 Python 环境

建议使用 Python 3.9 或 3.10。

```bash
conda create -n splice-rna python=3.10
conda activate splice-rna
```

### 3.2 基础依赖

```bash
pip install numpy pandas scikit-learn matplotlib seaborn tqdm biopython pyfaidx pyyaml
pip install torch torchvision torchaudio
pip install transformers datasets accelerate einops
```

### 3.3 模型依赖

```bash
# RNA-FM
pip install fair-esm

# RNABERT / multimolecule 模型
pip install multimolecule

# SpliceAI
pip install spliceai

# 可选：MMSplice / Kipoi 生态可能依赖较复杂，建议作为扩展
pip install mmsplice kipoi kipoi_veff || true
```

Pangolin、Borzoi、AlphaGenome 可能需要额外环境或官方接口。课程项目中可以采用两种方式：

1. **实跑方式**：安装官方代码或可用 API，对小规模 case study 跑预测。
2. **论文式 case study 方式**：将其作为长程调控模型扩展模块，描述输入、输出、预期结果和与 SpliceAI / RNA-FM 的差异，若无法运行则不作为主实验指标。

---

## 4. 数据准备

### 4.1 原始数据

| 文件 | 用途 | 示例来源 |
| --- | --- | --- |
| `genome.fa` | 提取候选位点和变异附近序列 | GRCh38 / hg38 |
| `gencode.gtf` | 提取 exon-intron 边界 | GENCODE |
| `clinvar.vcf` | 构造真实剪接相关变异 | ClinVar |
| `gtex_sqtl.tsv` | 构造组织特异性 sQTL / splicing QTL 扩展实验 | GTEx，可选 |
| `known_splice_events.tsv` | 组织特异性 alternative splicing case study | 可选 |

实验一和实验二只需要 `genome.fa` 与 `gencode.gtf`。实验三的最小版本可以使用人工构造突变，不依赖 ClinVar。

### 4.2 样本类别定义

| 标签 | 类别 | 生物学含义 |
| --- | --- | --- |
| `0` | donor | 5' splice site，通常位于 exon → intron 边界 |
| `1` | acceptor | 3' splice site，通常位于 intron → exon 边界 |
| `2` | non_splice | 含 GT / AG motif 但不是注释剪接位点的负样本 |

负样本分为两类：

| 负样本类型 | 构造方式 | 作用 |
| --- | --- | --- |
| easy negative | 从非剪接区域随机采样 | sanity check |
| hard negative | 从同一基因或相邻区域采样 GT/AG motif，但排除注释剪接位点 | 检验模型是否学习上下文，而不是只识别 GT/AG |

### 4.3 序列方向

构造样本时必须统一转录方向：

1. 正链基因按基因组序列直接提取。
2. 负链基因取反向互补序列。
3. 所有样本最终都应以转录方向保存。
4. 输入 RNA-FM / RNABERT 时，将 DNA 序列中的 `T` 替换为 `U`。
5. 输入 SpliceAI / Pangolin / Borzoi / AlphaGenome 时，通常保留 DNA 字母 `A/C/G/T`。

### 4.4 数据切分

推荐按染色体或 gene_id 切分，避免同一基因的相似窗口同时出现在训练集和测试集。

```text
train: chr1-chr16
valid: chr17-chr18
test : chr19-chr22, chrX
```

另外建议构造一个更严格的测试集：

```text
cross_gene_test.csv
```

该测试集要求测试基因与训练基因完全不重叠，用于检验模型是否真正泛化到新基因。

---

## 5. 实验一：剪接位点三分类

### 5.1 实验目的

验证 RNA foundation model 的 embedding 是否包含剪接位点识别所需的序列上下文信息。

核心问题：

> RNA-FM / RNABERT + MLP 是否优于普通 CNN baseline？它们与任务专用模型 SpliceAI 之间还有多大差距？

### 5.2 输入与输出

默认输入窗口：

```text
中心位点 ±200 nt，总长度 401 nt
```

输出：

```text
P(donor), P(acceptor), P(non_splice)
```

### 5.3 数据构造步骤

1. 从 `gencode.gtf` 中读取转录本和 exon 坐标。
2. 对每个转录本按转录方向排序 exon。
3. 对相邻 exon 之间的 intron 边界提取 donor 和 acceptor。
4. 从 `genome.fa` 提取中心 ±200 nt 的窗口。
5. 从同一基因附近采样 GT / AG motif 作为 non-splice hard negative。
6. 过滤 `N` 比例大于 5% 的序列。
7. 保证 donor:acceptor:non_splice 大致为 1:1:1。
8. 按染色体或 gene_id 切分 train / valid / test。

输出 CSV：

```text
sample_id,chrom,start,end,strand,center,label,label_name,negative_type,sequence,gene_id,transcript_id
```

### 5.4 模型设置

#### CNN baseline

```text
one-hot sequence
→ Conv1D kernel sizes 3/5/9
→ ReLU
→ MaxPool
→ Dropout
→ Linear
→ 3-class softmax
```

#### RNA-FM + MLP

```text
RNA sequence
→ RNA-FM frozen encoder
→ CLS embedding 或 mean pooling
→ MLP classifier
→ 3-class softmax
```

优先冻结 RNA-FM，只训练 MLP head。若时间充足，可以加一个 fine-tuning 版本：只微调最后 1-2 层。

#### RNABERT + MLP

与 RNA-FM 保持同样数据、同样窗口长度、同样 pooling、同样训练轮数，保证比较公平。

#### SpliceAI baseline

SpliceAI 输出每个位置作为 donor / acceptor 的概率。可将中心附近的最大 donor / acceptor 分数转换为三分类：

```text
donor_score = max donor probability near center
acceptor_score = max acceptor probability near center
non_splice_score = 1 - max(donor_score, acceptor_score)
```

### 5.5 训练配置

```text
batch_size: 16 / 32
epochs: 10-20
optimizer: AdamW
learning_rate:
  CNN: 1e-3
  MLP head: 1e-4 或 5e-4
weight_decay: 1e-4
early_stopping: valid macro-F1 连续 3 轮不提升则停止
loss: CrossEntropyLoss
```

### 5.6 评价指标

```text
Accuracy
Macro Precision
Macro Recall
Macro F1-score
AUROC
AUPRC
Donor-F1
Acceptor-F1
Non-splice-F1
Confusion Matrix
```

推荐结果表：

```text
results/tables/experiment_1_metrics.csv
```

| Model | Accuracy | Macro-F1 | AUROC | AUPRC | Donor-F1 | Acceptor-F1 | Non-splice-F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CNN |  |  |  |  |  |  |  |
| RNABERT + MLP |  |  |  |  |  |  |  |
| RNA-FM + MLP |  |  |  |  |  |  |  |
| SpliceAI |  |  |  |  |  |  |  |

---

## 6. 实验二：多尺度上下文与 hard negative 剪接建模

实验二是本项目升级后的重点。它不再只是普通“窗口长度消融”，而是从**局部 motif → 中程 exon/intron context → 长程 splice regulation → 组织特异性 splice usage**逐层分析模型能力。

### 6.1 实验目的

回答：

> 剪接位点识别究竟主要依赖局部 GT/AG motif，还是需要更长程的 intron/exon context、剪接增强子 / 沉默子、组织特异性调控和 splice junction topology？

### 6.2 实验 2A：Multi-scale context ablation

#### 输入长度设置

```text
Local context      : ±50 nt    -> 101 nt
Short-range context: ±100 nt   -> 201 nt
Regional context   : ±200 nt   -> 401 nt
Extended context   : ±400 nt   -> 801 nt
Gene-fragment      : ±1000 nt  -> 2001 nt，可选
Long-range context : 5 kb / 10 kb，可选，仅用于 SpliceAI / Pangolin / Borzoi / AlphaGenome
```

#### 模型分工

| 模型 | 适合上下文长度 | 作用 |
| --- | --- | --- |
| CNN | ±50 到 ±400 | 局部 motif 和普通深度学习基线 |
| RNA-FM | ±50 到 ±400，±1000 可选 | 检验 RNA foundation representation 的上下文利用能力 |
| RNABERT | ±50 到 ±400 | RNA foundation model 对照 |
| SpliceAI | ±200 到 10 kb | 剪接任务专用强基线 |
| Pangolin | 5 kb / 10 kb 或官方推荐输入 | 多组织 splice site strength / variant effect 扩展 |
| Borzoi / AlphaGenome | 长序列 case study | 长程 RNA-seq / splice junction / splice usage 分子后果预测 |

#### 数据要求

1. 不同窗口长度使用同一批中心位点。
2. train / valid / test 划分完全一致。
3. 长窗口越界时统一 padding 或统一过滤。
4. 不同模型的报告中要注明最大输入长度和截断策略。

#### 输出表

```text
results/tables/experiment_2A_multiscale_context.csv
```

| Model | ±50 | ±100 | ±200 | ±400 | ±1000 | 5kb/10kb |
| --- | --- | --- | --- | --- | --- | --- |
| CNN Macro-F1 |  |  |  |  |  | - |
| RNABERT Macro-F1 |  |  |  |  |  | - |
| RNA-FM Macro-F1 |  |  |  |  |  | - |
| SpliceAI Macro-F1 |  |  |  |  |  |  |
| Pangolin score | - | - | - | - |  |  |

#### 图表

```text
results/figures/exp2A_context_macro_f1.png
results/figures/exp2A_context_auprc.png
```

### 6.3 实验 2B：GT/AG hard negative benchmark

#### 设计动机

很多非剪接位置也包含 GT 或 AG motif。若负样本是普通随机序列，模型可能只需要识别 GT/AG 即可取得高分，无法证明其学到了真正的剪接上下文。

#### 测试集设计

| 测试集 | 负样本构造 | 难度 | 目的 |
| --- | --- | --- | --- |
| Test-Easy | 随机非剪接位置 | 低 | sanity check |
| Test-Hard-GT/AG | 含 GT/AG 但非注释剪接位点 | 高 | 检验上下文建模能力 |
| Test-Cross-Gene | 测试基因不出现在训练集 | 高 | 检验跨基因泛化能力 |
| Test-Rare-Motif | GC-AG 等少见剪接信号，可选 | 很高 | 检验模型是否过度依赖 canonical GT-AG |

#### 评价方式

对每个模型分别在多个测试集上报告：

```text
Macro-F1
AUPRC
Donor-F1
Acceptor-F1
Hard-negative false positive rate
```

Hard-negative false positive rate 定义为：

```text
含 GT/AG 的非剪接负样本中，被模型误判为 donor 或 acceptor 的比例
```

#### 输出表

```text
results/tables/experiment_2B_hard_negative.csv
```

| Model | Test-Easy Macro-F1 | Test-Hard Macro-F1 | Cross-Gene Macro-F1 | Hard-Neg FPR |
| --- | --- | --- | --- | --- |
| CNN |  |  |  |  |
| RNABERT + MLP |  |  |  |  |
| RNA-FM + MLP |  |  |  |  |
| SpliceAI |  |  |  |  |
| Pangolin |  |  |  |  |

### 6.4 实验 2C：组织特异性 splice usage case study

这是实验二的前沿扩展，不要求全量跑，但建议至少设计并完成小规模 case study。

#### 实验问题

> 同一个剪接位点或 splice junction 在不同组织中使用强度是否不同？长程调控模型是否能给出组织特异性预测？

#### 候选模型

```text
Pangolin
Borzoi
AlphaGenome
```

#### 候选组织

```text
brain
heart
liver
muscle
blood
```

#### 输入输出

输入：

```text
包含目标 exon-intron 结构的 genomic sequence，长度可为 10 kb 或模型推荐长度
```

输出：

```text
predicted splice site strength
predicted splice junction usage
predicted RNA-seq coverage
predicted tissue-specific Δusage
```

#### 结果展示

```text
results/figures/exp2C_tissue_splice_usage_heatmap.png
results/figures/exp2C_junction_usage_case_study.png
```

即使只做 3-5 个已知 alternative splicing event，也能作为“高级扩展实验”展示项目深度。

### 6.5 实验二预期解释

可能出现以下结果：

1. CNN 在 easy negative 上表现很好，但在 hard negative 上掉分明显，说明它主要依赖局部 motif。
2. RNA-FM / RNABERT 在 hard negative 上比 CNN 稳定，说明 foundation representation 学到了一定上下文规律。
3. SpliceAI / Pangolin 在 hard negative 和长窗口上表现更强，说明任务专用监督信号和长程上下文仍然重要。
4. Borzoi / AlphaGenome 的优势不一定体现在三分类 accuracy，而更体现在 splice junction、RNA-seq coverage 和组织特异性分子后果解释上。

---

## 7. 实验三：异常剪接变异效应预测

实验三将任务从“识别剪接位点”升级为“预测变异是否扰动剪接”。这是项目最能体现前沿性的部分。

### 7.1 实验目的

回答：

> 给定一个 SNV，模型能否预测它是否造成 donor loss、acceptor loss、donor gain 或 acceptor gain，并进一步对疾病相关 splice-altering variant 做优先级排序？

### 7.2 模型分组

| 类型 | 模型 | 作用 |
| --- | --- | --- |
| Foundation model zero-shot | RNA-FM、RNABERT、Nucleotide Transformer，可选 | 不重新训练或少量训练，考察预训练表示对变异扰动的敏感性 |
| 任务专用监督模型 | SpliceAI、Pangolin、MMSplice | 强基线，直接预测 splice variant effect |
| 传统统计模型 | MaxEntScan | 传统 splice site strength 基线 |
| 长程调控模型 | Borzoi、AlphaGenome | 预测 RNA-seq、splice junction、splice usage 等分子后果 |

### 7.3 数据方案 A：人工饱和突变

这是最建议完成的版本，标签清晰，可解释性强。

#### 构造方式

围绕真实 donor / acceptor 位点进行 saturation mutagenesis：

```text
对中心 ±20 nt 或 ±50 nt 的每个位置，分别替换为另外 3 种碱基。
```

同时构造三类重点变异：

```text
1. donor loss:
   将真实 donor 附近 GT 改为 AT / GC / TT

2. acceptor loss:
   将真实 acceptor 附近 AG 改为 AA / AC / GG

3. cryptic splice gain:
   在非剪接区域人为制造新的 GT / AG motif
```

#### 输出字段

```text
variant_id,chrom,pos,strand,ref,alt,variant_type,label,wt_sequence,mut_sequence,target_class,gene_id,transcript_id
```

### 7.4 数据方案 B：ClinVar 剪接相关变异

作为真实变异 benchmark 或扩展实验。

#### 阳性样本

ClinVar 中包含以下 consequence 或关键词的变异：

```text
splice donor variant
splice acceptor variant
splice region variant
splice-altering
pathogenic / likely pathogenic with splicing evidence
```

#### 阴性样本

```text
benign / likely benign variants
matched by chromosome / gene / distance to splice site / allele frequency if available
```

#### 任务

```text
positive: splice-altering variant
negative: likely neutral variant
```

### 7.5 数据方案 C：GTEx sQTL / tissue-specific splicing 扩展

该部分适合作为更高阶 case study。

任务：

```text
区分 sQTL variant 与 matched control variant，并分析其组织特异性 splice junction usage 变化。
```

输出：

```text
variant_id, tissue, target_gene, target_junction, observed_effect_direction, model_delta_score
```

### 7.6 打分方式

#### 7.6.1 Supervised classifier delta score

对实验一训练好的 CNN / RNA-FM / RNABERT 分类器：

```text
WT sequence  -> classifier -> P_wt(donor), P_wt(acceptor), P_wt(non_splice)
Mut sequence -> classifier -> P_mut(donor), P_mut(acceptor), P_mut(non_splice)
```

Donor loss：

```text
delta_donor_loss = P_wt(donor) - P_mut(donor)
```

Acceptor loss：

```text
delta_acceptor_loss = P_wt(acceptor) - P_mut(acceptor)
```

Cryptic gain：

```text
delta_splice_gain = max(P_mut(donor), P_mut(acceptor)) - max(P_wt(donor), P_wt(acceptor))
```

#### 7.6.2 Foundation model zero-shot scoring

不使用任务标签，直接用预训练模型衡量 WT 与 Mut 的表示差异。

可选分数：

```text
Delta embedding distance:
score = || embedding_wt - embedding_mut ||_2

Delta cosine distance:
score = 1 - cos(embedding_wt, embedding_mut)

Delta pseudo-likelihood:
score = PLL(mutant sequence) - PLL(wild-type sequence)
```

其中 Delta pseudo-likelihood 更适合 masked language model。若实现困难，可以优先使用 embedding distance。

#### 7.6.3 SpliceAI / Pangolin / MMSplice / MaxEntScan score

SpliceAI：

```text
score = max(DS_AG, DS_AL, DS_DG, DS_DL)
```

其中：

```text
DS_AG: acceptor gain
DS_AL: acceptor loss
DS_DG: donor gain
DS_DL: donor loss
```

Pangolin：

```text
score = max tissue-specific splice effect score
```

MMSplice：

```text
score = predicted ΔlogitPSI 或 splice impact score
```

MaxEntScan：

```text
score = splice_site_strength_wt - splice_site_strength_mut
```

#### 7.6.4 Borzoi / AlphaGenome molecular consequence score

用于 case study：

```text
WT sequence  -> predicted RNA-seq / splice junction / splice site usage
Mut sequence -> predicted RNA-seq / splice junction / splice site usage
Delta prediction -> molecular consequence
```

展示重点不只是“是否有害”，而是：

```text
影响哪个 junction？
影响哪个 tissue？
影响 RNA-seq coverage 还是 splice site usage？
是否伴随 gene expression 改变？
```

### 7.7 评价指标

```text
AUROC
AUPRC
Top-k recall
Enrichment@K
Delta score distribution
Calibration curve，可选
```

Top-k recall：按照模型预测变异影响分数排序，前 k 个变异中覆盖多少真实 splice-altering variant。

Enrichment@K：前 k 个高分变异中阳性比例相对于随机抽样的富集倍数。

### 7.8 推荐输出

```text
results/tables/experiment_3A_artificial_variant_metrics.csv
results/tables/experiment_3B_clinvar_variant_metrics.csv
results/tables/experiment_3C_sqtl_case_study.csv
results/figures/exp3_delta_score_boxplot.png
results/figures/exp3_variant_auroc.png
results/figures/exp3_variant_auprc.png
results/figures/exp3_saturation_mutagenesis_heatmap.png
results/figures/exp3_alphagenome_junction_case_study.png
```

### 7.9 实验三预期解释

可能结论：

1. RNA-FM / RNABERT zero-shot embedding distance 对关键 splice motif 突变有一定响应，但不一定能准确区分所有真实 ClinVar splicing variants。
2. SpliceAI / Pangolin 在真实 splice-altering variant benchmark 上更稳定，说明任务监督和长程上下文仍然关键。
3. MaxEntScan 在 canonical splice site 附近可能有效，但对 cryptic gain、深内含子变异、组织特异性剪接解释能力有限。
4. Borzoi / AlphaGenome 的优势在于提供分子后果层面的解释，例如 splice junction usage、RNA-seq coverage 和 tissue-specific effect。

---

## 8. 可解释性分析

### 8.1 Attention 可视化

适用于 RNA-FM / RNABERT。

步骤：

1. 选择预测正确且置信度高的 donor、acceptor、hard negative 样本。
2. 提取最后几层 attention map。
3. 对 attention head 求平均。
4. 绘制中心 ±50 nt 的 attention heatmap。

重点观察：

```text
donor site: 中心 GT 附近是否高亮
acceptor site: 中心 AG 附近是否高亮
acceptor 上游 polypyrimidine tract 是否有响应
hard negative: 模型是否能抑制非功能性 GT/AG motif
```

### 8.2 In silico mutagenesis

更推荐，因为它直观、可解释。

步骤：

1. 选择一个 donor 或 acceptor 样本。
2. 对序列每个位置依次突变为 A / C / G / T。
3. 重新预测目标类别概率。
4. 记录概率下降幅度。
5. 绘制位置 × 碱基的 heatmap。

Donor importance：

```text
importance = P_original(donor) - P_mutated(donor)
```

Acceptor importance：

```text
importance = P_original(acceptor) - P_mutated(acceptor)
```

输出：

```text
results/figures/ism_donor_heatmap.png
results/figures/ism_acceptor_heatmap.png
results/figures/ism_hard_negative_heatmap.png
```

### 8.3 Variant delta profile

对一个 ClinVar 或人工变异 case，画出 WT 与 Mut 在中心附近每个位置的 donor / acceptor score 曲线。

输出：

```text
results/figures/variant_delta_profile_donor_loss.png
results/figures/variant_delta_profile_cryptic_gain.png
```

该图能展示变异是否导致原有 splice site score 下降，或在新位置产生 cryptic splice site peak。

---

## 9. 推荐运行顺序

### Step 1：构建剪接位点数据集

产物：

```text
data/processed/splice_sites_pm200.csv
data/splits/train.csv
data/splits/valid.csv
data/splits/test.csv
```

检查：

```text
三类样本数量是否接近均衡
序列长度是否全部为 401
正负链是否正确
负样本是否包含 hard negative GT/AG
是否避免 train/test gene leakage
```

### Step 2：跑 CNN baseline

产物：

```text
results/tables/cnn_pm200_metrics.csv
results/figures/cnn_pm200_confusion_matrix.png
```

CNN 是 sanity check。若 CNN 无法学习，优先检查数据标签、序列方向和中心位点定义。

### Step 3：跑 RNA-FM + MLP

产物：

```text
results/tables/rnafm_pm200_metrics.csv
results/checkpoints/rnafm_mlp_best.pt
```

### Step 4：跑 RNABERT + MLP

产物：

```text
results/tables/rnabert_pm200_metrics.csv
results/checkpoints/rnabert_mlp_best.pt
```

### Step 5：补 SpliceAI baseline

产物：

```text
results/tables/spliceai_pm200_metrics.csv
```

### Step 6：实验二 Multi-scale context + hard negative

产物：

```text
results/tables/experiment_2A_multiscale_context.csv
results/tables/experiment_2B_hard_negative.csv
results/figures/exp2A_context_macro_f1.png
results/figures/exp2B_hard_negative_fpr.png
```

### Step 7：实验三 variant effect prediction

优先做人工饱和突变，再做 ClinVar 扩展。

产物：

```text
data/processed/artificial_variant_effect.csv
results/tables/experiment_3A_artificial_variant_metrics.csv
results/figures/exp3_delta_score_boxplot.png
results/figures/exp3_saturation_mutagenesis_heatmap.png
```

### Step 8：可解释性分析

产物：

```text
results/figures/ism_donor_heatmap.png
results/figures/ism_acceptor_heatmap.png
results/figures/variant_delta_profile_cryptic_gain.png
```

### Step 9：前沿 case study，可选

选择 1-3 个 case，用 Pangolin / Borzoi / AlphaGenome 展示组织特异性或长程分子后果预测。

产物：

```text
reports/alphagenome_borzoi_case_study.md
results/figures/exp2C_tissue_splice_usage_heatmap.png
results/figures/exp3_alphagenome_junction_case_study.png
```

---

## 10. 实验记录模板

每次训练建议记录：

```text
experiment_id:
date:
model:
task:
window_size:
negative_type:
dataset_version:
train_samples:
valid_samples:
test_samples:
learning_rate:
batch_size:
epochs:
freeze_encoder:
best_valid_macro_f1:
test_accuracy:
test_macro_f1:
test_auroc:
test_auprc:
notes:
```

---

## 11. 最终论文 / 汇报结构建议

```text
1. 引言
   RNA 剪接、异常剪接变异、RNA foundation model、AI for Science 背景。

2. 研究问题
   通用 RNA foundation model 是否真正学到了可迁移剪接语法？
   任务专用剪接模型和长程基因组调控模型分别解决什么问题？

3. 数据与方法
   GENCODE 剪接位点构造、hard negative、ClinVar / 人工突变、模型设置、评价指标。

4. 实验一：剪接位点三分类
   模型总体性能、混淆矩阵、donor / acceptor / non-splice 区分能力。

5. 实验二：多尺度上下文与 hard negative
   Multi-scale context ablation、GT/AG hard negative、cross-gene generalization、组织特异性 case study。

6. 实验三：异常剪接变异效应预测
   Artificial saturation mutagenesis、zero-shot variant scoring、ClinVar benchmark、SpliceAI / Pangolin 对照。

7. 可解释性分析
   Attention、in silico mutagenesis、variant delta profile。

8. 讨论
   RNA-FM / RNABERT 的迁移潜力，SpliceAI / Pangolin 的优势，Borzoi / AlphaGenome 的分子后果解释价值。

9. 结论
   总结多尺度 RNA 剪接建模的能力边界和后续改进方向。
```

---

## 12. 建议展示图表

```text
图 1：项目整体实验流程图
图 2：剪接位点三分类任务示意图
图 3：RNA foundation model、SpliceAI、长程调控模型的层次对比图
图 4：实验一模型性能柱状图
图 5：donor / acceptor / non-splice 混淆矩阵
图 6：Multi-scale context ablation 折线图
图 7：GT/AG hard negative false positive rate 对比图
图 8：人工饱和突变 heatmap
图 9：WT vs Mut delta score 箱线图
图 10：ClinVar variant prioritization AUROC / AUPRC
图 11：tissue-specific splice usage case study heatmap
图 12：in silico mutagenesis 可解释性热图
```

---

## 13. 最小可交付版本与高级版本

### 13.1 最小可交付版本

```text
模型：CNN、RNA-FM、RNABERT、SpliceAI
实验：
1. 剪接位点三分类
2. ±50 / ±100 / ±200 / ±400 多尺度上下文消融
3. GT/AG hard negative 测试
4. 人工 donor loss / acceptor loss / cryptic gain 变异效应预测
5. in silico mutagenesis 可解释性
```

### 13.2 高级展示版本

```text
额外模型：Pangolin、MMSplice、MaxEntScan、Borzoi / AlphaGenome case study
额外实验：
1. ClinVar splice-altering variant prioritization
2. GTEx sQTL / tissue-specific splicing case study
3. Borzoi / AlphaGenome 的 splice junction / RNA-seq coverage 分子后果预测
4. RNA-FM / RNABERT zero-shot pseudo-likelihood variant scoring
```

建议实际推进时先完成最小版本，再把高级版本作为论文扩展和汇报亮点。

---

## 14. 建议核心结论模板

最终论文可以围绕下面这种结论组织：

> 本项目表明，RNA-FM 和 RNABERT 等 RNA foundation model 能够通过预训练表示捕捉一定的剪接相关序列语法，在剪接位点三分类和 hard negative 测试中相比普通 CNN 具有更好的泛化潜力。然而，在异常剪接变异效应预测中，SpliceAI、Pangolin 等任务专用模型仍然更加稳定，说明长程上下文、剪接监督信号和组织特异性建模对于疾病相关变异解释至关重要。Borzoi、AlphaGenome 等长程基因组调控模型的价值不只在于提高分类分数，而在于进一步提供 splice junction、RNA-seq coverage 和 tissue-specific splice usage 等分子后果解释。整体来看，RNA foundation model 是转录组学任务的重要通用表征基础，但要走向可靠的剪接变异解释，还需要与任务监督、长程调控建模和可解释性分析结合。
