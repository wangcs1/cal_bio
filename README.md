# 基于 RNA 序列基础模型的剪接位点识别与异常剪接变异效应预测

英文题目：**Splice Site Recognition and Aberrant Splicing Variant Effect Prediction Based on RNA Foundation Models**

本项目围绕 RNA 剪接位点预测任务，评估 RNA 序列基础模型是否真正学到了可迁移的剪接信号。核心思路是将 RNA-FM、RNABERT 等通用 RNA foundation model 与剪接任务专用模型 SpliceAI、轻量 CNN baseline 放在统一实验框架下比较，回答两个问题：

1. 模型能不能识别正常 donor / acceptor 剪接位点？
2. 模型能不能感知突变对剪接信号的破坏？

建议将实验组织为三个部分：

| 实验 | 目标 | 必做程度 |
| --- | --- | --- |
| 实验一：剪接位点三分类 | 判断中心位置是 donor site、acceptor site 还是 non-splice site | 必做 |
| 实验二：上下文长度消融 | 比较不同窗口长度对模型性能的影响 | 必做 |
| 实验三：异常剪接变异效应预测 | 比较 wild-type 与 mutant 序列的剪接分数变化 | 加分项 |

---

## 1. 推荐项目结构

当前仓库只有 `guidance.md`，建议后续按下面结构补充代码和结果，便于复现实验和撰写论文。

```text
cal_bio/
├── README.md
├── guidance.md
├── data/
│   ├── raw/
│   │   ├── genome.fa
│   │   ├── gencode.gtf
│   │   └── clinvar.vcf
│   ├── processed/
│   │   ├── splice_sites_pm200.csv
│   │   ├── splice_sites_pm50.csv
│   │   ├── splice_sites_pm100.csv
│   │   ├── splice_sites_pm400.csv
│   │   └── variant_effect.csv
│   └── splits/
│       ├── train.csv
│       ├── valid.csv
│       └── test.csv
├── src/
│   ├── build_splice_site_dataset.py
│   ├── build_variant_dataset.py
│   ├── models/
│   │   ├── cnn.py
│   │   ├── rnafm_mlp.py
│   │   ├── rnabert_mlp.py
│   │   └── spliceai_wrapper.py
│   ├── train.py
│   ├── evaluate.py
│   ├── run_ablation.py
│   └── interpret.py
├── configs/
│   ├── cnn_pm200.yaml
│   ├── rnafm_pm200.yaml
│   ├── rnabert_pm200.yaml
│   └── ablation.yaml
├── results/
│   ├── tables/
│   ├── figures/
│   └── checkpoints/
└── reports/
    ├── experiment_1.md
    ├── experiment_2.md
    └── experiment_3.md
```

如果时间有限，可以先只实现 `data/processed/`、`src/train.py`、`src/evaluate.py` 和 `results/`，其余文件按需要逐步补齐。

---

## 2. 环境准备

### 2.1 Python 环境

建议使用 Python 3.9 或 3.10，并创建独立虚拟环境。

```bash
conda create -n splice-rna python=3.10
conda activate splice-rna
```

### 2.2 基础依赖

```bash
pip install numpy pandas scikit-learn matplotlib seaborn tqdm biopython pyfaidx
pip install torch torchvision torchaudio
pip install transformers datasets accelerate
```

### 2.3 模型依赖

不同模型可能需要额外安装：

```bash
# RNA-FM
pip install fair-esm

# RNABERT / multimolecule 模型
pip install multimolecule

# SpliceAI
pip install spliceai
```

如果 SpliceAI 安装失败，可以先将它作为论文中的强基线说明，或者使用官方预测结果 / 第三方实现作为对照。课程项目优先保证 RNA-FM、RNABERT、CNN baseline 三条线能完整跑通。

---

## 3. 数据准备总览

### 3.1 原始数据

建议使用人类基因组和注释文件：

| 文件 | 用途 | 示例来源 |
| --- | --- | --- |
| `genome.fa` | 提取剪接位点附近序列窗口 | GRCh38 / hg38 |
| `gencode.gtf` | 提取 exon-intron 边界 | GENCODE |
| `clinvar.vcf` | 构造真实剪接相关变异样本 | ClinVar |

实验一和实验二只需要 `genome.fa` 与 `gencode.gtf`。实验三可以先用人工突变，不依赖 ClinVar。

### 3.2 样本类别定义

三分类标签统一定义如下：

| 标签 | 类别 | 生物学含义 |
| --- | --- | --- |
| `0` | donor | 5' splice site，通常位于 exon 到 intron 边界 |
| `1` | acceptor | 3' splice site，通常位于 intron 到 exon 边界 |
| `2` | non_splice | 含 GT / AG motif 但不是注释剪接位点的负样本 |

注意：负样本不要从全基因组随机采样普通位置，否则任务会太简单。建议从同一基因或相邻区域中采样 GT / AG motif，但排除已注释的 donor / acceptor 位点，让模型必须学习上下文而不只是识别 GT / AG。

### 3.3 序列方向

构造样本时必须处理正负链：

1. 正链基因按基因组序列直接提取。
2. 负链基因需要取反向互补序列。
3. 所有样本最终应统一为转录方向上的序列窗口。
4. 如果输入 RNA-FM / RNABERT，需要将 DNA 序列中的 `T` 替换为 `U`。

### 3.4 数据切分

建议按染色体或基因切分，而不是随机按样本切分，避免同一基因的高度相似窗口同时出现在训练集和测试集。

推荐切分方式：

```text
train: chr1-chr16
valid: chr17-chr18
test : chr19-chr22, chrX
```

如果样本量较小，也可以按基因 ID 分组后做 8:1:1 切分。

---

## 4. 实验一：剪接位点三分类实验

### 4.1 实验目的

验证 RNA 基础模型的 embedding 是否包含剪接位点识别所需的序列上下文信息。

需要回答：

> RNA-FM / RNABERT + MLP 是否优于普通 CNN baseline？它们与任务专用模型 SpliceAI 之间还有多大差距？

### 4.2 输入与输出

输入是一段以候选位点为中心的序列窗口，建议默认使用：

```text
中心位点 ±200 nt，总长度 401 nt
```

输出为三分类概率：

```text
P(donor), P(acceptor), P(non_splice)
```

最终预测类别取概率最大的类别。

### 4.3 数据构造步骤

1. 从 `gencode.gtf` 中读取所有 protein-coding 或全部转录本的 exon 坐标。
2. 对每个转录本按基因组坐标排序 exon。
3. 对相邻 exon 之间的 intron 边界提取：
   - donor site：exon 结束到 intron 开始的边界。
   - acceptor site：intron 结束到 exon 开始的边界。
4. 对每个 donor / acceptor 位点，从 `genome.fa` 提取中心 ±200 nt 的窗口。
5. 从同一基因附近采样 GT / AG motif 作为 non-splice 负样本，并排除注释边界。
6. 过滤包含过多 `N` 的序列，建议丢弃 `N` 比例大于 5% 的样本。
7. 保证三类样本数量尽量均衡，例如 donor:acceptor:non_splice = 1:1:1。
8. 按染色体或基因切分训练集、验证集、测试集。

建议输出 CSV 格式：

```text
sample_id,chrom,start,end,strand,center,label,label_name,sequence,gene_id,transcript_id
```

示例：

```text
ENST000001_donor_1,chr1,10000,10400,+,10200,0,donor,AUGG...,GENE1,ENST000001
```

### 4.4 模型设置

#### CNN baseline

CNN baseline 用来检验深度基础模型是否真的带来增益。推荐结构：

```text
one-hot sequence
→ Conv1D kernel sizes 3/5/9
→ ReLU
→ MaxPool
→ Dropout
→ Linear
→ 3-class softmax
```

输入编码：

```text
A: [1,0,0,0]
C: [0,1,0,0]
G: [0,0,1,0]
T/U: [0,0,0,1]
N: [0,0,0,0]
```

#### RNA-FM + MLP

推荐先冻结 RNA-FM，只训练后接 MLP：

```text
RNA sequence
→ RNA-FM frozen encoder
→ CLS embedding 或 mean pooling
→ MLP classifier
→ 3-class softmax
```

MLP 可使用：

```text
Linear(hidden_dim, 256)
ReLU
Dropout(0.2)
Linear(256, 3)
```

#### RNABERT + MLP

RNABERT 与 RNA-FM 保持同样训练策略：

```text
RNA sequence
→ RNABERT frozen encoder
→ CLS / mean pooling
→ MLP classifier
→ 3-class softmax
```

为了公平比较，RNA-FM 与 RNABERT 应使用同一份数据切分、同样窗口长度、同样训练轮数和同样评价指标。

#### SpliceAI baseline

SpliceAI 是任务专用强基线。它通常输出每个位置作为 donor / acceptor 的概率，使用时需要把中心位置附近的 donor / acceptor 分数转成三分类结果。

建议处理方式：

```text
donor_score = SpliceAI 对中心附近 donor 的最大预测分数
acceptor_score = SpliceAI 对中心附近 acceptor 的最大预测分数
non_splice_score = 1 - max(donor_score, acceptor_score)
```

再将三个分数组合成预测类别。

### 4.5 训练配置

推荐起始配置：

```text
batch_size: 16 或 32
epochs: 10-20
optimizer: AdamW
learning_rate:
  CNN: 1e-3
  MLP head: 1e-4 或 5e-4
weight_decay: 1e-4
early_stopping: valid macro-F1 连续 3 轮不提升则停止
loss: CrossEntropyLoss
```

如果类别不平衡，使用 class weight 或 weighted sampler。

### 4.6 评价指标

必须报告：

```text
Accuracy
Macro Precision
Macro Recall
Macro F1-score
AUROC
AUPRC
Confusion Matrix
```

建议额外报告每一类的指标：

```text
Donor F1
Acceptor F1
Non-splice F1
```

尤其要观察 donor 和 acceptor 是否相互混淆。如果 donor / acceptor 与 non-splice 区分较好，但 donor 与 acceptor 混淆严重，说明模型可能学到了 GT / AG motif，却没有充分学到剪接方向和上下文。

### 4.7 推荐结果表

在 `results/tables/experiment_1_metrics.csv` 中保存：

| Model | Accuracy | Macro-F1 | AUROC | AUPRC | Donor-F1 | Acceptor-F1 | Non-splice-F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CNN |  |  |  |  |  |  |  |
| RNABERT + MLP |  |  |  |  |  |  |  |
| RNA-FM + MLP |  |  |  |  |  |  |  |
| SpliceAI |  |  |  |  |  |  |  |

推荐图表：

```text
results/figures/exp1_confusion_matrix_cnn.png
results/figures/exp1_confusion_matrix_rnafm.png
results/figures/exp1_confusion_matrix_rnabert.png
results/figures/exp1_model_comparison_barplot.png
```

### 4.8 论文中如何解释

如果 RNA-FM / RNABERT 强于 CNN，可以说明：

> RNA 序列基础模型的预训练表示包含一定剪接相关语法，能够迁移到剪接位点识别任务。

如果它们弱于 SpliceAI，可以说明：

> 通用 RNA foundation model 虽然具有迁移潜力，但任务专用监督模型在剪接识别上仍更稳定，尤其可能受益于更长上下文和剪接任务监督信号。

---

## 5. 实验二：上下文长度消融实验

### 5.1 实验目的

分析剪接位点识别到底主要依赖局部 motif，还是需要更长程上下文。

该实验与 AlphaGenome / sQTL prediction 中强调的 long-range context、splice competition、junction topology 思路对应。虽然本项目规模较小，但可以用窗口长度消融做一个可解释的验证。

### 5.2 窗口长度设置

建议使用：

```text
±50 nt   -> 总长度 101 nt
±100 nt  -> 总长度 201 nt
±200 nt  -> 总长度 401 nt
±400 nt  -> 总长度 801 nt
±1000 nt -> 总长度 2001 nt，可选
```

课程项目中，`±50 / ±100 / ±200 / ±400` 已经足够。`±1000` 对显存和模型最大输入长度要求更高，可以作为可选扩展。

### 5.3 数据构造

保持实验一的中心位点、标签和数据切分不变，只改变提取窗口长度。

关键要求：

1. 不同窗口长度必须对应同一批中心位点。
2. train / valid / test 划分必须完全一致。
3. 不要让窗口长度变化引入新的样本过滤偏差。
4. 如果某些长窗口越界，可以统一丢弃这些中心位点，或用 `N` padding，但所有窗口长度要使用同一处理规则。

推荐输出：

```text
data/processed/splice_sites_pm50.csv
data/processed/splice_sites_pm100.csv
data/processed/splice_sites_pm200.csv
data/processed/splice_sites_pm400.csv
```

### 5.4 模型训练

对每个窗口长度，重复训练相同模型：

```text
CNN baseline
RNABERT + MLP
RNA-FM + MLP
SpliceAI
```

为了减少工作量，可以采用两阶段策略：

1. 先跑 CNN 与 RNA-FM，确认消融曲线。
2. 再补 RNABERT 与 SpliceAI。

如果显存有限，长窗口下优先冻结基础模型，只训练 MLP head。

### 5.5 推荐实验矩阵

| Model | ±50 | ±100 | ±200 | ±400 | ±1000 |
| --- | --- | --- | --- | --- | --- |
| CNN | 必做 | 必做 | 必做 | 必做 | 可选 |
| RNABERT + MLP | 必做 | 必做 | 必做 | 必做 | 可选 |
| RNA-FM + MLP | 必做 | 必做 | 必做 | 必做 | 可选 |
| SpliceAI | 可选 | 可选 | 必做 | 必做 | 可选 |

### 5.6 评价指标

重点关注：

```text
Macro-F1
AUROC
AUPRC
Donor-F1
Acceptor-F1
```

推荐画折线图：

```text
x-axis: window size
y-axis: Macro-F1 / AUPRC
line: different models
```

保存路径：

```text
results/tables/experiment_2_ablation.csv
results/figures/exp2_window_size_macro_f1.png
results/figures/exp2_window_size_auprc.png
```

### 5.7 可能结果与解释

可能观察到：

```text
CNN 在 ±50 或 ±100 下已经能识别部分 motif，但长窗口提升有限。
RNA-FM / RNABERT 在 ±200 或 ±400 下性能更好，说明其表示能利用更丰富上下文。
SpliceAI 在较长窗口下表现最好，说明剪接任务专用模型对长程依赖建模更充分。
```

如果长窗口没有提升，也可以解释：

> 本项目的三分类任务仍以局部剪接 motif 为主，长程上下文优势可能需要在更复杂的 junction-level 或 variant-level 任务中才能充分体现。

这不会导致实验失败，反而可以自然引出实验三。

---

## 6. 实验三：异常剪接变异效应预测实验

### 6.1 实验目的

从“识别已有剪接位点”扩展到“判断突变是否破坏剪接信号”。

需要回答：

> RNA 基础模型能不能感知单碱基突变对 donor / acceptor 剪接信号的扰动？

### 6.2 数据方案 A：人工构造突变

如果课程时间有限，优先使用人工突变。它实现简单、标签清晰、可解释性强。

人工构造三类变异：

```text
1. donor 损伤突变：
   将真实 donor site 附近的 GT 改为 AT / GC / TT

2. acceptor 损伤突变：
   将真实 acceptor site 附近的 AG 改为 AA / AC / GG

3. cryptic splice gain 突变：
   在非剪接区域人为制造新的 GT / AG motif
```

标签定义：

| 标签 | 含义 |
| --- | --- |
| `1` | splice-altering variant，预期影响剪接 |
| `0` | neutral variant，预期不明显影响剪接 |

neutral variant 可以从远离剪接位点的位置随机替换同类碱基构造，例如不改变 GT / AG motif 的普通单碱基突变。

推荐输出：

```text
variant_id,chrom,pos,strand,ref,alt,variant_type,label,wt_sequence,mut_sequence,target_class
```

其中 `target_class` 可以记录该变异主要影响 donor 还是 acceptor。

### 6.3 数据方案 B：真实 ClinVar 变异

如果时间充足，可以使用 ClinVar 中带 splicing consequence 或相关 clinical significance 的变异。

基本步骤：

1. 读取 `clinvar.vcf`。
2. 筛选注释中包含 splice donor、splice acceptor、splice region、splice-altering 等关键词的变异作为阳性。
3. 从同基因或同区域中选取未标注剪接影响的变异作为阴性。
4. 对每个变异提取 wild-type 序列窗口。
5. 将中心位置替换为 alt allele，得到 mutant 序列窗口。
6. 使用实验一训练好的模型分别预测 WT 与 Mut。

真实数据标签更复杂，可能存在噪声。因此课程项目建议将真实 ClinVar 作为扩展分析，而不是唯一主实验。

### 6.4 模型打分方式

对 RNA-FM / RNABERT / CNN：

```text
WT sequence  -> classifier -> score_wt
Mut sequence -> classifier -> score_mut
Delta score = score_mut - score_wt
```

对 donor 损伤突变，可以定义：

```text
delta_donor = P_mut(donor) - P_wt(donor)
```

如果 `delta_donor` 显著为负，说明模型认为突变降低 donor 信号。

对 acceptor 损伤突变：

```text
delta_acceptor = P_mut(acceptor) - P_wt(acceptor)
```

对 cryptic splice gain 突变：

```text
delta_splice = max(P_mut(donor), P_mut(acceptor)) - max(P_wt(donor), P_wt(acceptor))
```

如果 `delta_splice` 显著为正，说明模型认为突变制造了新的剪接信号。

对 SpliceAI：

```text
直接使用 SpliceAI 输出的 variant effect score
```

或者使用 donor gain、donor loss、acceptor gain、acceptor loss 四类分数中的最大值作为总体变异效应分数。

### 6.5 评价指标

推荐报告：

```text
AUROC
AUPRC
Top-k recall
Delta score distribution
```

Top-k recall 的含义是：按照模型预测的变异影响分数从高到低排序，前 k 个变异中覆盖了多少真实 splice-altering variant。

推荐保存：

```text
results/tables/experiment_3_variant_metrics.csv
results/figures/exp3_delta_score_boxplot.png
results/figures/exp3_variant_auroc.png
results/figures/exp3_variant_auprc.png
```

### 6.6 推荐分析

重点比较：

1. 人工 donor GT 破坏后，模型 donor score 是否下降。
2. 人工 acceptor AG 破坏后，模型 acceptor score 是否下降。
3. 非剪接区域制造 GT / AG 后，模型是否出现 donor / acceptor gain。
4. RNA-FM / RNABERT 的 delta score 是否比 CNN 更稳定。
5. SpliceAI 是否在变异效应预测上明显更强。

### 6.7 论文中如何解释

如果 RNA-FM / RNABERT 能正确响应人工突变：

> 说明 RNA 基础模型不仅能识别静态剪接位点，还能对关键 motif 的单碱基扰动产生合理响应。

如果它们弱于 SpliceAI：

> 说明基础模型的通用表示虽然包含剪接信号，但要用于疾病相关变异解释，仍需要更长上下文建模和剪接监督微调。

---

## 7. 可解释性分析

可解释性分析可以作为实验一或实验三的补充，不一定单独作为主实验。

### 7.1 Attention 可视化

适用于 RNA-FM / RNABERT。

步骤：

1. 选择若干 donor、acceptor、non-splice 测试样本。
2. 提取最后几层 attention map。
3. 对 attention head 求平均，得到每个位置的重要性。
4. 将中心剪接位点附近 ±50 nt 的 attention 画成热图。

重点观察：

```text
donor site: 中心 GT 附近是否高亮
acceptor site: 中心 AG 附近是否高亮
acceptor 上游 polypyrimidine tract 是否有响应
```

输出图：

```text
results/figures/interpret_attention_donor.png
results/figures/interpret_attention_acceptor.png
```

### 7.2 In silico mutagenesis

更推荐使用该方法，因为它更直观。

步骤：

1. 选择一个模型预测正确且置信度高的 donor 或 acceptor 样本。
2. 对序列每个位置依次突变为 A / C / G / T。
3. 每次突变后重新预测。
4. 记录目标类别概率变化。
5. 画出位置 × 碱基的 importance heatmap。

对于 donor 样本：

```text
importance = P_original(donor) - P_mutated(donor)
```

对于 acceptor 样本：

```text
importance = P_original(acceptor) - P_mutated(acceptor)
```

如果 GT / AG 附近突变造成最大分数下降，就能说明模型学到了合理的剪接信号。

输出图：

```text
results/figures/ism_donor_heatmap.png
results/figures/ism_acceptor_heatmap.png
```

---

## 8. 推荐运行顺序

建议按下面顺序推进，避免一开始就陷入复杂模型调试。

### 第一步：完成数据集

目标产物：

```text
data/processed/splice_sites_pm200.csv
data/splits/train.csv
data/splits/valid.csv
data/splits/test.csv
```

检查点：

```text
三类样本数量是否接近均衡
序列长度是否全部为 401
是否正确处理正负链
是否没有 train/test 基因泄漏
```

### 第二步：跑 CNN baseline

目标产物：

```text
results/tables/cnn_pm200_metrics.csv
results/figures/cnn_pm200_confusion_matrix.png
```

CNN baseline 是整个项目的 sanity check。如果 CNN 都无法学习，优先检查标签、序列方向和数据构造。

### 第三步：跑 RNA-FM + MLP

目标产物：

```text
results/tables/rnafm_pm200_metrics.csv
results/checkpoints/rnafm_mlp_best.pt
```

先冻结 RNA-FM，只训练 MLP。确认流程跑通后，再考虑是否微调最后几层。

### 第四步：跑 RNABERT + MLP

目标产物：

```text
results/tables/rnabert_pm200_metrics.csv
results/checkpoints/rnabert_mlp_best.pt
```

保持与 RNA-FM 完全相同的数据和评价流程。

### 第五步：补 SpliceAI baseline

目标产物：

```text
results/tables/spliceai_pm200_metrics.csv
```

如果安装或运行 SpliceAI 困难，可以在论文中将 SpliceAI 作为文献强基线，并说明本项目主要完成 RNA-FM、RNABERT 与 CNN 的实证比较。

### 第六步：做上下文长度消融

目标产物：

```text
results/tables/experiment_2_ablation.csv
results/figures/exp2_window_size_macro_f1.png
```

优先跑 `±50 / ±100 / ±200 / ±400`。

### 第七步：做变异效应预测

目标产物：

```text
data/processed/variant_effect.csv
results/tables/experiment_3_variant_metrics.csv
results/figures/exp3_delta_score_boxplot.png
```

优先使用人工构造突变，确保结论清晰。

---

## 9. 实验记录模板

每次训练建议记录以下信息：

```text
experiment_id:
date:
model:
window_size:
dataset_version:
train_samples:
valid_samples:
test_samples:
learning_rate:
batch_size:
epochs:
best_valid_macro_f1:
test_accuracy:
test_macro_f1:
test_auroc:
test_auprc:
notes:
```

可以把每次实验写入 `reports/experiment_log.md`，避免后期整理结果时找不到配置。

---

## 10. 最终论文 / 汇报建议

### 10.1 论文结构

推荐结构：

```text
1. 引言
   介绍 RNA 剪接、剪接变异、RNA foundation model、SpliceAI。

2. 研究问题
   RNA foundation model 是否学到了可迁移的剪接信号？

3. 方法
   数据构造、模型、训练设置、评价指标。

4. 实验一：剪接位点三分类
   展示模型总体性能与混淆矩阵。

5. 实验二：上下文长度消融
   展示窗口长度与性能关系。

6. 实验三：异常剪接变异效应预测
   展示 WT/Mut delta score 与变异识别指标。

7. 可解释性分析
   展示 attention 或 in silico mutagenesis。

8. 讨论
   分析 RNA-FM / RNABERT 与 SpliceAI 差距。

9. 结论
   总结基础模型的迁移潜力与局限。
```

### 10.2 建议核心结论

可以围绕下面这类结论组织：

> RNA-FM 和 RNABERT 等 RNA 基础模型能够通过预训练表示捕捉一定的剪接信号，在剪接位点识别任务上优于普通 CNN baseline；但在异常剪接变异效应预测上，任务专用模型 SpliceAI 仍然更稳定。说明 RNA foundation model 具备转录组学任务迁移潜力，但要真正用于疾病相关剪接变异解释，还需要引入更长程上下文、剪接位点竞争机制和监督微调。

### 10.3 建议展示图表

最终汇报中至少准备：

```text
图 1：项目整体实验流程图
图 2：剪接位点三分类任务示意图
图 3：不同模型在实验一上的性能柱状图
图 4：混淆矩阵，展示 donor / acceptor / non-splice 混淆情况
图 5：上下文长度消融折线图
图 6：突变前后 delta score 箱线图
图 7：in silico mutagenesis heatmap
```

---

---



