# TODO

本项目后续实验**不使用百万级完整数据集**。所有主实验默认基于当前小样本 split：

- `data/shared/splits/train.csv`：855 条样本
- `data/shared/splits/valid.csv`：120 条样本
- `data/shared/splits/test.csv`：285 条样本

`data/raw/` 中的大型 `genome.fa`、`gencode.gtf` 只作为可选资源准备来源，不作为主实验运行前提。TODO 中所有任务都应围绕“小样本可复现、结果可解释、报告清楚”推进。

## 当前已完成

- 数据
  - `data/shared/splits/train.csv`
  - `data/shared/splits/valid.csv`
  - `data/shared/splits/test.csv`
  - `data/shared/processed/splice_sites_pm50.csv`
  - `data/shared/processed/splice_sites_pm100.csv`
  - `data/shared/processed/splice_sites_pm200.csv`
  - `data/shared/processed/splice_sites_pm400.csv`
  - `data/experiment_3/artificial_variant_effect.csv`

- 数据代码
  - `src/data/build_synthetic_splice_dataset.py`
  - `src/data/build_splice_site_dataset.py`
  - `src/data/split_dataset.py`
  - `src/data/build_variant_dataset.py`

- 实验一
  - `src/experiments/exp1/train.py`
  - `src/experiments/exp1/evaluate.py`
  - `src/experiments/exp1/run_classification.py`
  - `reports/experiment_1.md`

- 实验二
  - `src/experiments/exp2/run_multiscale.py`
  - `src/experiments/exp2/run_full_context.py`
  - `reports/experiment_2.md`

- 实验三
  - `src/experiments/exp3/run_variant_effect.py`
  - `src/experiments/exp3/run_interpretability.py`
  - `reports/experiment_3.md`

## 总原则

- [x] 所有实验入口默认使用当前小样本 split。
  - 说明：统一检查三个实验的默认参数和路径解析，确保直接运行实验脚本时读取 `data/shared/splits/` 下的 855/120/285 小样本。
  - 修改：`src/experiments/exp1/common.py`
  - 修改：`src/experiments/exp2/run_multiscale.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 修改：`src/experiments/exp3/run_interpretability.py`
  - 要点：不要再默认尝试生成或读取百万级全量数据。
  - 验收：各实验 `--help` 和默认运行路径中不再出现必须依赖 `data/raw/` 或旧全量 split 的逻辑。

- [x] README 中明确数据规模选择。
  - 说明：把“本项目主实验使用小样本 split”写入项目入口文档和三份实验报告，避免读者误以为当前结果来自百万级完整数据集。
  - 修改：`README.md`
  - 修改：`reports/experiment_1.md`
  - 修改：`reports/experiment_2.md`
  - 修改：`reports/experiment_3.md`
  - 要点：说明主实验使用 `train/valid/test = 855/120/285` 的小样本，而不是完整数据集。
  - 验收：README 和报告均能清楚区分“小样本主实验”“可选 raw 资源”“未运行的全量实验”。

## 实验一：剪接位点三分类

- [x] 补齐 SpliceAI 三分类 baseline，但只在当前小样本 test 上运行。
  - 说明：封装 SpliceAI 推理接口，把 donor/acceptor 信号转换成实验一的三分类输出，用作与 CNN、RNA-FM proxy、RNABERT proxy 对比的真实工具 baseline。
  - 新增：`src/models/spliceai_wrapper.py`
  - 修改：`src/experiments/exp1/common.py`
  - 修改：`src/experiments/exp1/run_classification.py`
  - 输出：`results/experiment_1/tables/experiment_1_metrics.csv`
  - 输出：`results/experiment_1/figures/experiment_1_macro_f1.png`
  - 要点：将中心附近 donor / acceptor score 转成 `donor, acceptor, non_splice` 三分类概率。
  - 完成备注：当前环境缺少真实 `spliceai` 包，实验一主表中的该行以 `backend=spliceai_signal_proxy` 明确标注；不会把 proxy fallback 写成真实 SpliceAI 完整推理。
  - 验收：metrics 表中新增 SpliceAI 行，且该行只基于 `data/shared/splits/test.csv` 计算。

- [x] 修正实验一 hard-negative FPR。
  - 说明：当前 hard-negative 指标需要严格限定在 `negative_type` 标注的困难负例上，避免把普通 non-splice 样本混入 FPR 统计。
  - 修改：`src/experiments/exp1/evaluate.py`
  - 修改：`src/experiments/exp1/common.py`
  - 依赖：`data/shared/splits/test.csv` 中应保留 `negative_type` 字段。
  - 输出：`results/experiment_1/tables/experiment_1_metrics.csv`
  - 验收：输出表中 hard-negative FPR 的分母、样本数和筛选条件可追溯，报告中能解释该指标的含义。

- [x] 明确 proxy 与真实 frozen encoder 结果的区别。
  - 说明：把当前 k-mer/token 特征近似实验和真实 pretrained backbone 小样本推理分开呈现，避免把 proxy 结果描述成真实 foundation model 结果。
  - 修改：`reports/experiment_1.md`
  - 修改：`src/experiments/exp1/run_classification.py`
  - 输出：`results/experiment_1/tables/real_foundation/`
  - 要点：报告中明确哪些是 k-mer/token proxy，哪些是真实 pretrained backbone 小样本运行。
  - 验收：实验一报告中每个模型名称都有清晰标签，例如 `proxy`、`frozen encoder` 或 `optional real model`。

- [x] 明确 RNA-FM / RNABERT 本地权重准备方式。
  - 说明：提供一个资源准备入口和文档说明，用于把可选 pretrained 权重放到本地固定目录，主实验没有权重时仍可运行 proxy 链路。
  - 新增：`src/resources/prepare_foundation_models.py`
  - 修改：`README.md`
  - 修改：`reports/experiment_1.md`
  - 输出说明：`models/hf/rnafm/`、`models/hf/rnabert/`
  - 要点：只需要支持当前小样本推理，不要求全量训练。
  - 验收：缺少本地权重时脚本给出明确提示；存在权重时可以对当前小样本完成 frozen encoder 推理。

## 实验二：多尺度上下文与 hard negative

- [x] 保持实验二窗口范围为当前 `pm50/pm100/pm200/pm400`。
  - 说明：将实验二主线限定为已经准备好的四个窗口文件，保证多尺度比较可以在当前仓库数据上复现。
  - 修改：`src/experiments/exp2/run_multiscale.py`
  - 修改：`reports/experiment_2.md`
  - 要点：不再把 `pm1000`、5kb、10kb 作为主实验 TODO；如需展示，只作为可选 case study。
  - 验收：默认运行只读取 `splice_sites_pm50/100/200/400.csv`，报告不把未生成的大窗口写成已完成实验。

- [x] 接入 Pangolin 小样本 case study。
  - 说明：补一个 Pangolin 包装器，在 test 或 hard-negative 小集合上跑少量样本，用来观察长上下文工具在困难负例上的表现。
  - 新增：`src/models/pangolin_wrapper.py`
  - 修改：`src/experiments/exp2/run_multiscale.py`
  - 修改：`src/resources/run_real_model_smoke.py`
  - 输出：`results/experiment_2/tables/experiment_2B_hard_negative.csv`
  - 输出：`results/experiment_2/figures/exp2B_hard_negative_fpr.png`
  - 要点：只对当前 test / hard-negative 小样本运行，不做全量 benchmark。
  - 验收：Pangolin 结果单独标注为 case study，失败或缺依赖时不影响实验二 proxy 主链路。

- [x] 增加 rare motif 小样本补充测试。
  - 说明：构造少量非 canonical motif 样本，例如 GC-AG，观察模型是否只记住 GT-AG/AG 规则。
  - 修改：`src/data/build_synthetic_splice_dataset.py`
  - 修改：`src/experiments/exp2/run_multiscale.py`
  - 输出：`data/shared/processed/rare_motif_splice_sites.csv`
  - 输出：`results/experiment_2/tables/experiment_2B_rare_motif.csv`
  - 要点：样本量控制在几十到几百条，用于讨论 GC-AG 等非 canonical motif。
  - 验收：rare motif 表中包含 motif 类型、样本数、模型预测分数和按 motif 分组的指标。

- [x] 将 tissue-specific usage 明确为小样本 case study。
  - 说明：把 GTEx 相关内容收缩为少量组织/事件展示，只用于说明未来可扩展方向，不作为组织特异性建模结论。
  - 修改：`src/resources/fetch_gtex_sqtl_cases.py`
  - 修改：`src/experiments/exp2/run_full_context.py`
  - 修改：`reports/experiment_2.md`
  - 输出：`results/experiment_2/tables/experiment_2C_gtex_tissue_usage.csv`
  - 输出：`results/experiment_2/figures/exp2C_tissue_splice_usage_heatmap.png`
  - 要点：不做 GTEx 全量建模，只做少量组织/事件展示。
  - 验收：输出表和热图标注 case 数量，报告中不使用“大规模 GTEx 验证”之类表述。

- [x] Borzoi / AlphaGenome 只保留为可选长程 case study。
  - 说明：长程调控模型暂不纳入主结果，只保留接口草案、文档示例或极少量可运行样例。
  - 可新增：`src/models/borzoi_wrapper.py`
  - 可新增：`src/models/alphagenome_case_study.md`
  - 修改：`reports/experiment_2.md`
  - 输出：`results/experiment_2/tables/long_range_regulatory_case_study.csv`
  - 要点：不要求真实模型全量推理；可用文档式 case study 或极少数示例。
  - 验收：报告明确这是 optional case study，缺少外部模型或 API 时不会阻塞主实验。

- [x] 报告中明确 proxy / real case study 边界。
  - 说明：整理实验二报告措辞，把 synthetic/proxy 分析、真实模型 smoke、GTEx case study 放在不同小节。
  - 修改：`reports/experiment_2.md`
  - 要点：tissue usage、junction topology 当前主要是 synthetic/proxy，不写成大规模真实生物结论。
  - 验收：读者可以从小节标题和图表注释直接判断每个结果的数据来源和证据强度。

## 实验三：异常剪接变异效应预测

- [x] 拆分 `cryptic_gain` 为 `donor_gain` 与 `acceptor_gain`。
  - 说明：人工变异集需要区分新增 donor 信号和新增 acceptor 信号，便于分析不同异常剪接机制。
  - 修改：`src/data/build_variant_dataset.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 修改：`src/experiments/exp3/run_interpretability.py`
  - 输出：`data/experiment_3/artificial_variant_effect.csv`
  - 输出：`results/experiment_3/tables/variant_effect_stratified_by_type.csv`
  - 要点：继续使用小样本人工变异，不扩展到大规模变异集。
  - 验收：变异类型统计中不再只有笼统 `cryptic_gain`，分层指标能分别显示 donor_gain 和 acceptor_gain。

- [x] ClinVar 只做小样本 smoke / case study。
  - 说明：从 ClinVar 中选少量与剪接相关的变异和 matched control，验证流程能跑通，不追求完整 ClinVar benchmark。
  - 新增：`src/data/build_clinvar_variant_dataset.py`
  - 修改：`src/resources/make_clinvar_smoke.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 输出：`data/experiment_3/clinvar_splicing_variants_smoke.csv`
  - 输出：`results/experiment_3/tables/experiment_3B_clinvar_smoke_metrics.csv`
  - 输出：`results/experiment_3/figures/exp3_clinvar_smoke_scores.png`
  - 要点：不做完整 ClinVar benchmark，只选少量 splice-related 与 matched control 示例。
  - 验收：ClinVar 结果文件中包含样本来源、标签构造方式和模型分数，报告标注为 smoke/case study。

- [x] 将真实 SpliceAI / Pangolin / MMSplice / MaxEntScan 纳入小样本变异对照。
  - 说明：为实验三增加真实剪接变异工具的可选 wrapper，在人工变异和 ClinVar smoke 上形成小规模对照。
  - 新增：`src/models/spliceai_wrapper.py`
  - 新增：`src/models/pangolin_wrapper.py`
  - 新增：`src/models/mmsplice_wrapper.py`
  - 新增：`src/models/maxentscan_wrapper.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 输出：`results/experiment_3/tables/experiment_3A_artificial_variant_metrics.csv`
  - 输出：`results/experiment_3/tables/experiment_3B_clinvar_smoke_metrics.csv`
  - 要点：运行对象为当前人工变异小样本和 ClinVar smoke 小样本。
  - 验收：每个工具的依赖、输入格式、失败跳过逻辑清楚；结果表能区分 proxy 模型和真实工具。

- [x] GTEx sQTL 只做小样本 tissue-specific case study。
  - 说明：选取少量 sQTL 事件，记录组织、基因、junction 和模型 delta score，用于展示变异效应与组织剪接事件的连接方式。
  - 修改：`src/resources/fetch_gtex_sqtl_cases.py`
  - 新增：`src/data/build_sqtl_variant_dataset.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 输出：`data/experiment_3/gtex_sqtl_variants_smoke.csv`
  - 输出：`results/experiment_3/tables/experiment_3C_sqtl_case_study.csv`
  - 要点：输出字段保留 `variant_id,tissue,target_gene,target_junction,observed_effect_direction,model_delta_score`，但不做全量 GTEx 评估。
  - 验收：报告中只讨论具体 case 的方向一致性或不一致性，不给出总体 GTEx 性能结论。

- [x] 实现 RNA-FM / RNABERT pseudo-likelihood zero-shot scoring。
  - 说明：在不训练新分类头的前提下，用 masked/pseudo-likelihood 分数比较 reference 与 alternate 序列，形成 zero-shot 变异影响分数。
  - 修改：`src/models/foundation_backbones.py`
  - 修改：`src/models/simple_splice_models.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 输出：`results/experiment_3/tables/experiment_3A_artificial_variant_scores.csv`
  - 要点：只在当前人工变异小样本上计算。
  - 验收：分数表中包含 ref_score、alt_score、delta_score 和模型名称，缺少权重时能跳过并保留说明。

- [x] 增加 calibration curve 与多阈值 Top-k / Enrichment。
  - 说明：补充比单点 AUC/F1 更细的排序与校准分析，用来观察小样本下高分变异是否更集中命中正例。
  - 修改：`src/utils.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 输出：`results/experiment_3/tables/experiment_3A_topk_enrichment_curve.csv`
  - 输出：`results/experiment_3/figures/exp3_calibration_curve.png`
  - 要点：基于当前小样本，报告中说明样本量限制。
  - 验收：Top-k 表至少包含 k、positive_count、precision/enrichment；校准图标注样本数和 bin 数。

- [x] 增加 ClinVar smoke case 的 variant delta profile。
  - 说明：对 ClinVar smoke 中的代表变异绘制突变前后窗口内 donor/acceptor 信号变化，帮助解释模型判断来自哪个局部区域。
  - 修改：`src/experiments/exp3/run_interpretability.py`
  - 依赖：`data/experiment_3/clinvar_splicing_variants_smoke.csv`
  - 输出：`results/experiment_3/tables/variant_delta_profile_clinvar_smoke_case.csv`
  - 输出：`results/experiment_3/figures/variant_delta_profile_clinvar_smoke_case.png`
  - 验收：每个 case 至少包含 variant_id、坐标、ref/alt 序列窗口、delta 最大位置和对应图。

## 可解释性

- [x] 增加 RNA-FM / RNABERT attention 可视化。
  - 说明：从当前 test 小样本中挑选 donor、acceptor、hard-negative 代表序列，导出 attention 热图用于解释 foundation encoder 关注区域。
  - 修改：`src/models/foundation_backbones.py`
  - 修改：`src/experiments/exp3/run_interpretability.py`
  - 输出：`results/experiment_3/figures/attention_donor_heatmap.png`
  - 输出：`results/experiment_3/figures/attention_acceptor_heatmap.png`
  - 输出：`results/experiment_3/figures/attention_hard_negative_heatmap.png`
  - 要点：只选择当前 test 小样本中的少量代表样本。
  - 验收：图中标出中心位点和 motif 区域；没有真实权重时报告中明确该项未运行或仅有 proxy。

- [x] 将 ISM 扩展到 CNN、RNA-FM、RNABERT。
  - 说明：把 in-silico mutagenesis 从当前 signal proxy 扩展到其他模型，逐碱基替换并记录预测分数变化。
  - 修改：`src/experiments/exp3/run_interpretability.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 输出：`results/experiment_3/tables/*_ism_matrix.csv`
  - 输出：`results/experiment_3/figures/ism_*_heatmap.png`
  - 要点：当前 ISM 主要用 `SpliceAI signal proxy`；扩展时只选少量样本，避免计算量失控。
  - 验收：每个模型的 ISM 输出包含样本 ID、位置、替换碱基和 delta score，运行样本数可配置。

## 数据与资源

- [x] 确认当前小样本 split 字段完整。
  - 说明：检查三个 split 是否保留训练、评估、hard-negative 分析和变异构建需要的最小字段集合。
  - 修改：`src/data/build_splice_site_dataset.py`
  - 修改：`src/data/split_dataset.py`
  - 检查：`data/shared/splits/train.csv`
  - 检查：`data/shared/splits/valid.csv`
  - 检查：`data/shared/splits/test.csv`
  - 必要字段：`sample_id,chrom,start,end,strand,center,label,label_name,negative_type,sequence,gene_id,transcript_id`
  - 验收：新增检查脚本或构建流程在字段缺失时直接报错，并提示缺失列名。

- [x] 增加小样本数据质控报告。
  - 说明：为当前小样本生成独立 QC 文档，记录类别分布、长度分布、motif 情况和 gene leakage 风险。
  - 新增：`src/data/qc_splice_dataset.py`
  - 输出：`reports/data_qc.md`
  - 检查：类别数量、序列长度、N 比例、中心 motif、hard negative motif 命中率、gene leakage。
  - 要点：质控对象是当前小样本 split，不是完整数据集。
  - 验收：`reports/data_qc.md` 可复现生成，并明确列出 train/valid/test 的样本数和主要 QC 结论。

- [x] 明确 raw data 与模型权重为可选资源。
  - 说明：把 3GB 级 raw 文件和大型模型权重从主链路中剥离，只作为需要真实资源时的手动准备项。
  - 修改：`README.md`
  - 新增：`src/resources/prepare_raw_resources.py`
  - 新增：`src/resources/prepare_foundation_models.py`
  - 输出：`reports/resource_setup.md`
  - 要点：`data/raw/` 和 `models/hf/` 不进 git；主实验不依赖完整 raw 数据。
  - 验收：没有 raw 文件或模型权重时，实验一到三的默认小样本流程仍能运行；资源文档说明如何补齐可选文件。

## 工程与文档

- [x] 补齐轻量环境依赖。
  - 说明：拆分主链路依赖和真实模型可选依赖，让新环境能先运行小样本 synthetic/proxy 实验。
  - 新增：`requirements.txt`
  - 修改：`requirements-real.txt`
  - 要点：`requirements.txt` 支持当前小样本 synthetic/proxy 主链路；`requirements-real.txt` 放真实模型可选依赖。
  - 验收：只安装 `requirements.txt` 时可运行数据检查、proxy 模型和报告生成；真实模型依赖只在需要时安装。

- [x] 增加配置目录。
  - 说明：用配置文件集中管理数据路径、输出路径、模型列表和样本数量，减少各实验脚本中的硬编码。
  - 新增：`configs/exp1_classification.yaml`
  - 新增：`configs/exp2_multiscale.yaml`
  - 新增：`configs/exp3_variant_effect.yaml`
  - 修改：`src/experiments/exp1/run_classification.py`
  - 修改：`src/experiments/exp2/run_multiscale.py`
  - 修改：`src/experiments/exp3/run_variant_effect.py`
  - 要点：配置中的默认数据路径应指向当前小样本 split。
  - 验收：三个实验脚本均支持 `--config`，默认配置不读取全量数据。

- [x] 增加自动化 smoke tests。
  - 说明：增加轻量测试来保护当前仓库结构、路径解析和核心 metrics，避免后续整理目录时再次破坏运行入口。
  - 新增：`tests/test_data_builders.py`
  - 新增：`tests/test_metrics.py`
  - 新增：`tests/test_pipeline_smoke.py`
  - 覆盖：小样本数据读取、variant 构建、metrics 计算、各实验 `--help`。
  - 要点：测试不应下载 raw 数据，也不应生成百万级数据。
  - 验收：`pytest` 能在没有 `data/raw/`、没有外部模型权重的环境中完成 smoke tests。

- [x] 增加模型卡和实验日志。
  - 说明：记录每个模型的数据来源、训练方式、是否 proxy、主要限制和当前实验运行记录，方便后续报告引用。
  - 新增：`reports/model_cards.md`
  - 新增：`reports/experiment_log.md`
  - 修改：`reports/experiment_1.md`
  - 修改：`reports/experiment_2.md`
  - 修改：`reports/experiment_3.md`
  - 验收：模型卡覆盖实验一到三出现的模型；实验日志记录运行日期、命令、输入数据和输出文件。

- [x] 清理报告生成脚本中的过时内容。
  - 说明：更新报告生成脚本，让它只引用当前 `reports/`、`results/experiment_*` 和小样本数据结构中的文件。
  - 修改：`src/reports/write_c_part_report.py`
  - 输出：`reports/experiment_1.md`、`reports/experiment_2.md`、`reports/experiment_3.md` 或统一总报告。
  - 要点：不要再引用不存在的 `REAL_RESOURCE_STATUS.md`、`C_PART_EXECUTION_DETAILS.md`；不要描述全量数据实验。
  - 验收：重新生成报告后不会恢复旧 report 文件，也不会在正文中声称已完成全量数据实验。

## 可选高级项

- [x] 支持 Nucleotide Transformer 或 RNAErnie 的小样本对照。
  - 说明：作为可选扩展模型，只在当前小样本上做 frozen/proxy 对照，用来丰富实验一模型比较。
  - 新增：`src/models/nucleotide_transformer_mlp.py`
  - 新增：`src/models/rnaernie_mlp.py`
  - 修改：`src/experiments/exp1/common.py`
  - 要点：只跑当前小样本，不做全量 fine-tuning。
  - 验收：缺少模型权重时该项可跳过；存在权重时结果追加到实验一 metrics 表并标注 optional。

- [x] RNA-FM / RNABERT fine-tuning 仅作为可选小样本实验。
  - 说明：如果后续需要尝试微调，也只在 855/120/285 split 上做小样本探索，不作为主结论来源。
  - 修改：`src/models/foundation_backbones.py`
  - 修改：`src/experiments/exp1/train.py`
  - 输出：`results/experiment_1/tables/experiment_1_finetune_metrics.csv`
  - 要点：不要求百万级训练。
  - 验收：fine-tuning 默认关闭，需要显式参数启用；报告中说明该结果受小样本规模限制。
