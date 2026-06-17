# 数据质控结果

QC target: real GENCODE/GRCh38 splice-site benchmark plus real ClinVar variant benchmark.

Important sampling definitions:

- Splice-site benchmark: donor / acceptor / non-splice hard-negative = 1:1:1 globally; binary positive base rate is 2000/3000 = 0.667.
- Negative splice-site examples: all real benchmark negatives are motif-matched GT/AG hard negatives; easy negatives are kept only for supplementary synthetic controls.
- ClinVar benchmark: pathogenic or likely pathogenic splice-related/near-splice SNVs vs benign or likely benign near-splice SNVs, balanced 250/250.
- ClinVar variants are restricted to held-out test chromosomes chr19, chr20, chr21, chr22, chrX and exclude genes already sampled in experiment 1/2 splits.

## Checklist Status

| section                        | status   | evidence                                                                                                                                                                    |
|:-------------------------------|:---------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1 Split integrity              | PASS     | Four windows share identical sample_id/split assignments; train/valid/test chromosomes and IDs are disjoint.                                                                |
| 2 Class balance                | PASS     | Three-class sampling is donor:acceptor:hard-negative = 1:1:1 globally; easy negatives are intentionally not used in the real benchmark.                                     |
| 3 Sequence QC                  | PASS     | Window lengths, alphabet, N fraction, and canonical donor/acceptor motif checks passed.                                                                                     |
| 4 Hard negatives               | PASS     | Hard negatives carry GT/AG decoys; overlap with GTF annotated splice sites = 0.                                                                                             |
| 5 Coordinates and strand       | PASS     | Positive motif checks pass after transcript-oriented reverse complement, including negative-strand records.                                                                 |
| 6 ClinVar                      | PASS     | ClinVar is balanced 250/250 on held-out chromosomes; REF and one-SNV WT/Mut checks passed; no sampled split gene is reused.                                                 |
| 7 Paralog leakage              | PASS     | Cross-split full-window duplicates, center-window duplicates, and high 9-mer Jaccard near-duplicates were audited; alignment-level paralog clustering remains out of scope. |
| 8 Reproducibility              | PASS     | Seed 42, GRCh38/GENCODE metadata, ClinVar fileDate, commands, and local resource status are recorded.                                                                       |
| 9 Cross-experiment consistency | PASS     | Experiments use GRCh38; ClinVar variants are restricted to the same held-out test chromosomes used by experiment 1/2.                                                       |

## 1. Split 完整性与防泄漏

Window-level identity and assignment check:

| window   |   rows |   duplicate_sample_ids | same_ids_as_pm50   | same_assignment_as_pm50   |
|:---------|-------:|-----------------------:|:-------------------|:--------------------------|
| pm50     |   3000 |                      0 | True               | True                      |
| pm100    |   3000 |                      0 | True               | True                      |
| pm200    |   3000 |                      0 | True               | True                      |
| pm400    |   3000 |                      0 | True               | True                      |

Split-file ID consistency across windows:

| split   |   pm50_rows | pm100_same_ids   | pm200_same_ids   | pm400_same_ids   |
|:--------|------------:|:-----------------|:-----------------|:-----------------|
| train   |        2339 | True             | True             | True             |
| valid   |         230 | True             | True             | True             |
| test    |         431 | True             | True             | True             |

| pair        |   sample_id_overlap | chrom_overlap   |   gene_overlap |
|:------------|--------------------:|:----------------|---------------:|
| train/valid |                   0 |                 |              0 |
| train/test  |                   0 |                 |              0 |
| valid/test  |                   0 |                 |              0 |

Conclusion: 四个窗口共用同一批 `sample_id -> split`；split 不是按行随机切。train/valid/test 的 sample_id、chromosome、gene 均互斥。

## 2. 类别构成与平衡

| split   |   rows |   donor |   acceptor |   easy_negative |   hard_negative |   positive_rate |   hard_negative_rate_among_negatives | chromosomes                                                                            |
|:--------|-------:|--------:|-----------:|----------------:|----------------:|----------------:|-------------------------------------:|:---------------------------------------------------------------------------------------|
| train   |   2339 |     792 |        786 |               0 |             761 |          0.6746 |                                    1 | chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9 |
| valid   |    230 |      79 |         67 |               0 |              84 |          0.6348 |                                    1 | chr17,chr18                                                                            |
| test    |    431 |     129 |        147 |               0 |             155 |          0.6404 |                                    1 | chr19,chr20,chr21,chr22,chrX                                                           |

Note: 真实主 benchmark 没有 easy negative；non-splice 类全部是 GT/AG hard negative。三分类全量为 1000/1000/1000，但按染色体 holdout 后各 split 比例会随染色体分布略有变化。

## 3. 序列层面 QC

| window   |   rows |   expected_length |   bad_length_rows | alphabet   |   n_fraction |   donor_gt_rate |   acceptor_ag_rate |   hard_gt_or_ag_rate |
|:---------|-------:|------------------:|------------------:|:-----------|-------------:|----------------:|-------------------:|---------------------:|
| pm50     |   3000 |               101 |                 0 | ACGT       |            0 |               1 |                  1 |                    1 |
| pm100    |   3000 |               201 |                 0 | ACGT       |            0 |               1 |                  1 |                    1 |
| pm200    |   3000 |               401 |                 0 | ACGT       |            0 |               1 |                  1 |                    1 |
| pm400    |   3000 |               801 |                 0 | ACGT       |            0 |               1 |                  1 |                    1 |

N strategy: rows with N fraction above the builder threshold are filtered; generated benchmark has N fraction 0.

## 4. 负例 / Hard-Negative 有效性

- pm200 hard-negative rows: 1000
- hard-negative center GT/AG motif rate: 1.0000
- overlap with any GTF annotated splice site: 0

GC fraction by category on pm200:

| category      |   rows |   mean |    std |    min |    max |
|:--------------|-------:|-------:|-------:|-------:|-------:|
| acceptor      |   1000 | 0.4599 | 0.1051 | 0.2469 | 0.7905 |
| donor         |   1000 | 0.4835 | 0.1145 | 0.2519 | 0.808  |
| hard_negative |   1000 | 0.4492 | 0.0933 | 0.1571 | 0.7456 |

## 5. 坐标与链

Donor GT and acceptor AG rates are 1.0 after transcript-oriented extraction, which indirectly checks GTF 1-based coordinates and negative-strand reverse-complement handling.

## 6. ClinVar 变异表

| label_name      | variant_type            | target_class_name   |   rows |
|:----------------|:------------------------|:--------------------|-------:|
| neutral         | clinvar_benign_snv      | acceptor            |    124 |
| neutral         | clinvar_benign_snv      | donor               |    126 |
| splice_altering | acceptor_clinvar_splice | acceptor            |    141 |
| splice_altering | donor_clinvar_splice    | donor               |    109 |

ClinVar distance to nearest annotated splice site:

| label_name      |   count |   mean |     std |   min |   25% |   50% |   75% |   max |
|:----------------|--------:|-------:|--------:|------:|------:|------:|------:|------:|
| neutral         |     250 | 20.34  | 13.9948 |     0 |     9 |  17   |    29 |    50 |
| splice_altering |     250 | 14.412 | 14.5075 |     0 |     2 |   9.5 |    24 |    50 |

Distance confounding diagnostic:

The pathogenic/splice-altering variants are closer to annotated splice sites than benign variants on average. This is biologically plausible, but it is also a confounder: a model can partially solve the task by learning `closer to splice site = more likely pathogenic`. To make this visible, QC reports a distance-only baseline and exports an exact-distance-matched subset.

| subset                 |   rows |   positive_rows |   neutral_rows |   distance_only_auroc |   distance_only_auprc |
|:-----------------------|-------:|----------------:|---------------:|----------------------:|----------------------:|
| full_clinvar           |    500 |             250 |            250 |                0.6475 |                0.6917 |
| exact_distance_matched |    326 |             163 |            163 |                0.5    |                0.5    |

Exact-distance matching summary:

|   nearest_splice_distance |   positive_rows |   neutral_rows |   matched_pairs |
|--------------------------:|----------------:|---------------:|----------------:|
|                         0 |               3 |              3 |               3 |
|                         1 |              56 |              3 |               3 |
|                         2 |              27 |              5 |               5 |
|                         3 |               8 |              9 |               8 |
|                         4 |               3 |              8 |               3 |
|                         5 |               5 |              6 |               5 |
|                         6 |               8 |              8 |               8 |
|                         7 |               7 |             12 |               7 |
|                         8 |               1 |              8 |               1 |
|                         9 |               7 |             11 |               7 |
|                        10 |               6 |              6 |               6 |
|                        11 |               4 |              9 |               4 |
|                        12 |               6 |              6 |               6 |
|                        13 |               6 |              5 |               5 |
|                        14 |               5 |              7 |               5 |
|                        15 |               3 |              4 |               3 |
|                        16 |               4 |              8 |               4 |
|                        17 |               4 |              9 |               4 |
|                        18 |               2 |              6 |               2 |
|                        19 |               6 |              4 |               4 |
|                        20 |               5 |             11 |               5 |
|                        21 |               4 |              3 |               3 |
|                        22 |               3 |              4 |               3 |
|                        23 |               3 |              4 |               3 |
|                        24 |               4 |              4 |               4 |
|                        25 |               1 |              6 |               1 |
|                        26 |               2 |              2 |               2 |
|                        27 |               2 |              4 |               2 |
|                        28 |               6 |              6 |               6 |
|                        29 |               4 |              7 |               4 |
|                        31 |               3 |              4 |               3 |
|                        32 |               2 |              3 |               2 |
|                        33 |               1 |              3 |               1 |
|                        34 |               3 |              2 |               2 |
|                        35 |               3 |              5 |               3 |
|                        37 |               1 |              1 |               1 |
|                        38 |               3 |              3 |               3 |
|                        40 |               3 |              5 |               3 |
|                        42 |               3 |              5 |               3 |
|                        43 |               3 |              5 |               3 |
|                        44 |               1 |              3 |               1 |
|                        45 |               2 |              2 |               2 |
|                        46 |               3 |              4 |               3 |
|                        48 |               3 |              3 |               3 |
|                        49 |               4 |              3 |               3 |
|                        50 |               1 |              1 |               1 |

Distance-matched ClinVar subset: `D:\CAL_BIO\data\experiment_3\clinvar_splicing_variants_distance_matched.csv`.

Model metrics on the distance-matched subset, if experiment-3 scores have been generated:

| status    | model                         |   rows |   auroc_distance_matched |   auprc_distance_matched |
|:----------|:------------------------------|-------:|-------------------------:|-------------------------:|
| available | RNA-FM frozen encoder + MLP   |    326 |                   0.6684 |                   0.6559 |
| available | Pangolin real sequence model  |    326 |                   0.5902 |                   0.6416 |
| available | CNN baseline (PyTorch Conv1D) |    326 |                   0.558  |                   0.565  |
| available | RNABERT frozen encoder + MLP  |    326 |                   0.5155 |                   0.5442 |
| available | SpliceAI real sequence model  |    326 |                   0.5135 |                   0.5693 |
| available | MMSplice real sequence model  |    326 |                   0.5    |                   0.5    |
| available | MaxEntScan real local score   |    326 |                   0.4668 |                   0.5496 |

ClinVar overlap with experiment 1/2 split genes/chromosomes:

| split   | chrom_overlap   |   gene_overlap |
|:--------|:----------------|---------------:|
| train   |                 |              0 |
| valid   |                 |              0 |
| test    | chr19           |              0 |

Interpretation: `chrom_overlap=chr19` for the test split is expected because ClinVar is restricted to held-out test chromosomes. `gene_overlap=0` for train/valid/test confirms that no sampled experiment-1/2 gene is reused in the ClinVar variant benchmark.

ClinVar allele checks:

|   genomic_ref_mismatches |   hamming_failures |   position_failures |
|-------------------------:|-------------------:|--------------------:|
|                        0 |                  0 |                   0 |

ClinVar filtering rule: SNVs only; positive labels are pathogenic/likely pathogenic records with splice-related consequence or within the near-splice threshold; negative labels are benign/likely benign SNVs within the same near-splice threshold. REF is checked against GRCh38 and strand-oriented REF/ALT are checked in WT/Mut windows.

Gain classes: not modeled as separate donor_gain/acceptor_gain labels in the real ClinVar benchmark. The experiment-3 story should be reported as real ClinVar splice-altering vs benign ranking, stratified by donor/acceptor target, not as synthetic gain/loss recovery.

## 7. 同源 / 旁系泄漏

Cross-split near-duplicate leakage audit:

| pair        |   full_sequence_duplicate_hashes |   center_161bp_duplicate_hashes |   kmer_jaccard_ge_0.90_pairs |
|:------------|---------------------------------:|--------------------------------:|-----------------------------:|
| train/valid |                                0 |                               0 |                            0 |
| train/test  |                                0 |                               0 |                            0 |
| valid/test  |                                0 |                               0 |                            0 |

High-similarity examples above the audit threshold, if any:

_No rows._

Interpretation: this audit checks exact full-window duplicates, exact center 161 bp duplicates, and high 9-mer Jaccard near-duplicates across train/valid/test. It is a practical homology-leakage screen for the current benchmark, but it is not a replacement for alignment-level paralog clustering with tools such as BLAST/CD-HIT/MMseqs.

## 8. 复现与溯源

| resource    | field       | value                                                                            |
|:------------|:------------|:---------------------------------------------------------------------------------|
| gencode.gtf | description | evidence-based annotation of the human genome (GRCh38), version 49 (Ensembl 115) |
| gencode.gtf | date        | 2025-07-08                                                                       |
| clinvar.vcf | fileDate    | 2026-05-30                                                                       |
| clinvar.vcf | source      | ClinVar                                                                          |
| clinvar.vcf | reference   | GRCh38                                                                           |

| step        | command                                                                                    |
|:------------|:-------------------------------------------------------------------------------------------|
| splice_site | python -m src.data.build_splice_site_dataset --max-per-class 1000 --windows 50 100 200 400 |
| clinvar     | python -m src.data.build_clinvar_variant_dataset                                           |
| qc          | python -m src.data.qc_splice_dataset                                                       |

Random seed: 42.

## 9. 跨实验一致性

Experiment 1/2 split uses the local GRCh38 + GENCODE annotation and chromosome holdout. Experiment 3 ClinVar windows use the same genome/annotation coordinate system and are restricted to held-out test chromosomes.

