# Real Resource Setup

The main C Part experiments now use real local resources.
`data/raw/genome.fa`, `data/raw/gencode.gtf`, and `data/raw/clinvar.vcf` are required for the default real-data pipeline.

| Resource | Status | Purpose |
| --- | --- | --- |
| `genome.fa` | present (3,151,417,447 bytes) | default real benchmark |
| `gencode.gtf` | present (4,688,772,094 bytes) | default real benchmark |
| `clinvar.vcf` | present (1,925,472,072 bytes) | default real benchmark |
| `gtex_sqtl.tsv` | optional / missing | optional real-resource case study |
| `known_splice_events.tsv` | optional / missing | optional real-resource case study |
| `gencode.db` | optional / missing | optional real-resource case study |

Large raw files and `models/hf/` pretrained weights are intentionally excluded from git.
Default experiment commands fail fast when required real resources are missing.
