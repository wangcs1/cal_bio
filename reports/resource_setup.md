# Optional Resource Setup

The main C Part experiments use the small synthetic/proxy split and do not require raw genome files.
`data/raw/` is reserved for optional real-data smoke tests and case studies.

| Resource | Status | Purpose |
| --- | --- | --- |
| `genome.fa` | present | optional real-resource smoke or case study |
| `gencode.gtf` | present | optional real-resource smoke or case study |
| `clinvar.vcf` | present | optional real-resource smoke or case study |
| `gtex_sqtl.tsv` | present | optional real-resource smoke or case study |
| `known_splice_events.tsv` | present | optional real-resource smoke or case study |
| `gencode.db` | present | optional real-resource smoke or case study |

Large raw files and `models/hf/` pretrained weights are intentionally excluded from git.
Default experiment commands remain runnable without these resources.
