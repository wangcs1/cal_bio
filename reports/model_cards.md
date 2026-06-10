# Model Cards

## CNN baseline

- Type: supervised small-sample classifier.
- Data: `data/shared/splits/train.csv` and `valid.csv`.
- Role: local motif and short-context baseline.
- Limitation: may over-rely on local motif-like features.

## RNA-FM frozen encoder + MLP

- Type: optional frozen encoder with deterministic proxy fallback.
- Local weights: `models/hf/rnafm/`.
- Main-pipeline behavior: if weights/dependencies are absent, k-mer/signal proxy embeddings are used.
- Limitation: proxy fallback is not a real pretrained RNA-FM result.

## RNABERT frozen encoder + MLP

- Type: optional frozen encoder with deterministic proxy fallback.
- Local weights: `models/hf/rnabert/`.
- Main-pipeline behavior: if weights/dependencies are absent, token/position/signal proxy embeddings are used.
- Limitation: proxy fallback is not a real pretrained RNABERT result.

## SpliceAI optional real tool

- Type: optional task-specific splice model wrapper.
- Main-pipeline behavior: proxy fallback maps donor/acceptor signal features to three-class probabilities.
- Limitation: fallback is a case-study compatible proxy, not full SpliceAI inference.

## Pangolin optional tool

- Type: optional long-context splice tool wrapper.
- Main-pipeline behavior: small hard-negative/tissue case-study proxy.
- Limitation: not a full Pangolin benchmark unless external dependency and inputs are explicitly supplied.

## MMSplice / MaxEntScan

- Type: optional splice variant tools or consensus-strength proxies.
- Main-pipeline behavior: deterministic small-sample delta score fallback.
- Limitation: no full ClinVar benchmark is claimed.

## Borzoi / AlphaGenome

- Type: optional long-range regulatory case studies.
- Main-pipeline behavior: documented small example outputs.
- Limitation: external API/weights are not required and no full benchmark is claimed.
