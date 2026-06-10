# AlphaGenome Optional Case Study

This repository keeps AlphaGenome as an optional long-range regulatory case study.
It is not required for the small-sample C Part main pipeline.

Expected inputs:

- A centered genomic interval around a splice-relevant event.
- Reference and alternate alleles for a variant case.
- Optional tissue or assay context if an external API supports it.

Expected outputs:

- Delta RNA-seq coverage or splice-junction usage.
- A short interpretation of whether the model supports donor loss, acceptor loss,
  donor gain, acceptor gain, or no strong splice effect.

Current status: documented optional case study; no full benchmark is claimed.
