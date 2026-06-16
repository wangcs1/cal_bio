# Model Cards

## CNN baseline

- Type: supervised small-sample classifier.
- Data: `data/shared/splits/train.csv` and `valid.csv`.
- Role: local motif and short-context baseline.
- Limitation: may over-rely on local motif-like features.

## RNA-FM frozen encoder + MLP

- Type: real local pretrained frozen encoder plus LogisticRegression classifier.
- Local weights: `models/hf/rnafm/`.
- Main-pipeline behavior: loads real local weights and fails if they are unavailable.
- Limitation: evaluated only on the current small synthetic split.

## RNABERT frozen encoder + MLP

- Type: real local pretrained frozen encoder plus LogisticRegression classifier.
- Local weights: `models/hf/rnabert/`.
- Main-pipeline behavior: loads real local weights and fails if they are unavailable.
- Limitation: evaluated only on the current small synthetic split.

## External splice tools

- SpliceAI: real package weights from `spliceai` are executed in the Python 3.10 splice-tools environment on padded synthetic WT/Mut sequences.
- Pangolin: real package weights are executed through the installed Pangolin/PyTorch environment on padded synthetic WT/Mut sequences.
- MMSplice: real `MMSplice().predict_on_seq` module scores are executed in the Python 3.10 splice-tools environment.
- MaxEntScan: real `maxentpy` score5/score3 local splice-window scores are executed in the Python 3.10 splice-tools environment.
- Boundary: these are real tool runs on synthetic sequence inputs, not full ClinVar/GTEx coordinate benchmarks.
