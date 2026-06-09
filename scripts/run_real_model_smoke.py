from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW = PROJECT_ROOT / "data/raw"
OUT = PROJECT_ROOT / "results/real_smoke"


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def ensure_clinvar_smoke(limit: int = 5) -> None:
    subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "scripts/make_clinvar_smoke.py"), "--limit", str(limit)],
        cwd=PROJECT_ROOT,
        check=True,
    )


def run_spliceai() -> None:
    out_vcf = OUT / "spliceai_clinvar_smoke.vcf"
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = "-1"
    env["TF_CPP_MIN_LOG_LEVEL"] = "2"
    subprocess.run(
        [
            "spliceai",
            "-I",
            str(RAW / "clinvar_smoke.vcf"),
            "-O",
            str(out_vcf),
            "-R",
            str(RAW / "genome.fa"),
            "-A",
            "grch38",
            "-D",
            "50",
        ],
        cwd=PROJECT_ROOT,
        env=env,
        check=True,
    )
    rows: list[dict[str, object]] = []
    with out_vcf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = fields[7] if len(fields) > 7 else ""
            spliceai = ""
            for item in info.split(";"):
                if item.startswith("SpliceAI="):
                    spliceai = item.split("=", 1)[1]
                    break
            rows.append(
                {
                    "chrom": fields[0],
                    "pos": fields[1],
                    "ref": fields[3],
                    "alt": fields[4],
                    "spliceai_info": spliceai,
                }
            )
    write_rows(
        OUT / "spliceai_clinvar_smoke_summary.csv",
        rows,
        ["chrom", "pos", "ref", "alt", "spliceai_info"],
    )


def run_pangolin() -> None:
    subprocess.run(
        [
            "pangolin",
            str(RAW / "clinvar_smoke.csv"),
            str(RAW / "genome.fa"),
            str(RAW / "gencode.db"),
            str(OUT / "pangolin_clinvar_smoke"),
            "-c",
            "CHROM,POS,REF,ALT",
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )


def run_foundation_models() -> None:
    import torch
    from multimolecule import RnaBertModel, RnaFmModel, RnaTokenizer

    rows: list[dict[str, object]] = []
    sequence = "AUGCUAGCUAGCUAGCUAGCUAGCUAGC"
    for repo, model_cls in [
        ("multimolecule/rnafm", RnaFmModel),
        ("multimolecule/rnabert", RnaBertModel),
    ]:
        tokenizer = RnaTokenizer.from_pretrained(repo)
        model = model_cls.from_pretrained(repo)
        device = "cuda" if torch.cuda.is_available() else "cpu"
        model.to(device).eval()
        inputs = tokenizer([sequence], return_tensors="pt")
        inputs = {name: tensor.to(device) for name, tensor in inputs.items()}
        with torch.no_grad():
            output = model(**inputs)
        hidden = output.last_hidden_state
        rows.append(
            {
                "model": repo,
                "device": device,
                "input_tokens": int(hidden.shape[1]),
                "hidden_size": int(hidden.shape[2]),
                "embedding_mean": float(hidden.mean().detach().cpu()),
                "embedding_std": float(hidden.std().detach().cpu()),
            }
        )
    write_rows(
        OUT / "foundation_model_smoke.csv",
        rows,
        ["model", "device", "input_tokens", "hidden_size", "embedding_mean", "embedding_std"],
    )


def run_maxent_mmsplice_worker() -> None:
    from maxentpy import maxent
    from mmsplice import MMSplice

    donor_seq = "CAGGTAAGT"
    acceptor_seq = "TTTTTTTTTTTTTTTTTTTTAGG"
    mmsplice_seq = ("A" * 80) + acceptor_seq[:-1] + ("A" * 40) + donor_seq + ("A" * 80)
    mmsplice_scores = MMSplice().predict_on_seq(mmsplice_seq, overhang=(80, 80))
    rows = [
        {
            "model": "MaxEntScan",
            "input": donor_seq,
            "score_name": "score5_donor",
            "score": float(maxent.score5(donor_seq)),
        },
        {
            "model": "MaxEntScan",
            "input": acceptor_seq,
            "score_name": "score3_acceptor",
            "score": float(maxent.score3(acceptor_seq)),
        },
    ]
    for index, score in enumerate(mmsplice_scores):
        rows.append(
            {
                "model": "MMSplice",
                "input": f"synthetic_splice_context_overhang_80_{index}",
                "score_name": f"module_{index}",
                "score": float(score),
            }
    )
    write_rows(OUT / "maxentscan_mmsplice_smoke.csv", rows, ["model", "input", "score_name", "score"])


def run_maxent_mmsplice() -> None:
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = "-1"
    env["TF_CPP_MIN_LOG_LEVEL"] = "2"
    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--mmsplice-worker"], env=env, check=True)


def main() -> None:
    if "--mmsplice-worker" in sys.argv:
        OUT.mkdir(parents=True, exist_ok=True)
        run_maxent_mmsplice_worker()
        return

    OUT.mkdir(parents=True, exist_ok=True)
    ensure_clinvar_smoke(limit=5)
    run_spliceai()
    run_pangolin()
    run_foundation_models()
    run_maxent_mmsplice()
    print(f"Wrote real model smoke outputs under {OUT}")


if __name__ == "__main__":
    main()
