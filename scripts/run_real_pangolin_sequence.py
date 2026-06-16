from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data.build_variant_dataset import build_and_write_variants
from src.models.external_splice_tools import ToolScore, pangolin_pair_delta
from src.utils import EXP3_DATA_DIR, ensure_dirs, exp3_data_file, write_dataframe


def _target_attr(row: pd.Series) -> str:
    name = str(row["target_class_name"]).replace("-", "_")
    if name not in {"donor", "acceptor", "non_splice"}:
        raise ValueError(f"Unexpected target_class_name: {row['target_class_name']!r}")
    return name


def score_variants(variants: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in variants.iterrows():
        target = _target_attr(row)
        loss_signal, gain_signal, impact = pangolin_pair_delta(str(row["wt_sequence"]), str(row["mut_sequence"]))
        if str(row["variant_type"]).endswith("loss"):
            delta = loss_signal
        elif "gain" in str(row["variant_type"]):
            delta = gain_signal
        else:
            delta = impact
        wt = ToolScore(donor=0.0, acceptor=0.0, non_splice=1.0)
        mut = ToolScore(donor=gain_signal, acceptor=loss_signal, non_splice=max(0.0, 1.0 - impact))
        rows.append(
            {
                "variant_id": row["variant_id"],
                "variant_type": row["variant_type"],
                "label": int(row["label"]),
                "label_name": row["label_name"],
                "target_class": int(row["target_class"]),
                "target_class_name": row["target_class_name"],
                "tool": "Pangolin real sequence model",
                "model": "Pangolin real sequence model",
                "source": "real_external_tool",
                "ref_donor": wt.donor,
                "ref_acceptor": wt.acceptor,
                "ref_non_splice": wt.non_splice,
                "alt_donor": mut.donor,
                "alt_acceptor": mut.acceptor,
                "alt_non_splice": mut.non_splice,
                "impact_score": delta,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run real Pangolin sequence scoring on real ClinVar variant windows.")
    parser.add_argument("--variants", type=Path, default=exp3_data_file("clinvar_splicing_variants.csv"))
    parser.add_argument("--out", type=Path, default=EXP3_DATA_DIR / "external_tools")
    args = parser.parse_args()
    ensure_dirs(args.out)
    if not args.variants.exists():
        raise FileNotFoundError(
            f"Variant table not found: {args.variants}. "
            "Build it with python -m src.data.build_clinvar_variant_dataset."
        )
    variants = pd.read_csv(args.variants)
    scored = score_variants(variants)
    write_dataframe(args.out / "pangolin_sequence_variant_scores.csv", scored)
    print(scored.groupby("tool")[["impact_score"]].mean().to_string())


if __name__ == "__main__":
    main()
