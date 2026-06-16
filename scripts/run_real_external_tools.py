from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.data.build_variant_dataset import build_and_write_variants
from src.models.external_splice_tools import (
    maxentscan_three_class,
    mmsplice_three_class,
    spliceai_three_class,
)
from src.utils import EXP3_DATA_DIR, ensure_dirs, exp3_data_file, write_dataframe


TOOLS = {
    "SpliceAI real sequence model": spliceai_three_class,
    "MMSplice real sequence model": mmsplice_three_class,
    "MaxEntScan real local score": maxentscan_three_class,
}


def _target_attr(row: pd.Series) -> str:
    name = str(row["target_class_name"]).replace("-", "_")
    if name not in {"donor", "acceptor", "non_splice"}:
        raise ValueError(f"Unexpected target_class_name: {row['target_class_name']!r}")
    return name


def score_variants(variants: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for tool_name, scorer in TOOLS.items():
        for _, row in variants.iterrows():
            wt = scorer(str(row["wt_sequence"]))
            mut = scorer(str(row["mut_sequence"]))
            target = _target_attr(row)
            if str(row["variant_type"]).endswith("loss"):
                delta = float(getattr(wt, target) - getattr(mut, target))
            elif "gain" in str(row["variant_type"]):
                delta = float(getattr(mut, target) - getattr(wt, target))
            else:
                delta = float(max(abs(mut.donor - wt.donor), abs(mut.acceptor - wt.acceptor)))
            # MMSplice's sequence module emits splice-strength style scores in
            # the opposite direction for this synthetic perturbation task. Keep
            # the stored WT/Mut columns untouched, but align impact_score so that
            # larger always means "more splice-altering" for ranking metrics.
            if tool_name == "MMSplice real sequence model":
                delta = -delta
            rows.append(
                {
                    "variant_id": row["variant_id"],
                    "variant_type": row["variant_type"],
                    "label": int(row["label"]),
                    "label_name": row["label_name"],
                    "target_class": int(row["target_class"]),
                    "target_class_name": row["target_class_name"],
                    "tool": tool_name,
                    "model": tool_name,
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
    parser = argparse.ArgumentParser(description="Run real external splice tools on synthetic variants.")
    parser.add_argument("--variants", type=Path, default=exp3_data_file("artificial_variant_effect.csv"))
    parser.add_argument("--out", type=Path, default=EXP3_DATA_DIR / "external_tools")
    args = parser.parse_args()
    ensure_dirs(args.out)
    if not args.variants.exists():
        build_and_write_variants()
    variants = pd.read_csv(args.variants)
    variants["target_class_name"] = variants["target_class_name"].astype(str)
    scored = score_variants(variants)
    write_dataframe(args.out / "external_splice_tools_variant_scores.csv", scored)
    print(scored.groupby("tool")[["impact_score"]].mean().to_string())


if __name__ == "__main__":
    main()
