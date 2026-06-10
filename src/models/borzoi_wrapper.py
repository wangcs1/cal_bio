from __future__ import annotations

import pandas as pd


def build_long_range_case_study() -> pd.DataFrame:
    rows = [
        {
            "case_id": "LR_CASE_001",
            "model": "Borzoi optional long-range case study",
            "input_scope": "synthetic 10kb-centered locus sketch",
            "predicted_track": "splice_junction_usage",
            "wt_signal": 0.72,
            "mut_signal": 0.41,
            "delta_signal": -0.31,
            "status": "documented proxy case; real model not required for main pipeline",
        },
        {
            "case_id": "LR_CASE_002",
            "model": "AlphaGenome optional long-range case study",
            "input_scope": "synthetic regulatory-context sketch",
            "predicted_track": "RNA-seq coverage / splice usage",
            "wt_signal": 0.36,
            "mut_signal": 0.58,
            "delta_signal": 0.22,
            "status": "documented proxy case; API/weights optional",
        },
    ]
    return pd.DataFrame(rows)
