from __future__ import annotations

import numpy as np
import pandas as pd

from src.utils import hard_negative_details, multiclass_metrics, topk_enrichment_curve


def test_hard_negative_details_uses_only_marked_hard_rows():
    frame = pd.DataFrame(
        {
            "label": [2, 2, 2, 0],
            "negative_type": ["hard_gtag", "easy_random", "hard_ag", "positive"],
        }
    )
    proba = np.asarray(
        [
            [0.8, 0.1, 0.1],
            [0.8, 0.1, 0.1],
            [0.1, 0.1, 0.8],
            [0.9, 0.05, 0.05],
        ]
    )
    details = hard_negative_details(frame, proba)
    assert details["hard_negative_rows"] == 2
    assert details["hard_negative_false_positives"] == 1
    assert details["hard_negative_fpr"] == 0.5


def test_metrics_and_topk_are_finite_for_small_input():
    y = [0, 1, 2]
    proba = np.asarray([[0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.2, 0.7]])
    metrics = multiclass_metrics(y, proba)
    assert metrics["macro_f1"] == 1.0
    curve = topk_enrichment_curve([1, 0, 1, 0], [0.9, 0.8, 0.7, 0.1], [0.5])
    assert int(curve.iloc[0]["k"]) == 2
