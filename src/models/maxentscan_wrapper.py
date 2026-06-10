from __future__ import annotations

from src.utils import acceptor_consensus_score, donor_consensus_score


def maxentscan_variant_delta(wt_sequence: str, mut_sequence: str, target_class: int, variant_type: str) -> tuple[float, float, float]:
    if target_class == 0:
        wt = donor_consensus_score(wt_sequence)
        mut = donor_consensus_score(mut_sequence)
    elif target_class == 1:
        wt = acceptor_consensus_score(wt_sequence)
        mut = acceptor_consensus_score(mut_sequence)
    else:
        wt = max(donor_consensus_score(wt_sequence), acceptor_consensus_score(wt_sequence))
        mut = max(donor_consensus_score(mut_sequence), acceptor_consensus_score(mut_sequence))
    if variant_type.endswith("loss"):
        delta = wt - mut
    elif "gain" in variant_type:
        delta = mut - wt
    else:
        delta = abs(mut - wt)
    return float(wt), float(mut), float(delta)
