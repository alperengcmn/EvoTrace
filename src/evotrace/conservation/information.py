"""Explicit information-theoretic metrics."""

import math
from collections import Counter
from typing import Dict, Iterable, Optional, Set


def shannon_entropy(
    values: Iterable[str],
    gap_policy: str = "ignore",
    alphabet_size: int = 20,
    valid_symbols: Optional[Set[str]] = None,
) -> Dict[str, float]:
    if gap_policy not in ("ignore", "include"):
        raise ValueError("gap_policy must be 'ignore' or 'include'.")
    vals = list(values)
    if valid_symbols is not None:
        vals = [x for x in vals if x in valid_symbols or (gap_policy == "include" and x == "-")]
    elif gap_policy == "ignore":
        vals = [x for x in vals if x not in "-?."]
    if not vals:
        return {
            "entropy": 0.0,
            "normalized_entropy": 0.0,
            "information_content": 0.0,
            "effective_alphabet_size": 1.0,
        }
    counts = Counter(vals)
    n = len(vals)
    entropy = -sum((count / n) * math.log(count / n, 2) for count in counts.values())
    effective_alphabet = (
        alphabet_size + 1
        if gap_policy == "include" and valid_symbols is not None
        else alphabet_size
    )
    maximum = math.log(max(1, effective_alphabet), 2)
    return {
        "entropy": entropy,
        "normalized_entropy": entropy / maximum if maximum else 0.0,
        "information_content": max(0.0, maximum - entropy),
        "effective_alphabet_size": 2**entropy,
    }


def jensen_shannon_divergence(p: Dict[str, float], q: Dict[str, float]) -> float:
    """Base-2 Jensen-Shannon divergence between nonnegative distributions, in bits."""
    keys = set(p) | set(q)
    if any(p.get(k, 0) < 0 or q.get(k, 0) < 0 for k in keys):
        raise ValueError("Distribution weights must be nonnegative.")
    ps, qs = sum(p.values()), sum(q.values())
    if ps <= 0 or qs <= 0:
        raise ValueError("Both distributions must have positive mass.")
    pn = {k: p.get(k, 0) / ps for k in keys}
    qn = {k: q.get(k, 0) / qs for k in keys}
    m = {k: (pn[k] + qn[k]) / 2 for k in keys}

    def kl(a: Dict[str, float], b: Dict[str, float]) -> float:
        return sum(v * math.log(v / b[k], 2) for k, v in a.items() if v > 0)

    return (kl(pn, m) + kl(qn, m)) / 2
