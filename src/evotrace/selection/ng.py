"""Descriptive pairwise Nei-Gojobori dN/dS for ungapped codon alignments.

This is an exploratory pairwise estimate, not a site/branch selection test.
"""

from itertools import permutations
from typing import Dict, List
from Bio.Data import CodonTable

TABLE = CodonTable.unambiguous_dna_by_name["Standard"]
CODONS = set(TABLE.forward_table) | set(TABLE.stop_codons)


def _syn_sites(codon):
    aa = TABLE.forward_table.get(codon, "*")
    total = 0.0
    for i in range(3):
        same = 0
        valid = 0
        for base in "ACGT":
            if base == codon[i]:
                continue
            alt = codon[:i] + base + codon[i + 1 :]
            if alt in TABLE.stop_codons:
                continue
            valid += 1
            if TABLE.forward_table.get(alt) == aa:
                same += 1
        total += same / valid if valid else 0.0
    return total


def _pair_changes(a, b):
    changed = [i for i in range(3) if a[i] != b[i]]
    if not changed:
        return 0.0, 0.0
    paths = []
    for order in permutations(changed):
        cur = list(a)
        syn = non = 0.0
        valid = True
        for i in order:
            nxt = cur[:]
            nxt[i] = b[i]
            c1, c2 = "".join(cur), "".join(nxt)
            aa1, aa2 = TABLE.forward_table.get(c1), TABLE.forward_table.get(c2)
            if aa1 is None or aa2 is None:
                valid = False
                break
            if aa1 == aa2:
                syn += 1
            else:
                non += 1
            cur = nxt
        if valid:
            paths.append((syn, non))
    if not paths:
        return None
    return (sum(x[0] for x in paths) / len(paths), sum(x[1] for x in paths) / len(paths))


def _jc(p):
    import math

    if p >= 0.75:
        return None
    return -0.75 * math.log(1 - 4 * p / 3)


def pairwise_dnds(alignment) -> List[Dict[str, object]]:
    if any(len(r.sequence) % 3 for r in alignment):
        return [
            {
                "status": "not_run",
                "reason": "Alignment length is not divisible by three; provide a codon-aware alignment.",
            }
        ]
    results = []
    for i, a in enumerate(alignment):
        for b in alignment[i + 1 :]:
            syn_sites = non_sites = syn_changes = non_changes = 0.0
            used = 0
            for k in range(0, len(a.sequence), 3):
                ca, cb = a.sequence[k : k + 3], b.sequence[k : k + 3]
                if (
                    ca == "---"
                    or cb == "---"
                    or "-" in ca
                    or "-" in cb
                    or set(ca + cb) - set("ACGT")
                ):
                    continue
                if ca in TABLE.stop_codons or cb in TABLE.stop_codons:
                    continue
                sa, sb = _syn_sites(ca), _syn_sites(cb)
                syn_sites += (sa + sb) / 2
                non_sites += 3 - (sa + sb) / 2
                changes = _pair_changes(ca, cb)
                if changes is None:
                    continue
                syn_changes += changes[0]
                non_changes += changes[1]
                used += 1
            ps = syn_changes / syn_sites if syn_sites else None
            pn = non_changes / non_sites if non_sites else None
            ds, dn = (_jc(ps) if ps is not None else None), (_jc(pn) if pn is not None else None)
            ratio = dn / ds if ds is not None and ds > 0 and dn is not None else None
            results.append(
                {
                    "sequence_a": a.id,
                    "sequence_b": b.id,
                    "codons_compared": used,
                    "dN": dn,
                    "dS": ds,
                    "dN_dS": ratio,
                    "status": "descriptive_only",
                    "interpretation": "Pairwise descriptive estimate; not a statistical positive-selection test.",
                }
            )
    return results
