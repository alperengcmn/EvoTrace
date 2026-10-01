"""Alignment-level conservation and entropy metrics."""

from collections import Counter
from typing import Any, Dict, List
from evotrace.models import SequenceRecord
from evotrace.conservation.information import shannon_entropy


def columns(
    alignment: List[SequenceRecord], alphabet: str = "protein", gap_policy: str = "ignore"
) -> List[Dict[str, Any]]:
    if not alignment:
        return []
    if gap_policy not in ("ignore", "include"):
        raise ValueError("gap_policy must be 'ignore' or 'include'.")
    n = len(alignment)
    length = len(alignment[0].sequence)
    rows = []
    for i in range(length):
        col = [r.sequence[i] for r in alignment]
        canonical = set("ACGT" if alphabet.lower() == "dna" else "ACDEFGHIKLMNPQRSTVWY")
        residues = [x for x in col if x in canonical]
        observed = residues + (["-" for x in col if x == "-"] if gap_policy == "include" else [])
        counts = Counter(observed)
        total = len(residues)
        dominant, count = Counter(residues).most_common(1)[0] if total else ("-", 0)
        alphabet_size = 4 if alphabet.lower() == "dna" else 20
        info = shannon_entropy(
            col, gap_policy=gap_policy, alphabet_size=alphabet_size, valid_symbols=canonical
        )
        entropy = info["entropy"]
        rows.append(
            {
                "position": i + 1,
                "dominant_residue": dominant,
                "dominant_frequency": count / total if total else 0.0,
                "entropy": entropy,
                "normalized_entropy": info["normalized_entropy"],
                "information_content": info["information_content"],
                "effective_alphabet_size": info["effective_alphabet_size"],
                "gap_fraction": sum(x == "-" for x in col) / n,
                "ambiguous_fraction": sum(x not in canonical and x != "-" for x in col) / n,
                "occupancy": total / n,
                "conservation_score": count / total if total else 0.0,
                "residue_frequencies": dict(counts),
            }
        )
    return rows


def alignment_stats(
    alignment: List[SequenceRecord], alphabet: str = "protein", gap_policy: str = "ignore"
) -> Dict[str, object]:
    c = columns(alignment, alphabet, gap_policy)
    n = len(alignment)
    pairs = []
    similarities = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = alignment[i].sequence, alignment[j].sequence
            comparable = [(x, y) for x, y in zip(a, b) if x not in "-?" and y not in "-?"]
            ident = sum(x == y for x, y in comparable) / len(comparable) if comparable else 0.0
            pairs.append(ident)
            if alphabet.lower() == "protein":
                from Bio.Align import substitution_matrices

                matrix = substitution_matrices.load("BLOSUM62")
                valid = [
                    (x, y) for x, y in comparable if x in matrix.alphabet and y in matrix.alphabet
                ]
                similarities.append(
                    sum(matrix[x, y] >= 0 for x, y in valid) / len(valid) if valid else 0.0
                )
            else:
                similarities.append(ident)
    L = len(alignment[0].sequence) if alignment else 0
    return {
        "sequence_count": n,
        "alignment_length": L,
        "gap_fraction": sum(float(x["gap_fraction"]) for x in c) / max(1, L),
        "conserved_columns": sum(
            float(x["conservation_score"]) == 1 and float(x["occupancy"]) > 0 for x in c
        ),
        "variable_columns": sum(len(x["residue_frequencies"]) > 1 for x in c),
        "mean_pairwise_identity": sum(pairs) / len(pairs) if pairs else None,
        "mean_pairwise_similarity": sum(similarities) / len(similarities) if similarities else None,
        "gap_policy": gap_policy,
        "pairwise_identity_values": pairs,
    }


def identity_matrix(
    alignment: List[SequenceRecord], alphabet: str = "protein"
) -> Dict[str, Dict[str, float]]:
    matrix: Dict[str, Dict[str, float]] = {}
    canonical = set("ACGT" if alphabet.lower() == "dna" else "ACDEFGHIKLMNPQRSTVWY")
    for a in alignment:
        matrix[a.id] = {}
        for b in alignment:
            pairs = [
                (x, y) for x, y in zip(a.sequence, b.sequence) if x in canonical and y in canonical
            ]
            matrix[a.id][b.id] = sum(x == y for x, y in pairs) / len(pairs) if pairs else 0.0
    return matrix
