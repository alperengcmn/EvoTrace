"""Optional Neighbor Joining tree using Biopython pairwise distances."""

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def neighbor_joining(alignment):
    names = [r.id for r in alignment]
    matrix = []
    for i, a in enumerate(alignment):
        row = []
        for b in alignment[: i + 1]:
            valid = [
                (x, y) for x, y in zip(a.sequence, b.sequence) if x not in "-?" and y not in "-?"
            ]
            row.append(1 - sum(x == y for x, y in valid) / len(valid) if valid else 0.0)
        matrix.append(row)
    dm = DistanceMatrix(names, matrix)
    return DistanceTreeConstructor().nj(dm)
