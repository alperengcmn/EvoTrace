"""Coordinate translations between raw sequences, alignments, and a reference."""

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass(frozen=True)
class ResidueCoordinate:
    sequence_id: str
    sequence_position: Optional[int]
    alignment_position: int
    reference_position: Optional[int]
    residue: str


class AlignmentMap:
    """1-based, gap-aware mappings for an aligned set of sequences."""

    def __init__(self, alignment, reference_id=None):
        if not alignment:
            raise ValueError("Cannot map an empty alignment.")
        lengths = {len(record.sequence) for record in alignment}
        if len(lengths) != 1:
            raise ValueError("All alignment records must have equal lengths.")
        ids = [record.id for record in alignment]
        if len(ids) != len(set(ids)):
            raise ValueError("Sequence identifiers must be unique.")
        self.length = lengths.pop()
        self.reference_id = reference_id or ids[0]
        if self.reference_id not in ids:
            raise ValueError(
                "Reference sequence {!r} is not in the alignment.".format(self.reference_id)
            )
        self._records = {r.id: r for r in alignment}
        self._seq_to_aln: Dict[str, Dict[int, int]] = {}
        self._aln_to_seq: Dict[str, Dict[int, int]] = {}
        for record in alignment:
            seq_to_aln, aln_to_seq = {}, {}
            seq_pos = 0
            for aln_pos, residue in enumerate(record.sequence, 1):
                if residue not in "-?.":
                    seq_pos += 1
                    seq_to_aln[seq_pos] = aln_pos
                    aln_to_seq[aln_pos] = seq_pos
            self._seq_to_aln[record.id] = seq_to_aln
            self._aln_to_seq[record.id] = aln_to_seq
        self.reference_sequence = self._records[self.reference_id].sequence

    def sequence_to_alignment(self, sequence_id: str, position: int) -> int:
        try:
            return self._seq_to_aln[sequence_id][position]
        except KeyError:
            raise ValueError("No residue at {} position {}.".format(sequence_id, position))

    def alignment_to_sequence(self, sequence_id: str, position: int) -> Optional[int]:
        if position < 1 or position > self.length:
            raise ValueError("Alignment position outside 1..{}".format(self.length))
        return self._aln_to_seq.get(sequence_id, {}).get(position)

    def alignment_to_reference(self, position: int) -> Optional[int]:
        return self.alignment_to_sequence(self.reference_id, position)

    def coordinate(self, sequence_id: str, alignment_position: int) -> ResidueCoordinate:
        if sequence_id not in self._records:
            raise ValueError("Unknown sequence identifier: {}".format(sequence_id))
        if alignment_position < 1 or alignment_position > self.length:
            raise ValueError("Alignment position outside 1..{}".format(self.length))
        residue = self._records[sequence_id].sequence[alignment_position - 1]
        return ResidueCoordinate(
            sequence_id,
            self.alignment_to_sequence(sequence_id, alignment_position),
            alignment_position,
            self.alignment_to_reference(alignment_position),
            residue,
        )

    def rows(self) -> List[Dict[str, object]]:
        result = []
        for seq_id, record in self._records.items():
            for aln_pos in range(1, self.length + 1):
                c = self.coordinate(seq_id, aln_pos)
                result.append(
                    {
                        "sequence_id": seq_id,
                        "sequence_position": c.sequence_position,
                        "alignment_position": aln_pos,
                        "reference_id": self.reference_id,
                        "reference_position": c.reference_position,
                        "residue": c.residue,
                    }
                )
        return result
