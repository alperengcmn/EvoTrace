"""User-supplied regular-expression motif scanning; no motif claims are inferred."""

import re
from typing import Dict, Iterable, List
from evotrace.models import SequenceRecord
from evotrace.mapping.residues import AlignmentMap


def scan_motifs(
    records: Iterable[SequenceRecord], motifs: list, alignment=None
) -> List[Dict[str, object]]:
    hits = []
    mapper = AlignmentMap(alignment) if alignment is not None else None
    for motif in motifs or []:
        if not isinstance(motif, dict) or not motif.get("name") or not motif.get("pattern"):
            raise ValueError("Each motif requires a name and regex pattern.")
        try:
            pattern = re.compile(str(motif["pattern"]))
        except re.error as exc:
            raise ValueError("Invalid regex for motif {}: {}".format(motif["name"], exc))
        for record in records:
            for hit in pattern.finditer(record.sequence):
                start, end = hit.start() + 1, hit.end()
                row = {
                    "motif": motif["name"],
                    "sequence_id": record.id,
                    "start": start,
                    "end": end,
                    "matched": hit.group(),
                }
                if mapper is not None:
                    row["alignment_position"] = mapper.sequence_to_alignment(record.id, start)
                    row["alignment_end"] = mapper.sequence_to_alignment(record.id, end)
                hits.append(row)
    return hits
