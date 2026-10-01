"""Validated back-translation of a protein alignment into a codon alignment."""

from typing import Any, Dict, List, Tuple
from Bio.Seq import Seq
from Bio.Data import CodonTable
from evotrace.models import SequenceRecord


class CodonAlignmentError(ValueError):
    pass


def backtranslate(
    protein_alignment: List[SequenceRecord], coding_sequences: List[SequenceRecord], table: int = 1
) -> Tuple[List[SequenceRecord], Dict[str, object]]:
    """Project each protein gap pattern onto its matching CDS codons.

    Protein residues must exactly match the CDS translation (terminal stop may be
    present and is removed); internal stops, ambiguous codons, frameshifts and
    mismatched identifiers/residues are rejected instead of guessed.
    """
    cds_by_id = {r.id: r.sequence.upper().replace("U", "T") for r in coding_sequences}
    if not protein_alignment or len({len(r.sequence) for r in protein_alignment}) != 1:
        raise CodonAlignmentError("Protein alignment records must be non-empty and equal length.")
    if len(cds_by_id) != len(coding_sequences):
        raise CodonAlignmentError("CDS identifiers are not unique.")
    protein_ids = {record.id for record in protein_alignment}
    extra_ids = set(cds_by_id) - protein_ids
    if extra_ids:
        raise CodonAlignmentError(
            "CDS identifiers have no matching protein: {}.".format(", ".join(sorted(extra_ids)))
        )
    out: List[SequenceRecord] = []
    report: Dict[str, Any] = {"genetic_code": table, "sequences": {}}
    for protein in protein_alignment:
        if protein.id not in cds_by_id:
            raise CodonAlignmentError(
                "No coding sequence supplied for protein {!r}.".format(protein.id)
            )
        dna = cds_by_id[protein.id]
        if len(dna) % 3:
            raise CodonAlignmentError(
                "{}: CDS length {} is not divisible by 3.".format(protein.id, len(dna))
            )
        if set(dna) - set("ACGT"):
            raise CodonAlignmentError(
                "{}: CDS contains ambiguous/non-ACGT bases.".format(protein.id)
            )
        codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]
        translation = str(Seq(dna).translate(table=table, to_stop=False))
        if "*" in translation[:-1]:
            raise CodonAlignmentError("{}: CDS contains an internal stop codon.".format(protein.id))
        if translation.endswith("*"):
            codons, translation = codons[:-1], translation[:-1]
        ungapped = protein.sequence.replace("-", "").replace(".", "")
        if ungapped != translation:
            mismatch = next(
                (i for i, (a, b) in enumerate(zip(ungapped, translation), 1) if a != b), None
            )
            raise CodonAlignmentError(
                "{}: protein/CDS translation mismatch{} (protein {}, translation {}).".format(
                    protein.id,
                    " at residue {}".format(mismatch) if mismatch else " in length",
                    len(ungapped),
                    len(translation),
                )
            )
        cursor, projected = 0, []
        for residue in protein.sequence:
            if residue in "-.":
                projected.append("---")
            else:
                if cursor >= len(codons):
                    raise CodonAlignmentError(
                        "{}: alignment consumes more residues than CDS.".format(protein.id)
                    )
                projected.append(codons[cursor])
                cursor += 1
        out.append(SequenceRecord(protein.id, "".join(projected), protein.description))
        report["sequences"][protein.id] = {
            "codons_used": cursor,
            "terminal_stop_removed": len(dna) // 3 > len(codons),
        }
    return out, report


def validate_codon_alignment(codon_alignment, protein_alignment=None, table=1):
    """Validate codon frame, symbols, stop codons and optional protein-column correspondence."""
    if not codon_alignment:
        raise CodonAlignmentError("Codon alignment is empty.")
    lengths = {len(record.sequence) for record in codon_alignment}
    if len(lengths) != 1 or next(iter(lengths)) == 0:
        raise CodonAlignmentError("Codon alignment records must have equal non-zero lengths.")
    length = lengths.pop()
    if length % 3:
        raise CodonAlignmentError("Codon alignment length must be divisible by three.")
    ids = [record.id for record in codon_alignment]
    if len(ids) != len(set(ids)):
        raise CodonAlignmentError("Codon alignment identifiers must be unique.")
    protein_by_id = None
    if protein_alignment is not None:
        protein_by_id = {record.id: record.sequence.upper() for record in protein_alignment}
        if set(protein_by_id) != set(ids):
            raise CodonAlignmentError("Codon and protein alignment identifiers must match exactly.")
        if any(len(sequence) * 3 != length for sequence in protein_by_id.values()):
            raise CodonAlignmentError("Protein and codon alignment columns do not match.")
    table_info = CodonTable.unambiguous_dna_by_id.get(table)
    if table_info is None:
        raise CodonAlignmentError("Unknown unambiguous DNA genetic code table: {}.".format(table))
    gap_codons = 0
    for record in codon_alignment:
        dna = record.sequence.upper().replace("U", "T")
        protein = protein_by_id.get(record.id) if protein_by_id is not None else None
        for offset in range(0, length, 3):
            codon = dna[offset : offset + 3]
            aa = protein[offset // 3] if protein is not None else None
            if codon in ("---", "..."):
                gap_codons += 1
                if aa is not None and aa not in "-.?":
                    raise CodonAlignmentError(
                        "{}: codon gap does not match protein residue at alignment position {}.".format(
                            record.id, offset // 3 + 1
                        )
                    )
                continue
            if "-" in codon or "." in codon:
                raise CodonAlignmentError(
                    "{}: partial codon gap at alignment bases {}-{}.".format(
                        record.id, offset + 1, offset + 3
                    )
                )
            if set(codon) - set("ACGT"):
                raise CodonAlignmentError(
                    "{}: ambiguous/non-ACGT codon at alignment position {}.".format(
                        record.id, offset // 3 + 1
                    )
                )
            translated = str(Seq(codon).translate(table=table, to_stop=False))
            if translated == "*":
                raise CodonAlignmentError(
                    "{}: stop codon at alignment position {}.".format(record.id, offset // 3 + 1)
                )
            if aa is not None and aa not in "X?" and translated != aa:
                raise CodonAlignmentError(
                    "{}: codon translates to {} but protein alignment has {} at position {}.".format(
                        record.id, translated, aa, offset // 3 + 1
                    )
                )
    return {
        "status": "PASS",
        "sequence_count": len(codon_alignment),
        "codon_columns": length // 3,
        "gap_containing_codons": gap_codons,
        "genetic_code": table,
        "protein_correspondence_checked": protein_by_id is not None,
    }
