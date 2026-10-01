"""FASTA alphabet detection and validation."""

from typing import List, Tuple
from evotrace.models import SequenceRecord

DNA = set("ACGTUNRYKMSWBDHV-")
PROTEIN = set("ABCDEFGHIKLMNPQRSTVWXYZ*UO-")


class ValidationError(ValueError):
    pass


def validate_records(
    records: List[SequenceRecord], requested: str = "auto", max_ambiguity: float = 0.5
) -> Tuple[str, List[str]]:
    if not records:
        raise ValidationError("Input has no sequences.")
    requested = requested.lower()
    if requested not in ("auto", "dna", "protein"):
        raise ValidationError("Alphabet must be one of: auto, dna, protein.")
    all_chars = set("".join(r.sequence.upper() for r in records))
    dna_ok = all_chars <= DNA
    alphabet = requested if requested != "auto" else ("dna" if dna_ok else "protein")
    allowed = DNA if alphabet == "dna" else PROTEIN
    errors = []
    warnings = []
    for rec in records:
        invalid = sorted(set(rec.sequence) - allowed)
        if invalid:
            errors.append("{}: unsupported characters {}".format(rec.id, "".join(invalid)))
        if len(rec.sequence) < 3:
            warnings.append(
                "{}: very short sequence ({} residues)".format(rec.id, len(rec.sequence))
            )
        informative = [c for c in rec.sequence if c not in "-?"]
        ambiguity_chars = set("NRYKMSWBDHVX") if alphabet == "dna" else set("XBZJ?")
        fraction = sum(c in ambiguity_chars for c in informative) / max(1, len(informative))
        if fraction > max_ambiguity:
            warnings.append("{}: high ambiguity fraction {:.1%}".format(rec.id, fraction))
    if errors:
        raise ValidationError("Invalid {} FASTA:\n- {}".format(alphabet, "\n- ".join(errors)))
    if len(records) < 2:
        warnings.append("Only one sequence: comparative metrics and phylogeny are limited.")
    return alphabet, warnings
