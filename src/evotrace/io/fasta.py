"""Strict, dependency-light FASTA parsing and writing."""

from pathlib import Path
from typing import Iterable, List, Set, Union
from evotrace.models import SequenceRecord


class FastaError(ValueError):
    """Raised for malformed FASTA input."""


def read_fasta(path: Union[str, Path]) -> List[SequenceRecord]:
    records: List[SequenceRecord] = []
    seen: Set[str] = set()
    header = None
    chunks: List[str] = []
    try:
        stream = open(path, "r", encoding="utf-8-sig")
    except OSError as exc:
        raise FastaError("Cannot read FASTA '{}': {}".format(path, exc))
    with stream:
        for line_no, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    _append(records, seen, header, chunks)
                header = line[1:].strip()
                if not header:
                    raise FastaError("Empty FASTA identifier at line {}".format(line_no))
                chunks = []
            else:
                if header is None:
                    raise FastaError(
                        "Sequence data before first FASTA header at line {}".format(line_no)
                    )
                chunks.append("".join(line.split()))
    if header is not None:
        _append(records, seen, header, chunks)
    if not records:
        raise FastaError("No FASTA records found in '{}'".format(path))
    return records


def _append(records, seen, header, chunks):
    seq = "".join(chunks).upper().replace(".", "-")
    ident = header.split()[0]
    if ident in seen:
        raise FastaError("Duplicate sequence identifier: '{}'".format(ident))
    if not seq:
        raise FastaError("Sequence '{}' is empty".format(ident))
    seen.add(ident)
    records.append(SequenceRecord(ident, seq, header[len(ident) :].strip()))


def write_fasta(records: Iterable[SequenceRecord], path: Union[str, Path]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for rec in records:
            out.write(">{}{}\n".format(rec.id, (" " + rec.description) if rec.description else ""))
            for start in range(0, len(rec.sequence), 80):
                out.write(rec.sequence[start : start + 80] + "\n")
