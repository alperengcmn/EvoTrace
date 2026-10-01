"""Core data structures shared by analysis stages."""

from dataclasses import dataclass
from typing import Dict, List


@dataclass(frozen=True)
class SequenceRecord:
    id: str
    sequence: str
    description: str = ""


@dataclass
class AnalysisResult:
    input_type: str
    alphabet: str
    records: List[SequenceRecord]
    alignment: List[SequenceRecord]
    alignment_stats: Dict[str, object]
    columns: List[Dict[str, object]]
    identity: Dict[str, Dict[str, float]]
    scores: List[Dict[str, object]]
    warnings: List[str]
