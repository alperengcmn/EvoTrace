"""Typed evidence records; measurements remain separate from interpretation."""

from dataclasses import asdict, dataclass
from typing import Optional


@dataclass(frozen=True)
class Evidence:
    source: str
    position: Optional[int]
    evidence_type: str
    strength: str
    description: str
    p_value: Optional[float] = None
    adjusted_p_value: Optional[float] = None
    method: Optional[str] = None
    evidence_family: Optional[str] = None
    result_level: Optional[str] = None

    def to_dict(self):
        return asdict(self)
