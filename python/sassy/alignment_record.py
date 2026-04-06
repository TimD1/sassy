from dataclasses import dataclass
from typing import Literal

from fgpyo.util.metric import Metric


@dataclass(frozen=True)
class AlignmentRecord(Metric["AlignmentRecord"]):
    """One alignment returned by search_all_alignments(), serializable to TSV.

    Only complete pattern matches are recorded (pattern_start=0,
    pattern_end=len(pattern)), so those fields are omitted.
    """

    ref_name: str
    text_start: int
    text_end: int
    cost: int
    cigar: str
    strand: Literal["Fwd", "Rev"]
