"""Coordinate contracts: all internal intervals are zero-based, half-open."""

from __future__ import annotations
from dataclasses import dataclass, field
import hashlib
from numbers import Integral
from typing import Optional, Tuple

Interval = Tuple[int, int]


@dataclass(frozen=True)
class Transcript:
    contig: str
    strand: str
    cds: Tuple[Interval, ...]
    source: str
    identifier: str
    phase: int = 0
    exons: Tuple[Interval, ...] = ()
    gene: str = ""
    origin: str = "prediction"

    def __post_init__(self):
        if not isinstance(self.phase, Integral) or isinstance(self.phase, bool):
            raise ValueError("CDS phase must be an integer")
        for intervals in (self.cds, self.exons):
            for iv in intervals:
                if len(iv) != 2 or any(
                    not isinstance(v, Integral) or isinstance(v, bool) for v in iv
                ):
                    raise ValueError("Coordinates must be integer interval pairs")
        object.__setattr__(self, "phase", int(self.phase))
        object.__setattr__(self, "cds", tuple((int(a), int(b)) for a, b in self.cds))
        object.__setattr__(
            self, "exons", tuple((int(a), int(b)) for a, b in self.exons)
        )
        if not self.contig or self.strand not in ("+", "-"):
            raise ValueError("A CDS transcript needs a contig and + or - strand")
        if self.phase not in (0, 1, 2) or not self.cds:
            raise ValueError("A CDS transcript needs a valid initial phase and CDS")
        previous = -1
        for a, b in self.cds:
            if a < 0 or b <= a or a < previous:
                raise ValueError(
                    f"Invalid or overlapping CDS intervals: {self.identifier}"
                )
            previous = b
        for a, b in self.exons:
            if a < 0 or b <= a:
                raise ValueError(f"Invalid exon interval: {self.identifier}")

    @property
    def key(self):
        return self.contig, self.strand, self.cds, self.phase

    @property
    def uid(self):
        return hashlib.sha256(repr(self.key).encode()).hexdigest()[:24]

    @property
    def introns(self):
        return tuple(
            (a[1], b[0]) for a, b in zip(self.cds, self.cds[1:]) if a[1] < b[0]
        )

    @property
    def length(self):
        return sum(b - a for a, b in self.cds)


@dataclass
class Candidate:
    model: Transcript
    sources: set = field(default_factory=set)
    original_ids: dict = field(default_factory=dict)
    gene_scores: list = field(default_factory=list)
    risks: list = field(default_factory=list)
    peptides: set = field(default_factory=set)
    complete: bool = False
    canonical: bool = False
    graph_path: bool = False
    graph_evidence: dict = field(default_factory=dict)

    @property
    def uid(self):
        return self.model.uid

    @property
    def q(self):
        return (
            sum(self.gene_scores) / len(self.gene_scores) if self.gene_scores else None
        )
