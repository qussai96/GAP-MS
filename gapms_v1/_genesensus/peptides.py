"""BED12 peptide evidence: deduplicate sequences and verify translation/frame."""

from __future__ import annotations
from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
import bisect, csv, json, re


@dataclass(frozen=True)
class Peptide:
    sequence: str
    contig: str
    strand: str
    blocks: tuple


def read_bed(path, sequence_table=None, keep_ambiguous=False):
    sequences = {}
    if sequence_table:
        with Path(sequence_table).open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if "name" not in row or "sequence" not in row:
                    raise ValueError(
                        "Peptide sequence table needs name and sequence columns"
                    )
                if (
                    row["name"] in sequences
                    and sequences[row["name"]] != row["sequence"]
                ):
                    raise ValueError("Conflicting peptide sequence mapping")
                sequences[row["name"]] = row["sequence"]
    placements = defaultdict(set)
    with Path(path).open() as f:
        for number, line in enumerate(f, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            row = line.rstrip().split("\t")
            if len(row) < 12:
                raise ValueError(
                    f"{path}:{number}: BED12 required to preserve split peptides"
                )
            chrom, left, right, name, _, strand = row[:6]
            a, b = int(left), int(right)
            n = int(row[9])
            sizes = [int(x) for x in row[10].rstrip(",").split(",")]
            starts = [int(x) for x in row[11].rstrip(",").split(",")]
            if (
                a < 0
                or b <= a
                or strand not in ("+", "-")
                or n < 1
                or len(sizes) != n
                or len(starts) != n
            ):
                raise ValueError(f"{path}:{number}: invalid BED12 geometry")
            blocks = tuple((a + x, a + x + size) for x, size in zip(starts, sizes))
            if (
                blocks[0][0] != a
                or blocks[-1][1] != b
                or any(x < 0 or size <= 0 for x, size in zip(starts, sizes))
                or any(x[1] > y[0] for x, y in zip(blocks, blocks[1:]))
            ):
                raise ValueError(f"{path}:{number}: invalid/overlapping peptide blocks")
            seq = sequences.get(name, name).upper().replace("I", "L")
            if not re.fullmatch("[ACDEFGHKLMNPQRSTVWY]+", seq):
                raise ValueError(
                    f"{path}:{number}: BED name must be a peptide sequence, or supply a name/sequence table"
                )
            if sum(y - x for x, y in blocks) != 3 * len(seq):
                raise ValueError(
                    f"{path}:{number}: peptide length differs from BED coding length"
                )
            placements[seq].add((chrom, strand, blocks))
    peptides = {
        seq: Peptide(seq, *next(iter(locations)))
        for seq, locations in sorted(placements.items())
        if len(locations) == 1
    }
    if keep_ambiguous:
        from .peptide_context import PeptideEvidence
        peptides = PeptideEvidence(peptides, placements)
    return peptides, dict(
        sequences=len(placements),
        unique_placement_sequences=len(peptides),
        ambiguous_sequences_excluded=0 if keep_ambiguous else sum(len(v) > 1 for v in placements.values()),
        ambiguous_sequences_retained=sum(len(v) > 1 for v in placements.values()) if keep_ambiguous else 0,
    )


def peptide_offset(model, peptide):
    if (model.contig, model.strand) != (peptide.contig, peptide.strand):
        return None
    starts = [a for a, b in model.cds]
    indices = []
    prefix = [0]
    for a, b in model.cds:
        prefix.append(prefix[-1] + b - a)
    for a, b in peptide.blocks:
        i = bisect.bisect_right(starts, a) - 1
        if i < 0 or b > model.cds[i][1]:
            return None
        indices.append(i)
    for j, (a, b) in enumerate(zip(indices, indices[1:])):
        if (
            b != a + 1
            or peptide.blocks[j][1] != model.cds[a][1]
            or peptide.blocks[j + 1][0] != model.cds[b][0]
        ):
            return None
    if model.strand == "+":
        offset = (
            prefix[indices[0]]
            + peptide.blocks[0][0]
            - model.cds[indices[0]][0]
            - model.phase
        )
    else:
        i = indices[-1]
        offset = (
            model.length
            - prefix[i]
            - (peptide.blocks[-1][1] - model.cds[i][0])
            - model.phase
        )
    return offset if offset >= 0 and offset % 3 == 0 else None


def compatible(model, peptide, protein):
    offset=peptide_offset(model,peptide)
    return (offset is not None and
        protein[offset//3:offset//3+len(peptide.sequence)].replace("I","L")==peptide.sequence)


def attach_peptides(candidates, peptides, genome):
    bins = defaultdict(set)
    for seq, p in peptides.items():
        for a, b in p.blocks:
            for k in range(a // 10000, (b - 1) // 10000 + 1):
                bins[p.contig, p.strand, k].add(seq)
    for c in candidates:
        m = c.model
        ids = set()
        for a, b in m.cds:
            for k in range(a // 10000, (b - 1) // 10000 + 1):
                ids.update(bins[m.contig, m.strand, k])
        if not ids:
            continue
        _, _, protein = genome.assess(m)
        c.peptides = {seq for seq in ids if compatible(m, peptides[seq], protein)}
