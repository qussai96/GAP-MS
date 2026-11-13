"""Reusable exact, I/L-equivalent peptide-to-protein search.

Build one C-backed peptide automaton and scan independent protein FASTAs.
Return all overlapping occurrences with one-based protein coordinates.
"""
from collections import defaultdict, OrderedDict
import json
from pathlib import Path
import tempfile
import time

import ahocorasick
from Bio import SeqIO

from peplocus.io import normalize_peptide, read_peptides


class ProteinMatcher:
    """Share a peptide automaton between prediction and BAM-derived proteins."""

    def __init__(self, peptides):
        started = time.monotonic()
        self.peptides = tuple(sorted({normalize_peptide(p) for p in peptides}))
        if not self.peptides:
            raise ValueError("no peptides to search")
        grouped = defaultdict(list)
        for index, peptide in enumerate(self.peptides):
            grouped[peptide.replace("I", "L")].append(index)
        automaton = ahocorasick.Automaton()
        for pattern, indices in grouped.items():
            automaton.add_word(pattern, (len(pattern), tuple(indices)))
        automaton.make_automaton()
        self.automaton = automaton
        self.build_seconds = time.monotonic() - started
        self.searches = 0
        # Reuse sequence searches across identical source/consensus proteins.
        # Bound both sequence storage and highly repetitive hit lists.
        self._sequence_hits = OrderedDict()
        self._cached_residues = 0
        self._cached_hits = 0
        self._max_cached_residues = 32_000_000
        self._max_cached_hits = 1_000_000

    @classmethod
    def from_file(cls, path, column=None):
        """Read validated plain, FASTA, CSV or TSV peptide input."""
        peptides, _ = read_peptides(path, column)
        return cls(peptides)

    def map(self, proteins, output):
        """Write a new mapping TSV and timing/count summary for one FASTA."""
        output = Path(output)
        if output.exists():
            raise ValueError(f"output already exists: {output}")
        started = time.monotonic()
        hits = defaultdict(list)
        seen = set()
        residues = 0
        scanned = reused = 0
        for record in SeqIO.parse(proteins, "fasta"):
            if record.id in seen:
                raise ValueError(f"duplicate protein FASTA identifier: {record.id}")
            seen.add(record.id)
            # Preserve stops/unknown residues as barriers and retain coordinates.
            sequence = str(record.seq).upper().replace("I", "L")
            residues += len(sequence)
            matches = self._sequence_hits.get(sequence)
            if matches is None:
                scanned += 1
                matches = tuple(self.automaton.iter(sequence))
                if len(sequence) <= self._max_cached_residues and len(matches) <= self._max_cached_hits:
                    while self._sequence_hits and (self._cached_residues + len(sequence) > self._max_cached_residues or self._cached_hits + len(matches) > self._max_cached_hits):
                        old_sequence, old_matches = self._sequence_hits.popitem(last=False)
                        self._cached_residues -= len(old_sequence)
                        self._cached_hits -= len(old_matches)
                    self._sequence_hits[sequence] = matches
                    self._cached_residues += len(sequence)
                    self._cached_hits += len(matches)
            else:
                reused += 1
                self._sequence_hits.move_to_end(sequence)
            for end, (length, indices) in matches:
                for index in indices:
                    hits[index].append((record.id, end - length + 2))
        if not seen:
            raise ValueError("protein FASTA is empty")
        scan_seconds = time.monotonic() - started

        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", dir=output.parent, prefix=output.name + ".",
                suffix=".partial", delete=False,
            ) as handle:
                temporary = Path(handle.name)
                handle.write("peptide\tprotein\tlocation\n")
                for index, peptide in enumerate(self.peptides):
                    placements = hits.get(index)
                    if placements:
                        for protein, start in sorted(placements):
                            handle.write(f"{peptide}\t{protein}\t{start}\n")
                    else:
                        handle.write(f"{peptide}\tUNMAPPED\t0\n")
            temporary.replace(output)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

        self.searches += 1
        stats = dict(
            backend="aho-corasick", il_equivalent=True,
            peptides=len(self.peptides), proteins=len(seen),
            protein_residues=residues, mapped_peptides=len(hits),
            protein_sequences_scanned=scanned, protein_sequences_reused=reused,
            cached_protein_sequences=len(self._sequence_hits),
            cached_protein_residues=self._cached_residues, cached_pattern_hits=self._cached_hits,
            unmapped_peptides=len(self.peptides) - len(hits),
            mapping_rows=sum(map(len, hits.values())), automaton_builds=1,
            automaton_searches=self.searches, build_seconds=self.build_seconds,
            scan_seconds=scan_seconds,
            search_and_output_seconds=time.monotonic() - started,
        )
        output.with_suffix(".summary.json").write_text(json.dumps(stats, indent=2) + "\n")
        return output
