"""Vectorized, ambiguity-aware six-frame translation and coordinate projection."""
from functools import lru_cache
from itertools import product

import numpy as np
from Bio.Seq import translate as bio_translate

ALPHABET = "ACGTRYSWKMBDHVNN"
COMPLEMENT = bytes.maketrans(b"ACGTRYSWKMBDHVNUacgtryswkmbdhvnu", b"TGCAYRSWMKVHDBNAtgcayrswmkvhdbna")
ENCODE = np.full(256, 15, dtype=np.uint16)
for _i, _base in enumerate(ALPHABET[:15]):
    ENCODE[ord(_base)] = ENCODE[ord(_base.lower())] = _i
ENCODE[ord("U")] = ENCODE[ord("u")] = 3


@lru_cache(maxsize=16)
def translation_lut(table=1):
    codons = "".join("".join(c) for c in product(ALPHABET, repeat=3))
    protein = bio_translate(codons, table=table)
    return np.frombuffer(protein.encode("ascii"), dtype=np.uint8)


def translate(dna, table=1):
    """Translate complete codons; deterministic IUPAC codons are resolved."""
    if isinstance(dna, str):
        dna = dna.encode("ascii")
    size = len(dna) // 3 * 3
    if not size:
        return ""
    codes = ENCODE[np.frombuffer(dna[:size], dtype=np.uint8)]
    indexes = codes[0::3] * 256 + codes[1::3] * 16 + codes[2::3]
    return translation_lut(table)[indexes].tobytes().decode("ascii")


def reverse_complement(dna):
    if isinstance(dna, str):
        return dna.encode("ascii").translate(COMPLEMENT)[::-1].decode("ascii")
    return dna.translate(COMPLEMENT)[::-1]


def canonical(peptide, il_equivalent=True):
    return peptide.replace("L", "I") if il_equivalent else peptide


def project(start, end, segments, strand):
    """Project transcript-relative [start,end) into ascending genomic BED blocks.

    Segments are ordered in translation direction. Codons may cross boundaries.
    """
    blocks = []
    offset = 0
    for left, right in segments:
        length = right - left
        a, b = max(start, offset), min(end, offset + length)
        if a < b:
            if strand == "+":
                blocks.append((left + a - offset, left + b - offset))
            else:
                blocks.append((right - (b - offset), right - (a - offset)))
        offset += length
    return tuple(sorted(blocks))


def join_blocks(reference, chrom, strand, blocks):
    dna = "".join(reference.fetch(chrom, a, b) for a, b in sorted(blocks))
    return reverse_complement(dna) if strand == "-" else dna
