"""Validated peptide inputs and read-only indexed assembly access."""
import csv
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path

import pysam

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
HEADERS = {"sequence", "peptide", "peptides", "peptide_sequence", "peptide sequence", "stripped.sequence", "stripped_sequence", "modified sequence"}


def normalize_peptide(text):
    text = text.strip()
    if re.match(r"^[A-Za-z-]\.\S+\.[A-Za-z-]$", text):
        text = text[2:-2]
    text = re.sub(r"\[[^\]]*\]|\([^)]*\)", "", text).strip("_").upper()
    if not text or set(text) - AMINO_ACIDS:
        raise ValueError(f"invalid peptide sequence {text!r}; expected the 20 standard amino acids (modifications in [] or () are accepted)")
    return text


def read_peptides(path, column=None, min_length=1):
    """Read plain lists, FASTA, CSV, or TSV, preserving first occurrence order."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError("peptide input is empty")
    values = []
    if lines[0].startswith(">"):
        sequence = []
        for line in lines:
            if line.startswith(">"):
                if sequence:
                    values.append("".join(sequence))
                    sequence = []
            else:
                sequence.append(line.strip())
        if sequence:
            values.append("".join(sequence))
    else:
        delimiter = "\t" if "\t" in lines[0] else "," if "," in lines[0] else None
        rows = list(csv.reader(lines, delimiter=delimiter)) if delimiter else [[line.strip()] for line in lines]
        header = [c.strip().lower() for c in rows[0]]
        if column is not None:
            if column.lower() not in header:
                raise ValueError(f"peptide column {column!r} is absent from {rows[0]!r}")
            index = header.index(column.lower())
            rows = rows[1:]
        elif any(c in HEADERS and (delimiter or not (rows[0][i].isupper() and set(rows[0][i]) <= AMINO_ACIDS)) for i, c in enumerate(header)):
            # An uppercase plain-list PEPTIDE/PEPTIDES is a real amino-acid
            # sequence, while lower/mixed-case labels denote a header.
            index = next(i for i, c in enumerate(header) if c in HEADERS)
            rows = rows[1:]
        elif delimiter:
            raise ValueError("multi-column input requires a recognised peptide header or --peptide-column")
        else:
            index = 0
        for number, row in enumerate(rows, 1):
            if len(row) <= index:
                raise ValueError(f"missing peptide field in input row {number}")
            values.append(row[index])
    unique = {}
    skipped = 0
    for number, value in enumerate(values, 1):
        try:
            peptide = normalize_peptide(value)
        except ValueError as error:
            raise ValueError(f"peptide record {number}: {error}") from error
        if len(peptide) < min_length:
            skipped += 1
            continue
        unique.setdefault(peptide, None)
    if not unique:
        raise ValueError("no peptides remain after length filtering")
    return list(unique), {"input_records": len(values), "unique_peptides": len(unique), "below_min_length_records": skipped}


def prepare_reference(genome, cache):
    """Index in the output cache, never alongside or inside the user's assembly."""
    genome, cache = Path(genome).resolve(), Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    stat = genome.stat()
    signature = {"path": str(genome), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    token = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:20]
    if genome.suffix.lower() in {".gz", ".bgz", ".bgzf"}:
        expanded = cache / f"{token}.fa"
        if not expanded.exists():
            temporary = expanded.with_suffix(".partial")
            with gzip.open(genome, "rb") as source, temporary.open("wb") as destination:
                shutil.copyfileobj(source, destination, 8 * 1024 * 1024)
            temporary.replace(expanded)
        genome = expanded
    index = cache / f"{token}.fai"
    if not index.exists():
        seen = set()
        with open(genome, "rb") as source:
            for line in source:
                if line.startswith(b">"):
                    fields = line[1:].split()
                    if not fields:
                        raise ValueError("assembly contains an empty FASTA identifier")
                    name = fields[0]
                    if name in seen:
                        raise ValueError(f"assembly has duplicate FASTA sequence identifier: {name.decode()}")
                    seen.add(name)
        partial_index = index.with_suffix(".fai.partial")
        pysam.faidx("-o", str(partial_index), str(genome))
        partial_index.replace(index)
    # Reject duplicate FASTA IDs: HTSlib otherwise silently selects one copy.
    names = [line.split("\t", 1)[0] for line in index.read_text().splitlines()]
    if len(names) != len(set(names)):
        raise ValueError("assembly has duplicate FASTA sequence identifiers")
    with pysam.FastaFile(str(genome), filepath_index=str(index)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
    if not lengths:
        raise ValueError("assembly contains no sequences")
    return str(genome), str(index), lengths
