"""Read-only indexed FASTA access; indexes never modify source FASTA directories."""

from pathlib import Path
import mmap
from Bio.Seq import Seq


class Genome:
    def __init__(self, path):
        self.path = Path(path)
        self.index = {}
        with self.path.open("rb") as f:
            name = None
            length = offset = bases = width = 0
            short = False
            while True:
                start = f.tell()
                line = f.readline()
                if not line:
                    if name is not None:
                        self.index[name] = (length, offset, bases, width)
                    break
                if line.startswith(b">"):
                    if name is not None:
                        self.index[name] = (length, offset, bases, width)
                    name = line[1:].split()[0].decode()
                    if name in self.index:
                        raise ValueError(f"Duplicate FASTA contig: {name}")
                    length = bases = width = 0
                    offset = f.tell()
                    short = False
                else:
                    if name is None:
                        raise ValueError("FASTA sequence precedes its header")
                    seq = line.rstrip(b"\r\n")
                    if not seq:
                        raise ValueError("Blank FASTA sequence lines are unsupported")
                    if short:
                        raise ValueError(f"Irregular FASTA line wrapping: {name}")
                    if not bases:
                        bases = len(seq)
                        width = len(line)
                    if len(seq) > bases:
                        raise ValueError(f"Irregular FASTA line wrapping: {name}")
                    if len(seq) < bases:
                        short = True
                    if any(c in seq for c in (b" ", b"\t")):
                        raise ValueError("Whitespace inside FASTA sequence")
                    length += len(seq)
        if not self.index:
            raise ValueError("Empty genome FASTA")
        self._file = self.path.open("rb")
        self._map = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)

    @classmethod
    def from_index(cls, path, index):
        """Reuse a validated parent index in a worker without rescanning FASTA."""
        self = cls.__new__(cls)
        self.path = Path(path)
        self.index = index
        self._file = self.path.open("rb")
        self._map = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        return self

    @property
    def lengths(self):
        return {name: data[0] for name, data in self.index.items()}

    def fetch(self, contig, start, end):
        if contig not in self.index:
            raise ValueError(f"Unknown genome contig: {contig}")
        length, offset, bases, width = self.index[contig]
        if not 0 <= start <= end <= length:
            raise ValueError(f"Genome interval out of bounds: {contig}:{start}-{end}")
        if start == end:
            return ""
        a = offset + start // bases * width + start % bases
        b = offset + (end - 1) // bases * width + (end - 1) % bases + 1
        return self._map[a:b].replace(b"\n", b"").replace(b"\r", b"").decode().upper()

    def coding_sequence(self, model):
        dna = "".join(self.fetch(model.contig, a, b) for a, b in model.cds)
        if model.strand == "-":
            dna = str(Seq(dna).reverse_complement())
        return dna

    def assess(self, model):
        dna = self.coding_sequence(model)
        coding = dna[model.phase :]
        aa = str(Seq(coding[: len(coding) // 3 * 3]).translate())
        complete = (
            model.phase == 0
            and len(dna) % 3 == 0
            and dna.startswith("ATG")
            and aa.endswith("*")
            and "*" not in aa[:-1]
            and "X" not in aa
        )
        canonical = True
        for a, b in model.introns:
            if b - a < 4:
                canonical = False
                continue
            left, right = self.fetch(model.contig, a, a + 2), self.fetch(
                model.contig, b - 2, b
            )
            if model.strand == "-":
                left, right = str(Seq(right).reverse_complement()), str(
                    Seq(left).reverse_complement()
                )
            canonical &= (left, right) in {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
        return complete, bool(canonical), aa

    def close(self):
        self._map.close()
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
