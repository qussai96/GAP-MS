"""GFF3/GTF readers and deterministic, phase-consistent GFF3 output."""

from pathlib import Path
from collections import defaultdict
from urllib.parse import unquote, quote
import re
from .models import Transcript
from .intervals import merge


def attributes(text):
    if re.match(r"\s*[^\s;=]+=", text):
        return {
            k: unquote(v)
            for k, v in (
                part.strip().split("=", 1) for part in text.split(";") if "=" in part
            )
        }
    return dict(re.findall(r'([^\s;]+)\s+"([^"\\]*(?:\\.[^"\\]*)*)"\s*;?', text))


def read_annotation(
    path, source="prediction", allow_orphan_cds=False, origin="prediction"
):
    parents = {}
    cds = defaultdict(list)
    exons = defaultdict(list)
    with Path(path).open() as f:
        for number, line in enumerate(f, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) != 9:
                raise ValueError(f"{path}:{number}: expected nine columns")
            chrom, _, kind, left, right, _, strand, phase, text = row
            if kind not in ("mRNA", "transcript", "CDS", "exon"):
                continue
            if strand not in ("+", "-"):
                raise ValueError(
                    f"{path}:{number}: CDS/transcript strand must be + or -"
                )
            a, b = int(left) - 1, int(right)
            if not 0 <= a < b:
                raise ValueError(f"{path}:{number}: invalid coordinates")
            attr = attributes(text)
            if kind in ("mRNA", "transcript"):
                tid = attr.get("ID", attr.get("transcript_id"))
                if not tid:
                    raise ValueError(f"{path}:{number}: missing transcript identifier")
                key = (chrom, strand, tid)
                gene = attr.get("Parent", attr.get("gene_id", tid))
                if key in parents and parents[key] != gene:
                    raise ValueError(f"{path}:{number}: conflicting transcript parents")
                parents[key] = gene
            else:
                ids = attr.get("Parent", attr.get("transcript_id"))
                if not ids:
                    raise ValueError(f"{path}:{number}: missing parent")
                for tid in ids.split(","):
                    key = (chrom, strand, tid)
                    if kind == "CDS":
                        if phase not in ("0", "1", "2"):
                            raise ValueError(
                                f"{path}:{number}: missing/invalid CDS phase"
                            )
                        cds[key].append((a, b, int(phase)))
                    else:
                        exons[key].append((a, b))
    models = []
    for (chrom, strand, tid), segments in sorted(cds.items()):
        key = (chrom, strand, tid)
        if key not in parents and not allow_orphan_cds:
            raise ValueError(f"{path}: CDS has no transcript parent: {key}")
        ss = sorted(set(segments))
        first = ss[0 if strand == "+" else -1][2]
        coords = tuple((a, b) for a, b, p in ss)
        models.append(
            Transcript(
                chrom,
                strand,
                coords,
                source,
                tid,
                first,
                tuple(sorted(set(exons[key]))) or coords,
                parents.get(key, tid),
                origin,
            )
        )
    if not models:
        raise ValueError(f"No coding transcripts found: {path}")
    return models


def write_gff3(models, path):
    models = sorted(models, key=lambda m: (m.contig, m.cds[0][0], m.strand, m.uid))
    with Path(path).open("w") as f:
        f.write("##gff-version 3\n")
        for m in models:
            g = "GS.g" + m.uid
            t = g + ".t1"
            # CDS bases are exonic, including explicitly normalized stop codons.
            ee = merge((m.exons or m.cds) + m.cds)
            lo = min(a for a, b in ee)
            hi = max(b for a, b in ee)
            attr = (
                f'ID={t};Parent={g};origin={quote(m.origin,safe="")};candidate={m.uid}'
            )
            f.write(
                f"{m.contig}\tGeneSensus\tgene\t{lo+1}\t{hi}\t.\t{m.strand}\t.\tID={g}\n{m.contig}\tGeneSensus\tmRNA\t{lo+1}\t{hi}\t.\t{m.strand}\t.\t{attr}\n"
            )
            for i, (a, b) in enumerate(ee, 1):
                f.write(
                    f"{m.contig}\tGeneSensus\texon\t{a+1}\t{b}\t.\t{m.strand}\t.\tID={t}.e{i};Parent={t}\n"
                )
            phases = {}
            used = -m.phase
            for a, b in m.cds if m.strand == "+" else reversed(m.cds):
                phases[a, b] = (-used) % 3
                used += b - a
            for i, (a, b) in enumerate(m.cds, 1):
                f.write(
                    f"{m.contig}\tGeneSensus\tCDS\t{a+1}\t{b}\t.\t{m.strand}\t{phases[a,b]}\tID={t}.c{i};Parent={t}\n"
                )
