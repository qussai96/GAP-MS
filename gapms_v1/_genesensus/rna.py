"""Optional RNA-seq junction adapter. No alignment or gene inference subprocesses."""

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json


def alignment_index(bam, lengths):
    """Select CSI explicitly; automatic BAI preference can truncate long contigs."""
    bam = Path(bam)
    csi = [Path(str(bam) + '.csi'), bam.with_suffix('.csi')]
    bai = [Path(str(bam) + '.bai'), bam.with_suffix('.bai')]
    for index in csi:
        if index.is_file():
            return str(index)
    if any(length > 2**29 for length in lengths.values()):
        raise ValueError(
            'RNA BAM with contigs longer than 2^29 bases requires a CSI index; '
            'BAI retrieval may silently omit alignments. Create a CSI index '
            '(samtools index -c reads.bam), or supply a validated junction JSON.'
        )
    for index in bai:
        if index.is_file():
            return str(index)
    raise ValueError('RNA BAM index is missing; indexed BAM or junction JSON required')


def transcript_strand(read, library):
    if library == "unstranded":
        return "."
    if library not in ("firststrand", "secondstrand"):
        raise ValueError("Unknown RNA library strandedness")
    same = library == "secondstrand"
    if read.is_paired and read.is_read2:
        same = not same
    reverse = read.is_reverse if same else not read.is_reverse
    return "-" if reverse else "+"


def _extract_contig(args):
    import pysam

    bam, contig, settings, index = args
    counts = defaultdict(set)
    stats = Counter()
    with pysam.AlignmentFile(bam, "rb", index_filename=index) as reader:
        for read in reader.fetch(contig):
            stats["records"] += 1
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_qcfail
            ):
                stats["filtered_flags"] += 1
                continue
            if settings["skip_duplicates"] and read.is_duplicate:
                stats["filtered_duplicates"] += 1
                continue
            nh = read.get_tag("NH") if read.has_tag("NH") else None
            if nh is not None and nh != 1:
                stats["filtered_multimappers"] += 1
                continue
            if read.mapping_quality == 255:
                if nh != 1:
                    stats["filtered_unknown_mapq"] += 1
                    continue
            elif read.mapping_quality < settings["min_mapq"]:
                stats["filtered_mapq"] += 1
                continue
            cigar = read.cigartuples or []
            pos = read.reference_start
            strand = transcript_strand(read, settings["library"])
            rg = read.get_tag("RG") if read.has_tag("RG") else ""
            for i, (op, n) in enumerate(cigar):
                if op == 3:
                    left = (
                        cigar[i - 1][1] if i > 0 and cigar[i - 1][0] in (0, 7, 8) else 0
                    )
                    right = (
                        cigar[i + 1][1]
                        if i + 1 < len(cigar) and cigar[i + 1][0] in (0, 7, 8)
                        else 0
                    )
                    if (
                        settings["min_intron"] <= n <= settings["max_intron"]
                        and min(left, right) >= settings["min_overhang"]
                    ):
                        counts[pos, pos + n, strand].add((rg, read.query_name))
                    else:
                        stats["filtered_junction_geometry"] += 1
                if op in (0, 2, 3, 7, 8):
                    pos += n
            stats["accepted_records"] += 1
    return (
        contig,
        [
            (a, b, strand, len(names))
            for (a, b, strand), names in sorted(counts.items())
        ],
        dict(stats),
    )


def extract_junctions(
    bam,
    genome_lengths,
    threads=1,
    min_mapq=20,
    min_overhang=8,
    min_intron=20,
    max_intron=500000,
    library="unstranded",
    skip_duplicates=True,
):
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("RNA evidence requires pip install genesensus[rna]") from exc
    if (
        threads < 1
        or min_mapq < 0
        or min_overhang < 1
        or not 1 <= min_intron <= max_intron
    ):
        raise ValueError("Invalid RNA extraction thresholds")
    if library not in ("unstranded", "firststrand", "secondstrand"):
        raise ValueError("Invalid RNA library mode")
    with pysam.AlignmentFile(str(bam), "rb") as reader:
        lengths = dict(zip(reader.references, reader.lengths))
        bad = {
            name: (length, genome_lengths.get(name))
            for name, length in lengths.items()
            if genome_lengths.get(name) != length
        }
        if bad:
            raise ValueError(
                "BAM/genome reference mismatch: " + repr(list(bad.items())[:10])
            )
    index = alignment_index(bam, lengths)
    with pysam.AlignmentFile(str(bam), "rb", index_filename=index) as reader:
        if not reader.has_index():
            raise ValueError("RNA BAM must have a readable index")
    settings = dict(
        min_mapq=min_mapq,
        min_overhang=min_overhang,
        min_intron=min_intron,
        max_intron=max_intron,
        library=library,
        skip_duplicates=skip_duplicates,
    )
    args = [(str(bam), c, settings, index) for c in sorted(lengths)]
    if threads == 1:
        results = list(map(_extract_contig, args))
    else:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            results = list(executor.map(_extract_contig, args))
    return dict(
        format_version=1,
        coordinate_system="0-based-half-open",
        index_kind="CSI" if index.endswith('.csi') else "BAI",
        settings=settings,
        junctions=[
            [chrom, a, b, s, n] for chrom, rows, stats in results for a, b, s, n in rows
        ],
        stats=dict(sum((Counter(stats) for c, rows, stats in results), Counter())),
        references=lengths,
    )


def write_junctions(data, path):
    Path(path).write_text(
        json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n"
    )
