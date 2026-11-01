"""PepLocus command-line interface."""
import argparse
import csv
import json
import logging
from pathlib import Path
import sys
import tempfile

import pysam

from . import __version__
from .io import prepare_reference
from .junctions import extract_bam, write_junctions
from .annotation import read_annotation
from .pipeline import map_peptides
from .rescue import one_edit
from .sequence import canonical, join_blocks, translate


def positive(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def nonnegative(text):
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("must be zero or a positive integer")
    return value


def junction_options(parser):
    parser.add_argument("--min-mapq", type=nonnegative, default=10, help="minimum BAM mapping quality (default: 10)")
    parser.add_argument("--min-read-anchor", type=positive, default=8, help="aligned bases on each side of a BAM intron (default: 8)")
    parser.add_argument("--min-junction-support", type=positive, default=2, help="minimum primary, nonduplicate alignment observations per junction (default: 2)")
    parser.add_argument("--min-intron", type=positive, default=20)
    parser.add_argument("--max-intron", type=positive, default=500000)


def run_junctions(args):
    out = Path(args.out)
    if out.exists():
        raise ValueError(f"output already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="peplocus-junctions-") as cache:
        genome, index, lengths = prepare_reference(args.genome, cache)
        with pysam.FastaFile(genome, filepath_index=index) as fasta:
            junctions, stats = [], {}
            if args.bam:
                junctions, bam_stats = extract_bam(args.bam, fasta, lengths, args.min_mapq, args.min_read_anchor, args.min_junction_support, args.min_intron, args.max_intron)
                stats["bam"] = bam_stats
            if args.gtf:
                annotation_junctions, annotation_stats = read_annotation(args.gtf, fasta, lengths)
                junctions.extend(annotation_junctions)
                stats["annotation"] = annotation_stats
        write_junctions(out, junctions)
    print(json.dumps(stats, indent=2))


def run_verify(args):
    count = 0
    with tempfile.TemporaryDirectory(prefix="peplocus-verify-") as cache:
        genome, index, lengths = prepare_reference(args.genome, cache)
        with pysam.FastaFile(genome, filepath_index=index) as fasta, open(args.mappings) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"peptide", "chrom", "strand", "blocks", "matched_sequence"} <= set(reader.fieldnames or []):
                raise ValueError("expected a PepLocus mappings.tsv file")
            for row in reader:
                blocks = json.loads(row["blocks"])
                if row["strand"] not in {"+", "-"} or row["chrom"] not in lengths or not blocks:
                    raise ValueError(f"invalid mapping: {row}")
                previous = -1
                for a,b in blocks:
                    if not 0 <= a < b <= lengths[row["chrom"]] or a < previous:
                        raise ValueError(f"invalid genomic blocks: {blocks}")
                    previous = b
                observed = translate(join_blocks(fasta, row["chrom"], row["strand"], blocks), args.genetic_code)
                if row.get("match_type") == "one_edit":
                    edit = one_edit(row["peptide"], observed, args.il_equivalent)
                    valid = (observed == row["matched_sequence"] and edit is not None
                             and row.get("edit_distance") == "1"
                             and all(str(edit[key]) == row.get(key) for key in
                                     ("edit_type", "edit_query_position", "query_aa", "reference_aa", "alignment_cigar")))
                else:
                    valid = observed == row["matched_sequence"] and canonical(observed,args.il_equivalent) == canonical(row["peptide"],args.il_equivalent)
                if not valid:
                    raise ValueError(f"mapping does not reproduce peptide/alignment: {row['peptide']} at {row['chrom']}:{blocks}")
                count += 1
    print(json.dumps({"verified_loci": count, "valid": True}))


def parser():
    root = argparse.ArgumentParser(description="PepLocus: map peptides directly to eukaryotic assemblies and supported splice paths")
    root.add_argument("--version", action="version", version=f"PepLocus {__version__}")
    commands = root.add_subparsers(dest="command", required=True)
    mapping = commands.add_parser("map", help="map a peptide list to an assembly")
    mapping.add_argument("--genome", required=True, help="assembly FASTA (plain or gzip; inputs remain unchanged)")
    mapping.add_argument("--peptides", required=True, help="peptide list, FASTA, CSV or TSV")
    mapping.add_argument("--peptide-column", help="column name in CSV/TSV input")
    mapping.add_argument("--out", required=True, help="output directory without existing PepLocus results")
    mapping.add_argument("--gtf", "--gff", dest="gtf", help="optional transcript annotation (GTF or GFF3-style Helixer); supplies exon junctions without BAM")
    mapping.add_argument("--bam", help="optional RNA-seq BAM; no index required")
    mapping.add_argument("--junctions", help="optional reusable PepLocus junction TSV")
    mapping.add_argument("--threads", type=positive, default=1, help="worker processes for contiguous search (default: 1)")
    mapping.add_argument("--chunk-size", type=positive, default=2000000, help="genomic bases per worker task (default: 2000000)")
    mapping.add_argument("--min-length", type=positive, default=1, help="minimum peptide length; short matches can be highly ambiguous")
    mapping.add_argument("--il-equivalent", action="store_true", default=True, help="compatibility option: I/L equivalence is always enabled")
    mapping.add_argument("--genetic-code", type=positive, default=1, help="NCBI translation table; applies to every contig (default: 1)")
    mapping.add_argument("--rescue-one-edit", action="store_true", help="after exact matching, rescue still-unmapped peptides with one amino-acid substitution/insertion/deletion; separate outputs")
    mapping.add_argument("--rescue-min-length", type=positive, default=12, help="minimum length for one-edit rescue (default: 12); exact-stage --min-length is independent")
    mapping.add_argument("--igv-display-bed", action="store_true", help="also write a display-only BED with reversed negative-strand labels and GAPMS-style colors")
    mapping.add_argument("--no-verify", action="store_true", help="skip final extraction/translation verification (enabled by default)")
    junction_options(mapping)
    mapping.add_argument("--max-junctions", type=positive, default=3, help="maximum introns crossed by a supported peptide path (default: 3)")
    mapping.add_argument("--max-paths-per-junction", type=nonnegative, default=1000, help="splice graph path cap per first intron; 0 means unlimited")
    mapping.add_argument("--discover-junctions", action="store_true", help="enable bounded, canonical single-intron discovery without BAM support")
    mapping.add_argument("--discover-all", action="store_true", help="run discovery for mapped peptides too (requires --discover-junctions)")
    mapping.add_argument("--novel-anchor", type=positive, default=4, help="amino acids required at each terminal discovery anchor (default: 4)")
    mapping.add_argument("--novel-cache-mb", type=nonnegative, default=256, help="memory budget for caching small contigs during discovery; 0 disables it (default: 256 MiB)")
    mapping.add_argument("--max-seed-hits", type=nonnegative, default=2000, help="discovery seed-hit cap per contig/strand; 0 means unlimited")
    mapping.add_argument("--max-novel-pairs", type=nonnegative, default=10000, help="discovery candidate-junction cap per peptide/contig/strand; 0 means unlimited")
    mapping.set_defaults(func=map_peptides)
    junction = commands.add_parser("junctions", help="extract reusable junctions from RNA-seq BAM and/or annotation")
    junction.add_argument("--genome", required=True)
    junction.add_argument("--bam", help="RNA-seq BAM")
    junction.add_argument("--gtf", "--gff", dest="gtf", help="GTF or Helixer GFF3-style annotation")
    junction.add_argument("--out", required=True, help="new junction TSV path")
    junction_options(junction)
    junction.set_defaults(func=run_junctions)
    verify = commands.add_parser("verify", help="independently extract and translate every reported locus")
    verify.add_argument("--genome", required=True)
    verify.add_argument("--mappings", required=True)
    verify.add_argument("--il-equivalent", action="store_true", default=True, help="compatibility option: I/L equivalence is always enabled")
    verify.add_argument("--genetic-code", type=positive, default=1)
    verify.set_defaults(func=run_verify)
    return root


def main(argv=None):
    root = parser()
    args = root.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.command == "junctions" and not (args.bam or args.gtf):
        root.error("junctions requires --bam and/or --gtf")
    if hasattr(args,"min_intron") and args.min_intron > args.max_intron:
        root.error("--min-intron cannot exceed --max-intron")
    if getattr(args,"discover_all",False) and not args.discover_junctions:
        root.error("--discover-all requires --discover-junctions")
    if getattr(args,"rescue_min_length",2) < 2:
        root.error("--rescue-min-length must be at least 2")
    try:
        args.func(args)
    except (ValueError, OSError, KeyError) as error:
        logging.getLogger("peplocus").error("%s", error)
        return 2
    except KeyboardInterrupt:
        logging.getLogger("peplocus").error("Interrupted; no completed result summary was published")
        return 130
    return 0
