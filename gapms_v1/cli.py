"""Four biological input types; all other options control technical policy."""
import argparse
import logging
from pathlib import Path

from . import __version__


def positive(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return value


def parser():
    p = argparse.ArgumentParser(description='GAPMS v1: peptide-driven genome mapping and annotation evidence')
    p.add_argument('--version', action='version', version=f'GAPMS v1 {__version__}')
    p.add_argument('-a', '--assembly', type=Path, required=True)
    p.add_argument('-p', '--peptides', type=Path, required=True)
    p.add_argument('-g', '--gtf', type=Path, action='append', default=[], help='repeat for automatic peptide-guided consensus of all supplied sources')
    p.add_argument('-b', '--bam', type=Path)
    p.add_argument('-o', '--output', type=Path, required=True)
    p.add_argument('-t', '--threads', type=positive, default=8)
    p.add_argument('--peptide-column')
    p.add_argument('--min-length', type=positive, default=7, help='minimum exact-match peptide length (default: 7); excluded queries are reported')
    p.add_argument('--genetic-code', type=positive, default=1)
    p.add_argument('--rescue-one-edit', action='store_true', help='separate exploratory rescue; never used as exact model/event evidence')
    p.add_argument('--rescue-min-length', type=positive, default=12)
    p.add_argument('--min-junction-support', type=positive, default=2)
    p.add_argument('--min-mapq', type=int, default=10)
    p.add_argument('--min-read-anchor', type=positive, default=8)
    p.add_argument('--max-junctions', type=positive, default=3)
    p.add_argument('--max-paths-per-junction', type=int, default=1000)
    p.add_argument('--min-orf-aa', type=positive, default=100, help='TransDecoder minimum ORF size (default: 100 aa; partial/single-peptide ORFs remain reportable)')
    p.add_argument('--max-peptide-qvalue', type=float, default=0.01)
    p.add_argument('--refinement-window', type=positive, default=3000, help='maximum contiguous CDS extension searched in nucleotides (default: 3000)')
    p.add_argument('--refinement-splice-window', type=int, default=30, help='splice-site shift radius in nucleotides (default: 30; zero searches single-exon skips at existing sites)')
    p.add_argument('--max-refinement-paths', type=positive, default=64, help='maximum canonical phase-preserving paths per local exon pair')
    p.add_argument('--min-refinement-anchor', type=positive, default=9, help='minimum peptide nucleotides on each side of a discovered junction')
    p.add_argument('--no-refinement-splice', action='store_true', help='disable annotation-guided local splice discovery; terminal refinement remains enabled')
    p.add_argument('--consensus-beam', type=positive, default=8, help='bounded paths per graph start and phase')
    p.add_argument('--consensus-max-nodes', type=positive, default=2000)
    p.add_argument('--consensus-max-starts', type=positive, default=32)
    p.add_argument('--consensus-max-candidates', type=positive, default=32, help='new graph candidates per locus; original models are never pruned')
    p.add_argument('--rna-library', choices=['unstranded','firststrand','secondstrand'], default='unstranded', help='RNA library protocol for consensus evidence; default: unstranded')
    p.add_argument('--keep-work', action='store_true', help='retain temporary search files after successful completion; failed runs always retain them')
    p.add_argument('--resume', action='store_true', help='resume only when inputs, parameters and implementation fingerprint agree')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.refinement_splice_window < 0:parser().error('refinement splice window must be nonnegative')
    if args.min_mapq < 0 or args.max_paths_per_junction < 0 or not 0 <= args.max_peptide_qvalue <= 1:
        parser().error('MAPQ/path cap must be nonnegative and q-value must be in [0,1]')
    from .workflow import run
    try:
        run(args)
    except (ValueError, OSError, RuntimeError) as error:
        logging.getLogger('gapms_v1').error('%s', error)
        return 1
    return 0
