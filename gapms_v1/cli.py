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
    p.add_argument('-g', '--gtf', type=Path, action='append', default=[], help='repeat to assess independent annotation sources')
    p.add_argument('--search-gtf', type=Path, action='append', default=[], help='additional GTF/GFF coding paths for mapping and ambiguity checks only; not an annotation assessment target')
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
    p.add_argument('--resume', action='store_true', help='resume only when inputs, parameters and implementation fingerprint agree')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.min_mapq < 0 or args.max_paths_per_junction < 0 or not 0 <= args.max_peptide_qvalue <= 1:
        parser().error('MAPQ/path cap must be nonnegative and q-value must be in [0,1]')
    from .workflow import run
    try:
        run(args)
    except (ValueError, OSError, RuntimeError) as error:
        logging.getLogger('gapms_v1').error('%s', error)
        return 1
    return 0
