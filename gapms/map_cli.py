"""Standalone GAPMS assembly mapper, backed by bundled PepLocus."""
import argparse
from pathlib import Path


def main(argv=None):
    from ._peplocus.cli import parser as engine_parser
    # Reuse every engine search option and its validation, adding GAPMS aliases.
    parser = engine_parser()
    map_parser = parser._subparsers._group_actions[0].choices['map']
    aliases = {'--genome': ['-a', '--assembly'], '--peptides': ['-p'],
               '--gtf': ['-g'], '--bam': ['-b'], '--out': ['-o', '--output'],
               '--threads': ['-t']}
    for name, names in aliases.items():
        action = map_parser._option_string_actions[name]
        for alias in names:
            action.option_strings.append(alias)
            map_parser._option_string_actions[alias] = action
    map_parser._option_string_actions['--out'].required = False
    map_parser._option_string_actions['--out'].default = 'GAPMS_Map'
    args = parser.parse_args(['map'] + (list(argv) if argv is not None else __import__('sys').argv[1:]))
    for name in ['genome', 'peptides', 'gtf', 'bam', 'junctions']:
        value = getattr(args, name)
        if value is not None and not Path(value).is_file():
            map_parser.error(f'input file does not exist: {value}')
    if args.min_intron > args.max_intron:
        map_parser.error('--min-intron cannot exceed --max-intron')
    if args.discover_all and not args.discover_junctions:
        map_parser.error('--discover-all requires --discover-junctions')
    from .peptides_to_genome import map_assembly
    # Forward parsed engine options through its own CLI for one source of truth.
    excluded = {'command', 'genome', 'peptides', 'gtf', 'bam', 'out', 'threads',
                'peptide_column', 'func', 'igv_display_bed', 'il_equivalent'}
    extra = []
    for name, value in vars(args).items():
        if name in excluded or value is None or value is False:
            continue
        option = '--' + name.replace('_', '-')
        extra += [option] if value is True else [option, str(value)]
    try:
        map_assembly(args.genome, args.peptides, args.out, gtf=args.gtf, bam=args.bam,
                     peptide_column=args.peptide_column, threads=args.threads,
                     extra_args=extra)
    except (ValueError, OSError, RuntimeError) as error:
        map_parser.error(str(error))
    return 0
