"""Command line entry point. Evaluation references are never accepted by run."""

import argparse, json, sys, logging
from pathlib import Path
from . import __version__


def parser():
    p = argparse.ArgumentParser(
        prog="genesensus",
        description="Evidence-guided gene reconstruction and refinement",
    )
    p.add_argument("--version", action="version", version="GeneSensus " + __version__)
    p.add_argument("--debug", action="store_true", help="Show a traceback on error")
    commands = p.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="Integrate predictor annotations and optional evidence"
    )
    run.add_argument("config", help="Inference-only JSON configuration")
    run.add_argument("-o", "--output", required=True, help="New output directory")
    run.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help="Independent graph worker processes",
    )
    run.add_argument("--cache", help="Reusable source, graph and RNA feature cache directory")
    run.add_argument(
        "--resume",
        action="store_true",
        help="Verify and reuse a completed identical output",
    )
    validate = commands.add_parser(
        "validate", help="Validate configuration, input paths and genome coordinates"
    )
    validate.add_argument("config")
    prepare = commands.add_parser(
        "prepare-rna", help="Prepare observed splice and exon evidence from an indexed BAM"
    )
    prepare.add_argument("--bam", required=True)
    prepare.add_argument("--genome", required=True)
    prepare.add_argument("-o", "--output", required=True)
    prepare.add_argument("-t", "--threads", type=int, default=1)
    prepare.add_argument(
        "--library",
        choices=["unstranded", "firststrand", "secondstrand"],
        default="unstranded",
    )
    prepare.add_argument("--min-mapq", type=int, default=20)
    prepare.add_argument("--min-overhang", type=int, default=8)
    prepare.add_argument("--minimum-templates", type=int, default=3, help="Minimum templates for novel RNA exon discovery; must match the run configuration")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    try:
        if args.command == "run":
            from .pipeline import run

            result = run(
                args.config, args.output, args.threads, args.cache, args.resume
            )
            print(
                json.dumps(
                    dict(
                        output=str(result.output),
                        annotation=str(result.annotation),
                        reused=result.reused,
                        selected=result.report["selected_count"],
                        novel_graph_selected=result.report["selected_graph_candidates"],
                        seconds=result.report["timings"]["total_seconds"],
                    ),
                    sort_keys=True,
                )
            )
            return 0
        if args.command == "validate":
            from .config import load_config, input_paths
            from .annotation import read_annotation
            from .genome import Genome

            cfg = load_config(args.config)
            paths = input_paths(cfg)
            counts = {}
            with Genome(cfg["genome"]) as genome:
                lengths = genome.lengths
                for source, spec in cfg["predictions"].items():
                    models = read_annotation(spec["annotation"], source)
                    counts[source] = len(models)
                    for m in models:
                        if m.contig not in lengths or m.cds[-1][1] > lengths[m.contig]:
                            raise ValueError(
                                "Annotation/genome coordinate mismatch: " + m.identifier
                            )
            print(
                json.dumps(
                    dict(valid=True, predictors=counts, input_files=len(paths)),
                    sort_keys=True,
                )
            )
            return 0
        if args.command == "prepare-rna":
            from .genome import Genome
            from .rna import extract_junctions, write_junctions
            from .rna_geometry import augment_geometry

            out = Path(args.output)
            if out.exists():
                raise ValueError("RNA output already exists; refusing overwrite")
            out.parent.mkdir(parents=True, exist_ok=True)
            with Genome(args.genome) as genome:
                data = extract_junctions(
                    args.bam,
                    genome.lengths,
                    args.threads,
                    args.min_mapq,
                    args.min_overhang,
                    library=args.library,
                )
                data = augment_geometry(args.bam, genome, data, threads=args.threads, minimum_templates=args.minimum_templates)
            temp = out.with_name("." + out.name + ".partial")
            if temp.exists():
                raise ValueError("Partial RNA output already exists")
            write_junctions(data, temp)
            temp.rename(out)
            print(json.dumps(dict(output=str(out), junctions=len(data["junctions"]))))
            return 0
    except KeyboardInterrupt:
        print(
            "GeneSensus interrupted; incomplete run artifacts are preserved.",
            file=sys.stderr,
        )
        return 130
    except Exception as exc:
        if args.debug:
            raise
        print("GeneSensus: " + str(exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
