"""GeneSensus: reference-independent gene reconstruction and refinement."""

__version__ = "0.6.0a3"


def run(config_path, output, threads=None, cache=None, resume=False):
    """Run deterministic evidence integration; return a pipeline.RunResult."""
    from .pipeline import run as _run

    return _run(config_path, output, threads=threads, cache=cache, resume=resume)
