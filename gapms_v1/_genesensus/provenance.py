"""Content fingerprints and atomic JSON writes; no mutation of input files."""

from pathlib import Path
import hashlib, json, os, platform, sys, importlib.metadata, uuid
from . import __version__


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_object(obj):
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def verify_inputs(expected, context):
    """Do not publish derived cache data under fingerprints of changed inputs."""
    for path, digest in expected.items():
        if sha256(path) != digest:
            raise ValueError(f"Input changed before publishing {context}: {path}")


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(
        path.name + ".tmp." + str(os.getpid()) + "." + uuid.uuid4().hex
    )
    with tmp.open("w") as f:
        json.dump(value, f, sort_keys=True, indent=2, allow_nan=False)
        f.write("\n")
    os.replace(tmp, path)


def software_record():
    versions = {}
    for package in ["numpy", "biopython", "pysam", "h5py", "torch"]:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    return dict(
        genesensus=__version__,
        python=sys.version,
        platform=platform.platform(),
        dependencies=versions,
        package_code={p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")},
    )
