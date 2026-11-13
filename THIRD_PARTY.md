# Third-party software

## PepLocus 0.5.6

PepLocus is an external MIT-licensed Python dependency, installed by pip.
The GAPMS 1.3 wheel does not vendor the older `gapms._peplocus` implementation.
Source and license: https://github.com/qussai96/PepLocus.

## GeneSensus 0.6.0a3

The consensus engine in `gapms_v1/_genesensus` is vendored from the author's
GeneSensus 0.6.0a3 release under BSD-3-Clause. Its license and per-file original
SHA-256 hashes are included there. Absolute imports in rna_geometry.py were
changed to relative imports for the private package namespace. The original
inference algorithms and default biological coefficients are retained.

## Retained GAPMS algorithms

The protein matcher, confidence feature/seed rules and TransDecoder coordinate
projection are retained from GAPMS under this project's MIT license. They are
included directly in the v1 package, without a dependency on the legacy release.
