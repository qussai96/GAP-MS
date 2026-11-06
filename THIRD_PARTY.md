
## Bundled PepLocus

The gapms._peplocus package vendors PepLocus 0.3.1 under its MIT license.
It adds recognition of the FragPipe "Peptide Sequence" header. Source: the
PepLocus package developed alongside GAPMS. License: gapms/_peplocus/LICENSE.

## GeneSensus algorithm acknowledgment

GAPMS v1.2.0 consensus is independently implemented following the bounded
exon/junction graph and conserved peptide ambiguity principles documented
in GeneSensus 0.6.0a3 (2026, GeneSensus contributors; BSD-3-Clause). It does
not bundle the GeneSensus package or claim identical selection behavior.

## GeneSensus 0.6.0a3

The consensus engine in `gapms_v1/_genesensus` is vendored from the author's
GeneSensus 0.6.0a3 release under BSD-3-Clause. Its license and per-file original
SHA-256 hashes are included there. Only two absolute imports in rna_geometry.py
were changed to relative imports for the private package namespace. Inference
algorithms and default biological coefficients are unchanged.
