# Validation

## Legacy baseline

Before v1 changes, commit `40cabc4` passed all 105 legacy tests in Python 3.10.
The three expected warnings concern legacy protein-only projection without an
assembly. The original classifier, mapping namespace and model file are unchanged.

## Synthetic tests

The final local suite passes **127 tests**: 105 unchanged legacy regressions and
22 v1 tests. The installed v1 wheel contains only `gapms_v1`, not a second copy
of the legacy `gapms` namespace. The local dependency-wheel builder was also
exercised; its wheel includes the pretrained model and excludes the removed
ProteoMapper scripts.

The v1 suite exercises genome-only, annotation, RNA and integrated workflows.
Synthetic RNA workflow tests replace external transcript reconstruction with
known coding models; BAM parsing, junction extraction, genome mapping, CDS
association, evidence integration and output generation remain real. Real
StringTie/TransDecoder validation is conducted in crop Slurm runs, separately
from these controlled tests.

Fixtures cover positive/negative strands, reading-frame offsets, I/L aliases,
multiple placements, split codons, negative-strand microexons, multiple annotation
sources, no-observed-support categories, malformed inputs, alternative genetic
code, guarded resume and read-only input indexes. Coordinates are independently
re-extracted and translated with BioPython. Event tests cover two-peptide
intergenic candidates, shared single-peptide ambiguity, nested evidence,
overlapping alternative frames, N-terminal specificity, coding extensions and
direct splice-junction evidence. Missing C-terminal peptides never create a
shortening event. Original q-values, PEP fields and spectral metadata remain
available for inspection. Additional fixtures check that adjacent blocks form
one locus and that equally compatible ORF structures remain ambiguous.

## Nine-crop Casanovo evaluation

The requested integrated runs use assembly, Helixer annotation, coordinate-
sorted VARUS BAM and `outputs/casanovo/*.casanovo_peps.tsv`. Output directories
are `outputs/gapms/GAPMS_V1_casanovo`. The BAM filename in these datasets is
`sorted_VARUS.bam` (uppercase VARUS). Unrelated outputs are preserved. Maize's
`.fna.upper` is a redundant sidecar, not a second assembly.

Actual job identifiers, completion states, runtimes, mapping ambiguity and
model/event counts will be recorded after the jobs are inspected. Submission
alone is not validation. Casanovo inputs currently supply sequences without
calibrated identification q-values; high-confidence event status must not be
invented from mapping alone.
