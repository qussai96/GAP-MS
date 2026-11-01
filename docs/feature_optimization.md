# Feature extraction optimization in GAPMS 0.2.1

GAPMS now defaults to eight mapping workers in both the full workflow and
`gapms map`. Explicit `-t`/`--threads` values override the default. This changes
GAPMS defaults; the independently installed PepLocus package is unchanged.

CDS nucleotide offsets now use grouped cumulative sums instead of constructing
and iterating over one DataFrame per protein. The first CDS phase is still trimmed,
subsequent phases retain split codons, and minus-strand order is preserved. Inputs
with duplicate start positions or duplicate row indexes retain the original
reference implementation to preserve its tie-order semantics.

Adjacent CDS introns are computed with grouped shifts. Junction support across
repeated peptide occurrences is computed by counting distinct supporting
occurrences and requiring that count to equal the total occurrence count. Boolean
group reductions replace Python callbacks for unambiguous terminal/splice flags.
These changes preserve the original evidence rules and output ordering.

## Carrot validation

Slurm job 188270 ran the same full Helixer/FragPipe workflow as integration job
188267, with the previous protein mapping held fixed. The new run omitted `-t`;
its mapping summary confirms eight workers. Times below are log timestamp
intervals on the same cluster, measured once per configuration; they are not
replicated scaling estimates.

| Measurement | Previous | Optimized |
|---|---:|---:|
| Build CDS features | 198 s | 2 s |
| Splice-flag processing | 275 s | 9 s |
| Full workflow | 11.80 min | 4.13 min |
| PepLocus mapping | 133.941 s | 131.8 s |

The full workflow was approximately 2.86 times faster (65% less elapsed time).
These two feature stages are predominantly serial grouped operations; the speedup
comes from reducing Python/DataFrame overhead rather than the worker default.
Protein sequence coverage, peptide positioning, classifier behavior, and the
mapping algorithm were not changed.

The following outputs were compared byte-for-byte against job 188267:
feature scores, supported protein FASTA/GTF, all protein-selection lists, unmapped
peptides, both BED files, plot images, and detailed engine mapping/evidence tables.
Run paths, timestamps, timings, and version metadata naturally differ in logs and
summary JSON and are excluded from byte-equality checks.

The new tests compare optimized CDS/junction calculations with the original
implementation across randomized annotations, both strands, phases, adjacent
segments, unusual identifiers, empty inputs, and tie-order cases. Existing tests
cover repeated peptide occurrences and shared-junction support rules. The complete
GAPMS suite passes 82 tests.

Artifacts: `evaluation/feature_speed/run_188270/equality.json`, workflow logs and
Slurm logs. Reproduction: `scripts/validate_feature_speed.sbatch`. The script checks
all listed result files and fails if any differs from the saved baseline.
