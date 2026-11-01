# Fast model metadata attachment (GAPMS 0.2.2)

Assembly runs with reference comparison now call `model_metadata.attach_model_metadata`
after PepLocus, instead of the historical whole-protein projection routine.
The standard and IGV genomic BEDs and the global-ambiguity/report join are unchanged.

## Implementation

The historical exporter rebuilt and walked the complete CDS nucleotide-coordinate
list separately for each peptide occurrence, repeatedly copying DataFrames and
sorting CDSs. The new method constructs that coordinate interpretation once per
protein and compresses it into contiguous genomic runs. A cumulative run index
and binary search project each peptide onto only its intersecting runs. No assembly
search or translation is repeated for metadata attachment.

The same established first-CDS phase, strand, interval ordering, validated protein
positions, terminal classifications, preceding residue, and deduplication rules
are preserved. Overlapping/reordered pathological CDS walks retain an exact cached
coordinate-array fallback. Repeated Protein/Gene annotation rows are removed before
the peptide join, avoiding the redundant CDS expansion while retaining first-row
order. The existing BED serializer preserves field formatting and sort order.

The compatibility filename remains `protein_projection/peptides_mapped.bed` to
preserve report and external interfaces. Its name does not mean the previous slow
routine is still called. Reports still attach metadata only to corresponding
PepLocus genomic loci and retain all other genomic alternatives in ambiguity tests.
Protein-only runs without an assembly retain their historical exporter.

## Validation

Slurm job 188282 ran isolated full metadata attachment on tomato, spinach and
carrot. Existing nine-species workflow outputs were not changed. On this single
run per species:

| Species | Metadata records | Fast attachment elapsed | Legacy comparison |
|---|---:|---:|---|
| Tomato | 236,531 | 42.1 s | Byte-identical |
| Spinach | 198,171 | 35.5 s | Byte-identical |
| Carrot | 297,855 | 50.9 s | Legacy baseline not complete when checked |

Tomato's prior log spans about 382 seconds from PepLocus completion to the start
of reference comparison, including legacy metadata work. The new 42.1-second
measurement is approximately nine times shorter, although these single-run timings
include different overhead/cache effects and are not a controlled scaling study.

Job 188283 regenerated the complete tomato annotation report using the new metadata
and the same genomic BED, annotation, reference, scores and TMAP. Both TSV reports
and all six event-specific FASTA/GTF exports were byte-identical to the legacy run.
The complete suite passes 92 tests, including randomized nucleotide-projection
oracles over strands and phases, adjacent/overlapping CDSs, byte-identical metadata,
and identical reference-report DataFrames and files.

Reproduction scripts: `scripts/validate_metadata.sbatch` and
`scripts/validate_metadata_report.sbatch`. Detailed local results are in
`evaluation/model_metadata/`; compact archived checks are in
[model_metadata_results.json](model_metadata_results.json).

GAPMS 0.2.2 applies to newly launched workflows. Already running nine-species jobs
were not restarted or rewritten and retain their loaded adapter implementation.
