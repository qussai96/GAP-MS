# Validation

## Legacy baseline

Before v1 changes, commit `40cabc4` passed all 105 legacy tests in Python 3.10.
The three expected warnings concern legacy protein-only projection without an
assembly. The original classifier, mapping namespace and model file are unchanged.

## Synthetic tests

The final local suite passes **129 tests**: 105 unchanged legacy regressions and
24 v1 tests. The installed v1 wheel contains only `gapms_v1`, not a second copy
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

All nine runs and the final aggregate audit completed successfully. See the
[full evaluation](evaluation/EVALUATION.md), including raw historical coordinate
overlap, independent sequence checks, mapping multiplicity, event counts,
length-stratified ambiguity and single-run timings. The figure and its source
summary data are included alongside the report.

Across all nine datasets, 1,316,152 of 2,420,861 searched normalized query labels
mapped (minimum length seven AA). All 2,747,138 original identification records
were retained. A further 84,134 distinct short sequences were explicitly listed
as excluded. At the exact coordinate level, 707,306 of 707,670 independently
sequence-reproducing historical reference paths were recovered (99.9486%). No
verified contiguous location was missed. All 364 remaining verified paths use
junctions absent from the supplied annotation and retained RNA evidence.
Raw reference overlap and reference sequence inconsistency are reported
separately; no claim of biological peptide-identification precision is made.

Independent audit translated 18,000 sampled peptide paths and 18,000 whole ORFs,
checked all BED schemas/count identities, and inspected diagnostic event
associations. The final Casanovo metadata correction was separately proven to
leave normalized query bytes and every mapped peptide ID/quality unchanged.
Archived manifests and command records disclose compatible-stage reuse.

The biological runs took 15.74–39.54 minutes in these single measurements;
the genome/path engine took 1.65–5.04 minutes. TransDecoder and RNA processing
account for most remaining time. These are descriptive runs, not replicated
scaling estimates or a demonstrated end-to-end speedup against legacy GAPMS.

A further full audit compared every one of the 438,013 assembly-derived ORF
proteins against TransDecoder's transcript-space prediction. There were 438,011
literal protein matches. Two cases (eight residues total) use deterministic
IUPAC codons: TransDecoder masks them as X, whereas v1 resolves the amino acid.
For both, transcript and genome CDS nucleotide paths were identical; enumerating
every concrete codon gave exactly the reported residue at each difference.
No unexplained projection differences remain. The initial literal-only audit
flagged these two cases and is preserved in the local audit history.

Whole ORFs with unresolved residues cannot receive high-confidence candidate
status, even if their termini and peptide q-values pass other criteria. Such
uncertainty is explicit in candidate limitations. The final dataset update
changed only those limitation annotations; all placements and confidence counts
remained unchanged for these unknown-quality inputs.
