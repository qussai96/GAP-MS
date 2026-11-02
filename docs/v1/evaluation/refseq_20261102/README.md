# RefSeq-peptide rerun and annotation-removal evaluation

The nine production runs use the assembly, Helixer annotation, VARUS BAM and
technical mapping settings recorded in the preceding Casanovo jobs manifest.
Only the peptide file and output directory change. The classifier correction
separates same-frame refinements from novel ORFs; it does not change mapping.
The installed implementation SHA256 is recorded in each run manifest.

## Annotation-removal benchmark

Select up to 100 RefSeq genes per crop, with up to half single-exon and half
multi-exon targets. Eligibility requires at least two distinct, non-nested
I/L-equivalent peptide sequences in an assembly-translated reference CDS without
model issues. Select deterministically by SHA256 of `GAPMS-holdout-2026:gene`;
selection does not inspect RNA recovery. For a gene with several eligible
isoforms, select the isoform with the largest non-nested peptide set, breaking
ties by stable model ID. Report recovery against that selected target isoform.

Remove all isoforms of selected genes, plus every annotation model whose span
overlaps any selected gene span, on either strand. Retained annotations cannot
supply held-out junctions. RefSeq identifiers are namespaced by contig and strand
to handle reused organellar IDs. Models with internally overlapping exon/CDS
features are quarantined and listed explicitly, rather than silently repaired.

Run GAPMS with the masked RefSeq annotation, original assembly/BAM and all RefSeq
peptides. Reuse only the prior BAM-derived RNA reconstruction after checking
that StringTie had no annotation-guide options and TransDecoder used the
RNA-derived transcript FASTA. Mapping and the annotation/RNA splice search graph
are rebuilt. No held-out structures are supplied to GAPMS. Selection and scoring
are benchmark operations, not inputs to the mapper. Reference-database peptide
identification occurred before masking: this tests annotation rediscovery, not
blind identification of novel spectra.

Score exact strand-aware peptide block coordinates; exact full CDS paths;
start and stop boundaries; exact full coding splice chains; and individual
intron boundary-pair precision and recall. A peptide-supported model merely
overlapping a target is not an exact CDS recovery. Start, stop and splice-chain
recovery are separate measures and need not come from the same model; only the
full-CDS metric requires one model to recover everything. Single-exon targets
are excluded from the splice-chain denominator. Splice precision compares all
peptide-supported overlapping candidate junctions compatible with at least one
target peptide against the selected target isoform; other real isoforms can be
counted as unmatched. Annotation-derived target structures are not independent
experimental validation of complete CDS structures.

## Ambiguity figures

`scripts/plot_v1_ambiguity.py` reads completed `peptide_locus_summary.tsv` files. Each point
is one crop at one peptide length. Numerator: distinct I/L-equivalent query
sequences with more than one genomic placement. Denominator: mapped sequences
at that length, after collapsing I/L aliases. Boxes show the unweighted median
and interquartile range across available crops; whiskers show the range. No
mapped queries means an absent point, not zero ambiguity. Every observed length
is included; lengths below the unchanged minimum of seven were not evaluated.
The lower panel shows mapped-query sample sizes. TSV also supplies percentages
using all searched peptides as the denominator. PNG, PDF and SVG are exported.

## Reproduction scripts

The scripts take the recorded jobs manifest rather than hardcoded assembly
paths. The holdout runner expects the crop input directory layout and the
previous `GAPMS_V1_casanovo` RNA reconstruction beside the new production run.
Only annotation-independent reconstruction is reused; expected peptide paths
are used for scoring after the blinded run has completed.

```bash
python scripts/benchmark_v1_gene_holdout.py jobs.json 0 --genes 100
python scripts/audit_v1_gene_holdout.py jobs.json
python scripts/validate_v1_refseq.py jobs.json
python scripts/plot_v1_ambiguity.py jobs.json --out figures/refseq --label 'RefSeq FragPipe'
python scripts/summarize_v1_refseq.py /path/to/experiment_directory
```

Repeat the holdout runner with indices 0 through 8. Its output directory must
be new. Keep the jobs manifest, masked annotations, source reconstruction
provenance and per-target recovery tables alongside the final summaries.

## Interpretation and reproducibility

The production runs repeat all inputs and technical settings except the peptide
list and output directory. Fresh RNA exon-structure multisets were identical
to the Casanovo runs in all nine crops. Fresh projected TransDecoder CDS
structure multisets shared 96.0–99.7% with the earlier predictions, despite the
same assembly, BAM and settings. Therefore differences in ORF counts between
peptide datasets cannot all be attributed to the peptide list. This is an
observed repeatability limitation of the external ORF-prediction stage; its
specific cause has not been established here. The controlled holdout benchmark
uses fixed RNA predictions to remove this source of variation.

The benchmark runner calls the installed v1 workflow and substitutes only the
reconstruction function with audited cached output. This is a benchmark runner,
not an undocumented CLI cache option. Genome mapping, model association and
classification use the installed implementation. Each run records the reused
RNA source and commands. The independent audit checks source assembly/BAM
identities, identical cached RNA GTF bytes, absence of overlapping retained
annotations and all withheld peptide coordinate projections/translations.

A peptide-supported candidate ORF does not establish every CDS boundary. The
full-CDS and splice-chain recovery metrics test that stronger claim directly.
Neither candidate confidence labels nor annotation-concordance precision are
calibrated probabilities that every proposed boundary is biologically correct.
Historical Casanovo ORF proposal counts preceded the classification correction;
18,515 models occurred in both novel and refinement lists. Their mapping and
ambiguity counts are unaffected, but those proposal totals must not be quoted
as distinct novel-ORF counts.


Detailed per-target CDS recovery tables and the 31 missed peptide paths remain in the local evaluation directory. This repository report contains aggregate results and figure data.
