# Scientific methods

## Search space and coordinates

Peptides are normalized to the 20 standard amino acids; common bracketed or
parenthesized modifications and flanking residues are removed for searching,
while table records are preserved. Exact matching always treats I/L as
indistinguishable. Original normalized query labels remain in output. Evidence
counts deduplicate I/L-equivalent sequences and repeated identification records.
The default minimum length is seven amino acids, with an explicit excluded list.
This is a computational ambiguity policy, not a peptide quality filter.

PepLocus translates overlapping genomic chunks in six frames. A C-backed
Aho–Corasick automaton finds all exact peptide occurrences. The engine does not
materialize a whole-genome six-frame FASTA. Process workers scan bounded chunks;
a disk-backed store deduplicates chunk overlaps. The pipeline searches provided
annotation and RNA-supported splice paths, including split codons and microexons.
BAM junction observations use primary, nonduplicate, non-QC-failed alignments;
minimum MAPQ is 10, adjacent aligned anchors are eight bases, introns are
20–500,000 bases, and at least two observations are required. Counts are alignment
observations, not necessarily independent molecules. XS/ts tags and canonical
splice motifs resolve strand; unstranded noncanonical evidence can support both
strand hypotheses. This uncertainty must be considered during review.

The default local path cap is three introns per peptide and 1,000 paths per
first junction. Caps are recorded by the mapping engine, including
`mapping/mapping_summary.json`. Full provided and reconstructed CDS translations
are also searched and projected without the local path cap. This recovers
supported known full-model paths without claiming arbitrary unobserved splice
connectivity. Multiplicity is relative to the searched assembly, available
models and supported local paths. Absence from unsearched paths is not proof of
global biological uniqueness. De novo unsupported splice discovery is not
silently enabled.

Each input annotation has an independent source identity. Duplicate feature rows
are collapsed within a transcript; overlapping features, unknown contigs and
invalid coordinates fail explicitly. GTF, GFF3 Parent lists and CDS-only
annotations are supported. Exon-only models contribute junction context but are
not assigned invented coding frames. The first CDS phase is trimmed once, in
translation direction; internal phases are consistency checks, not bases to
remove at every exon. Unknown first phase, inconsistent internal phase, internal
stop or terminal partial codon is recorded. Models with such issues cannot receive
unambiguous peptide-supported status or high-confidence candidate status.

RNA mode uses coordinate-sorted BAM, StringTie transcript reconstruction,
Python extraction of spliced transcript DNA, and TransDecoder LongOrfs/Predict.
Unstranded StringTie transcripts retain their `.` orientation in a transcript
orientation table and are extracted in genomic-plus direction; TransDecoder
then determines the ORF direction. The minimum ORF length is 100 aa by default
and can be lowered. ORF coordinates
are projected back with the existing tested GAPMS adapter, then independently
translated against the assembly for the model search. Complete and partial
candidates remain distinguishable; partial ORFs are retained rather than being
silently removed. BAM reconstruction currently requires the standard genetic
code. Other NCBI codes are supported for mapping and annotation-only modes.
The conservative start-completeness flag recognizes ATG, not all possible
non-AUG initiation sites. ORFs absent from TransDecoder's search space, including
shorter ORFs below its configured minimum, are not discovered by this release.

The new evidence SQLite database contains the union of genomic and verified CDS
placements. Identical placements are reconciled across sources and models. Adjacent blocks
are merged before locus counting, so a segmentation boundary without an intron
does not create artificial ambiguity.
A locus is an exact contig, strand and ordered genomic block path; repeated
positions in one gene can therefore make a peptide non-unique. Protein searches
build one automaton over the union of source-aware models. Every projected hit
is re-extracted and translated. Coordinates are zero-based, half-open internally
and in TSV/BED; GTF is one-based, inclusive.

## Model support and peptide quality

A model is `peptide_supported` when it has a locus-discriminating peptide and no
recorded coding-path issues. Shared-only or structurally uncertain support is
`ambiguous`. No detected sequence is `no_observed_peptide_support`, which never
means that the prediction is biologically incorrect. `no_CDS_annotation` is
reported separately. Model and gene counts collapse repeated sequences and I/L
aliases. They do not claim isoform specificity merely from locus specificity:
compatible isoforms can share the same mapping ID.

Recognized q-value columns are retained and compared with 0.01 by default.
Any available passing identification of a normalized sequence establishes
`qvalue_pass`; all available recognized q-values failing establishes
`qvalue_fail`. Otherwise quality is `unknown`. PEPs, spectrum counts,
modifications and other fields are preserved in raw identification JSON but are
not automatically interpreted as q-values. The software does not estimate FDR,
calibrate Casanovo confidence, or infer independence of repeated spectra.
A sequence in a header-only Casanovo list has unknown identification quality.
Mapping hits are not validated identifications. The original XGBoost model was
trained for the legacy protein-assessment context; v1 does not extrapolate it
to novel ORFs or present its scores as calibrated probabilities.

## Novel ORFs and event specificity

Only RNA-derived candidate coding models are proposed as new ORFs in this
release. Models with identical coding paths to an input annotation are not
novel. Supported differing candidates are classified against overlapping
annotation models as intergenic, antisense, intronic, UTR-associated,
overlapping alternative-frame or new transcript isoform. A same-frame,
overlapping coding model is assessed as an annotation refinement, not counted
as a novel ORF. Antisense and alternative-frame ORFs can overlap an annotated
locus without being refinements of its existing coding path. The classification uses annotation CDS/exon intervals
and relative coding frame at shared bases; missing or incomplete annotations
can change that relationship. UTR-associated means exonic outside an available
CDS, not experimentally validated UTR translation. Annotation source conflicts
remain visible rather than being resolved into a presumed truth.

`orfs/candidate_classification.tsv` assigns each proposed RNA model one primary
category. A model with just one supported refinement is classified as
`alternative_start`, `alternative_stop`, `c_terminal_extension` or `alternative_splice`; a model with
several supported changes is `combined_refinement`. Its component changes
remain in the event tables. Thus event counts can overlap between refinement
types, but primary model counts are mutually exclusive and no refinement is
also counted as a novel ORF. An unsupported boundary difference is not promoted
to a proposal merely to give every candidate an event label.

Diagnostic placements are peptide paths absent from the input annotations'
compatible peptide placements. Single-peptide ORFs are retained with lower
confidence. Nested peptide sequences are not counted as two independent
non-nested sequences. `high_confidence_candidate` requires at least two
non-nested, locus-discriminating diagnostic sequences, coherent complete ATG-to-
stop translation with a fully resolved standard-amino-acid sequence, passing recognized peptide q-value evidence and, for spliced
candidates, sufficient support for every coding junction. A structurally different coding model compatible with all diagnostic sequences
forces an ambiguous assignment, even when the locus is unique. Identical coding
paths are treated as the same hypothesis. It is still a candidate, not an
automatically validated annotation. Unknown quality, partial
coding paths and shared-only evidence preclude that category. Other categories
are `rna_peptide_supported`, `peptide_supported` or `ambiguous`; confidence is
rule-based, not a calibrated probability.

Alternative starts require a changed coding start plus an ATG-compatible
candidate and a diagnostic peptide at protein residue zero, or residue one
with plausible initiator methionine removal (A/C/G/P/S/T/V at residue two),
extending outside the existing start. Previous amino acid and tryptic context
are recorded for every model occurrence: K/R followed by a non-proline is
tryptic. A distant internal tryptic peptide cannot nominate a new start.
N-terminal modifications remain available in the identification records but
are not given a custom initiation likelihood. Proteolytic processing remains
an explicit competing interpretation; non-tryptic status alone is insufficient.
The current conservative detector primarily supports coding extensions and
will miss internal alternative starts whose peptides are indistinguishable
from an existing protein.

Alternative CDS ends require a diagnostic peptide extending beyond an
annotation's coding end. An `alternative_stop` additionally requires a complete in-frame stop in the assembly-derived ORF; without it the event is `c_terminal_extension` and termination is unresolved. Inferred stop positions are labeled as such. Terminal changes accompanied by splice-chain changes are explicitly flagged.
Missing downstream detections never establish shortening. Some alternative
ends, especially shortening, cannot be resolved with these peptide inputs.

Alternative splice candidates require an exact peptide directly spanning a
candidate CDS junction absent from overlapping same-strand reference CDS paths.
A peptide elsewhere in a protein is not diagnostic. RNA support for every
candidate coding junction is recorded, and novel splice candidates cannot
reach high confidence without it. Local junction support does not establish
that all junctions coexist in the same full-length transcript. StringTie
connectivity is a reconstruction hypothesis; the uncertainty is explicit.

## Optional one-edit rescue

Exact engine mapping runs first. Optional one-edit rescue searches only engine-stage exact-unmapped
peptides of at least 12 aa by default, allowing one amino-acid substitution,
insertion or deletion. I/L equivalence does not consume the edit. Rescue
alignments, limits and ambiguity are reported separately in `mapping/rescue/`.
They never count as exact model support or diagnostic annotation evidence.
The nine-crop default evaluation uses exact matching; rescue must be requested
explicitly to avoid silently increasing ambiguity and runtime.

## Further scope limits

Assemblies are treated as linear contigs. Circular-origin coding paths are not
automatically wrapped. A single genetic code applies to every contig in a run;
mixed nuclear/organelle translation tables need separate configured runs. The
reference assembly is the searched sequence: optional one-edit placements do
not validate a variant allele. External reconstruction needs sufficient RNA
coverage and ORFs for the assembler/predictor; tool failures are reported with
logs rather than silently inventing candidate models.

Deterministic IUPAC codons are resolved when every possible concrete codon
encodes the same residue. This can differ from TransDecoder's conservative X
masking. Unresolved amino acids elsewhere in an ORF remain a separate uncertainty;
they preclude high-confidence candidate status even when mapped peptides avoid
that region. The whole-ORF projection audit checks nucleotide identity and all
possible codon expansions before accepting an IUPAC representation difference.

One-edit rescue is scoped to the bundled engine's exact-unmapped set. Full-CDS
exact projections are integrated afterward and can add exact paths beyond local
splice caps. Interpret exploratory rescue against the final exact mapping
table; it never contributes to model assessment or annotation proposals.


## Local evidence and full-model review in 1.0.1

Diagnostic peptide blocks and new junctions are explicit event fields. They establish local coding compatibility, not every feature of the RNA model. Full candidate GTFs remain context models. `whole_model_status` separates `local_evidence_only` from `complete_model_for_manual_review`; neither is an automatic replacement. Review reasons flag partial starts/stops, coding-path issues, ambiguous diagnostic support, junctions below the RNA threshold, fewer than two independent diagnostic peptides, multiple overlapping annotation gene identities, and loss of compared annotation CDS sequence. The last condition deliberately prevents local evidence from implicitly validating wider deletions. It can also flag valid changes, so the evidence remains reportable. Full-structure correctness still requires independent evaluation.

A coding-structure ID hashes contig, strand and phased translation blocks. It groups identical paths under different RNA transcript IDs within an assembly; it is not globally unique across assemblies. Counts of RNA model IDs, distinct CDS paths and event comparisons have different denominators.

The search uses the supplied assessment annotations and BAM-supported or reconstructed paths. An evaluation reference is not added to the mapping graph. Uniqueness remains relative to the searched paths; unknown splice structures can contain other locations. Version 1.0.2 removes the auxiliary search-annotation input without removing the stop-codon or local-evidence fixes. Historical reference-assisted results from 1.0.1 must be distinguished from this workflow.
