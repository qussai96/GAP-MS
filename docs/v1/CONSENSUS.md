# Automatic consensus in GAPMS 1.2.1

Supplying two or more distinct `-g` files always produces a consensus. There is
no consensus switch. A single annotation retains the existing assessment and
local-refinement workflow; BAM is optional.

```bash
gapms-v1 -a genome.fa -p peptides.tsv -g helixer.gtf -g annevo.gtf -o results -t 8
gapms-v1 -a genome.fa -p peptides.tsv -g helixer.gtf -g annevo.gtf \
  -b reads.bam --rna-library firststrand -o results_rna -t 8
```

## Which engine is used

GAPMS 1.2.0 independently implemented part of GeneSensus's method. It did not
reproduce its selection behavior. Version 1.2.1 replaces that implementation
with the actual frozen GeneSensus 0.6.0a3 engine, included privately in the wheel
under BSD-3-Clause. The only engine edits change two absolute imports to relative
imports. Original file hashes, the change description and license are shipped
in `gapms_v1/_genesensus/`; a regression test checks the original algorithm bytes.

The fixed policy is `state_molecular_guided`, `molecular_strategy=source_context`,
`peptide_strategy=graph_contrast`, with the engine's biological defaults.
The molecular policy name does not imply molecular data exist: GAPMS supplies
no protein-likeness scores, residue profiles or SpuriousORF classifications.
Those missing channels abstain. There is no RefSeq-based source ranking or
reference-assisted search. All supplied annotations are inference inputs.

The [frozen release method](GENESENSUS_0.6.0a3.md) describes the underlying
algorithm and its broader optional inputs. Only inputs exposed by GAPMS are
used here; that document's molecular examples are not implicit GAPMS settings.
Multi-annotation consensus currently requires genetic code 1.

## Inputs and mapping

GAPMS's annotation parser continues to accept its supported GTF/GFF variants,
including compressed files, orphan CDS rows and shared GFF3 Parent lists.
The adapter writes explicit per-source transcript parents, safe IDs, exons,
CDS coordinates and initial phases. It recomputes downstream phase fields from
the same initial phase used for peptide projection, preserving the coding path.
Existing phase issues are recorded, not hidden as a sequence-based repair.
An unspecified first phase retains GAPMS's existing zero-phase assumption,
which is reported in the source audit. Input gene families remain grouped and
source IDs remain distinct. Exon-only models remain in assessment but do not
become coding consensus candidates. Every supplied annotation must contain CDS.

Genome/splice mapping and source-protein association run first. The frozen
engine then builds structural and peptide-guided graph candidates. Their proteins
are searched with the same ProteinMatcher automaton already built for source
proteins. Every peptide projection is checked against the assembly before being
added to the shared mapping database. This automaton is distinct from the
assembly mapper and local-refinement searches.

All I/L-equivalent sequences are one evidence class. Each eligible sequence's
complete reported placement universe is passed to GeneSensus, including locations
outside the prediction set. Removing a competing gene model does not make a
peptide unique. Classes with only failed identification q-values are excluded
from consensus evidence; unknown confidence remains unknown. Standard mapping
outputs still retain their identification metadata. One-edit rescue is separate
and never enters exact consensus evidence.

New graph paths can reveal more placements. After local refinement, final
selection checks every graph candidate against the existing association search.
If additional placements change the evidence universe, selection reruns with
updated ambiguity. Successful runs stop only after this feedback stabilizes;
an eight-pass safety bound raises an error rather than publishing a completed
result with stale evidence. The number of passes and new placements are recorded.
This closes feedback for enumerated paths, not for every possible genomic isoform.

## Graph and selection rules

The frozen engine constructs bounded exon-segment/splice-junction graphs in
prediction regions. Structural paths are retained alongside peptide-guided paths;
guidance uses full observed peptide chains, not independent fragments that can
be freely recombined. Graph ORFs must translate as complete, canonical models.
Exon pairing and source-frame guards restrict new combinations. Optional RNA
provides observed exon/junction geometry rather than proof of full transcripts.

Defaults are beam 8, at most 2,000 nodes, 32 starts and 32 graph candidates per
search region. The existing `--consensus-*` options control those bounds.
Original source candidates are not pruned by the graph candidate limit.
Original partial predictions can be retained, including models crossing assembly
gaps. These are incomplete source hypotheses, not newly validated complete ORFs;
inspect `complete=False` in `candidates.tsv.gz` and `decisions.tsv`.
Search reports record bounds and skipped work; this is not exhaustive discovery.

Selection uses the engine's coding/noncoding state utilities, source agreement,
regional peptide evidence and optional RNA. The fixed noncoding coefficient is
0.75; the minimum source-candidate coding length is 150 nt and minimum intron
length is 20 nt. These are inherited defaults, not thresholds fitted to the nine
references. Missing peptides do not become a negative evidence channel.

Peptide coding coverage saturates locally; multiple overlapping sequences cannot
multiply a base's support indefinitely. A sequence with N placements contributes
1/N per placement. Split-peptide support is capped by observed coding coverage
and divided across its junctions, avoiding an intron-length bonus. Fractions
are evidence bookkeeping, not posterior probabilities.

Actual coding-conflict components are optimized exactly up to 28 candidates;
larger components use the frozen bounded solver with two greedy orders and
improving exchanges. Intronic nested models need not conflict. Coordinate-based
IDs resolve ties deterministically. The engine's fusion guard uses separated
multi-source anchors, replacing 1.2.0's broader single-source gene-name veto.

Conservative local repair can replace one retained model while preserving its
unique peptides and supported RNA geometry. Strong peptide evidence means one
unique split sequence of at least nine amino acids with at least three nt per
block, or two nonnested unique sequences of at least nine amino acids. These
are the frozen engine rules, not the stricter nine-nt block and unconditional
changed-junction/shared-phase rules introduced in GAPMS 1.2.0. Full native
selection and repair decisions remain available for inspection.

## RNA and secondary isoforms

BAM consensus evidence uses the frozen engine's defaults: MAPQ at least 20,
at least eight matched overhang bases, introns 20..500,000 nt, at least three
qualified templates for junction admission, and duplicate/secondary/supplementary
filtering. NH-aware filtering and mate/template deduplication are retained.
`--rna-library` specifies unstranded, firststrand or secondstrand. GAPMS assembly
mapping thresholds remain separate; a junction used for exploratory peptide
mapping is not automatically qualified consensus evidence.

RNA reconstruction still supplies separately reported peptide-supported ORF
candidates. These reconstructed ORFs do not become another independent gene
predictor vote in consensus.

The frozen `event_guarded` secondary-isoform layer preserves primary selection
and admits at most two supported source isoforms per primary. Its routes require
measured molecular corroboration or event-specific RNA. GAPMS does not currently
supply molecular scores, so additional isoforms require RNA; peptides alone do
not activate this secondary layer. This deliberately differs from 1.2.0's custom
peptide-only admission rule. Unselected original isoforms remain audited, and
GAPMS's annotation assessment still checks all supplied source isoforms before
proposing a novel or alternative event.

## Evaluation and limits

A matched comparison on nine crops used Helixer and Annevo with the same final
GAPMS 1.2.0 peptide BEDs, without BAM or molecular scores. Exact full-CDS F1
against RefSeq increased from 45.61% for the independent implementation to
48.22% for the frozen engine, averaged equally across crops. Annevo averaged
46.99%; the engine improved on the old implementation in all nine crops but
still fell below Annevo in melon, tomato and barley. This is a development
comparison, not blinded validation or a guarantee of improvement.

Corrected end-to-end runs were also completed for all nine crops and reproduced
the same 48.22% macro F1. Each matched the standalone engine
byte for byte on identical inference inputs. Standard, IGV and splice peptide
BEDs were byte-identical to the earlier cohort in all nine crops; selection
changed while the reported peptide placements remained the same.

Graph reconstruction can change peptide placements on other inputs. Biological output equivalence
is checked against the standalone frozen release on identical normalized
annotations and final peptide evidence. References enter only the independent
evaluator, after inference is complete. CDS/protein translations, parents,
phase fields, graph evidence and peptide projections are independently checked.
Both earlier adverse results and corrected results are retained locally.

More compatible peptides do not establish every start, stop, splice boundary
or whole transcript. RefSeq itself can contain omissions and errors. Bounded
mapping limits apparent uniqueness; no peptide absence establishes gene absence.
The final annotation is a reviewable integration of available evidence.
