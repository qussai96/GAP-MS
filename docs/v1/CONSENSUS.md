# Automatic consensus in GAPMS 1.2.0

With two or more distinct `-g` annotation files, GAPMS always constructs a
consensus. There is no enable/disable consensus option. Assembly + peptides
alone and single-annotation modes retain the 1.1.0 behavior. BAM is optional.

```bash
gapms-v1 -a genome.fa -p peptides.tsv -g helixer.gtf -g annevo.gtf -o results -t 8
# Optional RNA evidence; declare library protocol when known:
gapms-v1 -a genome.fa -p peptides.tsv -g helixer.gtf -g annevo.gtf \
  -b reads.bam --rna-library firststrand -o results_rna -t 8
```

## Evidence flow

All supplied annotations are inference inputs. There is no hidden reference or
auxiliary search annotation. References used to evaluate a result must be kept
outside inference. Sources are ordered by resolved path for repeatable naming;
reordering the same set of `-g` arguments preserves annotation output. IDs are
scoped to sources. Input CDS phases, assembly translations and existing model
issues are inherited from the existing annotation parser. Coding-only source
structures are collapsed by chromosome, strand and normalized coding blocks,
retaining each contributing model and predictor. Distinct sources count once
per coding base or source junction; repeated isoforms cannot multiply votes.
If identical structures have both problematic and valid source records, the
valid record is preferred. Noncoding/exon-only models remain in source
assessment but do not become coding consensus candidates.

The usual genome/splice peptide mapping and source-protein association run
first. The consensus graph then generates additional bounded candidate
structures; their proteins are searched with the **same existing peptide
ProteinMatcher automaton** used for source proteins. This is not a claim that
the genome mapper, local refinement and protein search share one automaton:
they have separate search roles. All projected new paths are genome verified
and entered into the shared mapping database. Local refinement and final
selection use the updated ambiguity counts. No mapping is discarded because
its location lies outside the supplied predictions, and no peptide becomes
unique merely because one supplied model explains it. Exact I/L-equivalent
matches are eligible; exploratory one-edit rescue never enters consensus.

## Bounded graph reconstruction

Overlapping coding-span components are constructed separately per contig and
strand. These are comparison regions, not assertions that all models are one
gene. Input exon boundaries subdivide coding coverage into segments. Supplied
junctions connect segments in transcript direction; when BAM reconstruction
is used, peptide-supported RNA ORF models can contribute their geometry too.
The graph does not jump to arbitrary distant exon pairs. Each traversed coding
segment and splice edge must have a codon offset supplied by at least one
valid source context. Structurally invalid sources remain audited but do not
seed graph reconstruction or structural source votes.
Completed graph combinations must have an ATG start, an in-frame terminal stop,
no internal stop or unknown coding bases, and GT–AG, GC–AG or AT–AC splice sites
with introns at least 20 nt. Existing partial source models can remain partial;
their unobserved ends are not silently extended.

Two searches run: structural and uniquely placed peptide-guided, retaining their
union. Peptide hints require the complete observed block chain and codon offset;
separate junction fragments cannot freely combine as independent observations.
The structural beam score is source coding nt/1000 plus 0.2 per source-supported
splice edge. Complete unique-peptide hints add observed coding nt/1000 plus 0.2
per crossed junction, once per sequence. It is a search heuristic.

Defaults per region: beam 8 per start/phase, at most 2000 nodes, 32 starts and
32 new candidates. `--consensus-beam`, `--consensus-max-nodes`,
`--consensus-max-starts`, `--consensus-max-candidates` expose these bounds.
Original models are not removed by graph candidate pruning. Skips, candidate
limits and paths pruned are written to `graph_report.tsv`. Bounded search is
not exhaustive isoform discovery.

## Model evidence and selection

For a canonical peptide sequence with N distinct genomic block paths, each
placement contributes 1/N. This is an ambiguity allocation, **not a posterior
probability**. I/L aliases, repeated spectra and duplicate paths are not extra
independent sequences. Coding evidence integrates weighted peptide-covered
bases, saturated at one per base. A peptide crossing J junctions divides its
splice allocation by J. Each intron's peptide credit is capped by the observed
weighted coding bases and its length, avoiding a full-intron bonus for a short
junction peptide. Identification q-value failures do not supply selection
utility or hard support. Unknown identification quality remains unknown.

Selection uses the coding-versus-noncoding state comparison documented in
GeneSensus. Phase-compatible coding bases and exact phase-compatible introns
receive one vote per predictor. A source intron vote is in genomic nucleotide
units; this is distinct from the **bounded peptide** splice credit above.

```
positive = source_coding_nt + source_splice_nt
         + weighted_peptide_coding_nt + bounded_peptide_splice_nt
         + qualified_RNA_coding_nt + qualified_RNA_splice_nt
utility = (positive - 0.75 * whole_span_intergenic_nt) / 1000
```

Whole-span intergenic opportunity cost is the candidate span not covered by
any gene span of each supplied predictor, summed per predictor. This is a
structural prior from competing predictions, not experimentally observed
noncoding sequence. Admission additionally requires positive support to reach
0.75 times the sum of source-absent coding bases, coding bases overlapping
source introns, opposite-strand intronic overlap and exact splice support.
The fixed 0.75 coefficient follows the documented GeneSensus/EVM convention;
it is not fitted to the crop references. Utilities are not probabilities.
Neutral valid source hypotheses retain a small positive value. No new minimum
length is imposed on existing source genes. Missing peptide observations add
no penalty. Structurally disfavored sources remain unresolved, not assertions
that the gene does not exist. Physical sequence/phase defects remain audited.

RNA coding coverage is the union of qualified matched bases. A junction with
N>=3 qualified templates contributes intron_length * N/(N+3), a bounded evidence
utility rather than a transcript probability. With no BAM, both RNA terms are
zero. A peptide can support only its observed translated bases and junction
context; sparse support is never inflated to full-gene coverage.

Same-strand coding overlap defines primary conflicts. An exact maximum-utility
compatible set is found for up to 18 candidates per region; larger regions use
a deterministic greedy bound and record that limitation. Equal preferences
are resolved by stable coordinate-derived IDs, while losing hypotheses remain
in the unresolved audit; a tie is not claimed as molecular resolution. Nested
or overlapping genes sharing coding bases can remain unresolved by this rule.
A candidate spanning two coding-disjoint, differently identified genes from
one source is vetoed to guard against silent fusions.

Every graph exon must retain a supplied first/internal/last/single-exon role
and codon offset (including peptide-supported reconstructed RNA models), or
be fully observed as an internal block of a unique peptide chain. This prevents
freely pairing unsupported exon halves. New graph combinations need strong local translation (two nonnested unique
sequences >=9 aa, or one unique split sequence >=9 aa with >=9 nt in each block)
or qualified RNA-junction evidence. They must also distinguish a splice change
from an overlapping source, keep all overlapping source unique peptides and
qualified source RNA junctions, and satisfy the fusion/ORF guards. Other graph
paths remain provisional/rejected, never emitted just because enumerated.

After structural optimization, a conservative one-for-one splice repair may
prefer a complete compatible candidate explaining a new unique split peptide
or two nonnested unique sequences. Evidence must lie inside the incumbent
locus, and the candidate must retain every incumbent unique sequence and
qualified RNA junction. Shared coding bases must retain phase. At least one
splice edge changes; this rule does not force a terminal extension. Changes
are recorded in `peptide_repairs.tsv`. A strict gain of unique sequence support
prevents repair cycles; two retained genes cannot be silently fused.

## Alternatives and RNA

At most two additional full source-predicted coding isoforms are emitted per
primary. Each must overlap exactly one primary, share at least 50% of each
CDS with every shared base in frame, share an exact exon or junction, have a
different intron chain, and be a complete canonical ORF. Terminal-only variants
and unphased new graph combinations are not emitted as secondary isoforms.
They remain available for review. Alternative transcript parents share one
consensus gene, whose span is extended to contain all its transcripts.

An alternative needs either a unique >=9-aa peptide that crosses a junction
distinguishing it from the primary, or RNA support for the changed event. The
current peptide admission route intentionally abstains on exon-internal-only
support: this does not independently establish the exon boundaries. Source
prediction provides full structure but is not experimental connectivity proof.

Consensus RNA evidence uses MAPQ >=max(20, `--min-mapq`), no secondary,
supplementary, duplicate, QC-failed or NH>1 alignments; MAPQ255 requires NH=1.
Junction overhangs require `--min-read-anchor` (default 8 nt), and introns must
be 20..500000 nt. Templates are deduplicated by read group and query name at
each junction. Library defaults to unstranded; `--rna-library` interprets
firststrand/secondstrand for the consensus channel. This option does not change
the legacy StringTie reconstruction call or local-refinement RNA policy.
Every new junction needs at least three templates, and at least 95% of newly
coding bases need qualified matched-base coverage. A retained primary intron
needs at least three templates with one matched block spanning the entire
intron plus 8 nt on each flank. Coverage comes from matches, not deletions or
skipped regions. If an input BAM lacks an index, a CSI is created inside the
run's work directory; input files are never indexed in place. Reads are fetched
only over supplied candidate regions. RNA can support these local events
without proving full-transcript phasing or translation initiation.

## Files

- `consensus/annotation.gff3`: selected primaries and supported source alternatives.
- `primary_annotation.gff3`, `additional_isoforms.gff3`: their separate subsets.
- `proteins.faa`: exact assembly translations of emitted CDSs; IDs match GFF3.
- `unresolved_alternatives.gff3`: competing hypotheses, not the final annotation.
- `decisions.tsv`: every coding candidate, utility components, support, issues,
  status, distinguishing observations, boundary status and competing primaries.
- `source_model_mapping.tsv`: original source transcript to consensus candidate.
- `graph_report.tsv`, `selection_report.tsv`, `peptide_repairs.tsv`: search and optimization limits.
- `summary.json`: counts, evidence policy and measured stage times.

The GFF3 is coding-only: it does not invent UTRs, transcript ends or transcription
start sites. Phases are recomputed from the normalized coding path. Start/stop
flags mean genomic ORF compatibility, not experimentally confirmed boundaries.
Local start/stop/splice review proposals and the standard/IGV peptide BEDs are
still exported by their existing stages. A feature already explained by an
input isoform is not relabeled as novel. Unsupported RNA ORFs are not added to
the consensus solely because they were predicted.

## Relation to GeneSensus

This implementation adopts the bounded graph, conserved peptide budget and
conservative evidence principles described in GeneSensus 0.6.0a3. It is an
independent GAPMS policy, not a promise of identical GeneSensus results or a
port of its ProtAudit/state-utility optimizer. No molecular embeddings, protein
scores or target-reference training are introduced. GeneSensus 0.6.0a3 itself
does not admit peptide-only secondary isoforms; the unique junction-peptide
route above is GAPMS's explicit extension. Its thresholds are conservative
heuristics requiring evaluation, not protein-level FDR estimates.

## Development evaluation scope

The initial positive-only selector (e888085) passed sequence integrity checks
but had lower exact-CDS F1 than either predictor on the nine-crop development
cohort. That unfavorable result is retained locally. The released policy adds
the documented structural-state competition and conservative repair after
that diagnostic. Each implementation is frozen before its own evaluation;
these are retrospective development comparisons, not blind transfer tests.
No crop-specific thresholds, target-reference features or numerical weight
search are used. Cohort performance must be reported against both predictors
and a matched structural-only consensus control, including losses.
