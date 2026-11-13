# GeneSensus 0.6.0a3 — evidence-guarded coding isoforms

This CPU prerelease integrates the selected **event-guarded GS-Isoform v0.1**
policy into GeneSensus. Primary graph reconstruction and selection retain the
0.5.0a3 method. A final, additive stage can retain supported alternative coding
transcripts from the input predictions, with shared gene parents in the GFF3.
The prior 0.5.0a3 release and ongoing cohort predictions remain unchanged.

The isoform policy is enabled by default in this version. It **abstains without
Scores or RNA evidence**. It does not promise higher accuracy for every species.
The 18-crop/Arabidopsis transfer campaign is still running; preliminary results
below cover only Phaseolus and carrot.

## Installation and running

```bash
python -m pip install '/home/students/q.abbas/Tools/EVM_modified/Genesensus/releases/0.6.0a3/dist/genesensus-0.6.0a3-py3-none-any.whl'
genesensus --version
genesensus run config.json -o results -t 4 --cache cache
```

Python >=3.9 on Linux/POSIX, NumPy and Biopython are required. Install `genesensus[molecular]`
for HDF5 residue profiles and `genesensus[rna]` for BAM input when installing
from this package directory; equivalently install `h5py` and `pysam`. Graph
reconstruction and isoform admission run on CPU. CUDA bulk scoring remains
experimental. ProtAudit/ProtT5 inference is external preprocessing and is not
silently run by GeneSensus.

A complete two-predictor Scores + RNA example:

```json
{
  "genome": "genome.fa",
  "predictions": {
    "braker3": {
      "annotation": "braker3.gtf",
      "protein_fasta": "braker3.faa",
      "scores": "braker3.scores.tsv",
      "residue_profiles": "braker3.profiles.h5",
      "diagnostics": "braker3.diagnostics.tsv"
    },
    "galba": {
      "annotation": "galba.gtf",
      "protein_fasta": "galba.faa",
      "scores": "galba.scores.tsv",
      "residue_profiles": "galba.profiles.h5",
      "diagnostics": "galba.diagnostics.tsv"
    }
  },
  "rna": {"bam": "sorted_RNA.bam", "library": "unstranded"},
  "parameters": {
    "selection_policy": "state_molecular_guided",
    "molecular_strategy": "splice_orf",
    "peptide_strategy": "graph_contrast",
    "confidence_gate": 0.437721461,
    "localizer_gate": 0.437721461,
    "isoform_policy": "event_guarded",
    "isoform_score_gate": 0.437721461,
    "compress_diagnostics": true
  }
}
```

Remove optional paths/sections for unavailable evidence. For a prediction-only
run, retain the genome and named prediction annotations, remove molecular
paths and RNA, and use `molecular_strategy: "source_context"`. One or multiple
predictors are supported. Existing peptide BED12 and protein-alignment inputs
retain their 0.5.0a3 primary reconstruction behavior. RNA can also be supplied
as the validated GeneSensus junction/geometry JSON instead of BAM.

The example thresholds belong to the frozen manuscript ProtAudit classifier.
Use the producer's declared positive threshold consistently for your own score
model, including `isoform_score_gate`; do not optimize it on target-reference
agreement. Values must be finite probabilities. Predictor names are arbitrary:
source **gene/transcript relationships**, not names such as BRAKER, establish
multi-isoform families. Retain correct `gene_id`/`transcript_id` in GTF or
`ID`/`Parent` in GFF3. Malformed AUGUSTUS parent rows and missing annotated stop
codons need validated normalization before running; the prerelease does not
silently guess those source repairs.

## What happens, in order

1. Validate inputs, contigs, supplied protein correspondence, CDS phases and
   measured molecular profiles using the inherited pipeline and caches.
2. Construct and score exon/junction graph candidates; select primary models
   with the existing structural, molecular, RNA, peptide and protein logic.
3. Save those primary transcripts unchanged in `primary_annotation.gff3`.
4. Examine **full source-predicted alternatives** from gene families containing
   at least two distinct coding structures in at least one predictor. An exact
   coding structure found in another predictor can corroborate that alternative.
   UTR-only differences do not create additional coding-isoform candidates.
5. Apply the physical, event-evidence and admission rules below. Append admitted
   alternatives to the main annotation and protein FASTA. Group each alternative
   with its primary gene, extending the gene span to contain its transcripts.

This stage does not invent new combinations of source isoforms or replace an
incorrect primary transcript. Primary corrections still happen in step 2.
Candidate graph paths are not eligible as new secondary isoforms unless their
complete coding structure is present among the eligible source predictions.
The GFF3 writer serializes exons as the union of supplied exons and coding
intervals so every CDS base (including an explicitly annotated terminal stop)
is contained by its exon, transcript and gene spans. Coding coordinates, phase
and protein translation are unchanged by this serialization correction.

Stable transcript IDs use each model's existing candidate hash; secondary IDs
may have a different hash from the gene parent they share. Every transcript has
one valid parent, and its own exon/CDS children retain that transcript ID.

## Fixed structural safeguards

Every admitted alternative must:

- have a complete ORF and canonical splice sites, CDS length >=150 nt and each
  intron >=20 nt;
- overlap exactly one selected primary on the same chromosome and strand;
- share >=50% of each model's CDS, with **every shared coding base in frame**;
- share an exact coding exon or junction with the primary;
- have a different intron chain: same-chain terminal-only variants are excluded.

These 150/20 nt isoform floors are fixed independently of the inherited primary
graph thresholds. The stage admits at most **two** additional isoforms per primary. Duplicate
alternative intron chains are skipped. Among qualifying alternatives the fixed
ranking is: RNA event support, weakest alternative-junction template count,
measured q, number of distinct exact-structure predictor sources, shared CDS
fraction, then a stable candidate identifier. Unavailable evidence contributes
nothing to ranking. Primary CDS geometries and reading frames are immutable.

## Evidence requirements and each value's role

An alternative needs at least one of these routes:

**Scores route:** its minimum measured original ProtAudit q across identical
coding records reaches `isoform_score_gate`, and the **complete exact coding
structure** occurs in at least two distinct predictors. Several isoforms from
one predictor never create multiple votes. A missing q is unobserved, not zero.
One predictor + Scores alone cannot pass this two-source route.

**RNA route:** every newly introduced junction has >=3 template observations,
and >=95% of newly coding bases are covered by qualified RNA alignment blocks.
If the alternative only retains primary introns, each retained intron must be
spanned by an observed internal candidate exon with >=3 whole-exon observations.
Internal-exon observations are alignment-geometry counts; they are not a claim
of independent molecules or full-transcript phasing. When there are no newly
coding bases, no extra coding-base coverage assertion is needed. Shared exons
alone do not support a new splice event. One multi-isoform predictor + suitable
RNA evidence can pass this route.

**Residue-error probabilities:** remain active in the inherited `splice_orf`
primary selection, along with q, splice/ORF context and available diagnostics.
The selected secondary-isoform policy **does not apply the experimental
`risk_guarded` veto**, including its additional SpuriousORF veto. An alternative
can pass the RNA event route even if its q is negative. A high-q protein with no localization profile has
unobserved residue risk, not a measured risk of zero. This release does not
score reconstructed proteins or claim new whole-protein q for them.

**Peptides and protein alignments:** retain their inherited primary-model
roles, including bounded peptide ambiguity/translation support. They are not
new secondary-isoform admission routes in 0.6.0a3. A peptide-only run therefore
gets the inherited peptide-guided primaries and no extra isoforms from this
layer. Absence of peptides is never a new gene-removal rule.

| Available evidence | New secondary isoforms? |
|---|---|
| Predictions only | No; exact abstention |
| Peptides and/or protein alignments only | No admission route in this release |
| Scores, one predictor | No two-predictor molecular route |
| Scores, >=2 predictors | Possible with corroborated positive full structures |
| RNA, one or multiple predictors | Possible with isoform-event evidence |
| RNA + Scores | Either route; fixed evidence-aware ranking |

Disable this layer with `"isoform_policy": "off"` to retain the previous primary
selection/output behavior. Both `off` and `event_guarded` are valid; policy
thresholds and source evidence must be chosen before reference evaluation.

## Outputs and audit trail

- `annotation.gff3`: final primary + admitted secondary transcripts, grouped by
  primary gene; coordinates are ordinary 1-based GFF3.
- `proteins.faa`: one translated protein per final transcript; identifiers match
  the GFF3 transcript IDs.
- `primary_annotation.gff3`: unchanged pre-isoform selection (enabled mode).
- `additional_isoforms.gff3`: admitted secondary structures alone; use the main
  annotation for shared-parent locus grouping.
- `isoforms.json`: summary, accepted/rejected event decisions and parent mapping.
  Residue-veto fields from the experimental rules are not an active admission
  veto; this integration does not compute a new residue audit for alternatives.
- `report.json`: `primary_selected_count`, final `selected_count`, `isoforms`
  summary, stage timings and inherited evidence diagnostics.
- `configuration.json`, `provenance.json`, candidate/decision diagnostics:
  resolved configuration, input/output/code checksums and decisions. Resume
  checks include the version, code and settings; use a fresh output for the new
  release. Do not precreate the output directory.

The generic input formats, existing score logic, caches and advanced selection
modes are documented in [the package guide](../../README.md) and
[the inherited molecular guide](MOLECULAR_STRATEGY.md). The shared-parent
isoform grouping also applies to full/consensus annotations when the associated
primary is present. A subset lacking that primary keeps a valid standalone gene
parent for its retained alternative. Provisional-only files remain standalone.

## Measured results so far

Transcript changes versus **matched 0.5.0a3 BRAKER3 + GALBA controls**, in
percentage points. These are not comparisons against a different predictor set.

| Evidence | Phaseolus Δ precision | Δ recall | Δ F1 | Carrot Δ precision | Δ recall | Δ F1 |
|---|---:|---:|---:|---:|---:|---:|
| None | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Scores | +0.191 | +0.556 | +0.338 | +0.114 | +0.412 | +0.229 |
| RNA | +0.304 | +1.028 | +0.485 | +0.272 | +0.872 | +0.457 |
| Scores + RNA | +0.497 | +1.500 | +0.898 | +0.354 | +1.279 | +0.706 |

Phaseolus is the development plant; carrot is the first completed transfer
assessment. The fixed policy admitted 0/237/442/650 and 0/212/423/649 alternatives,
respectively. Remaining 18-crop and Arabidopsis transfer results are pending.
Arabidopsis uses explicitly labeled BRAKER2 2.1.6 AUGUSTUS-hints + GALBA inputs.
No claim of universal precision/F1 improvement or 18-crop superiority is made.

These accuracy measurements came from the frozen experimental policy. The
installed 0.6.0a3 implementation subsequently reproduced **all four Phaseolus
final CDS structures and frames exactly**, with every original primary retained,
valid shared-parent GFF3 spans, and matching GFF3/protein FASTA transcript IDs.
The installed wheel passed **246 tests**. Completed records are in
[validation](validation/). Individual two-thread CPU timings, including input validation and output
writing, are retained in the validation records. External ProtAudit embedding/score
generation is excluded. These are not EVM speed comparisons. Evaluation uses strict CDS-projected gffcompare
(`--strict-match -e 0 --no-merge`); transcript/intron-chain percentages are
recomputed from integer counts. References are evaluation-only. The user
requested preliminary evaluation of already frozen species while other jobs
continue; neither policy nor pending predictions may be tuned from those results.

A source-predicted, score-positive transcript is a protein-plausible hypothesis,
not proof of expression. BRAKER and GALBA share AUGUSTUS ancestry and are not
independent observations. Short RNA reads support local events, not complete
isoform phasing. Incomplete reference catalogs can classify real alternatives
as unmatched. Conservative guards miss some genuine isoforms; additive
preservation cannot remove bad primaries. Frozen external ProtAudit models
retain their documented pretraining overlaps; this is not sequence-disjoint
validation. See [the prior evidence-scope audit](../0.5.0a3/validation/PRETRAINED_EVIDENCE_SCOPE.json).

## Reproducibility and release scope

The precursor is [the fixed transfer experiment](../../experiments/isoform_transfer_v1/README.md),
with separate results under its `results/` symlink. Its already-submitted jobs
continue with their frozen experimental code; they do not silently switch to
this package. The implementation generalizes family detection to arbitrary
source names; published preliminary accuracy currently covers BRAKER + GALBA.
The default numerical guards match that policy; `isoform_score_gate` explicitly
exposes the positive threshold for other predeclared score models.

The previous release is preserved at [0.5.0a3](../0.5.0a3/README.md). This is an
alpha CPU release adopting the user-selected policy while broader transfer
validation remains in progress, not a replacement claim for finished cohort
results or a new reference-trained model.
