# Current architecture (1.3.0)

The pip distribution is `gapms`, with the `gapms` console command and a
`gapms-v1` compatibility alias. The package contains only `gapms_v1` and the
frozen GeneSensus engine. It depends explicitly on PepLocus 0.5.6; it neither
vendors an older mapper nor needs a GAPMS 0.3.0 wheel.

The reusable protein matcher, confidence features/seed rules and TransDecoder
coordinate projection are included directly in v1. All peptide mapping uses the
standalone PepLocus API. Annotation junctions and thresholded RNA junctions are
supplied separately; RNA alignments are scanned once. One protein automaton is
reused across source/ORF/consensus associations. All exact placements remain in
the union evidence store, with a separately counted tryptic-compatible subset.
Approximate rescue stays outside model and annotation-change evidence.

The public help lists biological inputs and threads. `--help-all` exposes
technical controls. Pip installs CPU-only XGBoost on Linux; StringTie and
TransDecoder remain explicit external requirements for RNA reconstruction.

Public report generation deduplicates identical event IDs. `outputs.refresh_report`
can rebuild completed reports without rerunning mapping or classification;
scientific tables and core provenance are preserved and publication provenance
is recorded separately.

## Historical development notes

# Architecture and baseline

The baseline is GAPMS 0.3.0, commit `40cabc4`. Its 105 tests passed before
implementation (three expected warnings for legacy projection without an
assembly). Development is isolated on `feature/gapms-1.0.0` in the GAP-MS-v1
worktree. The original `gapms/` implementation and pretrained model remain
unchanged. Main remains at the 0.3.0 release.

The new `gapms-v1` distribution contains only `gapms_v1` and depends on
`gapms==0.3.0`. It adds a separate console command. Install the local 0.3.0 wheel
first when the dependency is unavailable from an index. The legacy checkout on
main is independently installable. Do not install a different legacy dependency
version without repeating compatibility validation.

Modules have concrete responsibilities:

- `cli`: four biological input types and technical policy.
- `peptides`: validated sequences, identification metadata and quality status.
- `models`: source-aware CDS/exon structures, phase checks and translation.
- `rnaseq`: sorted-BAM validation, StringTie, read-only transcript extraction,
  TransDecoder, and the existing genomic ORF projection adapter.
- `mapping`: a disk-backed union of genomic hits and verified CDS projections.
- `assessment`: model support, diagnostic evidence and conservative events.
- `workflow`: input provenance, stage orchestration, output and guarded resume.
- `common`: stable identifiers, deterministic TSVs and safe subprocess logging.

Reusable components are bundled PepLocus's chunked six-frame engine, splice-path
mapping, junction extraction, reference indexing, peptide normalization, and
translation; the reusable Aho–Corasick protein matcher; and the mature
TransDecoder-to-genome projection adapter. Full CDS paths supplement bounded
local junction searches. Original annotations are neither merged biologically
nor rewritten: source-specific model IDs coexist in the shared evidence space.

The evidence store reconciles identical peptide placements before counting
support. I/L-equivalent observed sequences constitute one evidence unit.
Multiplicity includes every exact locus found by the genome/path search and
full model projections, rather than just proteins within an individual branch.
The legacy XGBoost decision logic is not modified or applied to novel ORFs.
