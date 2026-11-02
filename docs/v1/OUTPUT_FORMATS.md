# Output formats

All TSV files have a header, including empty result tables. Lists, dictionaries
and genomic blocks within cells use JSON. Coordinates in TSV/BED are zero-based,
half-open. Protein `aa_start`/`aa_end` are also zero-based, half-open. GTF uses
one-based inclusive coordinates and preserves CDS phase. Ordering is deterministic
for a given input order, source list and configuration. Runtime, host and command
metadata naturally differ between executions.

## Run record

`run_manifest.json` records version, workflow mode, run key, parameter values,
input paths, sizes, modification times and checksums. Files up to 32 MiB get full
SHA256; larger files get a clearly labeled first/last-1-MiB sample hash, size and
mtime, not a falsely claimed full checksum. A source-code fingerprint guards
resume. External versions/paths, logged command arguments, Slurm job, host,
status and timestamps support audit. Failed runs retain their error and logs.
The source fingerprint also covers the installed legacy dependency's Python
modules. `reports/annotation_sources.tsv` links source identifiers to input files.

## Mapping

- `mapped_peptides.bed`: standard BED12, all exact loci. Name is the normalized
  original peptide. Score is 1,000 for one placement and zero for multiple.
- `mapped_peptides.igv.bed`: visualization BED12; minus-strand labels reversed.
- `splice_peptides.bed`: standard BED12 subset with more than one block.
- `mapped_peptides.tsv`: `mapping_id`, `peptide_id`, `peptide`, `chrom`, `start`,
  `end`, `strand`, `blocks`, `matched_sequence`, `match_type`, `mapping_type`,
  `evidence`, `loci_for_peptide`, `identification_quality`.
- `peptide_locus_summary.tsv`: peptide, length, locus count and unique/multiple/
  unmapped status. I/L aliases share multiplicity but retain separate labels.
- `unmapped_peptides.tsv`: peptide ID, sequence, length and
  `unmapped_in_searched_space` status. This is not a proof of absence.
- `excluded_peptides.tsv`: sequences below minimum length, length and reason.
- `peptide_identifications.tsv`: identification ID, peptide ID, normalized
  sequence, parsed data-row index and complete original table row as JSON.
  Plain lists have no invented identification metrics.
- `model_peptide_occurrences.tsv`: model and mapping IDs, peptide, protein
  coordinates, preceding residue, N-terminal context and methionine-processing
  compatibility. Mapping IDs collapse I/L-equivalent placements.
- `mapping_summary.json`: final counts, engine timings/search limits and policy.
- Optional `rescue/`: original explicitly labeled PepLocus one-edit outputs;
  they are not merged into exact evidence or proposals.

## Annotation and ORFs

`gene_model_evidence.tsv` and `orf_peptide_evidence.tsv` share these fields:
`model_id`, `source`, `original_transcript_id`, `gene_id`, `chrom`, `start`, `end`,
`strand`, `category`, `distinct_peptides`, `locus_discriminating_peptides`,
`non_nested_discriminating_peptides`, `peptides`, `mapping_ids`, `complete_start`,
`complete_stop`, `issues`, `quality`, `reasons`.

`peptides` contains canonical I/L-equivalent evidence sequences, not new
identifications. `source` disambiguates identical original transcript IDs from
different files. `gene_evidence.tsv` aggregates each source/gene without counting
a shared sequence twice. No-support and ambiguous tables are subsets;
`supported_gene_models.gtf` contains model paths with unambiguous locus support.
It does not infer isoform-specific support where isoforms share placements.

`transcripts/assembled_transcripts.gtf` is StringTie output. The accompanying
transcript FASTA is extracted from the read-only indexed assembly.
`rnaseq_junctions.tsv` contains `chrom`, `start`, `end`, `strand`, `support`,
`source`, `motif`; boundaries delimit the removed intron and support counts
alignment observations. `candidate_orfs.gtf` contains projected TransDecoder
CDS paths. `candidate_orfs.faa` contains their assembly-derived translations
keyed by stable model ID. GTF includes original transcript identity; original
TransDecoder outputs remain under transcripts for provenance.
`partial_or_ambiguous_orfs.tsv` retains incomplete/uncertain models.

## Annotation proposals

`refinement/annotation_change_proposals.tsv` is the unified proposal table.
Separate `alternative_start_candidates.tsv`, `alternative_stop_candidates.tsv`,
`alternative_splice_candidates.tsv`, and `orfs/novel_orf_candidates.tsv` use the
same schema:

`event_id`, `event_type`, `candidate_model_id`, `reference_model_ids`,
`annotation_sources`, `chrom`, `start`, `end`, `strand`, `candidate_cds_blocks`,
`reference_cds_blocks`, `relationship`, `peptides`, `mapping_ids`,
`distinct_diagnostic_peptides`, `mapping_multiplicity`, `rna_junctions`,
`minimum_rna_support`, `confidence`, `reasons`, `limitations`.

The candidate and reference blocks show the proposed coding structures.
Peptide/mapping identifiers trace back to placements and identifications.
Junction objects retain coordinates, strand and RNA observations. Reasons and
limitations are mandatory evidence interpretation, not decoration. GTF subsets
contain full candidate models suitable for review alongside annotations; these
are proposed paths, not minimal edits to apply automatically. Novel-candidate
FASTA and GTF retain single-peptide and ambiguous candidates as well as stronger
ones. Empty proposal files mean no event met the implemented diagnostic rules,
not proof that no biological change exists.

## Summary and intermediates

`reports/summary.json` gives input, mapping, model-category and event counts,
runtime and limitations. `logs/` contains the workflow log and each external
command's output. `work/` retains search junctions, normalized source models,
protein search results, indexed reference, the locus database and PepLocus
engine outputs. These support resume and audit. Do not delete them while a job
is running. Very large compressed assemblies require decompression space within
the output. Completed outputs can be archived after scientific review.
