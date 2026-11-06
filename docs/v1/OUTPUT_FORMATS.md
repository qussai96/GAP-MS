# Public result layout (1.2.2)

The main directory contains `report.html`, `report.log`, `peptides.bed`,
`peptides.igv.bed` and `unmapped_peptides.tsv`. Annotation assessment also produces
`supported_genes.gtf` and `all_proteins_scores.tsv`; multiple annotations add
`consensus.gff3` and `consensus_proteins.faa`. Only nonempty proposal categories
appear in `review/`. See [the user guide](README.md#finding-your-results).

All native subdirectory paths described below are now **relative to
`diagnostics/`**. For example, `mapping/mapped_peptides.tsv` is
`diagnostics/mapping/mapped_peptides.tsv`. The workflow log used by users is
`report.log` in the output root. Internal workflow logs remain diagnostic.
Temporary `work/` is removed on success unless `--keep-work` is supplied.
The frozen consensus algorithm and native scientific schemas are unchanged.

# Output formats

All TSV files have a header, including empty result tables. Lists, dictionaries
and genomic blocks within cells use JSON. Coordinates in TSV/BED are zero-based,
half-open. Protein `aa_start`/`aa_end` are also zero-based, half-open. GTF uses
one-based inclusive coordinates and preserves CDS phase. Ordering is deterministic
for a given input order, source list and configuration. Runtime, host and command
metadata naturally differ between executions.

IDs are deterministic in their reference/source context. Peptide IDs identify
normalized sequences; locus IDs include the I/L-equivalent sequence and genomic
path and are scoped to the run's assembly. Preserve assembly and run identities
when joining outputs from different references. Model IDs retain source identity.

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

`gene_model_evidence.tsv`, `orf_peptide_evidence.tsv` and `all_predicted_orf_evidence.tsv` share these fields:
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
alignment observations. `all_predicted_orfs.gtf` contains all projected TransDecoder CDS paths;
`all_predicted_orfs.faa` and `all_predicted_orf_evidence.tsv` retain their
translations and complete assessment for audit. These include unsupported ORFs.

`candidate_orfs.gtf`, `candidate_orfs.faa` and `orf_peptide_evidence.tsv`
contain only RNA-derived ORFs with at least one exact/I-L-equivalent peptide
placement. Shared and ambiguous placements are retained and labeled, so this
is not a uniqueness or identification-quality filter. Known annotated isoforms
can also be peptide-supported candidates; novelty is assessed separately.
`orfs_without_peptide_support.tsv` lists predictions without observed support.
`partial_or_ambiguous_orfs.tsv` contains only peptide-supported partial or
uncertain candidates. A header-only TSV and empty GTF/FASTA are valid when no
ORF has peptide support. Original TransDecoder files remain under transcripts.

In summary counts, `predicted_orfs` counts the complete RNA-derived set;
`candidate_orfs` counts the peptide-supported subset, and
`orfs_without_peptide_support` counts the remainder. Distinct CDS-path counts
are provided separately for predictions and supported candidates. Historical
versions through 1.0.3 used `candidate_orfs` filenames/counts for all predictions.

## Annotation proposals

`orfs/candidate_classification.tsv` contains one row per proposed candidate:
`candidate_model_id`, `category`, `component_changes`, `event_ids`. Primary
categories are mutually exclusive; `combined_refinement` records multiple
supported changes without counting the model as a novel ORF. Component event
tables can contain multiple rows for the same candidate or comparisons with
multiple annotation models. Use `candidate_categories` in the summary for
model counts and `annotation_events` for event counts.


`refinement/annotation_change_proposals.tsv` is the unified proposal table.
Separate `alternative_start_candidates.tsv`, `alternative_stop_candidates.tsv`,
`alternative_splice_candidates.tsv`, `c_terminal_extension_candidates.tsv`, and `orfs/novel_orf_candidates.tsv` use the
same schema:

`event_id`, `event_type`, `candidate_model_id`, `reference_model_ids`,
`annotation_sources`, `chrom`, `start`, `end`, `strand`, `candidate_cds_blocks`,
`reference_cds_blocks`, `relationship`, `peptides`, `mapping_ids`,
`distinct_diagnostic_peptides`, `mapping_multiplicity`, `rna_junctions`,
`minimum_rna_support`, `confidence`, `reasons`, `limitations`,
`coding_structure_id`, `local_peptide_blocks`, `diagnostic_junctions`,
`splice_chain_changed`, `whole_model_status`, `whole_model_review_reasons`,
`uniqueness_scope`, `evidence_origin`, `boundary_status`, `cleavage_context`,
`upstream_peptide_support`, `proposed_feature_blocks`.

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


## Local-evidence review outputs in 1.0.1

`refinement/local_evidence_only.tsv` and `complete_models_for_manual_review.tsv` partition the integrated event table by whole-model review status. Neither is an approved correction list. `local_peptide_evidence.bed` is BED6 with one row per observed peptide block per event, named `event_id.block_index`. It does not connect blocks into a reconstructed CDS. Event JSON fields retain diagnostic junctions and mapping IDs for that purpose. The directory README explains the full-context GTF policy.

`orfs/annotation_isoform_checks.tsv` contains `candidate_model_id`, `chrom`, `strand`, `coding_structure_id`, `exact_annotated_isoform_ids`, `comparable_annotated_isoform_ids`, `known_start_isoform_ids`, `known_end_isoform_ids`, `new_junctions`, `observed_peptide_sequences`, `diagnostic_peptide_sequences`, `interpretation`. It audits all peptide-supported RNA candidates, including exact annotated isoforms that correctly produce no event. Gene-level `gene_evidence.tsv` now includes `chrom` and `strand` to distinguish reused gene names across loci.


## Annotation-guided local proposals (1.1)

`refinement/annotation_guided_events.tsv` reports local events inferred from the supplied annotations and peptide evidence, with or without BAM. These rows also appear in the unified `annotation_change_proposals.tsv` and event-type candidate tables. `evidence_origin` distinguishes `annotation_guided` from `RNA_ORF`. RNA-model comparisons and local events can describe the same feature; their summed row counts are not counts of distinct genes or biological events.

For annotation-guided rows, `candidate_model_id`, `coding_structure_id` and `candidate_cds_blocks` are empty because no full transcript is inferred. `local_peptide_blocks` and `mapping_ids` trace the exact observed placements; `proposed_feature_blocks` records an inferred start/stop codon when justified. `boundary_status` distinguishes unresolved extension boundaries, genomic stop inference, ATG start candidates and sequence-verified junctions. `cleavage_context` contains reference-specific flanking amino acids and tryptic-context flags. `upstream_peptide_support` is descriptive for internal starts, not a filter requiring absent upstream translation.

`annotation_guided_features.gtf` contains only local start/stop codons or peptide-supported regions. Its transcript IDs are event grouping IDs, not reconstructed transcript models. These features are for review, not replacement annotations. `local_peptide_evidence.bed` is a BED6 view of the observed blocks; the standard mapping BED12 retains each real peptide path for visualization.

Full RNA-derived model contexts, when available, are now named `<event_type>_rna_model_context.gtf`. They are distinct from local feature GTFs; annotation-only runs do not create these full-model context files. This replaces the ambiguous refinement `<event_type>_candidates.gtf` names from earlier versions. Complete novel ORF reconstruction still requires BAM, and peptide-supported candidate ORF exports retain the v1.0.4 filtering rules.

`reports/summary.json` and `mapping/mapping_summary.json` record the annotation-guided search bounds, cap counts, path counts and runtime. Assessment summaries give `annotation_guided_event_rows` and `rna_event_rows` separately. Changing search bounds can change mapping coverage and apparent uniqueness.

## Automatic consensus (1.2.1)

Multiple annotation inputs produce `consensus/annotation.gff3`, `proteins.faa`,
`primary_annotation.gff3`, `additional_isoforms.gff3`, `summary.json`, and the
native GeneSensus reports: `candidates.tsv.gz`, `decisions.jsonl.gz`,
`graph_report.json`, `isoforms.json`, `report.json`, `configuration.json` and
`provenance.json`. The annotation can retain original exons/UTRs where supplied;
CDS structures and proteins carry coding evidence. GeneSensus coordinate-derived
`GS.g...` identifiers are preserved.

`source_model_mapping.tsv` links each original GAPMS source model to its engine
candidate and selection status. `decisions.tsv` is a compact candidate/status
index, not the full explanation; the native decision JSONL contains evidence,
exclusions, optimization and repair events. The former custom `graph_report.tsv`,
`selection_report.tsv` and `peptide_repairs.tsv` are superseded by native reports.

Unresolved excluded originals are written by the engine to `provisional.gff3`
when present, with an alias `unresolved_alternatives.gff3`. Other unselected
competitors remain in native candidates and decisions; not all competitors
are duplicated in the unresolved GFF3. Neither is part of the final annotation.
Inputs normalized for the engine, including an identity/phase audit and the
full eligible peptide placement universe, remain in `work/consensus_inputs/`.
See [CONSENSUS.md](CONSENSUS.md) for the policy and its limits.
