# GAPMS 1.3.0

GAPMS v1 maps mass-spectrometry peptide sequences to eukaryotic assemblies and
organizes evidence for reviewing gene annotations. Genomic peptide locations are
the common basis for annotation support, RNA-derived ORFs and proposed changes.
The software reports candidates and uncertainty; sequence compatibility alone
does not validate a peptide identification or a biological annotation change.

Version 1.3.0 is a self-contained pip distribution named `gapms`. Both `gapms`
and `gapms-v1` run this workflow. PepLocus 0.5.6 is an explicit dependency;
the old bundled mapper and GAPMS 0.3.0 wheel are not runtime dependencies.

## Installation

Python 3.10 or newer is required. Install the release wheel or this checkout:

```bash
python -m pip install gapms-1.3.0-py3-none-any.whl
# From the source checkout:
python -m pip install .
gapms --version
```

If you installed the earlier `gapms-v1` development distribution, uninstall it
first with `python -m pip uninstall gapms-v1`. It pins legacy GAPMS and shares the
Python namespace with this release. A fresh environment does not need this step.

BAM reconstruction additionally needs StringTie and TransDecoder on PATH:

```bash
conda install -c conda-forge -c bioconda stringtie transdecoder=5.5.0
```

Genome and annotation modes need no external mapping executable. The C-backed Python dependencies are installed automatically with PepLocus. V1 uses its own Python transcript
extraction and creates FASTA indexes inside the output, never beside inputs.

## Automatic multi-annotation consensus

Version 1.2.1 always produces a consensus coding annotation when multiple
`-g` files are supplied, guided by genome-verified peptides and optional RNA
evidence. Single-annotation runs retain assessment and local refinement.
See [the consensus method and limitations](CONSENSUS.md) for the graph,
selection policy, ambiguity budget, outputs and RNA requirements.
There is no consensus enable/disable switch. Since 1.2.1, this uses the frozen
GeneSensus 0.6.0a3 engine rather than the independent 1.2.0 approximation.
The engine is included in the Python wheel; a separate GeneSensus install is not needed.

## Run

```bash
# Genome-wide peptide coordinates
python -m gapms_v1 -a genome.fna -p peptides.tsv -o genome_results

# Add model assessment and annotated junction peptides
gapms-v1 -a genome.fna -p peptides.tsv -g helixer.gtf -o annotation_results

# Reconstruct RNA-supported ORFs
gapms-v1 -a genome.fna -p peptides.tsv -b sorted.bam -o rna_results

# Automatic consensus when two or more -g files are supplied
gapms-v1 -a genome.fna -p peptides.tsv -b sorted.bam \
  -g helixer.gtf -g another.gff3 -t 8 -o integrated_results
```

`gapms --help` shows common inputs; `gapms --help-all` shows technical controls.

Assembly and peptides are required. Annotation and coordinate-sorted BAM are
optional. GTF/GFF3 transcript relationships, compressed annotation/assembly,
and plain peptide lists, FASTA, CSV or TSV are supported. Table peptide columns
are detected by common names or set with `--peptide-column`. I/L equivalence is
always enabled. The default minimum peptide length is 7; shorter queries are
labeled in `unmapped_peptides.tsv` and recorded in
`diagnostics/mapping/excluded_peptides.tsv`. This threshold limits ambiguity,
not peptide identification error. Use a lower threshold explicitly for small
peptide research. No uniqueness claim can cover unsearched splice structures.

Optional `--rescue-one-edit` searches engine-stage exact-unmapped peptides of at least 9 aa
by default. Rescue outputs are separate and do not contribute to exact model
support or annotation proposals. All exact placements are retained.

A new output directory is required. `--resume` only accepts unchanged inputs,
parameters and implementation. An interrupted run can reuse completed
reconstruction and mapping stages; completed runs are returned unchanged.
Do not delete stage outputs selectively and expect resume to repair them.

## Results

Open `report.html` first. The main directory contains the peptide BEDs,
`unmapped_peptides.tsv`, the model assessment table and supported genes.
Multiple annotations also produce `consensus.gff3` and `consensus_proteins.faa`.
Nonempty proposal categories are under `review/`; technical detail is under
`diagnostics/`. See [Finding your results](#finding-your-results) for the file list.
Use `peptides.igv.bed` for display and `peptides.bed` for downstream analysis.

Candidates need direct diagnostic evidence. Missing downstream peptides never
justify a shorter coding sequence. A peptide elsewhere in a protein never
supports its changed splice junction. Non-tryptic termini alone do not establish
translation initiation. Unknown identification quality remains unknown.

The nine-crop Casanovo evaluation completed successfully; see the
[results, reference audit and figures](evaluation/EVALUATION.md).

See [methods](SCIENTIFIC_METHODS.md), [output schemas](OUTPUT_FORMATS.md),
[architecture](ARCHITECTURE.md), and [validation](VALIDATION.md).

## Cite GAPMS

[GAPMS bioRxiv manuscript](https://www.biorxiv.org/content/10.64898/2026.03.17.712294v3).

The [RefSeq follow-up evaluation](evaluation/refseq_20261102/RESULTS.md) includes the nine-crop rerun, a 900-gene annotation-removal benchmark, and ambiguity boxplots for every observed peptide length. Read its [protocol and repeatability notes](evaluation/refseq_20261102/README.md) when using these results.

In the diagnostic RNA outputs, `orfs/candidate_orfs.gtf` and `.faa` contain only ORFs with peptide support. The complete RNA prediction set is named `all_predicted_orfs.gtf`/`.faa` for audit; `orfs_without_peptide_support.tsv` lists unsupported predictions. Shared or ambiguous peptide support remains explicitly labeled. Peptide support alone does not establish novelty or validate a complete model.


With v1.1, `-a assembly.fna -p peptides.tsv -g models.gtf -o results` also produces annotation-guided local start, stop, extension and splice proposals without a BAM. Internal ATG-starting peptides with non-tryptic N-terminal context are considered even when their sequence is already contained in an annotated longer protein. Existing isoforms are checked before claiming a new feature. Uncertain boundaries remain unresolved; peptide-supported local changes do not validate a replacement transcript. Optional BAM junction evidence can strengthen a local proposal. See [technical methods](SCIENTIFIC_METHODS.md) for the default 3,000-nt terminal window, 30-nt splice-site window, single-exon skips and search caps. Full new ORF reconstruction continues to require BAM.


## Finding your results

Start with `report.html`. It links the result files and explains what the evidence
supports. `report.log` is beside it and records progress, warnings and errors.

With one annotation, the main files are:

```text
results/
  report.html
  report.log
  supported_genes.gtf
  all_proteins_scores.tsv
  peptides.bed
  peptides.igv.bed
  peptides.tryptic.bed
  unmapped_peptides.tsv
  review/                  # only when proposals exist
  diagnostics/
```

The input annotation is not duplicated. Its path and SHA-256 are recorded in the
report. `supported_genes.gtf` preserves the supplied transcript structures with
locus-discriminating peptide support; shared or uncertain support remains in the
score table. Support does not establish every exon or protein boundary.

With multiple annotations, `consensus.gff3` and `consensus_proteins.faa` are also
produced. They are newly generated results. Review proposals are kept separate
and are not automatically applied to the consensus.

`review/` contains `candidates.tsv` and nonempty `alternative_start.gff3`,
`alternative_stop.gff3`, `alternative_splice.gff3` and `novel_orfs.gff3` as relevant.
Start/stop files also contain explicitly labeled terminal-extension proposals.
Local evidence stays local; unresolved boundaries are not exported as invented
complete transcripts. Novel RNA-derived ORFs remain labeled candidate models.

`all_proteins_scores.tsv` uses the legacy 29-column feature and confidence format.
It describes the supplied annotation for a single GTF, or the final consensus for
multiple GTFs. Source-specific assessments, classifier scores, seed labels and
RNA support remain in `diagnostics/classification/`. The report embeds a simple
two-panel summary; exportable SVG/PNG files are in `diagnostics/reports/`.

`unmapped_peptides.tsv` distinguishes `unmapped_in_searched_space`, `rescue_only`
and `not_searched_below_min_length`. Rescue-only placements remain exploratory;
they are not added to the exact peptide BED or used for annotation decisions.

Normalized consensus annotation/peptide inputs are retained in
`diagnostics/consensus/inputs/`, with a relocation record for native provenance.
Detailed tables, native GeneSensus decisions, input provenance and reconstruction
records are under `diagnostics/`. The diagnostic subdirectory names retain their
previous schemas. Temporary `diagnostics/work/` is removed after successful
publication; use `--keep-work` to retain it for debugging or mapping-database work.
Failed runs retain work and diagnostics. `--resume` still requires identical
inputs, parameters and implementation. Existing 1.2.1 result folders are not
migrated or overwritten.

### Supported gene models (1.2.3)

Each annotation source is classified separately. Multiple-input runs also classify the final consensus, and `supported_genes.gtf` contains supported consensus transcripts. Single-input runs export supported supplied transcripts. RNA-only runs classify reconstructed ORFs. Source structures are preserved; classifier labels do not alter GeneSensus selection or peptide placement.

High-confidence seeds satisfy at least one legacy rule: two protein-specific peptides; at least 80% protein coverage; both termini plus a gene-specific peptide; or every coding splice junction supported plus a gene-specific peptide. Low-confidence seeds have fewer than two peptides, no gene-specific peptide, and length outside 200–500 aa. High confidence takes precedence if rules overlap. Specificity is assessed within each annotation pool, not equated with genomic placement uniqueness. I/L-equivalent peptides count once; identification-quality failures are excluded. Repeated peptide occurrences must all support a terminal or splice interpretation to receive that feature.

The legacy base and engineered features feed a 20-configuration randomized XGBoost search, with stratified cross-validation (up to five folds), random seed 42, and CPU use bounded by `-t`. Up to five self-training rounds label remaining models at scores ≥0.90 or ≤0.10 and refit. Final support includes all original high-confidence seeds and initially unlabeled models with final score ≥0.50; low-confidence seeds remain excluded. Fewer than two seeds in either class skips training: high seeds remain supported, other non-seeds unresolved.

`diagnostics/classification/scores.tsv` records `seed_label`, `iterative_label` (-1 means not recruited), `classifier_score`, `confidence`, `supported`, and the feature values. `diagnostics/classification/` retains per-pool scores, trained model JSON, seed counts, search parameters and recruitment counts. These scores and seed cross-validation measure agreement with heuristic labels, not calibrated biological probabilities or independent gene-model accuracy. Lack of peptide support does not establish that a gene is false. Novel/alternative proposals retain their separate molecular evidence requirements.

### Score table and report figure

`all_proteins_scores.tsv` retains the 29-column legacy GAPMS schema and order, including `supported` (+/−) and `confidence` (high/moderate/no-support). A multiple-GTF run scores the final consensus here; source annotations and classifier probabilities remain in `diagnostics/classification/`. Stable protein IDs join the exported supported GTF. `isoforms` counts coding isoforms per source/gene; transcript scores come from supplied transcript/mRNA records (0 when absent), and external scores are 0 because v1 has no external-score input.

The offline HTML embeds a two-panel figure, also available as SVG and PNG under `diagnostics/reports/`. Panel A shows mapped input peptide sequences and separate percentages of N-terminal, C-terminal, splice and internal peptide–protein pairs in the assessed annotation. Termini may overlap; pairs are used so shared peptides are not silently treated as isoform-specific evidence. Panel B counts distinct proposal IDs by event type, separating unresolved terminal extensions from proposed alternative boundaries. Proposals are not distinct genes or confirmed annotation corrections. Counts and denominators are recorded in `overview_counts.json`.
