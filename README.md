> This is the isolated GAPMS 1.0.0 development branch. The new `gapms-v1`
> command and its installation are documented in [docs/v1/README.md](docs/v1/README.md).
> The legacy GAPMS workflow below remains available through `gapms`.

# GAP-MS: Gene Model Assessment using Peptides from Mass Spectrometry

GAP-MS uses mass spectrometry peptide evidence to assess predicted gene models.
It maps peptides to proteins, measures support for each model, and separates
supported from unsupported predictions. You can also supply RNA-seq alignments to
search for additional models, or a reference annotation to compare gene structures.

The main command is `gapms`. For example, to assess Helixer predictions alongside
RNA-seq evidence and compare the results with a reference:

```bash
gapms -a genome.fna -g helixer.gtf -b rnaseq.bam \
  -p peptides.tsv -rg reference.gtf -rf reference.faa -o results
```

Current version: **0.3.0**. Python **3.10** is required.

## Installation

### Conda

Create an environment with the external tools, then install GAP-MS:

```bash
conda create -n gapms_env -c conda-forge -c bioconda \
  python=3.10 pip \
  gffread gffcompare samtools stringtie transdecoder
conda activate gapms_env

git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
python -m pip install .
gapms --version
```

The Python dependencies are installed with the package. The annotation workflow
uses gffread and the bundled Python Aho–Corasick peptide mapper. BAM runs
also use samtools, StringTie and TransDecoder. Reference comparisons use gffcompare.

To update an existing checkout:

```bash
git pull
python -m pip install --upgrade .
```

### Apptainer

Build the image from the current source:

```bash
git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
apptainer build --fakeroot gapms.sif gapms.def
```

Run the same command through the image:

```bash
apptainer run --bind /path/to/data:/data gapms.sif \
  -a /data/genome.fna -g /data/helixer.gtf -b /data/rnaseq.bam \
  -p /data/peptides.tsv -rg /data/reference.gtf -rf /data/reference.faa \
  -o /data/results
```

All input and output paths must be accessible inside the container. If you need a
writable sandbox, create one from the SIF:

```bash
apptainer build --sandbox gapms_sandbox gapms.sif
```

Updating a Conda installation does not update an existing image. Check the version
with `apptainer run gapms.sif --version` and rebuild when needed.
[Container validation details](docs/container_validation.md).

## Running GAP-MS

### Assess predicted models

Supply a genome assembly, prediction annotation and peptide list:

```bash
gapms -a genome.fna -g predictions.gtf -p peptides.tsv -o results
```

GAP-MS translates the predictions with gffread, maps peptides to those proteins,
extracts evidence features, and applies the bundled XGBoost classifier. It also
maps the full peptide list to the assembly and writes genomic BED tracks.

### Add RNA-seq evidence

```bash
gapms -a genome.fna -g predictions.gtf -b rnaseq.bam \
  -p peptides.tsv -o results
```

This runs two branches:

- `prediction_search` assesses the supplied predictions.
- `bam_search` assembles transcripts with StringTie, predicts ORFs with
  TransDecoder, and assesses those models using the same peptide evidence.

The branches are compared, and supported models are merged. You can also run the
BAM branch alone:

```bash
gapms -a genome.fna -b rnaseq.bam -p peptides.tsv -o results
```

### Compare with a reference

Add a reference annotation and its protein FASTA to either workflow:

```bash
gapms -a genome.fna -g predictions.gtf -b rnaseq.bam \
  -p peptides.tsv -rg reference.gtf -rf reference.faa -o results
```

The reports include candidate novel models and peptide-supported differences in
start, stop and splice sites. Reference-event analysis considers all
peptide-mapped models, including models not retained by the classifier.

### Use existing proteins or mappings

If you already have the translated prediction proteins, supply them with `-f`.
A precomputed peptide-to-protein mapping can be supplied with `-m`:

```bash
gapms -a genome.fna -g predictions.gtf -f predictions.faa \
  -p peptides.tsv -m mapping.tsv -o results
```

Protein identifiers must agree between the FASTA, annotation and mapping file.
The built-in mapper reports one-based positions. GAP-MS validates supplied positions against
the protein sequence and also accepts validated zero-based positions.

A prediction run can omit the assembly if proteins are supplied with `-f`, but
that uses the older protein-coordinate BED export. Include `-a` to use direct
assembly mapping and obtain both standard and IGV BED files.

### Protein peptide search

GAP-MS searches the branch's translated proteins with a C-backed Aho–Corasick
matcher (`pyahocorasick`). It builds the peptide automaton once and reuses it for
the prediction and BAM-derived protein sets. Isoleucine and leucine are always
equivalent. Every overlapping occurrence is retained, with no mismatch rescue in
this protein search; PepLocus handles genomic matching separately.

Each branch writes `protein_mapping.tsv` (`peptide`, `protein`, `location`) and
`protein_mapping.summary.json` with counts and timings. Locations are one-based;
unmapped queries have protein `UNMAPPED` and location `0`. Protein identifiers
must be unique in each FASTA. Recognized table headers and modification labels
are handled by the peptide input parser.

ProteoMapper and its Perl scripts are no longer required. Existing mappings can
still be supplied with `-m`. Across nine species and both branches, the new mapper
reproduced all tested peptide–protein locations. See the
[comparison and timing details](docs/protein_mapper_evaluation.md).

### Inputs and options

Peptides may be supplied as a plain sequence list or a peptide table, including
FragPipe's `combined_peptide.peptides.tsv`. Annotation examples use `.gtf`, but
supported GFF attribute formats are accepted too. Sequence names must match the
assembly and BAM.

| Option | Input or setting |
|---|---|
| `-a`, `--assembly` | Genome FASTA; required for BAM runs and direct genomic mapping |
| `-g`, `--gtf` | Prediction annotation |
| `-b`, `--bam` | RNA-seq alignments |
| `-p`, `--peptides` | Peptide input; required |
| `-f`, `--proteins` | Existing prediction protein FASTA |
| `-m`, `--mapping` | Existing peptide-to-protein mapping |
| `-rg`, `--reference_gtf` | Reference annotation |
| `-rf`, `--reference_fasta` | Reference protein FASTA |
| `-s`, `--scores` | External scores CSV: `Protein` and `external_score` columns |
| `-i`, `--iterative` | Train an iterative classifier instead of using the bundled model |
| `-t`, `--threads` | Genomic mapping workers; default **8** |
| `-o`, `--output` | Output directory |

Use a fresh output directory for each run. If `-o` is omitted, GAP-MS creates
`GAPMS_Output` two directory levels above the supplied GTF, or above the BAM when
no GTF is supplied. Set `-o` explicitly to choose where results go.

`-t` controls the parallel genomic mapping stages, not every external tool or
feature calculation. Run `gapms --help` for the full workflow options.

## Results

For a prediction-only or BAM-only run, branch results are written directly under
the output directory. When both `-g` and `-b` are supplied, they are separated:

```text
results/
├── log.txt
├── prediction_search/
│   ├── all_proteins_scores.tsv
│   ├── supported_proteins.gtf
│   ├── supported_proteins.fa
│   ├── peptides_mapped.bed
│   ├── peptides_mapped.igv.bed
│   ├── peplocus/
│   ├── Figures/
│   ├── Txt/
│   └── Compare_to_Reference/
├── bam_search/                    # Same result files, plus transcript/ORF intermediates
├── comparisons/                   # Comparison between the two branches
└── combined_gapms_summary.png
```

Start with these files:

| File | What it contains |
|---|---|
| `protein_mapping.tsv` | Peptide-to-protein matches with one-based positions |
| `protein_mapping.summary.json` | Protein-search counts and timings |
| `all_proteins_scores.tsv` | Protein evidence features and scores |
| `supported_proteins.gtf` / `.fa` | Models retained by the classifier |
| `Txt/` | Supported, unsupported, high-confidence and low-confidence protein lists |
| `Figures/` | Evidence and classification plots |
| `peptides_mapped.bed` | Standard genomic BED12 for downstream analysis |
| `peptides_mapped.igv.bed` | Matching track with IGV display labels and colors |
| `peplocus/mappings.tsv` | Genomic placements, translations, evidence and multiplicity |
| `peplocus/unmapped.tsv` | Peptides without a reported genomic placement |
| `peplocus/summary.json` | Mapping settings, counts, runtime and search limits |

`Unmpped_pepides.txt` retains its historical filename and lists peptides unmapped
to the branch's **proteins**. It is different from `peplocus/unmapped.tsv`, which
reports the result of **genomic** mapping.

With reference comparison enabled, `Compare_to_Reference/` contains:

- `annotation_comparison_report.tsv` and `annotation_comparison_summary.tsv`.
- `Novel/`: candidate novel models and their evidence.
- `Different/`: peptide-supported start, stop and splice differences.

In the feature table, `splice_peptides` counts distinct peptide sequences;
`supported_splice_sites` counts distinct junctions with unambiguous peptide support.
A repeated peptide is not credited to a junction unless its occurrences support
the same junction.

## Map peptides without running model assessment

Use `gapms map` when you only need genomic peptide locations. This command uses
the bundled **PepLocus** mapper and does not run protein classification, StringTie
or TransDecoder.

```bash
# Genome and peptides only
gapms map -a genome.fna -p peptides.tsv -o mapped

# Supply splice junctions through an annotation
gapms map -a genome.fna -p peptides.tsv -g predictions.gtf -o mapped_gtf

# Or use RNA-seq alignments; annotation and BAM can also be combined
gapms map -a genome.fna -p peptides.tsv -b rnaseq.bam -o mapped_bam
```

Only assembly and peptides are required. The default output directory is
`GAPMS_Map`, and the default worker count is eight; change it with `-t`.
GTF/GFF2/GFF3 attribute variants, CDS-only annotations and gzip annotations are
supported. Gene-only files cannot define transcript junctions.

The main GAP-MS workflow uses this same mapper for its genomic BED step when an
assembly is supplied. It searches the full peptide list, including peptides not
found in predicted proteins. Protein classification still uses the protein-based
evidence; an extra genome match is not automatically assigned to a gene model.

### Matching and splice peptides

PepLocus searches all six assembly reading frames and always treats I/L as
equivalent. GTF/GFF or BAM input supplies junctions for splice-crossing peptides,
including split codons and microexons. Without either input, ordinary mapping finds
contiguous matches; optional novel-junction discovery is available through
`gapms map --help`.

All compatible placements are retained. A genomic match establishes sequence
compatibility; it does not validate the spectrum or prove which locus was expressed.
Splice paths are bounded by the configured intron/path limits, and multi-junction
paths can combine edges not observed together in one transcript.

### Standard and IGV BED files

Both files describe the same loci with zero-based, half-open genomic coordinates:

- **`peptides_mapped.bed`** keeps N-to-C peptide labels on both strands. Use this
  file for downstream analysis and sequence comparisons.
- **`peptides_mapped.igv.bed`** reverses negative-strand labels for the historical
  GAP-MS display convention. Contiguous mappings are gray and spliced mappings red.

A BED score of 1,000 means the peptide has one compatible placement; 0 means it has
multiple placements. The score is not a probability. Detailed tables and the
original engine BEDs are saved under `peplocus/`.

For reference reports, protein IDs and terminal-position metadata are attached
separately using cached CDS indexes. They are saved under `protein_projection/`;
the standard BED stays BED12. Extra genomic placements remain part of ambiguity
checks. [Metadata method and validation](docs/model_metadata.md).

### Optional one-edit rescue

To search still-unmapped peptides for one amino-acid substitution, insertion or
deletion:

```bash
gapms map -a genome.fna -p peptides.tsv -g predictions.gtf \
  --rescue-one-edit --rescue-min-length 12 -o mapped_rescue
```

Exact matching runs first. Rescue only searches peptides with no exact placement,
and saves its results separately under `peplocus/rescue_*`. I/L differences do not
consume an edit. The default minimum of 12 amino acids can be edited.

## Questions and bug reports

Please [open an issue](https://github.com/qussai96/GAP-MS/issues) 

## Citation

If you use GAP-MS in your work, please cite the [GAP-MS preprint on bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.17.712294v3).


### Candidate interpretation fixes in 1.0.1

A C-terminal peptide supports translated sequence beyond an annotation, but it does not establish a stop codon. `alternative_stop` now requires a complete in-frame stop in the assembly-derived ORF. Otherwise the event is `c_terminal_extension`, with termination explicitly unresolved. Events also record whether the splice chain changes; a terminal-exon change should not be interpreted as an isolated termination mechanism.

The full RNA models in candidate GTFs are context for manual review, not replacement annotations. Event tables now give the directly observed `local_peptide_blocks` and `diagnostic_junctions`, a run-scoped `coding_structure_id`, and `whole_model_status` with explicit review reasons. Partial models, ambiguous diagnostic evidence, missing RNA junction support, multiple overlapping annotation genes, fewer than two independent diagnostic peptides, and loss of compared annotation coding sequence keep a model in `local_evidence_only`. Even `complete_model_for_manual_review` is not a validated whole-CDS correction. The refinement directory contains separate tables for these statuses and a BED6 track of individual observed peptide blocks. These blocks are not a reconstructed CDS. Duplicate RNA IDs can represent the same coding structure; use distinct CDS paths when counting candidates, and always join IDs with the run or crop.

Peptide uniqueness is relative to the searched genome and splice paths. An unsearched splice path can contain another exact location. If a second annotation is available, add its coding paths without making it a target of annotation assessment:

```bash
gapms-v1 -a assembly.fna -p peptides.tsv -g helixer.gtf -b reads.bam \
  --search-gtf refseq.gtf -o results -t 8
```

`--search-gtf` can be repeated and accepts the same GTF/GFF syntax as `-g`. Its junctions expand the mapping graph, and its full coding paths participate in the shared, genome-verified protein projection search. These placements contribute to multiplicity before candidate assessment. Search-only models do not suppress Helixer-relative proposals and do not count as annotations being assessed. Reused transcript IDs on different contigs or strands are namespaced for this mapping-only input. Omitting this option preserves the original search inputs, with the same explicit scope limitation. Adding it changes search coverage and ambiguity, so comparisons must disclose the additional annotation. Reference-assisted ambiguity checks are not an independent validation of novelty or peptide identification.

Orphan exon records without transcript parents in mapping-only annotations are skipped and listed in `reports/skipped_search_annotation_features.tsv`; they cannot define a coding path. Orphan CDS records and malformed assessment-target models remain errors.
