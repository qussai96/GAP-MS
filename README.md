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

Current version: **0.2.2**. Python **3.10** is required.

## Installation

### Conda

Create an environment with the external tools, then install GAP-MS:

```bash
conda create -n gapms_env -c conda-forge -c bioconda \
  python=3.10 pip perl-xml-parser perl-xml-twig \
  gffread gffcompare samtools stringtie transdecoder
conda activate gapms_env

git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
python -m pip install .
gapms --version
```

The Python dependencies are installed with the package. The annotation workflow
uses gffread and ProteoMapper; the ProteoMapper Perl scripts are bundled. BAM runs
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
ProteoMapper positions are one-based. GAP-MS validates supplied positions against
the protein sequence and also accepts validated zero-based positions.

A prediction run can omit the assembly if proteins are supplied with `-f`, but
that uses the older protein-coordinate BED export. Include `-a` to use direct
assembly mapping and obtain both standard and IGV BED files.

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
consume an edit. The default minimum of 12 amino acids comes from the carrot
controls; it is a starting point, not an FDR threshold. Shorter rescue queries can
produce many incidental matches.

Mapping coordinates and translations are verified by default. To check an existing
mapping table again:

```bash
python -m gapms._peplocus verify --genome genome.fna \
  --mappings mapped/peplocus/mappings.tsv
```

Run `gapms map --help` for length filters, reusable junction input, search bounds
and the other mapping options.

## Tests and validation

To run the source tests:

```bash
python -m pip install pytest
python -m pytest -q
```

The current suite has 92 passing tests. Real-data checks and implementation details
are documented separately:

- [PepLocus integration and carrot comparison](docs/peplocus_integration.md)
- [Feature-extraction speedups with unchanged outputs](docs/feature_optimization.md)
- [Fast metadata attachment and reference-report equivalence](docs/model_metadata.md)
- [SIF and sandbox validation](docs/container_validation.md)

## Questions and bug reports

Please [open an issue](https://github.com/qussai96/GAP-MS/issues) with the GAP-MS
version, command, relevant log section and a small example if possible.

## Citation

If you use GAP-MS in your work, please cite the [GAP-MS preprint on bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.17.712294v3).
