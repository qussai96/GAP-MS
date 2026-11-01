# GAP-MS: Gene Model Assessment using Peptides from Mass Spectrometry

GAP-MS evaluates predicted gene models by integrating peptide evidence from mass spectrometry with in‑silico gene predictions. It produces supported/unsupported predictions and optional reference comparisons using a reproducible, command‑line pipeline.

---

## Direct peptide-to-assembly mapping (PepLocus)

GAPMS 0.2.2 bundles the PepLocus 0.3.1 engine. Map observed peptides directly to a
six-frame assembly translation and optional splice-junction evidence:

```bash
# Annotation-supported mapping; GTF/GFF2/GFF3 syntax is detected from attributes.
gapms map -a assembly.fa -p peptides.txt -g predictions.gtf -t 8 -o mapped_gtf

# RNA-seq-supported mapping; no StringTie or TransDecoder step is needed here.
gapms map -a assembly.fa -p peptides.txt -b rnaseq.bam -t 8 -o mapped_bam

# Both sources can be combined; FragPipe's "Peptide Sequence" column is recognized.
gapms map -a assembly.fa -p combined_peptide.peptides.tsv \
  -g predictions.gff3 -b rnaseq.bam -t 8 -o mapped_both

# Optional exactly-one-amino-acid-edit rescue of still-unmapped queries >=12 aa.
gapms map -a assembly.fa -p peptides.txt -g predictions.gtf \
  --rescue-one-edit --rescue-min-length 12 -t 8 -o mapped_rescue
```

Only assembly and peptides are required for `gapms map`. The default worker count is **8**, overridable with `-t`/`--threads`. The default output is
`GAPMS_Map`; use a fresh output directory for each run. `--assembly`/`--genome`,
`--peptides`, `--gtf`/`--gff`, `--bam`, `--threads`, and `--output`/`--out` are accepted
long options. `gapms map --help` lists the complete PepLocus search options,
including reusable junction TSV input and optional bounded novel-junction search.
Mapping requires the installed NumPy, Biopython, pyahocorasick, and pysam libraries;
it does not call a protein mapper or external alignment executable.

Every assembly mapping produces:

```text
mapped_gtf/
├── peptides_mapped.bed           # Standard BED12 for downstream analysis
├── peptides_mapped.igv.bed       # Display BED for IGV
└── peplocus/
    ├── mappings.bed             # Same records as the standard BED above
    ├── mappings.igv_display.bed  # Same records as the IGV BED above
    ├── mappings.tsv             # Translations, evidence, multiplicity, coordinates
    ├── peptides.tsv
    ├── unmapped.tsv
    ├── junctions.tsv            # When splice evidence is supplied
    └── summary.json             # Engine version, inputs, settings, timing, limits
```

The standard BED uses biological N-to-C peptide names on both strands and exactly
12 columns. The IGV BED reverses negative-strand names for historical GAPMS display
compatibility and uses gray/red for contiguous/spliced mappings. Coordinates and
loci are identical; use the standard file for sequence comparisons. Both use
zero-based, half-open coordinates and ascending genomic blocks. A BED score of
1,000 means one compatible placement, 0 means multiple; it is not a confidence
probability. All placements are retained. I/L equivalence is always enabled.

Rescue outputs remain separate under `peplocus/rescue_*`; they never add approximate
alternatives for an exact-mapped query. Twelve amino acids is the empirical rescue
starting point from the carrot controls, not an FDR threshold. Genome-only matching
is exhaustive for contiguous exact matches; splice-path and novel-discovery searches
have documented bounds. Translation verification is enabled by default. Recheck:

```bash
python -m gapms._peplocus verify --genome assembly.fa \
  --mappings mapped_gtf/peplocus/mappings.tsv
```

### Integration with the annotation-assessment workflow

Existing `gapms -g ... -a ... -p ...` workflows default to eight mapping workers
and use PepLocus for the genomic
BED step, searching the full observed peptide list rather than only peptides that
mapped to predicted proteins. `-t`/`--threads` selects mapping workers. Protein
mapping and classifier features remain protein-based; extra genomic placements
are not automatically assigned to a protein or included as classifier features.
The prediction and BAM branches each write both BED versions when assembly is
supplied.

Reference-comparison reports also need protein/terminal metadata. These are
attached using cached CDS interval indexes and saved separately in
`protein_projection/peptides_mapped.bed` when reference
comparison is requested. Metadata is attached only where the protein projection
matches a PepLocus locus under I/L equivalence. All additional genomic placements
remain in the report's ambiguity checks; inconsistent protein projections cannot
become genomic event evidence. The attachment preserves the previous metadata BED exactly, while avoiding
repeated full-protein nucleotide projection. See
[metadata validation](docs/model_metadata.md).

Historical protein-only workflows with `-f` but no assembly retain the previous
protein-coordinate exporter, with a warning. Supply `-a` to obtain the new genomic
search and both BED versions. Existing containers are unchanged; rebuild them to
include GAPMS 0.2.2 and the new mapping dependencies.

See [integration and validation details](docs/peplocus_integration.md). The bundled
engine's [MIT license](gapms/_peplocus/LICENSE) is included in built distributions.


For output-preserving feature-extraction speedups and the carrot benchmark, see
[feature optimization](docs/feature_optimization.md).

## Installation

### Option 1: Conda Installation

```bash
conda create -n gapms_env -c conda-forge -c bioconda \
  python=3.10 \
  perl-xml-parser perl-xml-twig \
  gffread gffcompare samtools stringtie transdecoder
conda activate gapms_env
git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
pip install .
```

For full annotation-assessment BAM runs (`gapms -b`, rather than `gapms map -b`), GAP-MS now expects these external tools to be available:
- `samtools`
- `StringTie`
- `gffread`
- `TransDecoder`
- `gffcompare`

### Option 2: Using Containers

**With Apptainer**:
```bash
apptainer pull oras://docker.io/qussaiab96/gapms:latest
apptainer run gapms_latest.sif -g annotations.gtf -f proteins.fasta -p peptides.txt
```

**With Docker/Podman**:
```bash
docker pull qussaiab96/gapms:latest
docker run --rm -v $(pwd):/data qussaiab96/gapms:latest -g /data/annotations.gtf -f /data/proteins.fasta -p /data/peptides.txt
```

**Build Locally from Apptainer Definition**:
```bash
git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
apptainer build --fakeroot gapms.sif gapms.def
apptainer run gapms.sif -g annotations.gtf -f proteins.fasta -p peptides.txt
```

---

## Required inputs

1) Peptides list (`-p`) [required]  
2) At least one search source:
   - Prediction GTF (`-g`), or
   - RNA-seq BAM (`-b`)
3) Genome assembly FASTA (`-a`) when:
   - proteins must be generated from a GTF, or
   - BAM mode is used
4) Optional:
   - prediction protein FASTA (`-f`)
   - precomputed peptide mapping (`-m`)
   - external scores CSV (`-s`)
   - reference GTF / FASTA (`-rg`, `-rf`)

---

## Usage

### 1) With both prediction GTF and BAM input

```bash
gapms -g predictions.gtf -b rnaseq_alignments.bam -a assembly.fasta -p peptides.txt [-rg reference_proteins.gtf -rf reference_proteins.faa]
```
### 2) With prediction GTF + protein FASTA

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt
```

### 3) With prediction GTF + genome assembly (proteins will be generated)

```bash
gapms -g annotations.gtf -a assembly.fasta -p peptides.txt
```

### 4) With BAM input only

```bash
gapms -b rnaseq_alignments.bam -a assembly.fasta -p peptides.txt
```

If `-o` is not provided, outputs go to `GAPMS_Output/` in the parent directory of the GTF or BAM file.

---

## Command‑line options

```
-g, --gtf                 Path to the prediction GTF file (optional unless --bam is omitted)
-b, --bam                 Optional RNA-seq BAM file for StringTie → TransDecoder search
-p, --peptides            Path to peptides TXT file (required)
-f, --proteins            Path to prediction protein FASTA (optional if -a is provided)
-a, --assembly            Path to genome assembly (required if -f is not provided or when --bam is used)
-m, --mapping             Optional peptide-to-protein mapping file
-s, --scores              Optional CSV with columns: Protein, external_score
-rg, --reference_gtf      Optional reference GTF for comparison
-rf, --reference_fasta    Optional reference protein FASTA
-o, --output              Output directory
-i, --iterative           Train an iterative XGBoost model instead of using the pre-trained classifier
```

---

## Examples

Use a precomputed mapping file:

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt -m mapping.tsv
```

Compare with reference annotations:

```bash
gapms -g predictions.gtf -a genome.fasta -p peptides.txt -rg reference.gtf -rf reference.fasta
```

Use iterative training instead of the bundled pre-trained XGBoost model:

```bash
gapms -g predictions.gtf -f proteins.fasta -p peptides.txt -i
```

---

## Output overview

GAP‑MS writes results to the output directory (default: `GAPMS_Output/`). The structure depends on whether you run a single branch or both prediction + BAM branches together.

### A) Single-branch output (`-g` only or `-b` only)

```
GAPMS_Output/
├── log.txt
├── all_proteins_scores.tsv          # Feature/evidence table for all proteins
├── supported_proteins.gtf           # Supported proteins in genome coordinates
├── supported_proteins.fa
├── peptide_mapped_proteins.gtf      # all mapped models used for event discovery
├── peptides_mapped.bed
├── Unmpped_pepides.txt
├── Figures/
│   ├── Peptides_mapping_bars.png
│   ├── roc_curve_*.png
│   ├── shap_summary_plot.png
│   └── ...
├── Txt/
│   ├── all_proteins.txt
│   ├── supported_proteins.txt
│   ├── unsupported_proteins.txt
│   └── ...
└── Compare_to_Reference/            # Produced when -rg/-rf is used
    ├── gffcmp.*
    ├── annotation_comparison_report.tsv
    ├── annotation_comparison_summary.tsv
    ├── Novel/
    │   ├── new_predicted_genes.tsv
    │   ├── new_predicted_proteins.gtf
    │   ├── new_predicted_proteins.fa
    │   └── new_predicted_proteins_scores.tsv
    └── Different/
        ├── peptide_support_different_start.gtf
        ├── peptide_support_different_stop.gtf
        └── peptide_support_different_splice.gtf
```

`supported_proteins.gtf/.fa` remain the classifier-supported protein set.
Reference-event discovery is intentionally broader: the strict novel-gene and
alternative start/stop/splice reports evaluate every peptide-mapped model, so a
valid event is not discarded by the protein classifier. In
`all_proteins_scores.tsv`, `splice_peptides` counts distinct peptide sequences,
whereas `supported_splice_sites` counts distinct annotated junctions with
unambiguous peptide support.

The bundled ProteoMapper output uses one-based peptide locations. GAP-MS
validates those positions against the protein sequence and retains compatibility
with validated zero-based precomputed mappings; missing positions are expanded
to every matching sequence occurrence.

### B) Combined output when both `-g` and `-b` are supplied

```
GAPMS_Output/
├── log.txt
├── combined_gapms_summary.png       # Parent-level summary figure for presentations
├── prediction_search/
│   ├── all_proteins_scores.tsv
│   ├── supported_proteins.gtf
│   ├── supported_proteins.fa
│   ├── peptides_mapped.bed
│   ├── Figures/
│   ├── Txt/
│   └── Compare_to_Reference/
├── bam_search/
│   ├── transcripts.gtf                      # StringTie transcript models
│   ├── transcripts.fa                       # transcript FASTA
│   ├── transcripts.fa.transdecoder.pep      # TransDecoder peptides
│   ├── all_proteins_scores.tsv
│   ├── supported_proteins.gtf               # supported BAM ORFs in genome
│   ├── supported_proteins.fa
│   ├── peptides_mapped.bed
│   ├── Figures/
│   ├── Txt/
│   └── Compare_to_Reference/
└── comparisons/
    ├── bam_vs_input_gtf_summary.tsv
    ├── bam_supported_no_gene_overlap.tsv
    ├── bam_supported_no_gene_overlap.gtf
    ├── bam_supported_no_gene_overlap.faa
    └── bam_vs_input_gtf_summary.png
```
---

## Feedback & issues

Please open an issue on the [GitHub repository](https://github.com/qussai96/GAP-MS/issues).
