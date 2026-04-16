# GAP-MS: Gene Model Assessment using Peptides from Mass Spectrometry

GAP-MS evaluates predicted gene models by integrating peptide evidence from mass spectrometry with in‑silico gene predictions. It produces supported/unsupported predictions and optional reference comparisons using a reproducible, command‑line pipeline.

---

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

For BAM-enabled runs (`-b`), GAP-MS now expects these external tools to be available:
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
gapms -g predictions.gtf -b rnaseq_alignments.bam -a assembly.fasta -p peptides.txt [-rf reference_proteins.gtf -rg reference_proteins.gff]
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
-c, --compute_psauron     Compute PSAURON scores
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

Compute PSAURON scores:

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt -c
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
    │   ├── new_predicted_proteins.gtf
    │   ├── new_predicted_proteins.fa
    │   └── new_predicted_proteins_scores.tsv
    └── Different/
        ├── peptide_support_different_start.gtf
        ├── peptide_support_different_stop.gtf
        └── peptide_support_different_splice.gtf
```

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
