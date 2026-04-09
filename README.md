# GAP-MS: Gene Model Assessment using Peptides from Mass Spectrometry

GAP-MS evaluates predicted gene models by integrating peptide evidence from mass spectrometry with in‑silico gene predictions. It produces supported/unsupported predictions and optional reference comparisons using a reproducible, command‑line pipeline.

---

## Installation

### Option 1: Conda Installation

```bash
conda create -n gapms_env python=3.10 perl-xml-parser perl-xml-twig -c conda-forge
conda activate gapms_env
git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
pip install .
```

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

1) Gene models in GTF format (`-g`)  
2) Peptides list (one peptide per line) (`-p`)  
3) Either:
   - Protein FASTA for the predictions (`-f`), or
   - Genome assembly FASTA (`-a`) to generate proteins from the GTF

---

## Usage

### With protein FASTA

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt
```

### With genome assembly (proteins will be generated)

```bash
gapms -g annotations.gtf -a assembly.fasta -p peptides.txt
```

If `-o` is not provided, outputs go to `GAPMS_Output/` in the parent directory of the GTF file.

---

## Command‑line options

```
-g, --gtf                 Path to the prediction GTF file (required)
-p, --peptides            Path to peptides TXT file (required)
-f, --proteins            Path to protein FASTA (optional if -a is provided)
-a, --assembly            Path to genome assembly (required if -f is not provided)
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

GAP‑MS writes results to the output directory (default: `GAPMS_Output/`) and creates subfolders used by the pipeline:

```
GAPMS_Output/
├── log.txt                          # Pipeline execution log
├── all_proteins_scores.tsv          # Feature/evidence table for all proteins
├── supported_proteins.gtf           # GTF with only supported proteins
├── supported_proteins.fa            # FASTA with only supported proteins
├── peptides_mapped.bed              # BED track of peptide genome coordinates
├── Unmpped_pepides.txt              # Peptides reported by Proteomapper as UNMAPPED
├── Figures/                         # Visualization plots
│   ├── roc_curve_*.png
│   ├── shap_summary.png
│   ├── sequence_coverage_hist.png
│   └── ...
├── Txt/                             # Text lists of proteins
│   ├── all_proteins.txt
│   ├── high_confident_proteins.txt
│   ├── supported_proteins.txt
│   ├── unsupported_proteins.txt
│   └── ...
└── Compare_to_Reference/            # Produced when -rg/-rf is used
    ├── gffcmp.*                     # gffcompare outputs (.tmap/.refmap/.stats/.tracking/...)
    ├── annotation_comparison_report.tsv
    ├── annotation_comparison_summary.tsv
    ├── Novel/
    │   ├── new_predicted_proteins.gtf
    │   ├── new_predicted_proteins.fa
    │   └── new_predicted_proteins_scores.tsv
    ├── Different/
    │   ├── peptide_support_different_start.gtf
    │   ├── peptide_support_different_stop.gtf
    │   └── peptide_support_different_splice.gtf
    └── ...

```

---

## Feedback & issues

Please open an issue on the [GitHub repository](https://github.com/qussai96/GAP-MS/issues).
