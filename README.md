# GAP-MS: Gene Model Assessment using Peptides from Mass Spectrometry

GAP-MS evaluates predicted gene models by integrating peptide evidence from mass spectrometry with in‑silico gene predictions. It produces supported/unsupported predictions and optional reference comparisons using a reproducible, command‑line pipeline.

---

## Installation

```bash
conda create -n gapms_env python=3.10 perl-xml-parser perl-xml-twig -c conda-forge
conda activate gapms_env
git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
pip install .
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
-c, --compute_psauron      Compute PSAURON scores
-cm, --compute_embeddings  Compute embeddings-based scores
-rg, --reference_gtf       Optional reference GTF for comparison
-rf, --reference_fasta     Optional reference protein FASTA
-o, --output               Output directory
```

Notes:
- `-c` and `-cm` are mutually exclusive.
- Embeddings scoring uses the model path defined in [gapms/main.py](gapms/main.py).

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

---

## Output overview

GAP‑MS writes results to the output directory (default: `GAPMS_Output/`) and creates subfolders used by the pipeline:

```
GAPMS_Output/
├── log.txt                          # Pipeline execution log
├── all_proteins_scores.tsv          # Feature scores for all proteins
├── supported_proteins.gtf           # GTF with only supported proteins
├── supported_proteins.fa            # FASTA with only supported proteins
├── peptides_mapped.bed        # BED format for genome 
|
├── Figures/                         # Visualization plots
│   ├── roc_curve_*.png
│   ├── shap_summary.png
│   ├── protein_coverage_hist.png
│   └── ...
├── Txt/                            # Text lists of proteins
│   ├── all_proteins.txt
│   ├── supported_proteins.txt
│   ├── unsupported_proteins.txt
│   └── ...
└── Novel/                          # Novel protein analysis (if -rg/-rf used)
```

Additional files (scores, filtered predictions, peptide‑genome tracks) are produced by the pipeline modules during execution.

---

## Feedback & issues

Please open an issue on the [GitHub repository](https://github.com/qussai96/GAP-MS/issues).
