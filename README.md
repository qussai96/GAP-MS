# GAP-MS: **Gene Model Assessment using Peptides from Mass Spectrometry**

**GAP-MS** is a pipeline designed to evaluate and refine predicted gene models by integrating proteogenomic evidence from mass spectrometry with in-silico gene predictions. This approach leverages an iterative machine learning framework to classify predicted proteins as *supported* or *unsupported*, improving both specificity and overall accuracy of gene models.

---

## 📦 Installation

We recommend setting up a dedicated Conda environment to manage dependencies:

```bash
conda create -n gapms_env python=3.10 perl-xml-parser perl-xml-twig -c conda-forge
conda activate gapms_env
git clone https://github.com/qussai96/GAP-MS.git
cd GAP-MS
pip install -e .
```

### Optional: Embeddings Support

To use the embeddings-based scoring feature (`-cm` flag), install additional dependencies:

```bash
pip install torch fair-esm
```

**Note:** This requires PyTorch and the ESM2 model (~2.6GB download), only needed if you plan to use `-cm` for computing embeddings-based external scores.

---

## 📂 Required Inputs

To run GAP-MS, you need the following input files:

1. **Gene models** in GTF or GFF format  
2. **Peptide list** - plain text (.txt) file containing peptides identified from a mass spectrometry experiment [Each peptide should be on a separate line]
3. **Either**:
   - Translated protein sequences (FASTA format) of the gene models, **OR**
   - Genome assembly (FASTA format) to generate proteins automatically

---

## 🧪 Running the Pipeline

### Basic Usage

**With translated protein FASTA file:**

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt
```

**With genome assembly (proteins will be generated automatically):**

```bash
gapms -g annotations.gtf -a assembly.fasta -p peptides.txt
```

> ℹ️ **Note**: If the `-o` output directory is not specified, results will be written to `GAPMS_Output/` in the parent directory of the GTF file.

---

## 🔬 Advanced Options

### Use Custom Mapping File

Skip the proteomapper step by providing a precomputed peptide-to-protein mapping:

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt -m mapping.tsv
```

The mapping file must include at least: `peptide` and `protein` columns.

### Add External Protein Scores

**Option 1: Use precomputed scores**

Provide a CSV file with `Protein` and `external_score` columns:

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt -s scores.csv
```

**Option 2: Compute PSAURON scores**

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt -c
```

**Option 3: Compute embeddings-based scores**

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt -cm
```

> ⚠️ **Note**: Options `-c` and `-cm` cannot be used together.

### Identify Novel Proteins

Compare predictions with reference annotations to find novel supported proteins:

```bash
gapms -g annotations.gtf -f proteins.fasta -p peptides.txt \
      -rg reference.gtf -rf reference_proteins.fasta
```

---

## 🔄 Complete Example

```bash
gapms -g predictions.gtf -a genome.fasta -p peptides.txt \
      -c -rg reference.gtf -rf reference.fasta -o results/
```

---

## � Output Files

GAP-MS generates comprehensive output files for analysis and visualization:

### Main Output Directory (default: `GAPMS_Output/`)

```
GAPMS_Output/
├── log.txt                          # Pipeline execution log
├── all_proteins_scores.tsv          # Feature scores for all proteins
├── supported_proteins.gtf           # GTF with only supported proteins
├── supported_proteins.fa            # FASTA with only supported proteins
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
├── Peptide_Tracks/                 # 🆕 Genomic peptide mappings
│   ├── peptides_mapped.bed        # BED format for genome browsers
│   ├── peptides_mapped.gff3       # GFF3 format
│   └── peptides_mapped.tsv        # Detailed table
└── Novel/                          # Novel protein analysis (if -rg/-rf used)
```

### Peptide Genomic Tracks (NEW!)

GAP-MS automatically generates genome browser tracks for all mapped peptides:

- **BED file**: Color-coded tracks showing peptide locations
  - 🔴 Red = Spliced peptides (spanning multiple exons)
  - 🟢 Green = N-terminal peptides
  - 🔵 Blue = C-terminal peptides
  - ⚪ Gray = Internal peptides

- **GFF3 file**: Hierarchical annotation format

- **TSV file**: Detailed table with all peptide coordinates

**Load these files in IGV, UCSC Genome Browser, or JBrowse to visualize peptide evidence on the genome!**

---

## �📁 Example Data

You can find examples of input files in the `tutorials/` directory to help you get started.

---

## 📬 Feedback & Issues

If you encounter bugs or have suggestions, please open an issue on the [GitHub repository](https://github.com/qussai96/GAP-MS/issues).

---
