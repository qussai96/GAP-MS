# GAP-MS

GAP-MS assesses gene predictions using mass-spectrometry peptides and optional
RNA-seq evidence. It maps peptides to the assembly with PepLocus, measures
support for supplied gene models, and reports specific annotation changes for review.
Multiple annotations are combined automatically using the included GeneSensus engine.

## Install

Requires Python 3.10 or newer. Install the release wheel with pip:

```bash
python -m pip install gapms-1.3.0-py3-none-any.whl
```

Or install this checkout:

```bash
python -m pip install .
```

The package installs the `gapms` command and PepLocus 0.5.6 automatically. It does
not require a separate legacy GAPMS installation. `gapms-v1` remains an alias.
The 1.3.0 artifacts must be published before `pip install gapms==1.3.0` can be
used from PyPI.

RNA reconstruction also needs StringTie and TransDecoder on PATH:

```bash
conda install -c conda-forge -c bioconda stringtie transdecoder=5.5.0
```

## Run

```bash
gapms -a assembly.fna -p peptides.tsv -g helixer.gtf -o results

# Include RNA-seq evidence
gapms -a assembly.fna -p peptides.tsv -g helixer.gtf -b sorted.bam -o results

# Multiple predictions: consensus is automatic
gapms -a assembly.fna -p peptides.tsv -g helixer.gtf -g another.gff3 -o results
```

Assembly, peptides and an output directory are required; GTF/GFF and BAM are
optional. Eight worker processes are used by default. Use `-t` to change this.
Isoleucine and leucine are always equivalent. Peptides shorter than seven
residues are excluded by default and reported separately.

`gapms --help` shows the common inputs. `gapms --help-all` lists detailed controls,
including optional one-edit rescue, with a nine-residue default minimum. Rescue
candidates are kept separate from exact evidence and do not support gene calls.

## Results

Open `results/report.html` first. `report.log` records progress and errors.

- `peptides.bed`: all exact/I/L-equivalent placements in standard BED12.
- `peptides.tryptic.bed`: tryptic-compatible subset, with its own multiplicity scores.
- `peptides.igv.bed`: display track, with reverse-strand labels arranged for IGV.
- `unmapped_peptides.tsv`: unmapped, rescue-only and excluded short queries.
- `all_proteins_scores.tsv`: original GAPMS feature/confidence table format.
- `supported_genes.gtf`: high-confidence seeds and classifier-supported unlabeled models.
- `review/`: nonempty novel-ORF and alternative start, stop or splice proposals.
- `diagnostics/`: detailed evidence, classification, provenance and reconstruction records.

Multiple input annotations also produce `consensus.gff3` and `consensus_proteins.faa`.
Review proposals are not automatically applied to the annotation. Shared peptides
and existing isoforms are considered before proposing a change. A peptide match
establishes sequence compatibility; it does not validate an entire gene model or
measure peptide identification error.

See the [technical guide](docs/v1/README.md), [methods](docs/v1/SCIENTIFIC_METHODS.md)
and [output schemas](docs/v1/OUTPUT_FORMATS.md).

## Cite

[GAPMS manuscript on bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.17.712294v3).
