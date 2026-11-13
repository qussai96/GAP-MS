# GAP-MS

GAP-MS uses mass-spectrometry peptides to assess gene predictions, combine
annotations and identify peptide-supported changes for review.

## Install

Python 3.10 or newer:

```bash
pip install gapms
```

For runs with a BAM file, also install the RNA reconstruction tools:

```bash
conda install -c conda-forge -c bioconda stringtie transdecoder=5.5.0
```

Or use the Apptainer image, which includes the RNA reconstruction tools:

```bash
apptainer pull gapms.sif oras://docker.io/qussaiab96/gapms:latest
apptainer run gapms.sif -a assembly.fa -p peptides.tsv -g genes.gtf
```

See [container instructions](packaging/CONTAINERS.md) for Docker Hub.

## Run

- **Assembly + peptides:** map peptides to the six-frame genome translation.
  `gapms -a assembly.fa -p peptides.tsv`
- **Add BAM:** include RNA-supported splice paths and reconstruct peptide-supported ORF candidates.
  `gapms -a assembly.fa -p peptides.tsv -b sorted.bam`
- **Add annotation:** map across annotated exons, score gene models and propose supported changes.
  `gapms -a assembly.fa -p peptides.tsv -g genes.gtf`
- **Annotation + BAM:** combine gene-model assessment with RNA-supported reconstruction.
  `gapms -a assembly.fa -p peptides.tsv -g genes.gtf -b sorted.bam`
- **Multiple annotations:** automatically build a peptide-guided consensus. BAM is optional here too.
  `gapms -a assembly.fa -p peptides.tsv -g helixer.gtf -g another.gff3`

GTF and GFF3 are accepted. BAM files must be coordinate-sorted and use the
assembly's sequence names. FragPipe, Casanovo and plain peptide lists are supported.
Use `-o results` to choose an output directory (default: `./GAPMS_output`) and
`-t 16` to change the default of **8 threads**.

Analysis uses **trypsin-compatible exact matches** by default; add `--non-tryptic`
for other enzymes. Isoleucine and leucine are treated as equivalent. All-match and
tryptic-match BED tracks are always retained. Optional one-edit rescue is separate
from the exact evidence used to support gene models. See `gapms --help-all`.

## Results

Start with **report.html**; **report.log** records the run. Main outputs include
`peptides.bed`, `peptides.tryptic.bed`, `peptides.igv.bed`,
`unmapped_peptides.tsv`, and, when gene models are available,
`all_proteins_scores.tsv` and `supported_genes.gtf`. Multiple annotations also
produce `consensus.gff3` and `consensus_proteins.faa`.

`review/` contains supported novel/alternative proposals when available.
They require review and are not automatically applied. Existing annotated isoforms
are considered before proposing changes. Without an input annotation, RNA-derived
ORFs are candidates whose novelty cannot be assessed. Peptide support does not
validate every boundary of a gene model.

[Technical guide](docs/v1/README.md) ·
[Five-mode validation of 1.3.2](docs/v1/evaluation/five_modes_hop_20261113/REPORT.md)

## Cite

[GAP-MS manuscript on bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.17.712294v3)
