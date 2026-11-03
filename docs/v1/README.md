# GAPMS 1.0.2

GAPMS v1 maps mass-spectrometry peptide sequences to eukaryotic assemblies and
organizes evidence for reviewing gene annotations. Genomic peptide locations are
the common basis for annotation support, RNA-derived ORFs and proposed changes.
The software reports candidates and uncertainty; sequence compatibility alone
does not validate a peptide identification or a biological annotation change.

V1 is a separate Python distribution and command. The existing GAPMS 0.3.0
workflow remains available as `gapms`, including its original classifier.

## Installation

Python 3.10 is required. From the separate development checkout:

```bash
# Build the pinned dependency from the preserved legacy source in this branch.
python -m pip install setuptools wheel
python scripts/build_legacy_wheel.py --out dist
python -m pip install dist/gapms-0.3.0-py3-none-any.whl .
gapms-v1 --version
```

The dependency wheel can also be reused from an existing GAPMS 0.3.0 release.
The build helper copies only source and model files to a temporary directory;
it does not alter or replace the legacy command.

BAM reconstruction additionally needs StringTie and TransDecoder on PATH:

```bash
conda install -c conda-forge -c bioconda stringtie transdecoder=5.5.0
```

Genome and annotation modes need no external mapping executable. The C-backed
Python dependencies come through GAPMS 0.3.0. V1 uses its own Python transcript
extraction and creates FASTA indexes inside the output, never beside inputs.

## Run

```bash
# Genome-wide peptide coordinates
python -m gapms_v1 -a genome.fna -p peptides.tsv -o genome_results

# Add model assessment and annotated junction peptides
gapms-v1 -a genome.fna -p peptides.tsv -g helixer.gtf -o annotation_results

# Reconstruct RNA-supported ORFs
gapms-v1 -a genome.fna -p peptides.tsv -b sorted.bam -o rna_results

# Unified annotation review; repeat -g for independent sources
gapms-v1 -a genome.fna -p peptides.tsv -b sorted.bam \
  -g helixer.gtf -g another.gff3 -t 8 -o integrated_results
```

Assembly and peptides are required. Annotation and coordinate-sorted BAM are
optional. GTF/GFF3 transcript relationships, compressed annotation/assembly,
and plain peptide lists, FASTA, CSV or TSV are supported. Table peptide columns
are detected by common names or set with `--peptide-column`. I/L equivalence is
always enabled. The default minimum peptide length is 7; shorter queries are
reported in `mapping/excluded_peptides.tsv`. This threshold limits ambiguity,
not peptide identification error. Use a lower threshold explicitly for small
peptide research. No uniqueness claim can cover unsearched splice structures.

Optional `--rescue-one-edit` searches engine-stage exact-unmapped peptides of at least 12 aa
by default. Rescue outputs are separate and do not contribute to exact model
support or annotation proposals. All exact placements are retained.

A new output directory is required. `--resume` only accepts unchanged inputs,
parameters and implementation. An interrupted run can reuse completed
reconstruction and mapping stages; completed runs are returned unchanged.
Do not delete stage outputs selectively and expect resume to repair them.

## Results

Start with `reports/summary.json`, `mapping/mapped_peptides.bed` and
`mapping/mapped_peptides.tsv`. The `.igv.bed` track reverses negative-strand
labels for display; use the standard BED for analysis.

Annotation mode adds `annotation/gene_model_evidence.tsv`, supported-model GTF,
and separate no-observed-support and ambiguous tables. RNA mode adds transcripts,
junction observations, ORF evidence and novel ORF candidates. Integrated mode
adds `refinement/annotation_change_proposals.tsv` and separate start, end and
splice candidate tables and GTFs. Sources and original model IDs are retained.

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
