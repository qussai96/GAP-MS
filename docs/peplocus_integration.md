# PepLocus integration in GAPMS 0.2.0

## Safety snapshot and provenance

Before integration, tracked GAPMS changes were committed as
`40d1f9a392e073e909be0cd3e93d5d0c7dd1682e` and tagged `pre-peplocus-20261101`.
The pre-existing untracked directory `tests/tutorial_image_output_20261001_141558/`
was preserved in place and was not added as generated data to Git.

The bundled engine is PepLocus 0.3.1, under `gapms/_peplocus/`, with its MIT license.
Its only integration-specific engine change is recognition of `Peptide Sequence`
as an input header. It requires NumPy, Biopython, pyahocorasick and pysam; the last
two have been added to GAPMS installation dependencies. No independently published
PepLocus distribution is required. Existing container images need rebuilding.

## Architecture and output contract

The console entry point is `gapms.cli:main`. A lightweight router dispatches
`gapms map` without importing the classifier. Existing workflow commands dispatch
to `gapms.main`. The module invocation `python -m gapms.main map` is also supported.

`gapms.map_cli` reuses the bundled engine's parser and adds short options `-a`,
`-p`, `-g`, `-b`, `-o`, and `-t`. Annotation and BAM inputs are optional and
combinable. I/L equivalence is mandatory; exact-first optional one-AA-edit rescue,
reusable junction TSVs, bounded novel junction discovery, and engine verification
remain available. Every assembly run enables the IGV display export.

`gapms.peptides_to_genome.map_assembly` invokes the engine and exposes two exact
copies of its BED outputs:

| GAPMS file | Engine file | Role |
|---|---|---|
| `peptides_mapped.bed` | `peplocus/mappings.bed` | Standard BED12, N-to-C query identity |
| `peptides_mapped.igv.bed` | `peplocus/mappings.igv_display.bed` | Historical reverse-strand display labels |

The nested engine TSV/JSON files retain sequence, evidence, mapping multiplicity,
search limits, and timing. Rescue output remains separate under `peplocus/`.
Existing products are not silently overwritten. Engine results are verified by
default before publication; BED aliases are copied only after successful mapping.

The full assessment workflow passes the assembly, original peptide input, branch
annotation, optional BAM, and worker count to this adapter. It still obtains
protein features and classifier decisions from protein mappings. Extra genomic
placements are not automatically linked to predicted proteins.

## Preserving reference-event semantics

The previous projection implementation is retained as `gapms.protein_projection`
for protein-only compatibility, CDS-coordinate helpers, and model metadata. It no
longer generates the authoritative genomic BED for assembly runs.

When a full workflow requests reference comparisons, a separate projection is
written to `protein_projection/peptides_mapped.bed`. The report joins its metadata
to genomic loci using I/L-equivalent peptide, contig, strand and all blocks.
Protein projections unsupported by the genomic mapping are excluded; additional
PepLocus placements remain in global ambiguity checks. Thus protein-specific
terminal claims cannot be inferred merely from a genome-only match, and a newly
found alternate locus can disqualify an otherwise unique event-support peptide.
As of GAPMS 0.2.2, cached CDS interval indexes generate these metadata with the
same records and ordering; see [metadata validation](model_metadata.md).

A historical workflow supplied only with protein FASTA and annotation, without
assembly, retains its previous export with a warning. Such a run cannot perform
six-frame mapping or produce the new pair of genomic BED files. Supply `-a` to
use the new mapping step.

## Validation protocol

Slurm script: `scripts/validate_peplocus.sbatch`; initial job: `188267`.
The script builds and installs the GAPMS 0.2.0 wheel, runs regression tests, changes
outside the checkout, and invokes the installed console command for:

1. Carrot assembly + all 234,952 FragPipe queries + Helixer annotation.
2. The same assembly and queries + RNA-seq BAM, independently of annotation.
3. The full GAPMS prediction/classifier workflow with the existing Helixer
   peptide-to-protein mapping held fixed, followed by the new genomic mapping.

The assembly is `GCF_001625215.2_DH1_v3.0_genomic.fna`, 441,119,532 bp, 11 contigs.
Eight CPUs and 24 GB were requested; one run per configuration is not a scaling
study. Original inputs and previous output directories were not modified.

The requested previous BED is:

```
/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/outputs/gapms/latest_w_bam/helixer/prediction_search/peptides_mapped.bed
```

Its SHA256 is
`f7e7496a8862e311c14965b1fa84c09fcf760fff331869ef5a8fa3edfc7a1bc8`, identical
to the original PepLocus evaluation reference. Its historical negative-strand
labels must be reversed before comparing biological peptide identity.

`scripts/compare_peplocus_beds.py` independently translates reference blocks with
Biopython, reports inconsistent/invalid records, and compares exact chromosome,
strand, block arrangement and I/L-equivalent peptide class. Recall uses only
sequence-consistent reference loci. Additional genomic placements are not
classified as biological false positives. It also checks the IGV-to-standard
round trip, equality to the previous standalone PepLocus GTF output, and equality
between the full workflow and standalone `gapms map`.

## Measured results

Validation jobs 188267 and 188268 completed successfully. The final wheel passed
70 GAPMS tests in an isolated install; the bundled engine passed all 153 PepLocus
tests. The first real-data job used the initial integration wheel; the final wheel
adds the separately tested reference-report ambiguity merge and documentation,
without changing the mapping engine or its genomic results.

| Mode | Mapped peptides | Loci | Reference loci recovered / 110,011 | Missed | Recall | Mapping seconds |
|---|---:|---:|---:|---:|---:|---:|
| gapms map, GTF | 234,952 | 416,611 | 109,954 | 57 | 99.9482% | 133.477 |
| gapms map, BAM | 233,591 | 415,977 | 109,532 | 479 | 99.5646% | 305.286 |
| Full workflow, GTF | 234,952 | 416,611 | 109,954 | 57 | 99.9482% | 133.941 |

The complete prediction/classifier workflow took **11.80 minutes**, including its
existing protein-feature/scoring stages. It used the prior protein mapping rather
than rerunning Proteomapper. BAM mapping time includes 186.5 seconds of junction
extraction. GTF mapping has 187,139 single-location and 47,813 multiple-location
peptides (20.35% ambiguity among mapped queries); BAM has 185,647 and 47,944.

The GTF run recovers at least one reference locus for 89,840/89,891 evaluable
reference peptide classes. Every one of the 57 missed reference loci contains an
intron boundary absent from the annotation-derived junction set. This observation
does not prove whether that reference placement is biologically correct.

All standard/IGV locus sets agree after undoing the IGV negative-strand display
labels. The GTF locus set is exactly equal to the previous standalone PepLocus GTF
run, and the full workflow locus set is exactly equal to standalone `gapms map`.
All mappings passed in-run coordinate/translation verification; a separate verify
command validated all 416,611 GTF loci.

The old reference audit reproduces the earlier result: 110,011 sequence-consistent
loci, 178,907 sequence-inconsistent records and 222 invalid block records. Reported
recall is conditional on the sequence-consistent subset. Additional compatible
placements are not automatically false positives, so these experiments do not
estimate biological precision.

Full results: `evaluation/peplocus_integration/run_188267/REPORT.md`,
`comparison.json`, `reference_audit.tsv`, per-mode missed-reference TSVs,
`gtf_miss_diagnosis.json`, and the three mapping directories. Generated data are
excluded from source control; a compact result record is in
[peplocus_integration_results.json](peplocus_integration_results.json).
