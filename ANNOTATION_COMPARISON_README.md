# GAP-MS Annotation Comparison Report

## Overview

The Annotation Comparison Report is a feature in GAP-MS that generates comprehensive TSV reports analyzing differences between peptide-supported predicted gene models and reference gene annotations. This helps identify potential annotation artifacts in reference genomes.

## Features

The script generates two TSV output files:

1. **annotation_comparison_report.tsv** - Detailed record for each supported gene model
   - Columns: `predicted_id`, `reference_id`, `gffcompare_class`, `annotation_category`, `mapped_peptides`, `sequence_coverage`
   - Contains information about how each predicted gene relates to the reference

2. **annotation_comparison_summary.tsv** - Summary statistics by annotation category
   - Columns: `Category`, `Total_Count`, `Peptide_Supported_Count`, `Gene_IDs`
   - Aggregated statistics showing how many genes fall into each category
   - Gene IDs listed for further analysis

## Annotation Categories

The script classifies genes into the following categories based on CDS coordinate comparison:

### Structure-Based Categories

- **Identical_matches**: Exact matches with reference genes (class code '=')
- **same_start_same_stop**: Same translation start and stop, identical exon structure
- **same_start_same_stop_different_exons**: Same start/stop but different splice sites
- **same_start_different_stop**: Translation start matches, but different stop position
- **different_start_same_stop**: Translation stop matches, but different start position
- **upstream_extension**: Predicted model extends upstream of reference
- **downstream_truncation**: Predicted model truncated downstream of reference
- **more_exons**: Predicted model has more exons than reference
- **fewer_exons**: Predicted model has fewer exons than reference

### Novel Gene Categories

- **Novel_genes**: Predicted genes with no reference match or different strand assignment
  - Includes class code `-` (unknown/unassigned)
  - Includes class code `u` (intergenic - in regions between reference genes)
  - Includes class code `x` (antisense - on opposite strand from reference)

### Data Quality Categories

- **missing_coordinates**: Unable to determine CDS coordinates for comparison
- **different_chromosome_or_strand**: Predicted and reference on different chromosomes/strands

## gffcompare Class Codes

The script uses gffcompare output to classify relationships:

| Code | Meaning |
|------|---------|
| `=` | Exact match - identical gene structure |
| `c` | Contained - predicted is within reference |
| `k` | Reverse complement |
| `m` | Matches opposite strand |
| `n` | Intron chain match |
| `j` | New isoform - same locus, different structure |
| `e` | Exon match |
| `o` | Other overlap on same strand |
| `s` | Intron chain share |
| `y` | Improperly contained |
| `p` | Possible polymerase run |
| `r` | Repeat |
| `u` | Intergenic - novel gene in intergenic region |
| `x` | Antisense - novel gene on opposite strand |
| `-` | Unknown - no clear reference relationship |

**Genes with class codes `u`, `x`, or `-` are classified as Novel_genes.**

## Usage

### In the GAP-MS Pipeline

When running GAP-MS with a reference GTF file, the annotation comparison report is generated automatically:

```bash
python -m gapms.main \
    --gtf predicted.gtf \
    --peptides peptides.txt \
    --assembly genome.fasta \
    --reference-gtf reference.gtf \
    --output output_directory
```

### Standalone Usage

To generate the report on existing GAP-MS output:

```bash
python gapms_annotation_report.py /path/to/gapms/output
```

**Optional arguments:**
- `--reference-gtf PATH` - Path to reference GTF (auto-detected if in standard location)
- `--supported-gtf PATH` - Path to supported proteins GTF (default: output_dir/supported_proteins.gtf)
- `--tmap PATH` - Path to gffcompare .tmap file (auto-detected if not provided)
- `--scores PATH` - Path to all_proteins_scores.tsv (auto-detected if not provided)

## Interpreting Results

### Key Insights

1. **same_start_different_stop** - Predicted gene has extended or truncated 3' end
   - May indicate missing stop codon in reference or real C-terminal extension in prediction
   - Cross-check with peptide coverage and number of mapped peptides

2. **same_start_different_stop with high peptide coverage** - Strong evidence for C-terminal difference
   - Could indicate annotation error in reference gene

3. **different_start_same_stop** - Predicted gene has different translation start
   - May indicate alternative start codon or missing upstream exons in reference
   - Check N-terminal peptide support to validate

4. **upstream_extension/downstream_truncation** - Predicted gene differs in boundary positions
   - Could be due to UTR differences or erroneous boundaries in reference

5. **same_start_same_stop_different_exons** - Same coding region but different splice structure
   - Usually indicates alternate isoforms
   - Check if exon differences are supported by peptides

6. **Novel_genes with peptide support** - Strong candidates for new genes
   - Genes in intergenic regions (class code `u`) or antisense (class code `x`)
   - High peptide support validates these as real genes
   - Good candidates for functional annotation

## Example Output

```
Category                           Total_Count  Peptide_Supported_Count
Identical_matches                  9506         9506
same_start_same_stop              4502         4502
same_start_same_stop_different_exons  2117      2117
same_start_different_stop          1435         1435
different_start_same_stop          1309         1309
missing_coordinates               671          671
Novel_genes                       687          687
downstream_truncation             223          223
upstream_extension                178          178
different_chromosome_or_strand    23           23
```

## Applications

This annotation comparison report is useful for:

1. **Quality assessment** of reference annotations against peptide-supported predictions
2. **Identifying annotation errors** in reference genomes
3. **Finding novel genes** in intergenic regions or on antisense strand
4. **Gene boundary validation** - confirming gene start/stop positions
5. **Isoform detection** - identifying alternative splice variants
6. **Variant calling validation** - checking if predicted models support known variants

## Technical Details

### Algorithm

1. Parses gffcompare TMAP file to identify predicted-to-reference relationships
2. Extracts CDS regions from both predicted and reference GTF files
3. Compares CDS coordinates on a per-strand basis
4. Classifies differences based on start and stop position changes
5. Integrates peptide support metrics from all_proteins_scores.tsv

### Novel Gene Classification

Genes are classified as Novel when gffcompare assigns class codes:
- `u` (intergenic) - Predicted gene located between reference genes
- `x` (antisense) - Predicted gene on opposite strand with no clear reference match
- `-` (unknown) - No clear reference relationship mapping

These represent genuinely novel or previously unannotated genes that should be scrutinized during manual curation.

### Strand-Aware Comparisons

For **positive strand** genes:
- First CDS start = leftmost CDS start (translation initiation)
- Last CDS end = rightmost CDS end (translation termination)

For **negative strand** genes:
- First CDS start = rightmost CDS end (translation initiation)
- Last CDS end = leftmost CDS start (translation termination)

## Output Files Location

When integrated in the GAP-MS pipeline, files are written to:
- `{output_dir}/annotation_comparison_report.tsv`
- `{output_dir}/annotation_comparison_summary.tsv`

## Dependencies

- pandas >= 0.24
- BioPython (optional, for reference file handling)

## Citation

If you use the annotation comparison report in your research, please cite the GAP-MS paper and note the use of gffcompare for structural comparison.

## References

- gffcompare: https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
