# GAP-MS Annotation Comparison Report - Quick Reference

## What It Does

Generates TSV reports comparing peptide-supported gene predictions with reference annotations to identify:
- Different start/stop codons
- Different exon structures  
- Extra exons or missing exons
- Novel genes not in reference (intergenic or antisense)
- Other annotation discrepancies

## Output Files

### annotation_comparison_summary.tsv
- **Category**: Type of annotation difference
- **Total_Count**: Number of genes in this category
- **Peptide_Supported_Count**: How many have peptide support
- **Gene_IDs**: All gene IDs in this category (in curly braces)

Format: `TSV with 4 columns`

### annotation_comparison_report.tsv  
- **predicted_id**: Your predicted gene model ID
- **reference_id**: Matching reference gene ID
- **gffcompare_class**: Type of structural match (=, c, o, j, u, x, etc.)
- **annotation_category**: Detailed category (same_start_different_stop, etc.)
- **mapped_peptides**: Number of peptides mapping to this gene
- **sequence_coverage**: Proportion of the protein sequence covered by peptides

Format: `TSV with 6 columns, one row per gene`

## Key Categories to Focus On

| Category | Meaning | Why It Matters |
|----------|---------|-----------------|
| `same_start_different_stop` | Different C-terminal | Reference may have wrong stop codon |
| `different_start_same_stop` | Different N-terminal | Reference may have wrong start or missing exon |
| `same_start_same_stop_different_exons` | Same boundaries, different splicing | Likely alternative isoforms |
| `upstream_extension` | Extends further upstream | Extended translation start region |
| `downstream_truncation` | Truncated at 3' end | Missing C-terminal region |
| `Novel_genes` | Intergenic (u), antisense (x), or unknown (-) | Completely new genes not in reference |

## Running the Analysis

### Automatic (In GAP-MS Pipeline)
```bash
python -m gapms.main \
    --gtf predictions.gtf \
    --peptides peptides.txt \
    --assembly genome.fasta \
    --reference-gtf reference.gtf \
    --output output_directory
```

### Manual (Standalone Tool)
```bash
python gapms_annotation_report.py /path/to/gapms/output
```

**Optional arguments:**
- `--reference-gtf PATH` - Path to reference GTF (auto-detected if in standard location)
- `--supported-gtf PATH` - Path to supported proteins GTF (default: output_dir/supported_proteins.gtf)
- `--tmap PATH` - Path to gffcompare .tmap file (auto-detected if not provided)
- `--scores PATH` - Path to all_proteins_scores.tsv (auto-detected if not provided)

## Interpreting Results

### High Peptide-Supported Genes with Different Start/Stop
- Strong evidence for annotation error in reference
- Consider manual curation
- Good candidates for genome re-annotation

### Same Start/Stop but Different Exons  
- Likely alternative isoforms
- Check peptide coverage of different exon combinations

### Upstream Extensions
- Could be alternative translation start sites
- Check for N-terminal peptide support

### Novel Genes (Class Codes u, x, -)
- Potentially new genes missed by reference annotation
- `u` = intergenic (between reference genes)
- `x` = antisense (opposite strand from reference)
- `-` = unmapped (no clear reference relationship)
- Check peptide quality and coverage before publication

## File Locations

In GAP-MS output directory:
```
output_dir/
├── annotation_comparison_report.tsv
├── annotation_comparison_summary.tsv
├── all_proteins_scores.tsv
├── supported_proteins.gtf
└── gffcmp.supported_proteins.gtf.tmap
```

## Tips for Analysis

1. **Filter by peptide support**: Focus on genes with mapped_peptides > 0
2. **Check coverage**: Higher sequence_coverage = more reliable predictions
3. **Prioritize categories**: Start with "same_start_different_stop" (likely real differences)
4. **Cross-validate**: Use IGV or other viewers to visualize differences
5. **Check splice sites**: Use peptide mapping to confirm exon boundaries

## Common Questions

**Q: What's the difference between annotation_comparison_report and annotation_comparison_summary?**
A: Report = detailed row per gene. Summary = statistics per category.

**Q: Can I filter the TSV files?**
A: Yes! Use standard tools like:
```bash
# Filter for genes with >50 peptides
awk -F'\t' '$5 > 50' annotation_comparison_report.tsv

# Get all genes with different stop codons
grep "same_start_different_stop" annotation_comparison_report.tsv
```

**Q: How do I find specific genes?**
A: Look in annotation_comparison_summary.tsv Gene_IDs column - genes are comma-separated in curly braces.

**Q: What does class code 'u' or 'x' mean?**
A: `u` = intergenic (predicted gene is in intergenic region), `x` = antisense (opposite strand). Both are now classified as Novel_genes.

## Contact & Support

For issues or questions about the annotation comparison report, refer to the GAP-MS documentation.
