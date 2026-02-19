#!/usr/bin/env python3
"""
Standalone tool to generate annotation comparison report for GAP-MS output.

Usage:
    python gapms_annotation_report.py <output_directory> [options]

This script analyzes gffcompare output to identify annotation differences between
predicted gene models (with peptide support) and reference annotations.

Output generates:
    - annotation_comparison_report.tsv: Detailed comparison for each gene
    - annotation_comparison_summary.tsv: Summary statistics by category
"""

import sys
from pathlib import Path

# Add parent directory to path to allow imports from gapms module
sys.path.insert(0, str(Path(__file__).parent))

from generate_annotation_comparison_report import generate_annotation_report


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nRequired Arguments:")
        print("  output_directory    Path to GAP-MS output directory")
        print("\nOptional Arguments:")
        print("  --supported-gtf     Path to supported proteins GTF (default: output_dir/supported_proteins.gtf)")
        print("  --reference-gtf     Path to reference GTF (searches output_dir/../../reference/*.ref.gtf)")
        print("  --tmap              Path to .tmap file (auto-detected if not provided)")
        print("  --scores            Path to all_proteins_scores.tsv (auto-detected if not provided)")
        print("\nExample:")
        print("  gapms_annotation_report.py /path/to/output_dir")
        sys.exit(1)
    
    output_dir = Path(sys.argv[1])
    
    if not output_dir.exists():
        print(f"❌ Error: Output directory not found: {output_dir}")
        sys.exit(1)
    
    # Parse optional arguments
    supported_gtf = output_dir / "supported_proteins.gtf"
    reference_gtf = None
    tmap_file = None
    scores_file = None
    
    i = 2
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--supported-gtf" and i + 1 < len(sys.argv):
            supported_gtf = Path(sys.argv[i + 1])
            i += 2
        elif arg == "--reference-gtf" and i + 1 < len(sys.argv):
            reference_gtf = Path(sys.argv[i + 1])
            i += 2
        elif arg == "--tmap" and i + 1 < len(sys.argv):
            tmap_file = Path(sys.argv[i + 1])
            i += 2
        elif arg == "--scores" and i + 1 < len(sys.argv):
            scores_file = Path(sys.argv[i + 1])
            i += 2
        else:
            i += 1
    
    # Auto-detect reference GTF if not provided
    if reference_gtf is None:
        # Try to find in standard GAP-MS locations
        possible_paths = [
            output_dir / "../../reference" / "*.ref.gtf",
            output_dir / "../../../inputs/reference" / "*.ref.gtf",
            output_dir.parent.parent / "inputs" / "reference" / "*.ref.gtf",
        ]
        
        for pattern in possible_paths:
            matches = list(Path(pattern.parent).glob(pattern.name))
            if matches:
                reference_gtf = matches[0]
                print(f"Auto-detected reference GTF: {reference_gtf}")
                break
        
        if reference_gtf is None:
            print("❌ Error: Could not find reference GTF file")
            print("Please specify with --reference-gtf")
            sys.exit(1)
    
    # Verify files exist
    if not supported_gtf.exists():
        print(f"❌ Error: Supported GTF not found: {supported_gtf}")
        sys.exit(1)
    
    if not reference_gtf.exists():
        print(f"❌ Error: Reference GTF not found: {reference_gtf}")
        sys.exit(1)
    
    print(f"📋 GAP-MS Annotation Comparison Report")
    print(f"=" * 80)
    print(f"Output directory:   {output_dir}")
    print(f"Supported GTF:      {supported_gtf}")
    print(f"Reference GTF:      {reference_gtf}")
    print(f"=" * 80)
    
    # Generate report
    try:
        report_df, summary_df = generate_annotation_report(
            output_dir,
            supported_gtf,
            reference_gtf,
            tmap_file=tmap_file,
            all_proteins_scores=scores_file
        )
        
        print(f"\n✅ Report generation completed successfully!")
        print(f"\nGenerated files:")
        print(f"  - annotation_comparison_report.tsv ({len(report_df)} rows)")
        print(f"  - annotation_comparison_summary.tsv ({len(summary_df)} rows)")
        
    except Exception as e:
        print(f"\n❌ Error generating report: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
