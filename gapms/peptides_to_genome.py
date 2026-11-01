"""PepLocus assembly mapping and explicit legacy protein-report compatibility.

The bundled engine is PepLocus 0.3.1 (MIT); no external executable is required.
"""
from pathlib import Path
import shutil
import warnings

# Preserve the public projection helpers used by external callers and old tests.
def __getattr__(name):
    if name in {'export_to_bed', 'get_protein_to_cds_mapping', 'map_codon_to_genome', 'annotate_cds_coordinates'}:
        from . import protein_projection
        return getattr(protein_projection, name)
    raise AttributeError(name)


def map_assembly(assembly, peptides, output_dir, *, gtf=None, bam=None,
                 peptide_column=None, threads=8, extra_args=()):
    """Run bundled PepLocus, always producing analysis and IGV BED files."""
    from ._peplocus.cli import main
    output_dir = Path(output_dir)
    # Refuse overwriting either alias before doing any mapping work.
    products = ['peptides_mapped.bed', 'peptides_mapped.igv.bed']
    for name in products:
        if (output_dir / name).exists():
            raise ValueError(f'output already exists: {output_dir / name}')
    argv = ['map', '--genome', str(assembly), '--peptides', str(peptides),
            '--out', str(output_dir / 'peplocus'), '--threads', str(threads),
            '--igv-display-bed']
    if gtf is not None:
        argv += ['--gtf', str(gtf)]
    if bam is not None:
        argv += ['--bam', str(bam)]
    if peptide_column is not None:
        argv += ['--peptide-column', str(peptide_column)]
    argv += list(extra_args)
    result = main(argv)
    if result:
        raise RuntimeError(f'PepLocus mapping failed (exit {result})')
    for source, target in [('mappings.bed', products[0]),
                           ('mappings.igv_display.bed', products[1])]:
        shutil.copyfile(output_dir / 'peplocus' / source, output_dir / target)
    return output_dir / products[0]


def map_peptides_to_genome(gtf_file, protein_fasta, mapping_file, output_dir,
                           *, assembly=None, peptides=None, bam=None, threads=8,
                           peptide_column=None, model_evidence=False):
    """Workflow adapter; assembly runs use the full observed peptide input.

    Protein-only historical workflows cannot run six-frame mapping and retain
    their old projection behavior with an explicit warning.
    """
    output_dir = Path(output_dir)
    if assembly is None:
        from .protein_projection import map_peptides_to_genome as project_proteins
        warnings.warn('No assembly supplied: using protein projection; add -a to '
                      'enable PepLocus genomic mapping.', UserWarning)
        project_proteins(gtf_file, protein_fasta, mapping_file, output_dir)
        return output_dir / 'peptides_mapped.bed'
    if peptides is None:
        raise ValueError('assembly mapping requires the original peptide input')
    bed = map_assembly(assembly, peptides, output_dir, gtf=gtf_file, bam=bam,
                       threads=threads, peptide_column=peptide_column)
    if model_evidence:
        from .protein_projection import map_peptides_to_genome as project_proteins
        # Reports need protein assignment/terminal metadata, not inferred labels
        # from a genome-only match. Keep that evidence separate from genomic BED.
        compatibility = output_dir / 'protein_projection'
        project_proteins(gtf_file, protein_fasta, mapping_file, compatibility)
    return bed


def main():
    from .map_cli import main as map_main
    return map_main()


if __name__ == '__main__':
    raise SystemExit(main())
