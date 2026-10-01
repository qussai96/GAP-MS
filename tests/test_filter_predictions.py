import pandas as pd
from pathlib import Path
import textwrap

import gapms.filter as filter_mod


def test_filter_predictions_does_not_force_gene_specific_evidence(tmp_path, monkeypatch):
    outdir = tmp_path / "out"
    outdir.mkdir()

    # Create minimal all_proteins_scores.tsv
    df = pd.DataFrame([
        {
            'Protein': 'P1',
            'mapped_peptides': 5,
            'gene_specific_peptides': 0,
            'protein_specific_peptides': 2,
            'sequence_coverage': 0.5,
            'N_terminal_peptides': 0,
            'C_terminal_peptides': 0,
            'splice_peptides': 0,
            'splice_sites': 0,
            'protein_length': 300,
        },
        {
            'Protein': 'P2',
            'mapped_peptides': 3,
            'gene_specific_peptides': 1,
            'protein_specific_peptides': 0,
            'sequence_coverage': 0.1,
            'N_terminal_peptides': 0,
            'C_terminal_peptides': 0,
            'splice_peptides': 0,
            'splice_sites': 0,
            'protein_length': 200,
        },
        {
            'Protein': 'P3',
            'mapped_peptides': 0,
            'gene_specific_peptides': 2,
            'protein_specific_peptides': 0,
            'sequence_coverage': 0.0,
            'N_terminal_peptides': 0,
            'C_terminal_peptides': 0,
            'splice_peptides': 0,
            'splice_sites': 0,
            'protein_length': 150,
        },
    ])
    scores_path = outdir / "all_proteins_scores.tsv"
    df.to_csv(scores_path, sep='\t', index=False)

    # Minimal fasta with all proteins
    fasta = outdir / "proteins.fa"
    fasta.write_text(textwrap.dedent(
        ">P1\nAAAA\n>P2\nBBBB\n>P3\nCCCC\n"
    ))

    # Stub the classifier to mark P2 as positive
    def fake_run_classifier(df, high_confident_df, low_confident_df, **kwargs):
        labeled = pd.DataFrame({'Protein': ['P2', 'P1', 'P3'], 'final_label': [1, 0, 0]})
        return labeled, 'stub'

    monkeypatch.setattr(filter_mod, 'run_protein_classifier', fake_run_classifier)

    # Stub GTF functions used by filter_predictions
    monkeypatch.setattr(filter_mod, 'gtf_to_df_with_genes', lambda p: pd.DataFrame({'Protein': ['P1','P2','P3'], 'Gene': ['G1','G2','G3'], 'Type': ['CDS','CDS','CDS']}))
    monkeypatch.setattr(filter_mod, 'save_gtf_subset_all_features', lambda path, proteins, outdir, name: None)

    # Run
    filter_mod.filter_predictions(gtf_file="dummy.gtf", protein_fasta=str(fasta), output_dir=outdir)

    # P1 is a high-confidence seed and remains supported even though it lacks a
    # gene-specific peptide; P2 is supported by the classifier.
    supported = (outdir / "Txt" / "supported_proteins.txt").read_text().splitlines()
    assert supported == ['P1', 'P2']

    # Check FASTA contains every supported protein.
    fasta_out = (outdir / "supported_proteins.fa").read_text()
    assert ">P1" in fasta_out and ">P2" in fasta_out and ">P3" not in fasta_out

    # The scores table retains gene_specific_peptides for optional downstream filtering.
    updated = pd.read_csv(scores_path, sep='\t')
    supported_marks = dict(zip(updated['Protein'], updated['supported']))
    assert supported_marks['P1'] == '+' and supported_marks['P2'] == '+'
    assert supported_marks['P3'] == '-'
    assert dict(zip(updated['Protein'], updated['gene_specific_peptides']))['P1'] == 0
