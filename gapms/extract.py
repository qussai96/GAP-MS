
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
from Bio import SeqIO

from .gtf_utils import gtf_to_df_with_genes, read_scores_csv
from .peptide_utils import (
    mapping_file_to_df,
    get_gene_protein_specific_peps,
    calculate_sequence_coverage,
)
from .plotting import (
    plot_mapped_percentage_bars
)


def _ts():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def extract_features(gtf_file, prediction_fasta, mapping_file, output_dir, external_scores_csv=None):
    gtf_file = Path(gtf_file)
    protein_fasta = Path(prediction_fasta)
    mapping_file = Path(mapping_file)
    external_scores_csv = Path(external_scores_csv) if external_scores_csv else None
    output_dir=Path(output_dir)

    print(f"{_ts()} [extract] Parsing GTF/GFF...")
    gtf_df = gtf_to_df_with_genes(gtf_file)
    print(f"{_ts()} [extract] GTF parsed: {len(gtf_df)} rows. Building CDS features...")
    print(f"Sample protein IDs from GTF: {gtf_df['Protein'].dropna().head().tolist()}")
    transcripts_df = gtf_df[gtf_df['Type'].isin(['transcript', 'mRNA'])][['Protein', 'Prediction_score']]
    transcripts_df = transcripts_df.rename(columns={"Prediction_score": "transcript_score"})


    # Filter for CDS features
    cds_df = gtf_df.loc[gtf_df['Type'] == 'CDS'].copy()

    # Compute isoforms per gene
    isoforms = (
        cds_df.groupby('Gene')['Protein']
        .nunique()
        .reset_index(name='isoforms')
    )

    # Compute protein start and end positions
    protein_coords = (
        cds_df.groupby('Protein')
        .agg(Protein_Start=('Start', 'min'),
            Protein_End=('End', 'max'))
        .reset_index()
    )

    # Merge precomputed info
    cds_df = cds_df.merge(protein_coords, on='Protein').merge(isoforms, on='Gene')
    cds_df = cds_df.merge(transcripts_df, on='Protein', how='left')

    # Prepare CDS dataframe for codon positions
    df_CDSs = (
        cds_df.sort_values(['Protein', 'Start'])
            .drop_duplicates(subset=['Seqid', 'Start', 'End', 'Protein'])
            .copy()
    )
    df_CDSs['cds_start'] = 0
    df_CDSs['cds_end'] = 0

    # Calculate CDS lengths in codons
    df_CDSs['cds_len'] = ((df_CDSs['End'] - df_CDSs['Start']) / 3).round().astype(int)

    # Forward strand (+)
    forward = df_CDSs['Strand'] == '+'
    df_CDSs.loc[forward, 'cds_start'] = (
        df_CDSs[forward].groupby('Protein')['cds_len'].cumsum() - df_CDSs[forward]['cds_len'] + 1
    )
    df_CDSs.loc[forward, 'cds_end'] = df_CDSs[forward].groupby('Protein')['cds_len'].cumsum()

    # Reverse strand (-)
    reverse = ~forward
    reverse_group = df_CDSs[reverse].groupby('Protein')['cds_len']
    cumsum_reverse = reverse_group.cumsum()
    max_len = reverse_group.transform('sum')
    df_CDSs.loc[reverse, 'cds_end'] = max_len - cumsum_reverse + df_CDSs.loc[reverse, 'cds_len']
    df_CDSs.loc[reverse, 'cds_start'] = df_CDSs.loc[reverse, 'cds_end'] - df_CDSs.loc[reverse, 'cds_len'] + 1

    # Compute number of splice sites per protein
    df_CDSs['splice_sites'] = df_CDSs.groupby('Protein')['Protein'].transform('count') - 1

    # Merge codon positions back into original CDS dataframe
    cds_df = cds_df.merge(
        df_CDSs[['Protein', 'Start', 'cds_start', 'cds_end', 'splice_sites']],
        on=['Protein', 'Start'],
        how='left'
    )

    # Ensure correct datatype
    cds_df = cds_df.astype({'cds_start': 'Int64', 'cds_end': 'Int64'})


    # 1. Load mapping
    print(f"{_ts()} [extract] Loading mapping file...")
    mapping_df = mapping_file_to_df(mapping_file)

    # 2. Get gene-protein specific peptides
    print(f"{_ts()} [extract] Resolving gene/protein-specific peptides...")
    pep_df = get_gene_protein_specific_peps(cds_df, mapping_df)

    # 3. Create protein sequence dictionary
    print(f"{_ts()} [extract] Loading protein FASTA sequences...")
    with open(prediction_fasta) as fasta_handle:
        seq_dict = {
            record.id.split()[0]: str(record.seq).replace('I', 'L').replace('*', '')
            for record in SeqIO.parse(fasta_handle, "fasta")
        }

    # 4. Map protein sequences
    pep_df['Prot_seq'] = pep_df['Protein'].map(seq_dict)
    # print(f"Peptides after mapping protein sequences: {len(pep_df)}")
    pep_df.dropna(subset=['Prot_seq'], inplace=True)
    # print(f"Peptides after filtering out rows with missing protein sequences: {len(pep_df)}")

    # 5. Compute peptide positions
    print(f"{_ts()} [extract] Computing peptide positions for {len(pep_df)} pep-protein pairs...")
    pep_starts = []
    pep_ends = []

    for prot_seq, peptide in zip(pep_df['Prot_seq'], pep_df['Peptide']):
        normalized_peptide = peptide.replace('I', 'L')  # Normalize for I/L isobaric amino acids
        start = prot_seq.find(normalized_peptide) + 1  # 1-based
        end = start + len(peptide) - 1 if start > 0 else 0
        pep_starts.append(start)
        pep_ends.append(end)

    pep_df['pep_start'] = pep_starts
    pep_df['pep_end'] = pep_ends
    pep_df['prot_len'] = pep_df['Prot_seq'].str.len()
    pep_df['N_terminal_peptide'] = np.where(
        (pep_df['pep_start'] == 1) |
        ((pep_df['pep_start'] == 2) & pep_df['Prot_seq'].str.startswith('M')),
        '+',
        '-'
    )

    # 6. Drop temporary column
    pep_df.drop(columns='Prot_seq', inplace=True)

    # 7. Compute sequence coverage (vectorized)
    if pep_df.empty:
        print("\n=== DEBUG INFO ===")
        print("Sample protein IDs from FASTA (seq_dict):")
        print(list(seq_dict.keys())[:5])
        print("\nSample protein IDs from mapping file:")
        print(mapping_df['Protein'].head().tolist())
        print("==================\n")
        raise ValueError("No peptides could be mapped to proteins. Check that protein IDs in the mapping file match the FASTA headers. We use only first part of the FASTA header up to the first space.")

    print(f"{_ts()} [extract] Computing sequence coverage...")
    coverage_map = calculate_sequence_coverage(pep_df)
    pep_df['sequence_coverage'] = pep_df['Protein'].map(coverage_map).fillna(0)

    print(f"{_ts()} [extract] Merging CDS and peptide features ({len(cds_df)} CDS rows × {len(pep_df)} pep-protein pairs)...")
    peptide_features_df = pd.merge(cds_df, pep_df, on=['Protein', 'Gene'])
    peptide_features_df = peptide_features_df[[
        'Peptide', 'Protein', 'Gene', 'isoforms', 'splice_sites', 'Protein_specific',
        'Gene_specific', 'cds_start', 'cds_end', 'peptide_length',
        'pep_start', 'pep_end', 'N_terminal_peptide', 'transcript_score', 'prot_len', 'sequence_coverage'
    ]]

    print(f"{_ts()} [extract] Computing splice-site flags on {len(peptide_features_df)} rows...")
    cds_intervals = peptide_features_df[['Protein', 'cds_start', 'cds_end']].drop_duplicates()
    pep_unique = peptide_features_df[['Protein', 'pep_start', 'pep_end']].drop_duplicates()

    check = pep_unique.merge(cds_intervals, on='Protein', how='left')
    check['contained'] = (
        (check['cds_start'] <= check['pep_start']) &
        (check['cds_end'] >= check['pep_end'])
    )
    internal_peps = (
        check[check['contained']]
        .drop_duplicates(subset=['Protein', 'pep_start', 'pep_end'])
        [['Protein', 'pep_start', 'pep_end']]
    )
    internal_peps['_internal'] = True
    peptide_features_df = peptide_features_df.merge(
        internal_peps, on=['Protein', 'pep_start', 'pep_end'], how='left'
    )
    peptide_features_df['Splice_peptide'] = np.where(
        peptide_features_df['_internal'].fillna(False), '-', '+'
    )
    peptide_features_df.drop(columns=['_internal'], inplace=True)

    peptide_features_df['C_terminal_peptide'] = np.where(peptide_features_df['pep_end'] == peptide_features_df['prot_len'], '+', '-')
    peptide_features_df = peptide_features_df.drop(columns=['cds_start', 'cds_end']).drop_duplicates()


    print(f"{_ts()} [extract] Aggregating protein-level scores...")
    proteins_scores_df = peptide_features_df.groupby("Protein").agg(
        protein_length=("prot_len", "mean"),
        isoforms=("isoforms", "mean"),
        splice_sites=("splice_sites", "mean"),
        transcript_score=("transcript_score", "mean"),
        sequence_coverage=("sequence_coverage", lambda x: round(x.mean(), 2)),
        mapped_peptides=("Peptide", "nunique"),
        N_terminal_peptides=("N_terminal_peptide", lambda x: (x == "+").sum()),
        C_terminal_peptides=("C_terminal_peptide", lambda x: (x == "+").sum()),
        protein_specific_peptides=("Protein_specific", lambda x: (x == "+").sum()),
        gene_specific_peptides=("Gene_specific", lambda x: (x == "+").sum()),
        splice_peptides=("Splice_peptide", lambda x: (x == "+").sum()),
        internal_peptides=("Splice_peptide", lambda x: (x == "-").sum())
    ).reset_index()


    other_proteins = cds_df[~cds_df['Protein'].isin(proteins_scores_df['Protein'])]
    other_proteins = other_proteins[['Protein', 'transcript_score', 'isoforms', 'splice_sites']].drop_duplicates().copy()


    other_proteins['Prot_seq'] = other_proteins['Protein'].map(seq_dict)
    other_proteins.dropna(subset=['Prot_seq'], inplace=True)
    other_proteins['protein_length'] = other_proteins['Prot_seq'].str.len()
    other_proteins.drop(columns='Prot_seq', inplace=True)

    for col in proteins_scores_df.columns.difference(['Protein', 'protein_length', 'isoforms', 'splice_sites', 'transcript_score']):
        other_proteins[col] = 0
    all_proteins_scores_df = pd.concat([proteins_scores_df, other_proteins], ignore_index=True)

    # external_scores
    if external_scores_csv:
        external_scores_df = read_scores_csv(external_scores_csv)
    else:
        external_scores_df = pd.DataFrame(columns=['Protein', 'external_score'])
        external_scores_df['Protein'] = all_proteins_scores_df['Protein']
        external_scores_df['external_score'] = 0
    
    external_scores_df['external_score'] = pd.to_numeric(external_scores_df['external_score'], errors='coerce')
    all_proteins_scores_df = pd.merge(all_proteins_scores_df, external_scores_df, on='Protein', how='left').fillna(0)


    all_proteins_scores_df.sort_values(by='mapped_peptides', ascending=False, inplace=True)
    all_proteins_scores_df = all_proteins_scores_df.astype({'isoforms': 'Int64', 'protein_length': 'Int64', 'splice_sites': 'Int64'})
    
    all_proteins_scores_df = engineer_features(all_proteins_scores_df)
    all_proteins_scores_df.to_csv(output_dir / 'all_proteins_scores.tsv', sep='\t', index=False)
    print(f"{_ts()} [extract] Saved all_proteins_scores.tsv")

    mapped_proteins = all_proteins_scores_df[all_proteins_scores_df['mapped_peptides'] > 0]['Protein']
    unmapped_proteins = all_proteins_scores_df[all_proteins_scores_df['mapped_peptides'] == 0]['Protein']
    Nterm_proteins = all_proteins_scores_df[all_proteins_scores_df['N_terminal_peptides'] > 0]['Protein']
    no_Nterm_proteins = all_proteins_scores_df[all_proteins_scores_df['N_terminal_peptides'] == 0]['Protein']
    Cterm_proteins = all_proteins_scores_df[all_proteins_scores_df['C_terminal_peptides'] > 0]['Protein']
    no_Cterm_proteins = all_proteins_scores_df[all_proteins_scores_df['C_terminal_peptides'] == 0]['Protein']

    splice_all_supported = all_proteins_scores_df[
        (all_proteins_scores_df['splice_sites'] > 0) &
        (all_proteins_scores_df['splice_peptides'] >= all_proteins_scores_df['splice_sites'])
    ]['Protein']
    splice_some_supported = all_proteins_scores_df[
        (all_proteins_scores_df['splice_sites'] > 0) &
        (all_proteins_scores_df['splice_peptides'] > 0) &
        (all_proteins_scores_df['splice_peptides'] < all_proteins_scores_df['splice_sites'])
    ]['Protein']
    splice_no_supported = all_proteins_scores_df[
        ~(
            ((all_proteins_scores_df['splice_sites'] > 0) &
             (all_proteins_scores_df['splice_peptides'] >= all_proteins_scores_df['splice_sites'])) |
            ((all_proteins_scores_df['splice_sites'] > 0) &
             (all_proteins_scores_df['splice_peptides'] > 0) &
             (all_proteins_scores_df['splice_peptides'] < all_proteins_scores_df['splice_sites']))
        )
    ]['Protein']
    
    print(f'Number of proteins with mapped peptide: {len(mapped_proteins)}')
    print(f'Number of proteins without mapped peptide: {len(unmapped_proteins)}')
    print(f'Number of proteins with mapped Nterm peptide: {len(Nterm_proteins)}')
    print(f'Number of proteins with mapped Cterm peptide: {len(Cterm_proteins)}')
    print(f'Number of proteins with all splice sites supported: {len(splice_all_supported)}')
    print(f'Number of proteins with some splice sites supported: {len(splice_some_supported)}')
    print(f'Number of proteins with no splice support: {len(splice_no_supported)}')
    
    stacked_bars_df = pd.DataFrame({
        'bar_name': ['Mapped\npeptides', 'N-terminal\nsupport', 'C-terminal\nsupport', 'Splice-site\nsupport'],
        'supported': [
            len(mapped_proteins),
            len(Nterm_proteins),
            len(Cterm_proteins),
            len(splice_all_supported)
        ],
        'partial': [0, 0, 0, len(splice_some_supported)],
        'unsupported': [
            len(unmapped_proteins),
            len(no_Nterm_proteins),
            len(no_Cterm_proteins),
            len(splice_no_supported)
        ]
    })
    plot_mapped_percentage_bars(stacked_bars_df, output_dir)


def engineer_features(df):
    df = df.copy()
    
    # Ensure required columns exist
    required_cols = [
        "sequence_coverage", "mapped_peptides", "protein_length",
        "protein_specific_peptides", "gene_specific_peptides",
        "splice_peptides", "internal_peptides",
        "N_terminal_peptides", "C_terminal_peptides"
    ]
    
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for feature engineering: {missing}")
    
    # Normalization features (rates)
    df['mapped_peptides_per_aa'] = (
        df['mapped_peptides'] / (df['protein_length'].clip(lower=1))
    ).fillna(0)
    
    df['protein_specific_rate'] = (
        df['protein_specific_peptides'] / df['mapped_peptides'].clip(lower=1)
    ).fillna(0)
    
    df['gene_specific_rate'] = (
        df['gene_specific_peptides'] / df['mapped_peptides'].clip(lower=1)
    ).fillna(0)
    
    df['splice_rate'] = (
        df['splice_peptides'] / df['mapped_peptides'].clip(lower=1)
    ).fillna(0)
    
    df['internal_rate'] = (
        df['internal_peptides'] / df['mapped_peptides'].clip(lower=1)
    ).fillna(0)
    
    df['terminal_rate'] = (
        (df['N_terminal_peptides'] + df['C_terminal_peptides']) / 
        df['mapped_peptides'].clip(lower=1)
    ).fillna(0)
    
    # Discriminative contrasts
    df['isoform_specificity'] = (
        df['protein_specific_rate'] - df['gene_specific_rate']
    ).fillna(0)
    
    # junction_evidence: splice-heavy support
    df['junction_evidence'] = (
        df['splice_peptides'] / df['internal_peptides'].clip(lower=1)
    ).fillna(0)
    
    # Binary flags (highly predictive and stable)
    df['has_N_term'] = (df['N_terminal_peptides'] > 0).astype(int)
    df['has_C_term'] = (df['C_terminal_peptides'] > 0).astype(int)
    df['has_both_term'] = ((df['N_terminal_peptides'] > 0) & 
                            (df['C_terminal_peptides'] > 0)).astype(int)
    df['has_any_protein_specific'] = (df['protein_specific_peptides'] > 0).astype(int)
    df['has_any_splice'] = (df['splice_peptides'] > 0).astype(int)
    
    return df