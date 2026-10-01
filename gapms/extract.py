
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
    resolve_peptide_starts,
)
from .plotting import (
    plot_mapped_percentage_bars
)
from .peptides_to_genome import annotate_cds_coordinates


def _ts():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def extract_features(gtf_file, prediction_fasta, mapping_file, output_dir, external_scores_csv=None):
    gtf_file = Path(gtf_file)
    protein_fasta = Path(prediction_fasta)
    mapping_file = Path(mapping_file)
    external_scores_csv = Path(external_scores_csv) if external_scores_csv else None
    output_dir=Path(output_dir)

    print(f"{_ts()} [extract] Loading protein FASTA sequences...")
    with open(prediction_fasta) as fasta_handle:
        seq_dict = {
            record.id.split()[0]: str(record.seq).replace('I', 'L').replace('*', '')
            for record in SeqIO.parse(fasta_handle, "fasta")
        }

    print(f"{_ts()} [extract] Parsing GTF/GFF...")
    gtf_df = gtf_to_df_with_genes(gtf_file, valid_protein_ids=seq_dict)
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

    # Prepare CDS dataframe for translation-space positions. This uses inclusive
    # coordinates and the first CDS phase, including codons split across exons.
    df_CDSs = (
        cds_df.drop_duplicates(subset=['Seqid', 'Start', 'End', 'Protein'])
            .copy()
    )
    df_CDSs = annotate_cds_coordinates(df_CDSs)

    # Record true introns in translation coordinates. Adjacent CDS records are
    # not splice junctions and must not inflate the junction count.
    splice_junction_records = []
    splice_site_counts = {}
    for protein, group in df_CDSs.groupby('Protein', sort=False):
        ordered = group.sort_values('cds_nt_start')
        junction_count = 0
        rows = list(ordered.itertuples())
        for current_cds, next_cds in zip(rows, rows[1:]):
            if current_cds.Seqid != next_cds.Seqid:
                continue
            if current_cds.Strand == '+':
                has_intron = int(next_cds.Start) > int(current_cds.End) + 1
            elif current_cds.Strand == '-':
                has_intron = int(current_cds.Start) > int(next_cds.End) + 1
            else:
                has_intron = False
            if has_intron:
                junction_count += 1
                splice_junction_records.append({
                    'Protein': protein,
                    'junction_nt': int(current_cds.cds_nt_end),
                })
        splice_site_counts[protein] = junction_count

    splice_junctions_df = pd.DataFrame(
        splice_junction_records,
        columns=['Protein', 'junction_nt'],
    )
    df_CDSs['splice_sites'] = (
        df_CDSs['Protein'].map(splice_site_counts).fillna(0).astype(int)
    )

    # Merge codon positions back into original CDS dataframe
    cds_df = cds_df.merge(
        df_CDSs[['Protein', 'Seqid', 'Start', 'End', 'cds_start', 'cds_end',
                 'cds_nt_start', 'cds_nt_end', 'splice_sites']],
        on=['Protein', 'Seqid', 'Start', 'End'],
        how='left'
    )

    # Ensure correct datatype
    cds_df = cds_df.astype({'cds_start': 'Int64', 'cds_end': 'Int64'})


    # 1. Load mapping
    print(f"{_ts()} [extract] Loading mapping file...")
    try:
        raw_mapping_df = pd.read_csv(mapping_file, sep='\t', header=0, names=['Peptide', 'Protein', 'location'])
        unmapped_peptides = (
            raw_mapping_df.loc[
                raw_mapping_df['Protein'].astype(str).str.upper().eq('UNMAPPED'),
                'Peptide'
            ]
            .dropna()
            .drop_duplicates()
            .sort_values()
        )
        unmapped_peptides.to_csv(output_dir / 'Unmpped_pepides.txt', index=False, header=False)
        print(f"{_ts()} [extract] Saved {len(unmapped_peptides)} unmapped peptides to {output_dir / 'Unmpped_pepides.txt'}")
    except Exception as exc:
        print(f"{_ts()} [extract] Warning: could not export unmapped peptides from mapping file: {exc}")

    mapping_df = mapping_file_to_df(mapping_file)

    # 2. Get gene-protein specific peptides
    print(f"{_ts()} [extract] Resolving gene/protein-specific peptides...")
    pep_df = get_gene_protein_specific_peps(cds_df, mapping_df)

    # 3. Map protein sequences
    pep_df['Prot_seq'] = pep_df['Protein'].map(seq_dict)
    # print(f"Peptides after mapping protein sequences: {len(pep_df)}")
    pep_df.dropna(subset=['Prot_seq'], inplace=True)
    # print(f"Peptides after filtering out rows with missing protein sequences: {len(pep_df)}")

    # 5. Compute peptide positions
    print(f"{_ts()} [extract] Computing peptide positions for {len(pep_df)} pep-protein pairs...")
    positioned_rows = []
    for row_index, row in pep_df.iterrows():
        starts = resolve_peptide_starts(
            row['Prot_seq'],
            row['Peptide'],
            row.get('location'),
        )
        for start in starts:
            positioned_rows.append((row_index, start))

    if positioned_rows:
        source_indices, pep_starts = zip(*positioned_rows)
        pep_df = pep_df.loc[list(source_indices)].copy().reset_index(drop=True)
        pep_df['pep_start'] = list(pep_starts)
        pep_df['pep_end'] = pep_df['pep_start'] + pep_df['Peptide'].str.len() - 1
    else:
        pep_df = pep_df.iloc[0:0].copy()
        pep_df['pep_start'] = pd.Series(dtype='Int64')
        pep_df['pep_end'] = pd.Series(dtype='Int64')
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
        'Gene_specific', 'cds_start', 'cds_end', 'cds_nt_start', 'cds_nt_end', 'peptide_length',
        'pep_start', 'pep_end', 'N_terminal_peptide', 'transcript_score', 'prot_len', 'sequence_coverage'
    ]]

    print(f"{_ts()} [extract] Computing splice-site flags on {len(peptide_features_df)} rows...")
    pep_unique = peptide_features_df[
        ['Protein', 'Peptide', 'pep_start', 'pep_end']
    ].drop_duplicates()
    pep_unique['pep_nt_start'] = (pep_unique['pep_start'] - 1) * 3 + 1
    pep_unique['pep_nt_end'] = pep_unique['pep_end'] * 3

    junction_check = pep_unique.merge(splice_junctions_df, on='Protein', how='left')
    junction_position = pd.to_numeric(
        junction_check['junction_nt'], errors='coerce'
    )
    junction_check['spans_junction'] = (
        junction_position.notna()
        & (junction_check['pep_nt_start'] <= junction_position.fillna(-1))
        & (junction_check['pep_nt_end'] > junction_position.fillna(-1))
    )
    spliced_occurrences = (
        junction_check[junction_check['spans_junction']]
        .drop_duplicates(subset=['Protein', 'Peptide', 'pep_start', 'pep_end'])
        [['Protein', 'Peptide', 'pep_start', 'pep_end']]
    )
    spliced_occurrences['_spliced'] = True
    peptide_features_df = peptide_features_df.merge(
        spliced_occurrences,
        on=['Protein', 'Peptide', 'pep_start', 'pep_end'],
        how='left',
    )
    peptide_features_df['Splice_peptide'] = np.where(
        peptide_features_df['_spliced'].eq(True), '+', '-'
    )
    peptide_features_df.drop(columns=['_spliced'], inplace=True)

    # Credit a junction only when the peptide supports that same junction at
    # every occurrence within the protein. This prevents several peptides at
    # one junction from masquerading as support for every junction.
    occurrence_junctions = {}
    for key, group in junction_check[junction_check['spans_junction']].groupby(
        ['Protein', 'Peptide', 'pep_start', 'pep_end'], sort=False
    ):
        occurrence_junctions[key] = set(group['junction_nt'].dropna().astype(int))

    supported_junctions_by_protein = {}
    for (protein, peptide), occurrences in pep_unique.groupby(
        ['Protein', 'Peptide'], sort=False
    ):
        junction_sets = [
            occurrence_junctions.get(
                (protein, peptide, row.pep_start, row.pep_end),
                set(),
            )
            for row in occurrences.itertuples()
        ]
        common_junctions = set.intersection(*junction_sets) if junction_sets else set()
        supported_junctions_by_protein.setdefault(protein, set()).update(common_junctions)

    peptide_features_df['C_terminal_peptide'] = np.where(peptide_features_df['pep_end'] == peptide_features_df['prot_len'], '+', '-')
    peptide_features_df = peptide_features_df.drop(columns=['cds_start', 'cds_end']).drop_duplicates()

    # The CDS merge above intentionally creates one row per CDS so a peptide
    # can be tested against every coding interval.  Protein evidence, however,
    # is a count of distinct peptide sequences, not peptide × CDS rows.  Collapse
    # repeated locations and CDS rows before summing the evidence features.
    # Terminal/junction support is only retained when every occurrence of that
    # peptide in the protein supports the same interpretation.
    def _unambiguous_support(values):
        return '+' if values.eq('+').all() else '-'

    peptide_features_df = (
        peptide_features_df
        .groupby(['Protein', 'Peptide'], as_index=False, sort=False)
        .agg(
            Gene=('Gene', 'first'),
            isoforms=('isoforms', 'mean'),
            splice_sites=('splice_sites', 'mean'),
            Protein_specific=('Protein_specific', 'first'),
            Gene_specific=('Gene_specific', 'first'),
            N_terminal_peptide=('N_terminal_peptide', _unambiguous_support),
            C_terminal_peptide=('C_terminal_peptide', _unambiguous_support),
            Splice_peptide=('Splice_peptide', _unambiguous_support),
            transcript_score=('transcript_score', 'mean'),
            prot_len=('prot_len', 'mean'),
            sequence_coverage=('sequence_coverage', 'mean'),
        )
    )


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
    proteins_scores_df['supported_splice_sites'] = (
        proteins_scores_df['Protein']
        .map(lambda protein: len(supported_junctions_by_protein.get(protein, set())))
        .astype(int)
    )

    evidence_count_columns = [
        'N_terminal_peptides', 'C_terminal_peptides',
        'protein_specific_peptides', 'gene_specific_peptides',
        'splice_peptides', 'internal_peptides',
    ]
    if proteins_scores_df[evidence_count_columns].gt(
        proteins_scores_df['mapped_peptides'], axis=0
    ).any().any():
        raise RuntimeError("Protein evidence counts exceeded the number of distinct mapped peptides")
    if proteins_scores_df['supported_splice_sites'].gt(proteins_scores_df['splice_sites']).any():
        raise RuntimeError("Supported splice-site counts exceeded annotated splice sites")


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
        (all_proteins_scores_df['supported_splice_sites'] == all_proteins_scores_df['splice_sites'])
    ]['Protein']
    splice_some_supported = all_proteins_scores_df[
        (all_proteins_scores_df['splice_sites'] > 0) &
        (all_proteins_scores_df['supported_splice_sites'] > 0) &
        (all_proteins_scores_df['supported_splice_sites'] < all_proteins_scores_df['splice_sites'])
    ]['Protein']
    splice_no_supported = all_proteins_scores_df[
        (all_proteins_scores_df['splice_sites'] == 0)
        | (all_proteins_scores_df['supported_splice_sites'] == 0)
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
