"""Original GAPMS feature engineering and high/low seed rules (MIT)."""
import pandas as pd

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

BASE_FEATURE_COLS = [
    "sequence_coverage",
    "protein_specific_peptides",
    "gene_specific_peptides",
    "splice_peptides",
    "internal_peptides",
    "mapped_peptides",
    "N_terminal_peptides",
    "C_terminal_peptides",
    "protein_length",
]

ENGINEERED_FEATURE_COLS = [
    "mapped_peptides_per_aa",
    "protein_specific_rate",
    "gene_specific_rate",
    "splice_rate",
    "internal_rate",
    "terminal_rate",
    "isoform_specificity",
    "junction_evidence",
    "has_N_term",
    "has_C_term",
    "has_both_term",
    "has_any_protein_specific",
    "has_any_splice",
]

ITERATIVE_FEATURE_COLS = BASE_FEATURE_COLS + ENGINEERED_FEATURE_COLS

def get_high_confident_proteins(df):
        high_confident = []
        for _, group in df.groupby("Protein"):
            for _, row in group.iterrows():
                cond1 = row["protein_specific_peptides"] >= 2
                cond2 = row["sequence_coverage"] >= 0.8
                cond3 = row["N_terminal_peptides"] >= 1 and row["C_terminal_peptides"] >= 1 and row["gene_specific_peptides"] > 0
                supported_splice_sites = row.get("supported_splice_sites", 0)
                cond4 = (
                    supported_splice_sites == row["splice_sites"]
                    and supported_splice_sites > 0
                    and row["gene_specific_peptides"] > 0
                )
                if cond1 or cond2 or cond3 or cond4:
                    high_confident.append(row)
                    break
        print(f"Number of high confident proteins: {len(high_confident)}")
        if not high_confident:
            return df.iloc[0:0].copy()
        return pd.DataFrame(high_confident)

def get_low_confident_proteins(df):
        low_confident = []
        for _, group in df.groupby("Protein"):
            for _, row in group.iterrows():
                cond1 = row["mapped_peptides"] < 2
                cond2 = row["gene_specific_peptides"] == 0
                cond3 = row["protein_length"] < 200 or row["protein_length"] > 500
                if cond1 and cond2 and cond3:
                    low_confident.append(row)
                    break
        print(f"Number of low confident proteins: {len(low_confident)}")
        if not low_confident:
            return df.iloc[0:0].copy()
        return pd.DataFrame(low_confident)
