from pathlib import Path
import pandas as pd
import numpy as np

from .gtf_utils import gtf_to_df_with_genes, save_gtf_subset
from .model_utils import (
    get_high_confident_proteins,
    get_low_confident_proteins,
    run_protein_classifier,
)
from .score_filter import find_apply_score_filter
from .plotting import (
    plot_roc_curve,
    plot_shap_summary,
    plot_mapped_percentage_bars,
    plot_sequence_coverage_groups_hist,
    plot_external_scores
)


def filter_predictions(
    gtf_file,
    protein_fasta,
    output_dir,
    plot_external_scores_enabled=False,
    use_iterative_training=False,
):
    all_scores_path = output_dir / "all_proteins_scores.tsv"
    all_scores_df = pd.read_csv(all_scores_path, sep='\t')

    high_confident_df = get_high_confident_proteins(all_scores_df)
    high_confident_proteins = set(high_confident_df['Protein'].unique())

    low_confident_df = get_low_confident_proteins(all_scores_df)
    low_confident_proteins = set(low_confident_df['Protein'].unique())

    # Step 2: Run protein classifier
    labeled_df, model_mode = run_protein_classifier(
        df=all_scores_df,
        high_confident_df=high_confident_df,
        low_confident_df=low_confident_df,
        use_iterative_training=use_iterative_training,
        pos_thr=0.90,
        neg_thr=0.10,
        n_iter=5,
        shap_output_dir=output_dir,
        plot_shap=plot_shap_summary,
    )
    print(f"Protein classifier mode used: {model_mode}")

    xgb_pos_proteins = set(labeled_df[labeled_df["final_label"] == 1]["Protein"])
    
    
    # Step 3: Apply ROC-based score filtering only when user provided external scores
    score_supported = set()
    if plot_external_scores_enabled:
        external_scores_cutoff, external_scores_supported = find_apply_score_filter(
            all_scores_df,
            high_confident_proteins,
            'external_score',
            output_dir=output_dir,
            plot_roc=plot_roc_curve,
        )
        score_supported = external_scores_supported

    print(f"Number of supported proteins by scores: {len(score_supported)}")


    # Collect stats
    all_proteins = set(all_scores_df["Protein"])
    original_supported = high_confident_proteins | xgb_pos_proteins | score_supported

    # Require at least one gene-specific peptide for a protein to be considered supported
    proteins_with_gene_specific_pep = set(
        all_scores_df.loc[all_scores_df.get('gene_specific_peptides', 0) > 0, 'Protein'].unique()
    )

    # Apply the requirement by intersecting with proteins that have ≥1 gene-specific peptide
    supported_proteins = original_supported & proteins_with_gene_specific_pep
    removed_no_gene_specific = original_supported - supported_proteins

    unsupported_proteins = all_proteins - supported_proteins

    print(f"\nPipeline completed\nOutputs written to: {output_dir}")
    print(f"\nNumber of all proteins = {len(set(all_scores_df['Protein']))}")
    print(f"Number of high confident proteins = {len(high_confident_proteins)}")
    print(f"Number of low confident proteins = {len(low_confident_proteins)}")
    print(f"Number of all supported proteins = {len(supported_proteins)}")
    print(f"Number of all un-supported proteins = {len(unsupported_proteins)}")
    
    # Step 7: Save protein sets
    txt_dir = output_dir / "Txt"
    txt_dir.mkdir(parents=True, exist_ok=True)

    def save_protein_list(protein_set, filename):
        with open(txt_dir / filename, "w") as f:
            for p in sorted(protein_set):
                f.write(f"{p}\n")
    save_protein_list(all_proteins, "all_proteins.txt")
    save_protein_list(high_confident_proteins, "high_confident_proteins.txt")
    save_protein_list(low_confident_proteins, "low_confident_proteins.txt")
    save_protein_list(supported_proteins, "supported_proteins.txt")
    save_protein_list(unsupported_proteins, "unsupported_proteins.txt")
    
    # Step 8: Save subsets 
    gtf_file = Path(gtf_file)
    gtf_df = gtf_to_df_with_genes(gtf_file)
    save_gtf_subset(gtf_df, supported_proteins, output_dir, 'supported_proteins.gtf')


    all_scores_df['supported'] = np.where(all_scores_df['Protein'].isin(supported_proteins), '+', '-')
    all_scores_df['confidence'] = np.select(
        [
            all_scores_df['Protein'].isin(high_confident_proteins),
            all_scores_df['Protein'].isin(supported_proteins),
        ],
        [
            'high',
            'moderate',
        ],
        default='no-support'
    )

    ordered_cols = ['Protein', 'supported', 'confidence']
    remaining_cols = [col for col in all_scores_df.columns if col not in ordered_cols]
    all_scores_df = all_scores_df[ordered_cols + remaining_cols]
    all_scores_df['_supported_sort'] = all_scores_df['supported'].map({'+': 0, '-': 1}).fillna(2)
    all_scores_df = all_scores_df.sort_values(
        by=['_supported_sort', 'confidence', 'mapped_peptides', 'Protein'],
        ascending=[True, True, False, True]
    ).drop(columns=['_supported_sort']).reset_index(drop=True)

     # Save updated all_scores_df
    all_scores_df.to_csv(all_scores_path, sep='\t', index=False)
    plot_mapped_percentage_bars(all_scores_df, output_dir)
    plot_sequence_coverage_groups_hist(all_scores_df, output_dir)
    if plot_external_scores_enabled:
        plot_external_scores(all_scores_df, output_dir)


    # Step 8: Filter FASTA
    protein_fasta = Path(protein_fasta)
    output_fasta = output_dir / "supported_proteins.fa"
    with open(protein_fasta) as fasta, open(output_fasta, "w") as out:
        header, seq = "", ""
        for line in fasta:
            if line.startswith(">"):
                if header and header[1:] in supported_proteins:
                    out.write(header + "\n" + seq + "\n")
                header, seq = line.strip(), ""
            else:
                seq += line.strip()
        if header and header[1:] in supported_proteins:
            out.write(header + "\n" + seq + "\n")
    if removed_no_gene_specific:
        print(f"Removed {len(removed_no_gene_specific)} proteins lacking gene-specific peptides")