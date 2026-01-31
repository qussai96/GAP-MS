from pathlib import Path
import pandas as pd
import numpy as np

from .gtf_utils import gtf_to_df_with_genes, save_gtf_subset
from .model_utils import (
    get_high_confident_proteins,
    get_low_confident_proteins,
    find_apply_score_filter,
    train_iterative_model
)
from .plotting import (
    plot_roc_curve,
    plot_shap_summary,
    plot_supported_percentage_bar,
    plot_protein_coverage_groups_hist,
    plot_external_scores
)


def filter_predictions(gtf_file, protein_fasta, output_dir,):    
    all_scores_path = output_dir / "all_proteins_scores.tsv"
    all_scores_df = pd.read_csv(all_scores_path, sep='\t')
    
    # Step 4: Find highly supported proteins
    high_confident_df = get_high_confident_proteins(all_scores_df)
    high_confident_proteins = set(high_confident_df['Protein'].unique())


    low_confident_df = get_low_confident_proteins(all_scores_df)
    low_confident_proteins = set(low_confident_df['Protein'].unique())



    # Step 6: Train iterative self-labeling model
    labeled_df = train_iterative_model(
        df=all_scores_df,
        high_confident_df=high_confident_df,
        low_confident_df=low_confident_df,
        pos_thr=0.90,
        neg_thr=0.10,
        n_iter=5,
        shap_output_dir=output_dir,
        plot_shap=plot_shap_summary
    )

    xgb_pos_proteins = set(labeled_df[labeled_df["final_label"] == 1]["Protein"])
    
    
    # Step 5: Apply ROC-based score filtering
    external_scores_cutoff, external_scores_supported = find_apply_score_filter(
        all_scores_df, high_confident_proteins, 'external_score',
        output_dir=output_dir,
        plot_roc=plot_roc_curve)
    
    
    score_supported = external_scores_supported
    # score_supported = set()


    print(f"Number of supported proteins by scores: {len(score_supported)}")


    # Collect stats
    all_proteins = set(all_scores_df["Protein"])
    supported_proteins = high_confident_proteins | xgb_pos_proteins | score_supported
    
    unsupported_proteins = all_proteins - supported_proteins

    plot_supported_percentage_bar(len(supported_proteins), len(unsupported_proteins), output_dir)

    print(f"\nPipeline completed\nOutputs written to: {output_dir}")
    print(f"\nNumber of all proteins = {len(set(all_scores_df['Protein']))}")
    print(f"Number of high confident proteins = {len(high_confident_proteins)}")
    print(f"Number of low confident proteins = {len(low_confident_proteins)}")
    print(f"Number of all supported proteins = {len(supported_proteins)}")
    print(f"Number of all un-supported proteins = {len(unsupported_proteins)}")

    # print(f"Number of supported N-terminal proteins = {len(set(all_scores_df['Protein'][all_scores_df['N_terminal_peptides'] > 0]))}")
    # print(f"Number of supported C-terminal proteins = {len(set(all_scores_df['Protein'][all_scores_df['C_terminal_peptides'] > 0]))}")
    # print(f"Number of supported splice junctions = {len(set(all_scores_df['Protein'][all_scores_df['splice_peptides'] > 0]))}")
    

    
    # Step 7: Save protein sets
    def save_protein_list(protein_set, filename):
        with open(output_dir / "Txt" / filename, "w") as f:
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


    all_scores_df['supported'] = np.where(all_scores_df['Protein'].isin(supported_proteins), 'Yes', 'No')
    all_scores_df['high_conf'] = np.where(all_scores_df['Protein'].isin(high_confident_proteins), 'Yes', 'No')
    all_scores_df['low_conf'] = np.where(all_scores_df['Protein'].isin(low_confident_proteins), 'Yes', 'No')
    all_scores_df['unsupported'] = np.where(all_scores_df['Protein'].isin(unsupported_proteins), 'Yes', 'No')

     # Save updated all_scores_df
    plot_protein_coverage_groups_hist(all_scores_df, output_dir)
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