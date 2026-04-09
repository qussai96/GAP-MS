from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def _get_figures_dir(output_dir):
    output_dir = Path(output_dir)
    figures_dir = output_dir / 'Figures'
    figures_dir.mkdir(parents=True, exist_ok=True)
    return output_dir, figures_dir


def plot_mapped_percentage_bars(df, output_dir):
    """
    Create a combined figure with:
      1) peptide evidence distribution across all proteins
      2) GAP-MS supported proteins within the same evidence categories

    Accepted inputs:
      - a summary DataFrame with columns: `bar_name`, `supported`, `partial`, `unsupported`
      - or the full `all_proteins_scores.tsv`-style DataFrame with peptide feature columns and
        an optional `supported` column.
    """
    df = df.copy()

    def build_summary_from_scores(scores_df):
        splice_all_mask = (
            (scores_df['splice_sites'] > 0) &
            (scores_df['splice_peptides'] >= scores_df['splice_sites'])
        )
        splice_some_mask = (
            (scores_df['splice_sites'] > 0) &
            (scores_df['splice_peptides'] > 0) &
            (scores_df['splice_peptides'] < scores_df['splice_sites'])
        )
        splice_none_mask = ~(splice_all_mask | splice_some_mask)

        evidence_df = pd.DataFrame({
            'bar_name': ['Mapped\npeptides', 'N-terminal\nsupport', 'C-terminal\nsupport', 'Splice-site\nsupport'],
            'supported': [
                int((scores_df['mapped_peptides'] > 0).sum()),
                int((scores_df['N_terminal_peptides'] > 0).sum()),
                int((scores_df['C_terminal_peptides'] > 0).sum()),
                int(splice_all_mask.sum())
            ],
            'partial': [0, 0, 0, int(splice_some_mask.sum())],
            'unsupported': [
                int((scores_df['mapped_peptides'] == 0).sum()),
                int((scores_df['N_terminal_peptides'] == 0).sum()),
                int((scores_df['C_terminal_peptides'] == 0).sum()),
                int(splice_none_mask.sum())
            ]
        })

        supported_summary_df = None
        if 'supported' in scores_df.columns:
            gapms_mask = scores_df['supported'] == 'Yes'
            total_proteins = len(scores_df)
            supported_summary_df = pd.DataFrame({
                'bar_name': evidence_df['bar_name'],
                'count': [
                    int((gapms_mask & (scores_df['mapped_peptides'] > 0)).sum()),
                    int((gapms_mask & (scores_df['N_terminal_peptides'] > 0)).sum()),
                    int((gapms_mask & (scores_df['C_terminal_peptides'] > 0)).sum()),
                    int((gapms_mask & (splice_all_mask | splice_some_mask)).sum())
                ],
                'total': [total_proteins] * len(evidence_df)
            })
        return evidence_df, supported_summary_df

    summary_cols = {'bar_name', 'supported', 'partial', 'unsupported'}
    score_cols = {'mapped_peptides', 'N_terminal_peptides', 'C_terminal_peptides', 'splice_sites', 'splice_peptides'}

    if summary_cols.issubset(df.columns):
        evidence_df = df.copy()
        supported_summary_df = None
    elif score_cols.issubset(df.columns):
        evidence_df, supported_summary_df = build_summary_from_scores(df)
    else:
        missing_cols = sorted(summary_cols - set(df.columns))
        raise ValueError(
            'Input to plot_mapped_percentage_bars must be either a summary DataFrame '
            f'or a GAP-MS scores DataFrame. Missing summary columns: {missing_cols}'
        )

    if {'mapped', 'unmapped'}.issubset(evidence_df.columns):
        evidence_df = evidence_df.rename(columns={'mapped': 'supported', 'unmapped': 'unsupported'})
    if 'partial' not in evidence_df.columns:
        evidence_df['partial'] = 0

    value_cols = ['supported', 'partial', 'unsupported']
    evidence_df[value_cols] = evidence_df[value_cols].fillna(0)
    totals = evidence_df[value_cols].sum(axis=1).replace(0, 1)

    colors = {
        'supported': '#2E8B57',
        'partial': '#F4A261',
        'unsupported': '#B0B0B0'
    }
    legend_labels = {
        'supported': 'Mapped',
        'partial': 'Partially mapped',
        'unsupported': 'No mapped peptides'
    }

    has_supported_panel = supported_summary_df is not None
    if has_supported_panel:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5.8), gridspec_kw={'width_ratios': [1.2, 1.0], 'wspace': 0.32})
        ax = axes[0]
        ax_supported = axes[1]
    else:
        fig, ax = plt.subplots(figsize=(9.5, 5.8))
        ax_supported = None

    fig.suptitle('Peptide support summary across proteins', fontsize=13, fontweight='bold', y=0.97)

    x = np.arange(len(evidence_df))
    bottoms = np.zeros(len(evidence_df))
    for col in value_cols:
        heights = evidence_df[col].to_numpy()
        bars = ax.bar(
            x,
            heights,
            bottom=bottoms,
            color=colors[col],
            edgecolor='black',
            width=0.7,
            label=legend_labels[col]
        )

        for i, (bar, height, total) in enumerate(zip(bars, heights, totals)):
            if height <= 0:
                continue
            pct = (height / total) * 100 if total else 0
            label = f"{int(height)} ({pct:.1f}%)"
            text_color = 'white' if height >= max(total * 0.08, 4) else 'black'
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bottoms[i] + height / 2,
                label,
                ha='center',
                va='center',
                fontsize=9,
                color=text_color,
                fontweight='bold'
            )

        bottoms += heights

    ax.set_xticks(x)
    ax.set_xticklabels(evidence_df['bar_name'])
    ax.set_ylabel('Number of proteins')
    ax.set_title('Peptide evidence across all proteins', fontweight='bold', fontsize=11)
    ax.grid(axis='y', linestyle='--', alpha=0.35)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper left', frameon=False, fontsize=8, handlelength=1.4)

    if has_supported_panel:
        supported_summary_df = supported_summary_df.copy()
        supported_summary_df['percent'] = np.where(
            supported_summary_df['total'] > 0,
            (supported_summary_df['count'] / supported_summary_df['total']) * 100,
            0
        )
        bar_colors = ['#4C78A8', '#72B7B2', '#E45756', '#B279A2']
        bars = ax_supported.bar(
            np.arange(len(supported_summary_df)),
            supported_summary_df['count'],
            color=bar_colors,
            edgecolor='black',
            width=0.7
        )
        for bar, count, pct in zip(bars, supported_summary_df['count'], supported_summary_df['percent']):
            ax_supported.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(supported_summary_df['count'].max() * 0.015, 0.5),
                f"{int(count)} ({pct:.1f}%)",
                ha='center',
                va='bottom',
                fontsize=9,
                fontweight='bold'
            )

        ax_supported.set_xticks(np.arange(len(supported_summary_df)))
        ax_supported.set_xticklabels(supported_summary_df['bar_name'])
        ax_supported.set_ylabel('Number of GAP-MS supported proteins')
        ax_supported.set_title('GAP-MS supported proteins by category', fontweight='bold', fontsize=11)
        ax_supported.grid(axis='y', linestyle='--', alpha=0.35)
        ax_supported.spines['top'].set_visible(False)
        ax_supported.spines['right'].set_visible(False)


    output_dir, figures_dir = _get_figures_dir(output_dir)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.16, top=0.86, wspace=0.32)
    plt.savefig(figures_dir / 'Peptides_mapping_bars.png', dpi=300, bbox_inches='tight')
    plt.close()




def plot_sequence_coverage_groups_hist(df, output_dir):
    _, figures_dir = _get_figures_dir(output_dir)
    plt.figure(figsize=(6, 4))
    plt.hist(
        df[df["supported"] == 'Yes']["sequence_coverage"],
        bins=20, density=True, color="lightgray", alpha=0.8, label="Supported"
    )
    plt.hist(
        df[df["unsupported"] == 'Yes']["sequence_coverage"],
        bins=20, density=True, color="dimgray", alpha=0.8, label="Not supported"
    )
    plt.hist(
        df[df["high_conf"] == 'Yes']["sequence_coverage"],
        bins=20, density=True, color="#ADD8E6", alpha=0.8, label="High confident"
    )
    plt.hist(
        df[df["low_conf"] == 'Yes']["sequence_coverage"],
        bins=20, density=True, color="#D8BFD8", alpha=0.8, label="Low confident"
    )
    plt.yscale("log")
    plt.xlabel("Sequence coverage")
    plt.ylabel("Density")
    plt.title("Sequence Coverage by Confidence Groups")
    plt.legend()
    plt.savefig(figures_dir / 'Sequence_coverage_groups.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_shap_summary(model, X_labeled, output_dir):
    _, figures_dir = _get_figures_dir(output_dir)
    import shap
    explainer = shap.Explainer(model)
    shap_values = explainer(X_labeled)
    shap.summary_plot(shap_values, X_labeled, show=False)
    plt.savefig(figures_dir / 'Shap_summary_plot.png', bbox_inches='tight')
    plt.close()


def plot_external_scores(df, output_dir):
    if 'external_score' not in df.columns or df['external_score'].dropna().empty:
        return

    _, figures_dir = _get_figures_dir(output_dir)
    plt.figure(figsize=(6, 4))
    plt.hist(
        df["external_score"],
        bins=20, density=True, color="navy", alpha=0.7, label="All proteins"
    )
    plt.hist(
        df[df["supported"] == 'Yes']["external_score"],
        bins=20, density=True, color="mediumseagreen", alpha=0.8, label="Supported proteins"
    )
    plt.yscale("log")
    plt.xlabel("External score")
    plt.ylabel("Density")
    plt.title("Distribution of user-provided external scores")
    plt.legend()
    plt.savefig(figures_dir / 'External_scores.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_roc_curve(fpr, tpr, cutoff, fpr_cut, tpr_cut, output_dir, score_column):
    _, figures_dir = _get_figures_dir(output_dir)
    plt.figure(figsize=(6, 4))
    plt.plot(fpr, tpr, marker='o', markersize=3, linestyle='-', color='purple', alpha=0.9, linewidth=2)
    plt.axvline(x=fpr_cut, color='red', linestyle='-', linewidth=2)
    plt.scatter([fpr_cut], [tpr_cut], color='red', s=100, zorder=5, marker='o')
    plt.text(fpr_cut + 0.02, tpr_cut - 0.05, f"score = {cutoff:.3f}", fontsize=11, color="black")
    plt.xlabel("False positive rate", fontsize=12)
    plt.ylabel("True positive rate", fontsize=12)
    plt.title("Score Cutoff", fontsize=13)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.savefig(figures_dir / f'ROC_curve_{score_column}.png', dpi=300, bbox_inches='tight')
    plt.close()

