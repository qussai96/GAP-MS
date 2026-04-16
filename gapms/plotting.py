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
    """Create a single stacked panel for peptide mapping support categories."""
    df = df.copy()

    def build_summary_from_scores(scores_df):
        total_proteins = len(scores_df)
        gapms_mask = scores_df['supported'].isin(['Yes', '+']) if 'supported' in scores_df.columns else pd.Series(False, index=scores_df.index)

        evidence_masks = {
            'Mapped\npeptides': scores_df['mapped_peptides'] > 0,
            'N-terminal\nsupport': scores_df['N_terminal_peptides'] > 0,
            'C-terminal\nsupport': scores_df['C_terminal_peptides'] > 0,
            # Partial splice support is intentionally collapsed into one evidence category.
            'Splice-site\nsupport': scores_df['splice_peptides'] > 0,
        }

        summary_rows = []
        for bar_name, evidence_mask in evidence_masks.items():
            evidence_count = int(evidence_mask.sum())
            gapms_supported = int((evidence_mask & gapms_mask).sum())
            peptide_only = max(evidence_count - gapms_supported, 0)
            no_mapped = max(total_proteins - evidence_count, 0)
            summary_rows.append({
                'bar_name': bar_name,
                'gapms_supported': gapms_supported,
                'peptide_evidence': peptide_only,
                'no_mapped_peptides': no_mapped,
                'total': total_proteins,
            })

        return pd.DataFrame(summary_rows)

    def build_summary_from_legacy(legacy_df):
        # Legacy summary input does not include GAP-MS support counts, so blue stays at zero.
        legacy_df = legacy_df.copy()
        if {'mapped', 'unmapped'}.issubset(legacy_df.columns):
            legacy_df = legacy_df.rename(columns={'mapped': 'supported', 'unmapped': 'unsupported'})
        if 'partial' not in legacy_df.columns:
            legacy_df['partial'] = 0

        summary_rows = []
        for _, row in legacy_df.iterrows():
            evidence_count = int(row.get('supported', 0)) + int(row.get('partial', 0))
            no_mapped = int(row.get('unsupported', 0))
            total = evidence_count + no_mapped
            summary_rows.append({
                'bar_name': row['bar_name'],
                'gapms_supported': 0,
                'peptide_evidence': evidence_count,
                'no_mapped_peptides': no_mapped,
                'total': total,
            })

        return pd.DataFrame(summary_rows)

    summary_cols = {'bar_name', 'supported', 'partial', 'unsupported'}
    score_cols = {'mapped_peptides', 'N_terminal_peptides', 'C_terminal_peptides', 'splice_sites', 'splice_peptides'}

    if summary_cols.issubset(df.columns) or {'bar_name', 'supported', 'unsupported'}.issubset(df.columns):
        evidence_df = build_summary_from_legacy(df)
    elif score_cols.issubset(df.columns):
        evidence_df = build_summary_from_scores(df)
    else:
        missing_cols = sorted(summary_cols - set(df.columns))
        raise ValueError(
            'Input to plot_mapped_percentage_bars must be either a summary DataFrame '
            f'or a GAP-MS scores DataFrame. Missing summary columns: {missing_cols}'
        )

    value_cols = ['gapms_supported', 'peptide_evidence', 'no_mapped_peptides']
    evidence_df[value_cols] = evidence_df[value_cols].fillna(0)
    totals = evidence_df['total'].replace(0, 1)

    colors = {
        'gapms_supported': '#4C78A8',
        'peptide_evidence': '#2E8B57',
        'no_mapped_peptides': '#B0B0B0'
    }
    text_colors = {
        'gapms_supported': 'white',
        'peptide_evidence': 'black',
        'no_mapped_peptides': 'red'
    }
    legend_labels = {
        'gapms_supported': 'Supported by GAP-MS',
        'peptide_evidence': 'Mapped peptides',
        'no_mapped_peptides': 'No mapped peptides'
    }

    fig, ax = plt.subplots(figsize=(10.2, 6.0))
    fig.suptitle('Peptide mapping summary across proteins', fontsize=13, fontweight='bold', y=0.97)

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
            text_color = text_colors.get(col, 'black')
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
    ax.grid(axis='y', linestyle='--', alpha=0.35)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(
        loc='lower center',
        bbox_to_anchor=(0.5, 1.02),  # center above plot
        ncol=max(1, len(legend_labels)),  # keep legend horizontal; avoid zero columns
        frameon=False,
        fontsize=8.5,
        handlelength=1.4
    )


    output_dir, figures_dir = _get_figures_dir(output_dir)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.86)
    plt.savefig(figures_dir / 'Peptides_mapping_bars.png', dpi=300, bbox_inches='tight')
    plt.close()


def _load_annotation_summary_counts(compare_dir):
    """Load counts for novel and peptide-supported difference categories."""
    compare_dir = Path(compare_dir)
    counts = {
        'novel': 0,
        'peptide_support_different_splice': 0,
        'peptide_support_different_start': 0,
        'peptide_support_different_stop': 0,
    }

    novel_scores = compare_dir / 'Novel' / 'new_predicted_proteins_scores.tsv'
    if novel_scores.exists() and novel_scores.stat().st_size > 0:
        try:
            novel_df = pd.read_csv(novel_scores, sep='\t')
            counts['novel'] = int(novel_df['Protein'].nunique()) if 'Protein' in novel_df.columns else int(len(novel_df))
        except pd.errors.EmptyDataError:
            pass

    summary_path = compare_dir / 'annotation_comparison_summary.tsv'
    if summary_path.exists() and summary_path.stat().st_size > 0:
        try:
            summary_df = pd.read_csv(summary_path, sep='\t')
            if {'Category', 'Count'}.issubset(summary_df.columns):
                for _, row in summary_df.iterrows():
                    category = str(row['Category']).strip()
                    if category in counts:
                        counts[category] = int(row['Count'])
        except pd.errors.EmptyDataError:
            pass

    return counts
def _plot_reference_count_bars(ax, compare_dir, title):
    """Plot reference comparison category counts as a compact labeled bar chart."""
    counts = _load_annotation_summary_counts(compare_dir)
    categories = [
        'novel',
        'peptide_support_different_splice',
        'peptide_support_different_start',
        'peptide_support_different_stop',
    ]
    labels = ['Novel', 'Diff splice', 'Diff start', 'Diff stop']
    values = [counts[key] for key in categories]
    colors = ['#4C78A8', '#F58518', '#54A24B', '#E45756']

    bars = ax.bar(labels, values, color=colors, edgecolor='black', width=0.68)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(max(values) * 0.02 if max(values) else 0.2, 0.2),
            str(int(value)),
            ha='center',
            va='bottom',
            fontsize=10,
            fontweight='bold'
        )

    ax.set_title(title, fontweight='bold', fontsize=11)
    ax.set_ylabel('Protein count')
    ax.grid(axis='y', linestyle='--', alpha=0.35)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, max(values) * 1.18 + 1 if any(values) else 1)


def _draw_existing_figure(ax, image_path, title):
    """Render an existing PNG into a panel, or show a placeholder if missing."""
    image_path = Path(image_path)
    ax.set_title(title, fontweight='bold', fontsize=11)
    ax.axis('off')

    if not image_path.exists():
        ax.text(0.5, 0.5, f'Missing figure:\n{image_path.name}', ha='center', va='center', fontsize=11, transform=ax.transAxes)
        return

    image = plt.imread(image_path)
    ax.imshow(image)
    ax.set_aspect('auto')


def plot_parent_run_summary(output_dir):
    """Create a single parent-level summary figure for combined prediction+BAM runs."""
    output_dir = Path(output_dir)
    prediction_dir = output_dir / 'prediction_search'
    bam_dir = output_dir / 'bam_search'
    comparisons_dir = output_dir / 'comparisons'

    fig, axes = plt.subplots(
        5,
        1,
        figsize=(14, 30),
        gridspec_kw={'height_ratios': [1.45, 1.45, 1.0, 1.0, 1.45]}
    )
    fig.suptitle('GAP-MS combined summary', fontsize=16, fontweight='bold', y=0.995)

    _draw_existing_figure(
        axes[0],
        prediction_dir / 'Figures' / 'Peptides_mapping_bars.png',
        '1. prediction_search peptide mapping summary'
    )
    _draw_existing_figure(
        axes[1],
        bam_dir / 'Figures' / 'Peptides_mapping_bars.png',
        '2. bam_search peptide mapping summary'
    )
    _plot_reference_count_bars(
        axes[2],
        prediction_dir / 'Compare_to_Reference',
        '3. prediction_search vs reference (novel and peptide-supported differences)'
    )
    _plot_reference_count_bars(
        axes[3],
        bam_dir / 'Compare_to_Reference',
        '4. bam_search vs reference (novel and peptide-supported differences)'
    )
    _draw_existing_figure(
        axes[4],
        comparisons_dir / 'bam_vs_prediction_supported_summary.png',
        '5. supported and novel overlap between branches'
    )

    fig.tight_layout(rect=[0, 0, 1, 0.985])
    output_path = output_dir / 'combined_gapms_summary.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_sequence_coverage_groups_hist(df, output_dir):
    _, figures_dir = _get_figures_dir(output_dir)
    plt.figure(figsize=(6, 4))

    supported_mask = df['supported'].isin(['Yes', '+']) if 'supported' in df.columns else pd.Series(False, index=df.index)
    unsupported_mask = ~supported_mask

    if 'confidence' in df.columns:
        high_mask = df['confidence'] == 'high'
        moderate_mask = df['confidence'] == 'moderate'
    else:
        high_mask = df['high_conf'] == 'Yes'
        moderate_mask = supported_mask & ~high_mask

    plt.hist(
        df[supported_mask]["sequence_coverage"],
        bins=20, density=True, color="lightgray", alpha=0.8, label="Supported"
    )
    plt.hist(
        df[unsupported_mask]["sequence_coverage"],
        bins=20, density=True, color="dimgray", alpha=0.8, label="Not supported"
    )
    plt.hist(
        df[high_mask]["sequence_coverage"],
        bins=20, density=True, color="#ADD8E6", alpha=0.8, label="High confidence"
    )
    plt.hist(
        df[moderate_mask]["sequence_coverage"],
        bins=20, density=True, color="#D8BFD8", alpha=0.8, label="Moderate confidence"
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
        df[df["supported"].isin(['Yes', '+'])]["external_score"],
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

