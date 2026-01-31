import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_mapped_percentage_bars(df, output_dir):
    """
    Plots stacked percentage bars from a DataFrame.
    Expected columns: ['bar_name', 'mapped', 'unmapped']
    """

    # Compute percentages
    totals = df['mapped'] + df['unmapped']
    df['percent_mapped'] = df['mapped'] / totals * 100
    df['percent_unmapped'] = df['unmapped'] / totals * 100

    # Plot
    fig, ax = plt.subplots(figsize=(len(df) * 2, 6))  # width scales with number of bars

    ax.bar(df['bar_name'], df['percent_mapped'],
           color='green', edgecolor='black', label='Proteins with mapped peptides')
    ax.bar(df['bar_name'], df['percent_unmapped'],
           bottom=df['percent_mapped'],
           color='red', edgecolor='black', label='Proteins with no mapped peptides')

    # Add text inside bars
    for i, row in df.iterrows():
        ax.text(i, row['percent_mapped'] / 2,
                f"{row['mapped']}",
                ha='center', va='center', fontsize=10, color='white')
        ax.text(i, row['percent_mapped'] + row['percent_unmapped'] / 2,
                f"{row['unmapped']}",
                ha='center', va='center', fontsize=10, color='white')

    # Formatting
    ax.set_ylim(0, 100)
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Percentage of mapped vs un-mapped proteins')
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    
    ax.set_xticklabels(df['bar_name'], rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / 'Figures' / 'Peptides_mapping_bars.png', bbox_inches='tight')
    plt.close()




def plot_protein_coverage_groups_hist(df, output_dir):
    plt.figure(figsize=(6, 4))
    plt.hist(
        df[df["supported"] == 'Yes']["protein_coverage"],
        bins=20, density=True, color="lightgray", alpha=0.8, label="Supported"
    )
    plt.hist(
        df[df["unsupported"] == 'Yes']["protein_coverage"],
        bins=20, density=True, color="dimgray", alpha=0.8, label="Not supported"
    )
    plt.hist(
        df[df["high_conf"] == 'Yes']["protein_coverage"],
        bins=20, density=True, color="#ADD8E6", alpha=0.8, label="High confident"
    )
    plt.hist(
        df[df["low_conf"] == 'Yes']["protein_coverage"],
        bins=20, density=True, color="#D8BFD8", alpha=0.8, label="Low confident"
    )
    plt.yscale("log")
    plt.xlabel("Protein coverage")
    plt.ylabel("Density")
    plt.title("Protein Coverage by Confidence Groups")
    plt.legend()
    plt.savefig(output_dir / 'Figures'/ 'Protein_coverage_groups.png', dpi=300, bbox_inches='tight')
    plt.close()



def plot_external_scores(df, output_dir):
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
    plt.xlabel("external score")
    plt.ylabel("Density")
    plt.title("external Scores")
    plt.legend()
    plt.savefig(output_dir / 'Figures'/ 'External_scores.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_roc_curve(fpr, tpr, cutoff, fpr_cut, tpr_cut, output_dir, score_column):
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
    plt.savefig(output_dir / 'Figures'/ f'ROC_curve_{score_column}.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_shap_summary(model, X_labeled, output_dir):
    import shap
    explainer = shap.Explainer(model)
    shap_values = explainer(X_labeled)
    shap.summary_plot(shap_values, X_labeled, show=False)
    plt.savefig(output_dir / 'Figures'/ 'Shap_summary_plot.png', bbox_inches='tight')
    plt.close()


def plot_supported_percentage_bar(value1, value2, output_dir):
    # Compute total and percentages
    total = value1 + value2
    percent1 = value1 / total * 100
    percent2 = value2 / total * 100

    # Plot
    fig, ax = plt.subplots(figsize=(2.5, 6))  # Tall format
    bars1 = ax.bar(0, percent1, color='green', edgecolor='black')
    bars2 = ax.bar(0, percent2, bottom=percent1, color='red', edgecolor='black')

    # Add text inside bars
    ax.text(0, percent1 / 2, f'{value1}', ha='center', va='center', fontsize=12, color='white')
    ax.text(0, percent1 + percent2 / 2, f'{value2}', ha='center', va='center', fontsize=12, color='white')

    # Legend and title
    ax.legend(['supported proteins', 'un-supported proteins'], loc='upper left', bbox_to_anchor=(1.6, 1))
    ax.set_title('Percentage of supported and un-supported proteins', fontsize=12)

    # Formatting
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, 100)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / 'Figures'/ 'Supported_proteins_bars.png', bbox_inches='tight')
    plt.close()




# def plot_protein_metrics(df, output_dir):
    
#     # 1. Protein length distribution
#     plt.figure(figsize=(8, 5))
#     plt.hist(df['protein_length'], bins=50, color='skyblue', edgecolor='black')
#     plt.title('Protein Length Distribution')
#     plt.xlabel('Protein Length')
#     plt.ylabel('Count')
#     plt.xlim(0, 2000)
#     plt.savefig(output_dir / 'Figures'/ 'Protein_length_distribution.png', bbox_inches='tight')
#     plt.close()

#     # 2. Splice site distribution
#     bins = list(range(0, 12)) + [float('inf')]  # [0,1,2,...,11,inf]
#     labels = ['0'] + [str(i) for i in range(1, 11)] + ['>10']
#     mapped_ranges = pd.cut(df['splice_sites'], bins=bins, labels=labels, right=False)
#     mapped_counts = mapped_ranges.value_counts().sort_index()
#     plt.figure(figsize=(10, 5))
#     ax = mapped_counts.plot(kind='bar', color='green', edgecolor='black')
#     plt.title('Splice Sites Distribution')
#     plt.xlabel('Number of Splice Sites')
#     plt.ylabel('Protein Count')
#     plt.xticks(rotation=90)
#     # Add data labels
#     for i, count in enumerate(mapped_counts):
#         ax.text(i, count + 0.5, str(count), ha='center', va='bottom', fontsize=9)

#     plt.tight_layout()
#     plt.savefig(output_dir / 'Figures'/ 'Splice_sites_distribution.png', bbox_inches='tight')
#     plt.close()


#     # 3. Mapped peptides bar plot (0–30, >30)
#     bins = [0, 1, 2] + list(range(11, 101, 10)) + [np.inf]
#     labels = ['0', '1-10'] + [f"{i}-{i+9}" for i in range(11, 100, 10)] + ['>100']

#     # Bin the data
#     mapped_ranges = pd.cut(df['mapped_peptides'], bins=bins, labels=labels, right=False)
#     mapped_counts = mapped_ranges.value_counts().sort_index()

#     # Total number of proteins
#     total_proteins = len(df)

#     # Plotting
#     plt.figure(figsize=(10, 5))
#     ax = mapped_counts.plot(kind='bar', color='green', edgecolor='black')
#     plt.title('Mapped Peptides Distribution (in ranges)')
#     plt.xlabel('Number of Mapped Peptides')
#     plt.ylabel('Protein Count')
#     plt.xticks(rotation=45)

#     # Add value labels and percentage for the first bar
#     for i, count in enumerate(mapped_counts):
#         # Add the count labels above the bar
#         ax.text(i, count + 0.5, str(count), ha='center', va='bottom', fontsize=9)

#         # Add percentage inside the first bar
#         if i == 0:
#             percentage = (count / total_proteins) * 100
#             ax.text(i, count / 2, f'{percentage:.2f}%', ha='center', fontsize=10, rotation=90,color='white')

#     plt.tight_layout()
#     plt.savefig(output_dir / 'Figures'/ 'Mapped_peptides_distribution.png', bbox_inches='tight')
#     plt.close()



#     # 4. protein isoforms bar plot (0–10, >10)
#     mapped_counts_10 = df['protein_isoforms'].apply(lambda x: str(int(x)) if x <= 10 else '>10')
#     mapped_counts_10 = mapped_counts_10.value_counts().sort_index(key=lambda x: [int(i.replace('>','')) if '>' not in i else 999 for i in x])
#     plt.figure(figsize=(8, 5))
#     mapped_counts_10.plot(kind='bar', color='purple', edgecolor='black')
#     plt.title('Mapped Peptides Distribution')
#     plt.xlabel('Number of Protein isoforms')
#     plt.ylabel('isoforms Count')
#     plt.xticks(rotation=0)
#     plt.savefig(output_dir / 'Figures'/ 'Protein_isoforms_distribution.png', bbox_inches='tight')
#     plt.close()