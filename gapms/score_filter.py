import numpy as np


def calculate_sensitivity(scores, cutoff):
    if not scores:
        return 0.0
    return sum(score >= cutoff for score in scores) / len(scores)


def calculate_fpr(unsupported_scores, cutoff):
    # Calculate FPR directly: unsupported proteins passing cutoff / all unsupported proteins
    if not unsupported_scores:
        return 0.0  # All proteins are high-confidence, FPR = 0
    unsupported_passing = sum(s >= cutoff for s in unsupported_scores)
    return unsupported_passing / len(unsupported_scores)


def calculate_specificity_tnr(unsupported_scores, cutoff):
    # TNR (True Negative Rate) = 1 - FPR
    return 1 - calculate_fpr(unsupported_scores, cutoff)


def find_cutoff(proteins_scores_df, high_confident_proteins, score_column, output_dir=None, plot_roc=None):
    scores = proteins_scores_df[score_column].tolist()

    # Case: all scores are zero
    if all(score == 0 for score in scores):
        return 0

    # Subset high-confidence protein scores
    high_confident_proteins_scores = proteins_scores_df[
        proteins_scores_df['Protein'].isin(high_confident_proteins)
    ][score_column].tolist()

    if not high_confident_proteins_scores:
        return 0  # fallback if no high-confidence proteins

    # Pre-calculate unsupported scores once (expensive operation)
    high_conf_set = set(high_confident_proteins_scores)
    unsupported_scores = [s for s in scores if s not in high_conf_set]

    # Generate cutoff thresholds from 1.0 → 0.0 with finer granularity
    cutoffs = np.linspace(1.0, 0.0, 1000)

    tpr = [calculate_sensitivity(high_confident_proteins_scores, c) for c in cutoffs]
    fpr = [calculate_fpr(unsupported_scores, c) for c in cutoffs]
    tnr = [calculate_specificity_tnr(unsupported_scores, c) for c in cutoffs]

    # Sort by FPR to create a proper ROC curve
    roc_points = sorted(zip(fpr, tpr, cutoffs, tnr), key=lambda x: x[0])
    fpr_sorted, tpr_sorted, cutoffs_sorted, tnr_sorted = zip(*roc_points)
    fpr_sorted = list(fpr_sorted)
    tpr_sorted = list(tpr_sorted)
    cutoffs_sorted = list(cutoffs_sorted)
    tnr_sorted = list(tnr_sorted)

    # Primary method: Youden's Index (J = TPR + TNR - 1)
    youden_indices = [tpr_sorted[i] + tnr_sorted[i] - 1 for i in range(len(tpr_sorted))]
    best_idx_youden = np.argmax(youden_indices)
    cutoff_youden = cutoffs_sorted[best_idx_youden]
    j_max = youden_indices[best_idx_youden]

    # Secondary method: elbow method on ROC curve slope
    cutoff_elbow = None
    for i in range(2, len(tpr_sorted)):
        delta_fpr = fpr_sorted[i] - fpr_sorted[i - 1]
        slope = (tpr_sorted[i] - tpr_sorted[i - 1]) / delta_fpr if delta_fpr != 0 else float("inf")

        if slope < 1:
            cutoff_elbow = cutoffs_sorted[i]
            break

    if j_max > 0.2:
        cutoff = cutoff_youden
        method = "Youden's Index"
    elif cutoff_elbow is not None:
        cutoff = cutoff_elbow
        method = "Elbow Method (fallback)"
    else:
        cutoff = 0
        method = "No valid cutoff found"

    print(f"  {score_column} - {method}: cutoff={cutoff:.3f}, Youden's J={j_max:.3f}")

    if plot_roc and output_dir:
        plot_roc(
            fpr_sorted,
            tpr_sorted,
            cutoff,
            fpr_sorted[best_idx_youden],
            tpr_sorted[best_idx_youden],
            output_dir,
            score_column,
        )

    return cutoff


def find_apply_score_filter(proteins_scores_df, high_confident_proteins, score_column, output_dir=None, plot_roc=None):
    cutoff = find_cutoff(proteins_scores_df, high_confident_proteins, score_column, output_dir, plot_roc)
    cutoff = max(cutoff, 0.5)
    high_scoring = set(proteins_scores_df[proteins_scores_df[score_column] > cutoff]['Protein'])
    print(f"Determined cutoff of {score_column}: {cutoff}")
    print(f"Number of {score_column} supported proteins: {len(high_scoring)}")
    return cutoff, high_scoring
