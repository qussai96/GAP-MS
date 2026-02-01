import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV
from xgboost import XGBClassifier

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
    # Finds the point with maximum distance from diagonal (most balanced)
    youden_indices = [tpr_sorted[i] + tnr_sorted[i] - 1 for i in range(len(tpr_sorted))]
    best_idx_youden = np.argmax(youden_indices)
    cutoff_youden = cutoffs_sorted[best_idx_youden]
    j_max = youden_indices[best_idx_youden]

    # Secondary method: Find cutoff where slope of ROC curve < 1 (elbow method)
    # Only use if Youden's index is too low (indicates poor classifier)
    cutoff_elbow = None
    for i in range(2, len(tpr_sorted)):
        delta_fpr = fpr_sorted[i] - fpr_sorted[i - 1]
        slope = (tpr_sorted[i] - tpr_sorted[i - 1]) / delta_fpr if delta_fpr != 0 else float("inf")

        if slope < 1:
            cutoff_elbow = cutoffs_sorted[i]
            break

    # Use Youden's method as primary; use elbow only if Youden gives very poor result
    if j_max > 0.2:  # Youden's index > 0.2 indicates decent classifier
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
        plot_roc(fpr_sorted, tpr_sorted, cutoff, fpr_sorted[best_idx_youden], tpr_sorted[best_idx_youden], output_dir, score_column)

    return cutoff



def find_apply_score_filter(proteins_scores_df, high_confident_proteins, score_column, output_dir=None, plot_roc=None):
    cutoff = find_cutoff(proteins_scores_df, high_confident_proteins, score_column, output_dir, plot_roc)
    cutoff = max(cutoff, 0.5)
    high_scoring = set(proteins_scores_df[proteins_scores_df[score_column] > cutoff]['Protein'])
    print(f"Determined cutoff of {score_column}: {cutoff}")
    print(f"Number of {score_column} supported proteins: {len(high_scoring)}")
    return cutoff, high_scoring 


def get_high_confident_proteins(df):
        high_confident = []
        for _, group in df.groupby("Protein"):
            for _, row in group.iterrows():
                cond1 = row["protein_specific_peptides"] >= 2
                cond2 = row["protein_coverage"] >= 0.8
                cond3 = row["N_terminal_peptides"] >= 1 and row["C_terminal_peptides"] >= 1
                cond4 = row["splice_peptides"] == row["splice_sites"] and row["splice_peptides"] > 0
                if cond1 or cond2 or cond3 or cond4:
                    high_confident.append(row)
                    break
        print(f"Number of high confident proteins: {len(high_confident)}")
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
        return pd.DataFrame(low_confident)

def train_iterative_model(df, high_confident_df, low_confident_df, pos_thr=0.90, neg_thr=0.10, n_iter=5, shap_output_dir=None, plot_shap=None):
    feature_cols = [
    "protein_coverage",
    "protein_specific_peptides",
    "gene_specific_peptides",
    "splice_peptides",
    "internal_peptides",
    "mapped_peptides",
    "N_terminal_peptides",
    "C_terminal_peptides",
    "protein_length"
    ]
    
    # Initialize label column
    df["label"] = np.nan

    df.loc[high_confident_df.index, "label"] = 1
    df.loc[low_confident_df.index, "label"] = 0

    # Prepare labeled data
    df_labeled = df.dropna(subset=["label"])
    X_labeled = df_labeled[feature_cols]
    y_labeled = df_labeled["label"]

    if (y_labeled == 0).sum() == 0:
        raise ValueError("No negative labels found. Cannot train model.")

    # Stratified CV for imbalance handling
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # Initialize model + hyperparameter search
    base_model = XGBClassifier(eval_metric="logloss", random_state=42)

    param_dist = {
        "n_estimators": [100, 200],
        "max_depth": [3, 5, 7],
        "learning_rate": [0.01, 0.1],
        "subsample": [0.8, 1.0],
        "colsample_bytree": [0.8, 1.0],
        "scale_pos_weight": [1, 2, 5],
    }

    random_search = RandomizedSearchCV(
        estimator=base_model,
        param_distributions=param_dist,
        n_iter=20,
        scoring="roc_auc",
        cv=skf,
        verbose=1,
        random_state=42,
        n_jobs=-1,
    )

    random_search.fit(X_labeled, y_labeled)
    model = random_search.best_estimator_

    # Optional SHAP analysis
    # if plot_shap and shap_output_dir:
    #     plot_shap(model, X_labeled, shap_output_dir)

    # Track new positives added
    total_new_pos = set()

    # Iterative labeling of unlabeled data
    for it in range(1, n_iter + 1):
        df_unlabeled = df[df["label"].isna()]
        if df_unlabeled.empty:
            break

        X_unlabeled = df_unlabeled[feature_cols]
        prob_pred = model.predict_proba(X_unlabeled)[:, 1]

        new_pos = df_unlabeled.index[prob_pred >= pos_thr]
        new_neg = df_unlabeled.index[prob_pred <= neg_thr]

        df.loc[new_pos, "label"] = 1
        df.loc[new_neg, "label"] = 0

        # Keep track of newly added positives
        total_new_pos.update(new_pos)

        print(f"Iteration {it}: added {len(new_pos)} new positives, {len(new_neg)} new negatives.")

        df_labeled = df.dropna(subset=["label"])
        X_labeled = df_labeled[feature_cols]
        y_labeled = df_labeled["label"]

        if y_labeled.nunique() == 2:
            model = XGBClassifier(eval_metric="logloss", random_state=42)
            model.fit(X_labeled, y_labeled)

        if len(new_pos) == 0 and len(new_neg) == 0:
            break

    # Final predictions for all proteins
    df["prob_pos"] = model.predict_proba(df[feature_cols])[:, 1]
    df["final_label"] = (df["prob_pos"] >= 0.5).astype(int)

    # Print XGBoost evaluation metrics on labeled data
    from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score
    y_true = df_labeled["label"]
    y_pred = model.predict(X_labeled)
    y_proba = model.predict_proba(X_labeled)[:, 1]
    # print("\nXGBoost Evaluation Metrics (on labeled data):")
    # print(f"Accuracy: {accuracy_score(y_true, y_pred):.4f}")
    # print(f"Precision: {precision_score(y_true, y_pred):.4f}")
    # print(f"Recall: {recall_score(y_true, y_pred):.4f}")
    # print(f"F1-score: {f1_score(y_true, y_pred):.4f}")
    # print(f"AUROC: {roc_auc_score(y_true, y_proba):.4f}")
    
    # Feature importance reporting using SHAP
    print("\nFeature importance (SHAP summary):")
    try:
        import shap
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_labeled)
        mean_abs_shap = np.abs(shap_values).mean(axis=0)
        feature_importance = dict(zip(feature_cols, mean_abs_shap))
        sorted_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)
        for feat, val in sorted_features:
            print(f"{feat}: {val:.4f}")
        # Overfitting check: warn if one feature dominates
        # if sorted_features[0][1] > 2 * sorted_features[1][1]:
            # print(f"WARNING: Model may be overfitting to '{sorted_features[0][0]}' (SHAP value much higher than others)")
    except Exception as e:
        print(f"Could not compute SHAP feature importance: {e}")

    # Report total number of new positives added beyond initial high-confidence ones
    print(f"Total new positives added during iterations (beyond high-confidence set): {len(total_new_pos)}")

    return df
