import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV, cross_val_predict
from xgboost import XGBClassifier

from .score_filter import (
    calculate_fpr,
    calculate_sensitivity,
    calculate_specificity_tnr,
    find_apply_score_filter,
    find_cutoff,
)

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
DEFAULT_PRETRAINED_MODEL_NAME = "xgb_protein_classifier.pkl"


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

def train_iterative_model(df, high_confident_df, low_confident_df, pos_thr=0.90, neg_thr=0.10, n_iter=5, shap_output_dir=None, plot_shap=None):
    feature_cols = ITERATIVE_FEATURE_COLS
    
    # Initialize label column
    df["label"] = np.nan

    df.loc[high_confident_df.index, "label"] = 1
    df.loc[low_confident_df.index, "label"] = 0

    # Prepare labeled data
    df_labeled = df.dropna(subset=["label"])
    X_labeled = df_labeled[feature_cols]
    y_labeled = df_labeled["label"]

    if y_labeled.nunique() < 2:
        raise ValueError("Both positive and negative seed labels are required to train the iterative model.")

    cv_folds = min(5, int(y_labeled.value_counts().min()))
    if cv_folds < 2:
        raise ValueError(
            "At least two positive and two negative seed labels are required "
            "to train the iterative model."
        )

    # Stratified CV for imbalance handling
    skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)

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

    # Report out-of-fold diagnostics on the original rule-derived seed labels.
    # These are not independent biological validation, but unlike training-set
    # metrics they do not score the fitted model on the same rows it saw.
    seed_y_pred = cross_val_predict(model, X_labeled, y_labeled, cv=skf, method="predict")
    seed_y_proba = cross_val_predict(
        model, X_labeled, y_labeled, cv=skf, method="predict_proba"
    )[:, 1]
    seed_metrics = {
        "accuracy": accuracy_score(y_labeled, seed_y_pred),
        "precision": precision_score(y_labeled, seed_y_pred, zero_division=0),
        "recall": recall_score(y_labeled, seed_y_pred, zero_division=0),
        "f1": f1_score(y_labeled, seed_y_pred, zero_division=0),
        "auroc": roc_auc_score(y_labeled, seed_y_proba),
    }

    # Optional SHAP analysis
    if plot_shap and shap_output_dir:
        plot_shap(model, X_labeled, shap_output_dir)

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
            model = clone(model)
            model.fit(X_labeled, y_labeled)

        if len(new_pos) == 0 and len(new_neg) == 0:
            break

    # Final predictions for all proteins
    df["prob_pos"] = model.predict_proba(df[feature_cols])[:, 1]
    df["final_label"] = (df["prob_pos"] >= 0.5).astype(int)

    print("\nCross-validated diagnostics on initial seed labels:")
    print(f"Accuracy: {seed_metrics['accuracy']:.4f}")
    print(f"Precision: {seed_metrics['precision']:.4f}")
    print(f"Recall: {seed_metrics['recall']:.4f}")
    print(f"F1-score: {seed_metrics['f1']:.4f}")
    print(f"AUROC: {seed_metrics['auroc']:.4f}")
    
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


def resolve_pretrained_model_path(model_path=None):
    """Resolve the best available pretrained model path."""
    if model_path is not None:
        candidate_paths = [Path(model_path)]
    else:
        module_dir = Path(__file__).resolve().parent
        repo_dir = module_dir.parent
        workspace_dir = repo_dir.parent
        candidate_paths = [
            module_dir / DEFAULT_PRETRAINED_MODEL_NAME,
            repo_dir / DEFAULT_PRETRAINED_MODEL_NAME,
            workspace_dir / "GAP-MS-ONE-MODEL" / DEFAULT_PRETRAINED_MODEL_NAME,
        ]

    for candidate in candidate_paths:
        if candidate.exists():
            return candidate

    searched = "\n - ".join(str(path) for path in candidate_paths)
    raise FileNotFoundError(
        "Pre-trained model not found. Checked:\n - " + searched
    )


def _get_pretrained_feature_columns(model, df):
    """Infer which features the pretrained model expects."""
    if hasattr(model, "feature_names_in_"):
        feature_cols = list(model.feature_names_in_)
    else:
        n_features = getattr(model, "n_features_in_", None)
        if n_features == len(ITERATIVE_FEATURE_COLS):
            feature_cols = ITERATIVE_FEATURE_COLS
        elif n_features == len(BASE_FEATURE_COLS):
            feature_cols = BASE_FEATURE_COLS
        else:
            raise ValueError(
                "Could not infer the pretrained model feature set. "
                f"Model reports {n_features} features."
            )

    # Backward/forward compatibility for renamed coverage term.
    legacy_coverage_col = "protein" + "_coverage"
    current_coverage_col = "sequence_coverage"
    if legacy_coverage_col in feature_cols and legacy_coverage_col not in df.columns and current_coverage_col in df.columns:
        df[legacy_coverage_col] = df[current_coverage_col]
    if current_coverage_col in feature_cols and current_coverage_col not in df.columns and legacy_coverage_col in df.columns:
        df[current_coverage_col] = df[legacy_coverage_col]

    missing_features = [col for col in feature_cols if col not in df.columns]
    if missing_features:
        raise ValueError(f"Missing required features for pretrained model: {missing_features}")

    return feature_cols


def apply_pretrained_model(df, model_path=None):
    """Apply a pre-trained XGBoost classifier to the protein feature table."""
    resolved_model_path = resolve_pretrained_model_path(model_path)
    print(f"Using pre-trained XGBoost model: {resolved_model_path}")

    with open(resolved_model_path, 'rb') as handle:
        model = pickle.load(handle)

    df_copy = df.copy()
    feature_cols = _get_pretrained_feature_columns(model, df_copy)
    X = df_copy[feature_cols]

    df_copy["prob_pos"] = model.predict_proba(X)[:, 1]
    df_copy["final_label"] = (df_copy["prob_pos"] >= 0.5).astype(int)

    print(
        f"Pre-trained model predictions complete: "
        f"{(df_copy['final_label'] == 1).sum()} positive, "
        f"{(df_copy['final_label'] == 0).sum()} negative"
    )
    return df_copy


def run_protein_classifier(
    df,
    high_confident_df,
    low_confident_df,
    use_iterative_training=False,
    model_path=None,
    pos_thr=0.90,
    neg_thr=0.10,
    n_iter=5,
    shap_output_dir=None,
    plot_shap=None,
):
    """
    Run protein classification in two modes:
      1) default: use the pre-trained XGBoost model
      2) iterative: train from high/low confidence sets when `use_iterative_training=True`

    Iterative training is only used when explicitly requested. A missing or
    incompatible pretrained model is reported instead of silently changing the
    analysis method.
    """
    if use_iterative_training:
        print("Using iterative XGBoost training mode (-i/--iterative).")
        labeled_df = train_iterative_model(
            df=df,
            high_confident_df=high_confident_df,
            low_confident_df=low_confident_df,
            pos_thr=pos_thr,
            neg_thr=neg_thr,
            n_iter=n_iter,
            shap_output_dir=shap_output_dir,
            plot_shap=plot_shap,
        )
        return labeled_df, "iterative"

    labeled_df = apply_pretrained_model(df=df, model_path=model_path)
    return labeled_df, "pretrained"
