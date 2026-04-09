#!/usr/bin/env python3
"""
Functions for converting protein sequences to ESM2 embeddings and scoring them.
"""
import gc
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import esm


# ============================================================================
# EMBEDDINGS EXTRACTION
# ============================================================================

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta_sequences(fasta_path: str) -> tuple[List[str], List[str]]:
    """
    Read sequences from FASTA file.
    Returns: (ids, sequences)
    """
    ids = []
    sequences = []
    
    with open(fasta_path, "r") as f:
        current_id = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith(">"):
                # Save previous sequence
                if current_id is not None:
                    ids.append(current_id)
                    sequences.append("".join(current_seq))
                
                # Start new sequence
                current_id = line[1:].split()[0]  # Get first part of header
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id is not None:
            ids.append(current_id)
            sequences.append("".join(current_seq))
    
    return ids, sequences


def filter_standard_aa(ids: List[str], sequences: List[str]) -> tuple[List[str], List[str]]:
    """Keep only sequences made of the 20 standard amino acids."""
    filtered_ids = []
    filtered_seqs = []
    
    for seq_id, seq in zip(ids, sequences):
        seq = seq.strip().upper()
        if set(seq) <= ALLOWED_AA and len(seq) > 0:
            filtered_ids.append(seq_id)
            filtered_seqs.append(seq)
    
    return filtered_ids, filtered_seqs


@torch.no_grad()
def embed_sequences(
    sequences: List[str],
    batch_size: int = 64,
    device: torch.device = torch.device("cpu"),
) -> np.ndarray:
    """Return a (N, D) numpy array of per-sequence mean embeddings."""
    print(f"Loading ESM2 model...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    embeddings = []
    i = 0
    total = len(sequences)
    
    while i < total:
        # Simple batching
        batch_seqs = sequences[i : i + batch_size]
        batch_data = [(str(j), s) for j, s in enumerate(batch_seqs)]
        
        try:
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(device)

            results = model(batch_tokens, repr_layers=[33], return_contacts=False)
            token_reps = results["representations"][33]  # (B, L, D)

            # token 0 is BOS; last token is EOS → average 1..(len-2)
            lens = (batch_tokens != alphabet.padding_idx).sum(1).tolist()
            for b, L in enumerate(lens):
                # exclude BOS (0) and EOS (L-1)
                rep = token_reps[b, 1 : L - 1].mean(0)  # (D,)
                embeddings.append(rep.detach().cpu().float().numpy())

            if (i // batch_size) % 10 == 0:
                print(f"Processed {i}/{total} sequences...")
            
            i += batch_size

        except torch.cuda.OutOfMemoryError:
            # fallback: reduce batch size or skip if already at 1
            if batch_size > 1:
                torch.cuda.empty_cache()
                batch_size = max(1, batch_size // 2)
                print(f"[WARN] CUDA OOM. Reducing batch size to {batch_size} and retrying...", file=sys.stderr)
            else:
                # extremely long sequence – skip to make progress
                print(f"[WARN] Skipping very long sequence at index {i} (len={len(batch_seqs[0])}).", file=sys.stderr)
                i += 1
                torch.cuda.empty_cache()

    return np.vstack(embeddings)


def convert_to_embeddings(fasta_file: Path, output_file: Path, batch_size: int = 64) -> tuple[Path, List[str]]:
    """
    Convert protein FASTA file to ESM2 embeddings.
    
    Args:
        fasta_file: Path to input protein FASTA file
        output_file: Path to output .npy embeddings file
        batch_size: Batch size for embedding computation
    
    Returns:
        (embeddings_file_path, list_of_protein_ids)
    """
    gc.collect()
    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        torch.cuda.empty_cache()
    else:
        device = torch.device("cpu")

    print(f"Device: {device}")
    print(f"Reading sequences from {fasta_file}...")
    
    ids, sequences = read_fasta_sequences(str(fasta_file))
    ids, sequences = filter_standard_aa(ids, sequences)
    
    print(f"Sequences after filtering: {len(sequences)}")
    
    if not sequences:
        raise ValueError("No valid sequences found after filtering.")
    
    print("Computing embeddings...")
    embs = embed_sequences(sequences, batch_size=batch_size, device=device)
    
    # Save embeddings
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    np.save(output_file, embs)
    print(f"Saved embeddings to: {output_file}")
    print(f"Shape: {embs.shape}")
    
    # Save IDs file
    ids_file = output_file.with_suffix('.ids.txt')
    with open(ids_file, 'w') as f:
        for protein_id in ids:
            f.write(f"{protein_id}\n")
    print(f"Saved protein IDs to: {ids_file}")
    
    return output_file, ids


# ============================================================================
# SCORING FROM EMBEDDINGS
# ============================================================================

class ScoringModel(nn.Module):
    """MLP model for scoring protein embeddings."""
    def __init__(self, dim: int, n_layers: int):
        super().__init__()
        layers = []
        for i in range(n_layers):
            c_in = dim if i == 0 else int(dim / (i * 2))
            c_out = int(dim / ((i + 1) * 2))
            layers.append(nn.Sequential(
                nn.Linear(c_in, c_out),
                nn.LayerNorm(c_out),
                nn.ReLU(),
            ))
        self.layers = nn.Sequential(*layers)
        self.classifier = nn.Linear(int(dim / (n_layers * 2)), 1)

    def forward(self, x):
        # x shape during training was (batch, 1, dim)
        out = self.layers(x)
        logits = self.classifier(out).squeeze(-1).squeeze(-1)  # -> (batch,)
        return logits


def load_scoring_model(pth_path: str, device: torch.device, dim: int = 1280, n_layers: int = 3) -> nn.Module:
    """Load the scoring model from checkpoint."""
    obj = torch.load(pth_path, map_location=device)

    # Case 1: checkpoint dict with 'model_state_dict'
    if isinstance(obj, dict) and 'model_state_dict' in obj:
        model = ScoringModel(dim, n_layers)
        model.load_state_dict(obj['model_state_dict'])
        print(f"Loaded model_state_dict from checkpoint (epoch {obj.get('epoch', '?')})")

    # Case 2: checkpoint dict with full model under 'model'
    elif isinstance(obj, dict) and 'model' in obj and isinstance(obj['model'], nn.Module):
        model = obj['model']
        print("Loaded full nn.Module from checkpoint dict['model']")

    # Case 3: obj is directly an nn.Module
    elif isinstance(obj, nn.Module):
        model = obj
        print("Loaded full nn.Module directly")

    else:
        raise ValueError(
            "Unrecognized checkpoint format. "
            "Expected a dict with 'model_state_dict' or a full nn.Module."
        )

    model.to(device)
    model.eval()
    return model


def batched_predict(model: nn.Module, embs: np.ndarray, device: torch.device,
                    batch_size: int = 1024) -> np.ndarray:
    """Run batched prediction on embeddings."""
    # Ensure (N, 1, dim) like in training
    if embs.ndim == 1:
        embs = embs.reshape(1, -1)
    if embs.ndim == 2:
        embs = embs[:, None, :]  # add dummy channel dim -> (N, 1, dim)
    if embs.ndim != 3 or embs.shape[1] != 1:
        raise ValueError(f"Unexpected embedding shape {embs.shape}; expected (N, 1280) or (N, 1, 1280).")

    ds = TensorDataset(torch.from_numpy(embs).float())
    dl = DataLoader(ds, batch_size=batch_size, shuffle=False, num_workers=0)
    probs_list = []
    
    with torch.no_grad():
        for (x,) in dl:
            x = x.to(device, non_blocking=True)
            logits = model(x)
            probs = torch.sigmoid(logits).float().cpu().numpy()
            probs_list.append(probs)
    
    return np.concatenate(probs_list, axis=0)


def get_external_scores_from_embeddings(
    embeddings_file: Path,
    model_file: Path,
    output_file: Path,
    dim: int = 1280,
    n_layers: int = 3,
    batch_size: int = 2048
) -> Path:
    """
    Score protein embeddings using a trained MLP model.
    
    Args:
        embeddings_file: Path to .npy embeddings file
        model_file: Path to trained model .pth file
        output_file: Path to output CSV file
        dim: Embedding dimension (default: 1280 for ESM2)
        n_layers: Number of MLP layers (default: 3)
        batch_size: Inference batch size
    
    Returns:
        Path to output CSV file with columns: Protein, external_score
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print("Loading embeddings...")
    embs = np.load(embeddings_file)
    
    # Load IDs file
    ids_file = embeddings_file.with_suffix('.ids.txt')
    if not ids_file.exists():
        raise FileNotFoundError(f"IDs file not found: {ids_file}")
    
    with open(ids_file) as f:
        ids = [line.strip() for line in f]

    if embs.shape[0] != len(ids):
        print(f"Warning: embeddings count ({embs.shape[0]}) != ids count ({len(ids)}). "
              f"Proceeding with min length.")
    
    n = min(embs.shape[0], len(ids))
    embs = embs[:n]
    ids = ids[:n]

    print("Loading scoring model...")
    model = load_scoring_model(str(model_file), device, dim, n_layers)

    print("Running inference...")
    probs = batched_predict(model, embs, device, batch_size=batch_size)

    print("Writing output CSV...")
    rows = []
    for protein_id, score in zip(ids, probs):
        rows.append([protein_id, round(float(score), 5)])
    
    df = pd.DataFrame(rows, columns=["Protein", "external_score"])
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)
    print(f"Saved scores to: {output_file}")
    
    gc.collect()
    return output_file
