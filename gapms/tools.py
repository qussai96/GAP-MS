import subprocess
import os
from pathlib import Path

def run_proteomapper(protdb_path: str, peptides_path: str, output_dir):
    # Get the directory of this Python file (assuming Perl scripts are in the same dir)
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # Paths to the local Perl scripts
    clips_path = os.path.join(base_dir, "clips.pl")
    promast_path = os.path.join(base_dir, "promast.pl")

    # Output file in the same directory as peptides
    output_file = output_dir / "proteomapper.tsv"
    # output_file = peptides_path.with_suffix(".proteomapper.tsv")

    # Run clips.pl
    clips_cmd = ["perl", clips_path, "-f", protdb_path]
    clips_result = subprocess.run(clips_cmd, capture_output=True, text=True)

    if clips_result.returncode != 0:
        raise RuntimeError(f"clips.pl failed:\n{clips_result.stderr}")

    # Run promast.pl
    promast_cmd = f"perl {promast_path} {protdb_path} {peptides_path} > {output_file}"
    promast_result = subprocess.run(promast_cmd, shell=True, capture_output=True, text=True)

    if promast_result.returncode != 0:
        raise RuntimeError(f"promast.pl failed:\n{promast_result.stderr}")
    
    pep_idx_file = str(protdb_path) +'.pep.idx'
    if os.path.exists(pep_idx_file):
        os.remove(pep_idx_file)
    return output_file


def run_gffread(genome_path: str, gtf_path: str, output_dir, name):
    output_file = output_dir / f"translated_{name}_proteins.fa"

    command = [
        "gffread",
        "-C",
        "-g", str(genome_path),
        "-y", str(output_file),
        str(gtf_path)
    ]
    print("Running gffread:", " ".join(command))

    gffread_result = subprocess.run(command, capture_output=True, text=True)

    if gffread_result.returncode != 0:
        raise RuntimeError(f"gffread failed:\n{gffread_result.stderr}")
    return output_file


def run_psauron(protdb_path: Path, output_dir):
    output_file = output_dir / "psauron.csv"

    command = [
        "psauron",
        "-i", str(protdb_path),
        "-o", str(output_file),
        "-p"
    ]

    psauron_result = subprocess.run(command, capture_output=True, text=True)

    if psauron_result.returncode != 0:
        raise RuntimeError(f"psauron failed:\n{psauron_result.stderr}")

    print("STDOUT:\n", psauron_result.stdout)
    print("STDERR:\n", psauron_result.stderr)

    return output_file


def run_embeddings(protdb_path: Path, model_file: Path, output_dir: Path):
    """
    Generate embeddings from protein FASTA and score them using a trained model.
    
    Args:
        protdb_path: Path to protein FASTA file
        model_file: Path to trained scoring model (.pth)
        output_dir: Directory to save output files
    
    Returns:
        Path to external_scores.csv file
    """
    # Lazy import to avoid loading heavy dependencies unless needed
    try:
        from gapms.embeddings import convert_to_embeddings, get_external_scores_from_embeddings
    except ImportError as e:
        raise RuntimeError(
            "Embeddings functionality requires additional dependencies. "
            "Please install: pip install torch fair-esm"
        ) from e
    
    embeddings_file = output_dir / "proteins_embeddings.npy"
    external_scores_file = output_dir / "external_scores.csv"

    try:
        # Step 1: Convert proteins to embeddings
        print("Converting proteins to embeddings...")
        convert_to_embeddings(protdb_path, embeddings_file, batch_size=64)
        
        # Step 2: Score embeddings using the trained model
        print("Scoring embeddings...")
        get_external_scores_from_embeddings(
            embeddings_file,
            model_file,
            external_scores_file,
            dim=1280,
            n_layers=3,
            batch_size=2048
        )
        
        return external_scores_file
        
    except Exception as e:
        raise RuntimeError(f"Embeddings processing failed: {str(e)}")


def run_gffcompare(reference_gtf_path: str, supported_gtf_path: str):

    command = [
        "gffcompare",
        "-r", str(reference_gtf_path),
        str(supported_gtf_path)
    ]
    print("Running gffcompare:", " ".join(command))

    gffcompare_result = subprocess.run(command, capture_output=True, text=True)

    if gffcompare_result.returncode != 0:
        raise RuntimeError(f"gffcompare failed:\n{gffcompare_result.stderr}")
    return gffcompare_result
