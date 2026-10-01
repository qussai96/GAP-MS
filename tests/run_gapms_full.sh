#!/bin/bash
#SBATCH --job-name=gapms_full
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=/home/students/q.abbas/Git/GAP-MS/tests/gapms_full_%j.out
#SBATCH --error=/home/students/q.abbas/Git/GAP-MS/tests/gapms_full_%j.err

set -euo pipefail

REPO=/home/students/q.abbas/Git/GAP-MS
OUTPUT_DIR="$REPO/tests/test_braker_carrot"
BAM=/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/inputs/varus/sorted_VARUS.bam
GENOME=/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/inputs/genomic/GCF_001625215.2_DH1_v3.0_genomic.fna
PEPTIDES="/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/outputs/fragpipe/braker/combined_peptide.peptides.tsv"
PRED_GTF="/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/inputs/braker/GCF_001625215.2_Daucus_carota.braker.gtf"
REF_GTF="/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/inputs/reference/GCF_001625215.2_Daucus_carota.ref.gtf"
REF_FASTA="/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/inputs/reference/GCF_001625215.2_Daucus_carota.ref.faa"

mkdir -p "$OUTPUT_DIR"
rm -rf "$OUTPUT_DIR"/*

copy_slurm_logs() {
  local job_id="${SLURM_JOB_ID:-}"
  local extension slurm_log

  [[ -n "$job_id" ]] || return 0
  for extension in out err; do
    slurm_log="$REPO/tests/gapms_full_${job_id}.${extension}"
    if [[ -f "$slurm_log" ]]; then
      cp "$slurm_log" "$OUTPUT_DIR/"
    fi
  done
}

# Preserve the job-specific logs on both success and failure. Slurm owns the
# original open files in tests/, so copying is safer than renaming them.
trap copy_slurm_logs EXIT

cd "$REPO"
CONDA_BIN="/home/students/q.abbas/anaconda3/bin/conda"

echo "[$(date)] Starting full GAP-MS BAM pipeline"
echo "Output directory: $OUTPUT_DIR"

"$CONDA_BIN" run -n gapms_env --no-capture-output python -m gapms.main \
  -g "$PRED_GTF" \
  -a "$GENOME" \
  -b "$BAM" \
  -p "$PEPTIDES" \
  -rg "$REF_GTF" \
  -rf "$REF_FASTA" \
  -o "$OUTPUT_DIR"

echo "[$(date)] GAP-MS job finished"
