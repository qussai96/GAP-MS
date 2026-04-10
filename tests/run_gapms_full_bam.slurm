#!/bin/bash
#SBATCH --job-name=gapms_bam_full
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/home/students/q.abbas/Git/GAP-MS/tests/gapms_bam_full_%j.out
#SBATCH --error=/home/students/q.abbas/Git/GAP-MS/tests/gapms_bam_full_%j.err

set -euo pipefail

REPO=/home/students/q.abbas/Git/GAP-MS
OUTPUT_DIR="$REPO/tests/bam_full_pipeline_output"
BAM=/home/students/q.abbas/GAPMS/0-raw_files/from18/GCF_000499845.2_Phaseolus_vulgaris/inputs/varus/sorted_VARUS.bam
GENOME=/home/students/q.abbas/GAPMS/0-raw_files/from18/GCF_000499845.2_Phaseolus_vulgaris/inputs/genomic/GCF_000499845.2_P._vulgaris_v2.0_genomic.fna
PEPTIDES="$REPO/tutorials/Helixer_MS_identified_peptides.tsv"
PRED_GTF="$REPO/tutorials/Helixer_prediction.gtf"
REF_GTF="$REPO/tutorials/Reference_annotation.gtf"
REF_FASTA="$REPO/tutorials/Reference_annotation.faa"

mkdir -p "$OUTPUT_DIR"
rm -rf "$OUTPUT_DIR"/*

cd "$REPO"
CONDA_BIN="/home/students/q.abbas/anaconda3/bin/conda"

echo "[$(date)] Starting full GAP-MS BAM pipeline"
echo "Output directory: $OUTPUT_DIR"

"$CONDA_BIN" run -n gapms_env --no-capture-output python -m gapms.main \
  -g "$PRED_GTF" \
  -b "$BAM" \
  -a "$GENOME" \
  -p "$PEPTIDES" \
  -rg "$REF_GTF" \
  -rf "$REF_FASTA" \
  -o "$OUTPUT_DIR"

echo "[$(date)] GAP-MS job finished"

mv gapms_bam_full_*.out "$OUTPUT_DIR"/.
mv gapms_bam_full_*.err "$OUTPUT_DIR"/.
