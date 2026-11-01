#!/bin/bash
#SBATCH --job-name=gapms-peplocus
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
set -euo pipefail
cd /home/students/q.abbas/Git/GAP-MS
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg PYTHONNOUSERSITE=1
PYTHON=/home/students/q.abbas/anaconda3/envs/gapms_env/bin/python
GAPMS=/home/students/q.abbas/anaconda3/envs/gapms_env/bin/gapms
DATA=/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota
GENOME=/home/students/q.abbas/Git/PepLocus/test_carrot/GCF_001625215.2_DH1_v3.0_genomic.fna
PEPTIDES=/home/students/q.abbas/Git/PepLocus/evaluation/fragpipe_helixer/combined_peptide.peptides.tsv
GTF="$DATA/inputs/helixer/GCF_001625215.2_Daucus_carota.helixer.gtf"
BAM=/home/students/q.abbas/Git/PepLocus/test_carrot/sorted_VARUS.bam
RUN="evaluation/peplocus_integration/run_${SLURM_JOB_ID}"
mkdir -p "$RUN"
"$PYTHON" -c 'from setuptools.build_meta import build_wheel,build_sdist;build_wheel("dist");build_sdist("dist")'
"$PYTHON" -m pip install --no-deps --force-reinstall dist/gapms-0.2.0-py3-none-any.whl
"$PYTHON" -m pytest -q
# Run outside the source tree to exercise the installed package.
PROJECT="$PWD"
cd /tmp
"$GAPMS" map -a "$GENOME" -p "$PEPTIDES" -g "$GTF" -t 8 -o "$PROJECT/$RUN/gtf"
"$GAPMS" map -a "$GENOME" -p "$PEPTIDES" -b "$BAM" -t 8 -o "$PROJECT/$RUN/bam"
"$PYTHON" -m gapms._peplocus verify --genome "$GENOME" --mappings "$PROJECT/$RUN/gtf/peplocus/mappings.tsv" > "$PROJECT/$RUN/gtf_verification.json"
# Full prediction/classifier branch with the prior protein mapping held fixed.
"$GAPMS" -a "$GENOME" -p "$PEPTIDES" -g "$GTF" \
  -f "$DATA/inputs/helixer/GCF_001625215.2_Daucus_carota.helixer.faa" \
  -m "$DATA/outputs/gapms/latest_w_bam/helixer/prediction_search/proteomapper.tsv" \
  -t 8 -o "$PROJECT/$RUN/workflow"
