#!/bin/bash
#SBATCH --job-name=gapms-bed-compare
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
cd /home/students/q.abbas/Git/GAP-MS
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/students/q.abbas/anaconda3/envs/gapms_env/bin/python scripts/compare_peplocus_beds.py evaluation/peplocus_integration/run_188267 \
 --reference /home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/outputs/gapms/latest_w_bam/helixer/prediction_search/peptides_mapped.bed \
 --genome /home/students/q.abbas/Git/PepLocus/test_carrot/GCF_001625215.2_DH1_v3.0_genomic.fna \
 --peptides /home/students/q.abbas/Git/PepLocus/evaluation/fragpipe_helixer/combined_peptide.peptides.tsv \
 --prior-peplocus /home/students/q.abbas/Git/PepLocus/evaluation/gtf_validation/run_188265/mappings.bed
