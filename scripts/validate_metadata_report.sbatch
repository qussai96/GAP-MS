#!/bin/bash
#SBATCH --job-name=gapms-meta-report
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
set -euo pipefail
cd /home/students/q.abbas/Git/GAP-MS
export PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
/home/students/q.abbas/anaconda3/envs/gapms_env/bin/python - <<'PY'
from pathlib import Path
import json,hashlib
from gapms.compare_supp_ref import generate_annotation_report
species=Path('/home/students/q.abbas/GAPMS/0-raw_files/GCF_036512215.1_Solanum_lycopersicum')
branch=species/'outputs/gapms/updated_w_peplocus/prediction_search'
new=Path('evaluation/model_metadata/tomato_reference_report').resolve();new.mkdir(parents=True,exist_ok=True)
(new/'protein_projection').mkdir(exist_ok=True)
(new/'peptides_mapped.bed').symlink_to(branch/'peptides_mapped.bed')
(new/'protein_projection/peptides_mapped.bed').symlink_to(Path('evaluation/model_metadata/GCF_036512215.1_Solanum_lycopersicum/peptides_mapped.bed').resolve())
generate_annotation_report(new,next((species/'inputs/helixer').glob('*.helixer.gtf')),
                           next((species/'inputs/reference').glob('*.ref.gtf')),
                           tmap_file=branch/'Compare_to_Reference/gffcmp.peptide_mapped_proteins_cleaned.gtf.tmap',
                           all_proteins_scores=branch/'all_proteins_scores.tsv',
                           protein_fasta=branch/'translated_predictions_proteins.fa')
checks={}
for path in (new/'Compare_to_Reference').rglob('*'):
 if path.is_file():
  relative=path.relative_to(new/'Compare_to_Reference');old=branch/'Compare_to_Reference'/relative
  checks[str(relative)]=old.exists() and path.read_bytes()==old.read_bytes()
(new/'equality.json').write_text(json.dumps(checks,indent=2)+'\n')
assert checks and all(checks.values()),checks
print(json.dumps(checks,indent=2))
PY
