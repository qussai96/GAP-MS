#!/bin/bash
#SBATCH --job-name=gapms-metadata
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
set -euo pipefail
cd /home/students/q.abbas/Git/GAP-MS
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$PWD"
/home/students/q.abbas/anaconda3/envs/gapms_env/bin/python - <<'PY'
from pathlib import Path
import json,time,hashlib
from gapms.model_metadata import attach_model_metadata
root=Path('/home/students/q.abbas/GAPMS/0-raw_files')
report={}
for crop in ['GCF_036512215.1_Solanum_lycopersicum','GCF_020520425.1_Spinacia_oleracea','GCF_001625215.2_Daucus_carota']:
 species=root/crop;branch=species/'outputs/gapms/updated_w_peplocus/prediction_search'
 output=Path('evaluation/model_metadata')/crop
 started=time.monotonic()
 result=attach_model_metadata(next((species/'inputs/helixer').glob('*.helixer.gtf')),
                              branch/'translated_predictions_proteins.fa',branch/'proteomapper.tsv',output)
 prior=branch/'protein_projection/peptides_mapped.bed'
 row={'elapsed_seconds':time.monotonic()-started,'metadata_records':sum(1 for _ in result.open())-1,
      'sha256':hashlib.sha256(result.read_bytes()).hexdigest(),'old_metadata_path':str(prior)}
 if prior.exists():
  row['old_sha256']=hashlib.sha256(prior.read_bytes()).hexdigest()
  row['byte_identical']=row['sha256']==row['old_sha256']
  assert row['byte_identical'],(crop,row)
 report[crop]=row
 print(crop,json.dumps(row),flush=True)
 Path('evaluation/model_metadata/results.json').write_text(json.dumps(report,indent=2)+'\n')
PY
