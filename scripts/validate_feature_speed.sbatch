#!/bin/bash
#SBATCH --job-name=gapms-speed
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
set -euo pipefail
cd /home/students/q.abbas/Git/GAP-MS
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$PWD"
PYTHON=/home/students/q.abbas/anaconda3/envs/gapms_env/bin/python
DATA=/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota
RUN="evaluation/feature_speed/run_${SLURM_JOB_ID}"
mkdir -p "$RUN"
/usr/bin/time -v "$PYTHON" -m gapms.main \
 -a /home/students/q.abbas/Git/PepLocus/test_carrot/GCF_001625215.2_DH1_v3.0_genomic.fna \
 -p /home/students/q.abbas/Git/PepLocus/evaluation/fragpipe_helixer/combined_peptide.peptides.tsv \
 -g "$DATA/inputs/helixer/GCF_001625215.2_Daucus_carota.helixer.gtf" \
 -f "$DATA/inputs/helixer/GCF_001625215.2_Daucus_carota.helixer.faa" \
 -m "$DATA/outputs/gapms/latest_w_bam/helixer/prediction_search/proteomapper.tsv" \
 -o "$RUN/workflow"
"$PYTHON" - "$RUN" <<'PY'
from pathlib import Path
import json,sys,hashlib
run=Path(sys.argv[1]);old=Path('evaluation/peplocus_integration/run_188267/workflow');new=run/'workflow'
names=['all_proteins_scores.tsv','supported_proteins.gtf','supported_proteins.fa',
       'peptides_mapped.bed','peptides_mapped.igv.bed','Unmpped_pepides.txt']
names += [str(p.relative_to(old)) for folder in ['Txt','Figures'] for p in (old/folder).glob('*') if p.is_file()]
names += [str(p.relative_to(old)) for p in (old/'peplocus').glob('*') if p.is_file() and p.name != 'summary.json']
checks={}
for name in names:
 a=old/name;b=new/name
 checks[name]=hashlib.sha256(a.read_bytes()).hexdigest()==hashlib.sha256(b.read_bytes()).hexdigest()
summary=json.loads((new/'peplocus/summary.json').read_text())
checks['default_mapping_threads_8']=summary['parameters']['threads']==8
(run/'equality.json').write_text(json.dumps(checks,indent=2)+'\n')
assert all(checks.values()),checks
print(json.dumps(checks,indent=2))
PY
