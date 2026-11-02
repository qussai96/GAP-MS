import json,re,sys
from collections import defaultdict,Counter
from pathlib import Path
p=Path(sys.argv[1]);jobs=json.loads(p.read_text());results=[]
feature=sys.argv[2] if len(sys.argv)>2 else 'CDS'
def structures(path):
 models=defaultdict(list)
 with path.open() as f:
  for line in f:
   if '\t'+feature+'\t' not in line:continue
   r=line.split('\t');tid=re.search(r'transcript_id "([^"]+)"',r[8]).group(1)
   models[(r[0],r[6],tid)].append((int(r[3]),int(r[4]),r[7]))
 return Counter((chrom,strand,tuple(sorted(blocks))) for (chrom,strand,tid),blocks in models.items())
for j in jobs:
 suffix='orfs/candidate_orfs.gtf' if feature=='CDS' else 'transcripts/assembled_transcripts.gtf'
 n=Path(j['output'])/suffix;o=Path(j['output']).with_name('GAPMS_V1_casanovo')/suffix
 old,new=structures(o),structures(n);common=sum((old&new).values());record=dict(crop=j['crop'],old_models=sum(old.values()),new_models=sum(new.values()),shared_structures_with_multiplicity=common,identical_structure_multiset=old==new,shared_old_percent=100*common/sum(old.values()),shared_new_percent=100*common/sum(new.values()))
 results.append(record);print(record,flush=True)
(p.parent/('RNA_structure_comparison.json' if feature=='CDS' else 'RNA_exon_structure_comparison.json')).write_text(json.dumps(results,indent=2)+'\n')
