"""Independent nucleotide-level check of withheld peptide coordinate scoring."""
import csv,json,sqlite3,sys,hashlib
from gapms_v1.common import input_identity
from pathlib import Path
from itertools import groupby
import ahocorasick
import pysam
from Bio.Seq import Seq

jobs=json.loads(Path(sys.argv[1]).read_text());root=Path(sys.argv[1]).parent;audits=[];missed=[]
for j in jobs:
 d=root/'heldout'/j['crop'];selection=json.loads((d/'selection.json').read_text());targets=selection['target_models']
 old=Path(selection['rna_reconstruction_origin']);origin=json.loads((old/'run_manifest.json').read_text())
 identities={x['path']:x for x in origin['inputs']}
 for source in [j['assembly'],j['bam']]:assert input_identity(Path(source))==identities[str(Path(source).resolve())]
 assert hashlib.sha256((old/'orfs/candidate_orfs.gtf').read_bytes()).digest()==hashlib.sha256((d/'run/orfs/candidate_orfs.gtf').read_bytes()).digest()
 regions={}
 for t in targets:
  blocks=t['exons'] or t['cds'];regions.setdefault(t['chrom'],[]).append((min(x for x,y in blocks),max(y for x,y in blocks)))
 with (d/'masked_reference.gtf').open() as f:
  for line in f:
   if '\ttranscript\t' not in line:continue
   r=line.split('\t');left,right=int(r[3])-1,int(r[4])
   assert not any(max(left,x)<min(right,y) for x,y in regions.get(r[0],[])),('annotation_leak',j['crop'],r[:5])
 with (d/'target_recovery.tsv').open() as f:records={r['gene']:r for r in csv.DictReader(f,delimiter='\t')}
 queries={x.strip().replace('I','L') for x in (d/'run/mapping/queries.txt').read_text().splitlines() if x.strip()}
 automaton=ahocorasick.Automaton()
 for seq in sorted(queries):automaton.add_word(seq,seq)
 automaton.make_automaton()
 con=sqlite3.connect(f'file:{d}/run/work/evidence.sqlite?mode=ro',uri=True)
 observed={(p.replace('I','L'),c,s,tuple(map(tuple,json.loads(b)))) for p,c,s,b in con.execute('SELECT canonical,chrom,strand,blocks FROM hits')};con.close()
 with (d/'run/work/search_junctions.tsv').open() as f:available={(r['chrom'],int(r['start']),int(r['end']),r['strand']) for r in csv.DictReader(f,delimiter='\t')}
 index=next((d/'reference_index').glob('*.fai'));n=0;recovered=0
 with pysam.FastaFile(j['assembly'],filepath_index=str(index)) as fa:
  for t in targets:
   # Expand the truth CDS into genomic nucleotide coordinates in coding order.
   coordinates=[pos for a,b in t['coding'] for pos in (range(a,b) if t['strand']=='+' else range(b-1,a-1,-1))]
   expected=set()
   for end,seq in automaton.iter(t['protein'].replace('I','L')):
    bases=sorted(coordinates[3*(end-len(seq)+1):3*(end+1)])
    blocks=[]
    for _,group in groupby(enumerate(bases),lambda pair:pair[1]-pair[0]):
     run=[pos for _,pos in group];blocks.append((run[0],run[-1]+1))
    nt=''.join(fa.fetch(t['chrom'],a,b) for a,b in blocks)
    if t['strand']=='-':nt=str(Seq(nt).reverse_complement())
    assert str(Seq(nt).translate()).replace('I','L')==seq,(j['crop'],t['gene'],seq)
    expected.add((seq,t['chrom'],t['strand'],tuple(blocks)))
   row=records[t['gene']]
   assert len(expected)==int(row['expected_peptide_loci']),(t['gene'],'expected_count')
   assert len(expected&observed)==int(row['recovered_peptide_loci']),(t['gene'],'recovered_count')
   n+=len(expected);recovered+=len(expected&observed)
   for seq,chrom,strand,blocks in sorted(expected-observed):
    edges={(chrom,a[1],b[0],strand) for a,b in zip(blocks,blocks[1:]) if a[1]<b[0]}
    absent=sorted(edges-available)
    reason='junction_absent_from_search_evidence' if absent else 'exceeds_maximum_three_junctions' if len(edges)>3 else 'all_junctions_present_requires_graph_search_review' if edges else 'contiguous_mapping_requires_review'
    missed.append(dict(crop=j['crop'],gene=t['gene'],peptide=seq,chrom=chrom,strand=strand,blocks=blocks,missing_junctions=absent,reason=reason))
 result=dict(crop=j['crop'],targets=len(targets),independently_projected_and_translated_truth_loci=n,recovered=recovered,status='passed',heldout_annotation_overlap_check='passed',reused_RNA_and_original_assembly_BAM_identity_check='passed')
 audits.append(result);print(j['crop'],'passed',n,flush=True)
(root/'heldout_truth_audit.json').write_text(json.dumps(audits,indent=2)+'\n')

(root/'heldout_missed_peptide_loci.json').write_text(json.dumps(missed,indent=2)+'\n')
