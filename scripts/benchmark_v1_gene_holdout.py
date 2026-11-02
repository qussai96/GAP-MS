"""Deterministic gene holdout; withheld structures never enter the GAPMS run.
RNA reconstruction is reused only from an audited annotation/peptide-independent
StringTie/TransDecoder run. Mapping, splice graph and evidence are rebuilt.
"""
import argparse,csv,json,hashlib,shutil,sqlite3
from pathlib import Path
from collections import defaultdict,Counter
from dataclasses import asdict
from unittest.mock import patch
import pysam
from gapms._peplocus.io import read_peptides,prepare_reference
from gapms._peplocus.sequence import canonical
from gapms.protein_mapper import ProteinMatcher
from gapms_v1.models import load_models,write_models,attributes
from gapms._peplocus.annotation import annotation_parents
from gapms_v1.assessment import non_nested
from gapms_v1.cli import parser
from gapms_v1.workflow import run
from gapms_v1.common import write_json,write_table,input_identity

p=argparse.ArgumentParser();p.add_argument('jobs',type=Path);p.add_argument('index',type=int);p.add_argument('--genes',type=int,default=100);a=p.parse_args()
j=json.loads(a.jobs.read_text())[a.index];old=Path(j['output']).with_name('GAPMS_V1_casanovo');root=a.jobs.parent/'heldout'/j['crop'];root.mkdir(parents=True,exist_ok=True)
refgtf=Path(j['assembly']).parents[2]/'reference'/(j['crop']+'.ref.gtf')
# Input locations: inputs/genomic/assembly -> inputs/reference.
if not refgtf.exists():
 candidates=list((Path(j['assembly']).parents[1]/'reference').glob('*.ref.gtf'))
 assert len(candidates)==1,candidates
 refgtf=candidates[0]
manifest=json.loads((old/'run_manifest.json').read_text())
commands=[r['command'] for r in manifest['commands'] if 'command' in r]
assert len(commands)>=3 and Path(commands[0][0]).name=='stringtie'
assert '-G' not in commands[0] and '-g' not in commands[0]
assert commands[0][1]==j['bam']
assert all('peptide' not in ' '.join(c).lower() for c in commands[:3])
oldrna=old/'orfs/candidate_orfs.gtf'
genome,index,lengths=prepare_reference(Path(j['assembly']),root/'reference_index')
peptides,_=read_peptides(j['peptides']);peptides=[x for x in peptides if len(x)>=7];matcher=ProteinMatcher(peptides)
# RefSeq reuses some organellar transcript IDs on distinct contigs. Namespace
# identities by contig/strand and explicitly quarantine unsupported overlapping CDS.
raw=[];parents={}
for line in refgtf.read_text().splitlines():
 if not line or line.startswith('#'):continue
 r=line.split('\t');at=attributes(r[8]);scope=(r[0],r[6])
 if r[2].lower() in {'mrna','transcript'}:parents[scope+(at.get('id') or at.get('transcript_id'),)]=at.get('gene_id') or at.get('parent')
 if r[2].lower() in {'cds','exon'}:
  for tid in annotation_parents(r[8]):raw.append((scope,tid,r,at))
features=defaultdict(list)
for scope,tid,r,at in raw:features[scope+(tid,r[2].lower())].append((int(r[3])-1,int(r[4])))
bad=set()
for k,blocks in features.items():
 blocks=sorted(set(blocks))
 if any(x[1]>y[0] for x,y in zip(blocks,blocks[1:])):bad.add(k[:3])
normalized_reference=root/'namespaced_reference.gtf'
with normalized_reference.open('w') as f:
 for scope,tid,r,at in raw:
  if scope+(tid,) in bad:continue
  gene=at.get('gene_id') or parents.get(scope+(tid,)) or tid
  r[8]='gene_id '+json.dumps(':'.join(scope+(gene,)))+'; transcript_id '+json.dumps(':'.join(scope+(tid,)))+';'
  f.write('\t'.join(r)+'\n')
write_json(root/'reference_normalization.json',dict(rule='all transcript and gene identifiers prefixed with contig:strand; coordinates and phases preserved',excluded_overlapping_feature_models=sorted(bad),input=input_identity(refgtf)))
with pysam.FastaFile(genome,filepath_index=index) as fa:truth=load_models(normalized_reference,'HOLDOUT_REFERENCE',fa,lengths)
bygene=defaultdict(list);placements=defaultdict(list)
for m in truth:
 bygene[m.gene].append(m)
 if not m.coding or m.issues:continue
 for end,(length,indices) in matcher.automaton.iter(m.protein.replace('I','L')):
  for i in indices:
   seq=canonical(matcher.peptides[i]);start=end-length+1
   placements[m.id].append(dict(peptide=seq,blocks=m.project(3*start,3*(end+1)),aa_start=start))
# At least two distinct, non-nested peptide sequences; eligibility ignores RNA recovery.
eligible=[]
for gene,models in bygene.items():
 choices=[m for m in models if len(non_nested(h['peptide'] for h in placements[m.id]))>=2]
 if choices:
  target=sorted(choices,key=lambda m:(-len(non_nested(h['peptide'] for h in placements[m.id])),m.id))[0]
  eligible.append(target)
order=lambda m:hashlib.sha256(('GAPMS-holdout-2026:'+m.gene).encode()).hexdigest()
single=sorted([m for m in eligible if not m.junctions],key=order);multi=sorted([m for m in eligible if m.junctions],key=order)
targets=single[:a.genes//2]+multi[:a.genes-a.genes//2]
if len(targets)<a.genes:targets+=sorted([m for m in eligible if m not in targets],key=order)[:a.genes-len(targets)]
assert targets
heldgenes={m.gene for m in targets};regions=defaultdict(list)
for m in truth:
 if m.gene in heldgenes:regions[m.chrom].append((min(x for x,y in m.intervals),max(y for x,y in m.intervals)))
def withheld(m):
 lo=min(x for x,y in m.intervals);hi=max(y for x,y in m.intervals)
 return m.gene in heldgenes or any(max(lo,x)<min(hi,y) for x,y in regions[m.chrom])
retained=[m for m in truth if not withheld(m)]
masked=root/'masked_reference.gtf';write_models(retained,masked)
write_json(root/'selection.json',dict(crop=j['crop'],requested_genes=a.genes,selected=len(targets),eligible=len(eligible),eligible_single_exon=len(single),eligible_multi_exon=len(multi),removed_models=len(truth)-len(retained),retained_models=len(retained),selection='SHA256(GAPMS-holdout-2026:gene); up to half single-exon and half multi-exon; >=2 non-nested peptides; RNA-independent',reference=input_identity(refgtf),target_models=[dict(**{k:v for k,v in asdict(m).items() if k!='phases'},phases=[dict(block=list(k),phase=v) for k,v in m.phases.items()]) for m in targets],rna_reconstruction_origin=str(old),rna_reconstruction_commands=commands[:3]))
args=parser().parse_args(['-a',j['assembly'],'-p',j['peptides'],'-g',str(masked),'-b',j['bam'],'-o',str(root/'run'),'-t','8'])
def reuse(args,reference,lengths,commands):
 (args.output/'orfs').mkdir(exist_ok=True);(args.output/'transcripts').mkdir(exist_ok=True)
 path=args.output/'orfs/candidate_orfs.gtf';shutil.copyfile(oldrna,path)
 commands.append(dict(component='annotation_independent_RNA_reconstruction_reuse',origin=str(old),input_identity=input_identity(oldrna),audited_commands=manifest['commands'][:3],reason='holdout never used during unguided StringTie/TransDecoder reconstruction'))
 return path
with patch('gapms_v1.workflow.reconstruct',reuse):run(args)
# Scoring starts only after the blinded run has completed.
with pysam.FastaFile(genome,filepath_index=index) as fa:predictions=load_models(root/'run/orfs/candidate_orfs.gtf','RNA_ORF',fa,lengths)
occ=defaultdict(set)
with (root/'run/mapping/model_peptide_occurrences.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):occ[r['model_id']].add(canonical(r['peptide']))
db=sqlite3.connect(root/'run/work/evidence.sqlite')
observed={(p,c,s,tuple(map(tuple,json.loads(b)))) for p,c,s,b in db.execute('SELECT canonical,chrom,strand,blocks FROM hits')};db.close()
def key(m,h):return h['peptide'],m.chrom,m.strand,tuple(map(tuple,h['blocks']))
def same_region(m,t):return m.chrom==t.chrom and m.strand==t.strand and any(max(x,u)<min(y,v) for x,y in m.coding for u,v in t.coding)
rows=[]
for t in targets:
 expected={key(t,h) for h in placements[t.id]};sequences={h['peptide'] for h in placements[t.id]}
 supported=[m for m in predictions if same_region(m,t) and occ[m.id]&sequences]
 exact=[m for m in supported if m.signature()==t.signature()]
 tj=set(t.junctions);pj={edge for m in supported for edge in m.junctions}
 rows.append(dict(gene=t.gene,target_transcript=t.original_id,chrom=t.chrom,strand=t.strand,cds_blocks=t.coding,target_peptides=len(sequences),expected_peptide_loci=len(expected),recovered_peptide_loci=len(expected&observed),supported_candidate_models=len(supported),exact_cds_recovered=bool(exact),exact_cds_model_ids=[m.id for m in exact],exact_splice_chain_recovered=bool(tj) and any(m.junctions==t.junctions for m in supported),start_recovered=any(m.coding[0][0 if t.strand=='+' else 1]==t.coding[0][0 if t.strand=='+' else 1] for m in supported),stop_recovered=any(m.coding[-1][1 if t.strand=='+' else 0]==t.coding[-1][1 if t.strand=='+' else 0] for m in supported),junction_tp=len(tj&pj),junction_fn=len(tj-pj),junction_fp=len(pj-tj),truth_junctions=len(tj)))
write_table(root/'target_recovery.tsv',list(rows[0]),rows)
sums={k:sum(r[k] for r in rows) for k in ['expected_peptide_loci','recovered_peptide_loci','exact_cds_recovered','exact_splice_chain_recovered','start_recovered','stop_recovered','junction_tp','junction_fn','junction_fp']}
summary=dict(crop=j['crop'],heldout_genes=len(rows),multi_exon_targets=sum(bool(r['truth_junctions']) for r in rows),**sums)
summary['peptide_locus_recall']=sums['recovered_peptide_loci']/sums['expected_peptide_loci']
summary['exact_cds_recall']=sums['exact_cds_recovered']/len(rows)
summary['junction_recall']=sums['junction_tp']/(sums['junction_tp']+sums['junction_fn']) if sums['junction_tp']+sums['junction_fn'] else None
summary['junction_precision']=sums['junction_tp']/(sums['junction_tp']+sums['junction_fp']) if sums['junction_tp']+sums['junction_fp'] else None
summary['limitations']=['Annotation-derived reference is a target, not experimentally proven full-length truth.','Only peptide-supported, >=2 non-nested-peptide genes were eligible; stratified sample is not a genome-wide prevalence estimate.','Reference-database peptide identification predates masking and is not a de novo identification benchmark.','Precision compares supported overlapping RNA candidate junctions to the selected target isoform; valid alternative isoforms may count as unmatched.','RNA reconstruction reused from prior unguided BAM-only reconstruction; genome mapping and annotation splice graph rebuilt from masked annotation.']
write_json(root/'recovery_summary.json',summary);print(json.dumps(summary,indent=2))
