"""Consensus regression tests: positive cases and scientific failure modes."""
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
import pysam
from Bio.Seq import Seq

from gapms_v1.cli import main,parser
from gapms_v1.consensus import (make_candidate,merge_blocks,peptide_offset,overlap_nt,same_frame,
    Observations,source_features,select_primary,fusion_reason,reconstruct_graph,finish,Candidate)
from gapms_v1.consensus_rna import RNAEvidence,subtract
from gapms_v1.mapping import EvidenceStore
from gapms_v1.models import Model,load_models
from test_mapping import dna,table,independently_verify


def args(tmp_path,**kw):
 d=dict(output=tmp_path,genetic_code=1,bam=None,min_mapq=10,min_read_anchor=8,rna_library='unstranded',consensus_max_nodes=2000,consensus_max_starts=32,consensus_max_candidates=32,consensus_beam=8)
 return SimpleNamespace(**dict(d,**kw))


def ref(tmp_path,sequence):
 path=tmp_path/'g.fa';path.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(path));return pysam.FastaFile(path)


def source(reference,blocks,source='A',strand='+',gene='g'):
 c=make_candidate('chr1',strand,blocks,reference,1,'source');c.model.source=source;c.model.gene=gene;c.model.id=source+c.id;c.members=[c.model];c.sources=(source,);return c


def store(tmp_path,peptides=()):
 s=EvidenceStore(tmp_path/'e.sqlite')
 for p,st,bb,matched in peptides:s.add(p,'chr1',st,bb,matched,'genome')
 s.db.commit();return s


def metadata(peps):return {p:dict(quality='unknown',peptide_id=p) for p in peps}


@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('cut',[10,11,12])
def test_exact_splice_geometry_and_phase(strand,cut):
 protein=dna('MPEPTIDEK*');blocks=[(0,cut),(50,50+len(protein)-cut)];bb=[(3,cut),(50,50+21-cut)]
 if strand=='-':
  blocks=sorted((100-b,100-a) for a,b in blocks);bb=sorted((100-b,100-a) for a,b in bb)
 assert peptide_offset(blocks,strand,bb)==3
 bad=[(a+1,b+1) for a,b in bb]
 assert peptide_offset(blocks,strand,bad) is None


def test_conserved_fractional_budget_and_nested_saturation(tmp_path):
 r=ref(tmp_path,dna('MPEPTIDEK*')+'N'*30+dna('MPEPTIDEK*'))
 c=source(r,[(0,30)]);s=store(tmp_path,[('PEPTIDEK','+',[(3,27)],'PEPTIDEK'),('PEPTIDEK','+',[(63,87)],'PEPTIDEK'),('EPTIDEK','+',[(6,27)],'EPTIDEK')]);o=Observations(s,metadata(['PEPTIDEK','EPTIDEK']));f=o.features(c)
 assert f['peptide_coding_nt']==22.5 # 3 bases at half weight plus 21 saturated bases
 assert 'PEPTLDEK' not in f['unique'] and 'EPTLDEK' not in f['unique'] # shorter than 9 aa
 s.db.close();r.close()


def test_source_votes_not_isoform_count(tmp_path):
 r=ref(tmp_path,dna('MPEPTIDEK*'));c=source(r,[(0,30)]);value=source_features(c,[c]);c.members.append(c.model);assert source_features(c,[c])==value
 d=source(r,[(0,30)],'B');assert source_features(c,[c,d])==2*value;r.close()


def test_fusion_guard_separate_source_genes(tmp_path):
 r=ref(tmp_path,dna('MPEPTIDEK*')+'N'*30+dna('MPEPTIDEK*'))
 a=source(r,[(0,30)],gene='one');b=source(r,[(60,90)],gene='two');f=source(r,[(0,90)],'B',gene='fused')
 assert fusion_reason(f,[a,b,f])=='spans_separate_source_genes';r.close()


def test_small_locus_exact_selection_and_large_bound(tmp_path):
 r=ref(tmp_path,'ATG'*100)
 cc=[source(r,[(0,30)],str(i)) for i in range(20)]
 for i,c in enumerate(cc):c.id=str(i);c.features={'utility':float(i+1)}
 chosen,solver,_=select_primary(cc[:10],args(tmp_path));assert solver=='exact' and [x.id for x in chosen]==['9']
 chosen,solver,_=select_primary(cc,args(tmp_path));assert solver=='greedy_bound' and [x.id for x in chosen]==['19'];r.close()


@pytest.mark.parametrize('strand',['+','-'])
def test_graph_combines_source_exons_preserves_strand_phase_and_complete_orf(tmp_path,strand):
 # Two choices at each of two exon edges. The cross combination is an intact ORF.
 seq=list('C'*300);seq[0:30]=dna('MPEPTIDEKA');seq[60:90]=dna('ACDEFGHIKL');seq[120:150]=dna('ACDEFGHIKL');seq[180:210]=dna('ACDEFGHIK*')
 for a,b in [(30,60),(90,180),(30,120),(150,180)]:seq[a:a+2]='GT';seq[b-2:b]='AG'
 sequence=''.join(seq);paths=[[(0,30),(60,90),(180,210)],[(0,30),(120,150),(180,210)]]
 # Both alternate branches are source paths; cross combination through contiguous edges must never create a novel jump.
 if strand=='-':sequence=str(Seq(sequence).reverse_complement());paths=[sorted((300-b,300-a) for a,b in pp) for pp in paths]
 r=ref(tmp_path,sequence);cc=[source(r,p,str(i),strand) for i,p in enumerate(paths)];s=store(tmp_path);obs=Observations(s,{})
 found,report=reconstruct_graph(cc,r,obs,args(tmp_path))
 source_jj={iv for c in cc for iv in __import__('gapms_v1.consensus',fromlist=['junctions']).junctions(c.blocks)}
 for c in found:
  assert c.model.complete_start and c.model.complete_stop and not c.model.issues
  assert set(__import__('gapms_v1.consensus',fromlist=['junctions']).junctions(c.blocks))<=source_jj
 _,bounded=reconstruct_graph(cc,r,obs,args(tmp_path,consensus_max_nodes=1));assert bounded['skipped']=='node_or_start_limit'
 s.db.close();r.close()


@pytest.mark.parametrize('strand',['+','-'])
def test_multigtf_automatic_outputs_and_input_order_stability(tmp_path,strand):
 seq='C'*30+dna('MPEPTIDEKACDEFGHIK*')+'C'*30
 if strand=='-':seq=str(Seq(seq).reverse_complement())
 genome=tmp_path/'genome.fa';genome.write_text('>chr1\n'+seq+'\n')
 ann=[]
 for name in ['helixer','annevo']:
  p=tmp_path/(name+'.gtf');p.write_text(f'chr1\tx\tCDS\t31\t{len(seq)-30}\t.\t{strand}\t0\tgene_id "g"; transcript_id "t";\n');ann.append(p)
 peps=tmp_path/'peps.txt';peps.write_text('PEPTIDEK\nACDEFGHIK\n');out=tmp_path/'out'
 cli=['-a',str(genome),'-p',str(peps),'-o',str(out),'-g',str(ann[0]),'-g',str(ann[1]),'-t','1']
 assert main(cli)==0
 summary=json.loads((out/'reports/summary.json').read_text());assert summary['consensus']['primary_models']==1
 decisions=table(out/'consensus/decisions.tsv');assert len(decisions)==1 and len(json.loads(decisions[0]['sources']))==2
 independently_verify(genome,table(out/'mapping/mapped_peptides.tsv'))
 assert not any(x in parser().format_help() for x in ['--consensus ','--assess-only'])
 second=tmp_path/'second';cli[cli.index(str(out))]=str(second);cli[cli.index(str(ann[0]))]=str(ann[1]);cli[-3]=str(ann[0])
 assert main(cli)==0
 assert (out/'consensus/annotation.gff3').read_bytes()==(second/'consensus/annotation.gff3').read_bytes()
 assert main(cli+['--resume'])==0


def test_single_annotation_no_consensus(tmp_path):
 genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+dna('MPEPTIDEK*')+'\n');p=tmp_path/'p';p.write_text('PEPTIDEK\n');g=tmp_path/'a.gtf';g.write_text('chr1\tx\tCDS\t1\t30\t.\t+\t0\tgene_id "g"; transcript_id "t";\n');out=tmp_path/'out'
 assert main(['-a',str(genome),'-p',str(p),'-g',str(g),'-o',str(out),'-t','1'])==0
 assert not (out/'consensus').exists()


def test_invalid_source_is_audited_not_emitted(tmp_path):
 r=ref(tmp_path,dna('MPEP*IDEK*'));c=source(r,[(0,30)]);s=store(tmp_path);result=finish(args(tmp_path),[c],s,{},r,{},[],0)
 assert result['selected_models']==0;assert table(tmp_path/'consensus/decisions.tsv')[0]['status']=='rejected';s.db.close();r.close()


def test_rna_requires_new_base_coverage_and_all_junctions(tmp_path):
 r=ref(tmp_path,'ATG'*100);p=source(r,[(0,30),(120,150)]);c=source(r,[(0,60),(120,150)],'B');e=RNAEvidence()
 e.junction_templates['chr1',60,120,'.']={('','a'),('','b'),('','c')}
 assert not e.supports_change(c,p)
 e.coverage['chr1','.']=[(30,60)];assert e.supports_change(c,p)
 e.junction_templates['chr1',60,120,'.'].pop();assert not e.supports_change(c,p);r.close()


def test_rna_intron_retention_needs_whole_spanning_templates(tmp_path):
 r=ref(tmp_path,'ATG'*100);p=source(r,[(0,30),(60,90)]);c=source(r,[(0,90)],'B');e=RNAEvidence();e.coverage['chr1','.']=[(0,90)]
 assert not e.supports_change(c,p)
 e.blocks['chr1','.']=[(0,90,('',str(i))) for i in range(3)];assert e.supports_change(c,p);r.close()


def test_rna_bam_filters_and_mate_deduplication(tmp_path):
 path=tmp_path/'reads.bam';header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':300}]}
 with pysam.AlignmentFile(path,'wb',header=header) as f:
  for i in range(8):
   x=pysam.AlignedSegment();x.query_name='pair' if i<2 else str(i);x.query_sequence='A'*40;x.flag=0;x.reference_id=0;x.reference_start=10;x.mapping_quality=60;x.cigar=[(0,20),(3,30),(0,20)];x.set_tag('NH',1)
   if i==3:x.flag=1024
   if i==4:x.flag=256
   if i==5:x.set_tag('NH',2)
   if i==6:x.mapping_quality=255;x.set_tags([])
   if i==7:x.mapping_quality=10
   f.write(x)
 pysam.index(str(path));r=ref(tmp_path,'ATG'*100);c=source(r,[(0,90)])
 e=RNAEvidence.from_bam(path,[[c]]);assert e.count('chr1','+',30,60)==2
 assert e.coverage['chr1','.']==[(10,30),(60,80)];r.close()


@pytest.mark.parametrize('strand',['+','-'])
def test_positive_graph_recombination_and_peptide_selected_primary(tmp_path,strand):
 sequence=list('C'*240)
 for (a,b),protein in [((0,30),'MPEMIDEKAA'),((60,90),'ACDEFGHIKL'),((120,150),'QWERTYVFA*'),((180,210),'VVFGHIKLM*')]:sequence[a:b]=dna(protein)
 for a,b in [(30,60),(90,120),(90,180)]:sequence[a:a+2]='GT';sequence[b-2:b]='AG'
 sequence=''.join(sequence);paths=[[(0,30),(60,90),(120,150)],[(9,30),(60,90),(180,210)]];target=[(0,30),(60,90),(180,207)]
 pep='MPEMIDEKAAACDEFGHIKLVVFGHIKLM'
 if strand=='-':sequence=str(Seq(sequence).reverse_complement());paths=[sorted((240-b,240-a) for a,b in bb) for bb in paths];target=sorted((240-b,240-a) for a,b in target)
 r=ref(tmp_path,sequence);cc=[source(r,bb,str(i),strand) for i,bb in enumerate(paths)];s=store(tmp_path,[(pep,strand,target,pep)]);obs=Observations(s,metadata([pep]))
 graphs,report=reconstruct_graph(cc,r,obs,args(tmp_path));assert len(graphs)==2
 c=[g for g in graphs if peptide_offset(g.blocks,strand,target)==0][0]
 result=finish(args(tmp_path),cc+graphs,s,metadata([pep]),r,{},[report],0)
 assert c.status=='primary' and result['graph_selected']==1
 assert not c.reasons
 assert result['primary_models']==1
 s.db.close();r.close()


@pytest.mark.parametrize('evidence',['shared','unique_split','ambiguous_split'])
def test_isoform_admission_requires_distinguishing_unique_event(tmp_path,evidence):
 sequence=list('C'*240)
 for (a,b),protein in [((0,30),'MPEPTIDEKA'),((60,90),'ACDEFGHIKL'),((120,150),'VVFGHIKLMN'),((180,210),'QWERTYVFA*')]:sequence[a:b]=dna(protein)
 for a,b in [(30,60),(90,180),(30,120),(150,180)]:sequence[a:a+2]='GT';sequence[b-2:b]='AG'
 r=ref(tmp_path,''.join(sequence));a=source(r,[(0,30),(60,90),(180,210)]);b=source(r,[(0,30),(120,150),(180,210)],'B')
 # More peptides support A; B needs a distinguishing junction peptide, not shared first-exon evidence.
 pep='MPEPTIDEKA';entries=[(pep,'+',[(0,30)],pep),('PTIDEKAACDEFGHIKL','+',[(9,30),(60,90)],'PTIDEKAACDEFGHIKL')]
 split='PTIDEKAVVFGHIKLMN';bb=[(9,30),(120,150)]
 if evidence!='shared':entries.append((split,'+',bb,split))
 if evidence=='ambiguous_split':entries.append((split,'+',[(300,321),(400,430)],split))
 s=store(tmp_path,entries);mm=metadata([x[0] for x in entries]);finish(args(tmp_path),[a,b],s,mm,r,{},[],0)
 if evidence=='unique_split':assert sum(c.status=='supported_alternative' for c in [a,b])==1
 else:assert sum(c.status=='supported_alternative' for c in [a,b])==0
 s.db.close();r.close()


def test_rna_can_admit_source_alternative_without_unique_splice_peptide(tmp_path,monkeypatch):
 sequence=list('C'*240)
 for (a,b),protein in [((0,30),'MPEPTIDEKA'),((60,90),'ACDEFGHIKL'),((120,150),'VVFGHIKLMN'),((180,210),'QWERTYVFA*')]:sequence[a:b]=dna(protein)
 for a,b in [(30,60),(90,180),(30,120),(150,180)]:sequence[a:a+2]='GT';sequence[b-2:b]='AG'
 r=ref(tmp_path,''.join(sequence));a=source(r,[(0,30),(60,90),(180,210)]);b=source(r,[(0,30),(120,150),(180,210)],'B')
 e=RNAEvidence();e.coverage['chr1','.']=[(0,30),(60,90),(120,150),(180,210)]
 for iv in [(30,60),(90,180),(30,120),(150,180)]:e.junction_templates[('chr1',*iv,'.')]={('',str(i)) for i in range(3)}
 monkeypatch.setattr(RNAEvidence,'from_bam',classmethod(lambda cls,*a,**k:e))
 s=store(tmp_path,[('MPEPTIDEKA','+',[(0,30)],'MPEPTIDEKA')]);finish(args(tmp_path),[a,b],s,metadata(['MPEPTIDEKA']),r,{},[],0)
 assert sum(c.status=='supported_alternative' for c in [a,b])==1
 s.db.close();r.close()


def test_failed_identifications_do_not_supply_hard_support(tmp_path):
 r=ref(tmp_path,dna('MPEPTIDEKACDEFGHIK*'));c=source(r,[(0,57)]);s=store(tmp_path,[('ACDEFGHIK','+',[(27,54)],'ACDEFGHIK')]);o=Observations(s,{'ACDEFGHIK':dict(quality='qvalue_fail')});f=o.features(c)
 assert len(f['hits'])==1 and not f['unique'] and f['peptide_coding_nt']==0
 s.db.close();r.close()


def test_source_junction_votes_require_codon_phase_not_just_coordinates(tmp_path):
 r=ref(tmp_path,'ATG'*80);a=source(r,[(0,30),(60,90)]);b=source(r,[(1,30),(60,90)],'B')
 # Incorrect source phase cannot supply a splice vote for A; only A's own source can.
 b.model.issues=[]
 own=source_features(a,[a]);together=source_features(a,[a,b]);assert together==own
 r.close()
