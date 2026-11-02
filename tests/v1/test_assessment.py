from collections import defaultdict
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
import pysam

from gapms_v1.cli import main
from gapms_v1.models import Model
from gapms_v1.assessment import assess, relationship, non_nested
from gapms_v1.mapping import EvidenceStore
from gapms_v1.common import stable_id
from test_mapping import dna, table


def make_bam(path,sequence,intron=None,reads=3):
    with pysam.AlignmentFile(path,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':len(sequence)}]}) as stream:
        for index in range(reads):
            read=pysam.AlignedSegment();read.query_name=f'r{index}';read.flag=0;read.reference_id=0;read.reference_start=0;read.mapping_quality=60
            read.cigartuples=[(0,10),(3,intron),(0,20)] if intron else [(0,min(30,len(sequence)))]
            read.query_sequence='A'*sum(n for op,n in read.cigartuples if op==0);read.query_qualities=pysam.qualitystring_to_array('I'*len(read.query_sequence));read.set_tag('XS','+')
            stream.write(read)


@pytest.mark.parametrize('with_annotation',[False,True])
def test_rna_and_integrated_end_to_end(tmp_path,monkeypatch,with_annotation):
    sequence=dna('MQQQQQQQK*')+'N'*100+dna('MPEPTIDEKACDEFGHK*')
    assembly=tmp_path/'g.fa';assembly.write_text('>chr1\n'+sequence+'\n')
    peptides=tmp_path/'p.tsv';peptides.write_text('sequence\tq-value\tmodification\nPEPTIDE\t0.001\t\nACDEFGH\t0.002\t\n')
    bam=tmp_path/'r.bam';make_bam(bam,sequence)
    reference=tmp_path/'a.gtf';reference.write_text('chr1\tx\tCDS\t1\t30\t.\t+\t0\tParent=ref;gene_id=g\n')
    def reconstruction(args,reference,lengths,commands):
        (args.output/'orfs').mkdir(exist_ok=True);(args.output/'transcripts').mkdir(exist_ok=True)
        path=args.output/'orfs/candidate_orfs.gtf';path.write_text(f'chr1\tx\tCDS\t131\t{len(sequence)}\t.\t+\t0\tParent=candidate;gene_id=n\n')
        return path
    monkeypatch.setattr('gapms_v1.workflow.reconstruct',reconstruction)
    monkeypatch.setattr('gapms_v1.workflow.tool_versions',lambda bam:{'synthetic_external_reconstruction':'mocked; BAM extraction and all mapping real'})
    output=tmp_path/'out';argv=['-a',str(assembly),'-p',str(peptides),'-b',str(bam),'-o',str(output),'-t','1']
    if with_annotation:argv+=['-g',str(reference)]
    assert main(argv)==0
    novel=table(output/'orfs/novel_orf_candidates.tsv');assert len(novel)==1
    assert novel[0]['confidence']=='high_confidence_candidate'
    assert novel[0]['relationship']=='intergenic'
    assert novel[0]['distinct_diagnostic_peptides']=='2'
    assert Path(output/'transcripts/rnaseq_junctions.tsv').exists()
    if with_annotation:
        assert len(table(output/'refinement/annotation_change_proposals.tsv'))==1
        assert len(table(output/'annotation/no_peptide_support_gene_models.tsv'))==1
    assert not Path(str(assembly)+'.fai').exists()
    assert not Path(str(bam)+'.bai').exists()


def test_model_relationships_and_non_nested_evidence():
    ref=Model('ref','ann','r','g','chr1','+',exons=[(0,30),(60,120)],cds=[(9,30),(60,90)],coding=[(9,30),(60,90)])
    assert relationship(Model('c','RNA_ORF','c','g','chr1','-',cds=[(9,30)],coding=[(9,30)]),[ref])=='antisense'
    assert relationship(Model('c','RNA_ORF','c','g','chr1','+',cds=[(31,59)],coding=[(31,59)]),[ref])=='intronic'
    assert relationship(Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,9)],coding=[(0,9)]),[ref])=='UTR_associated'
    assert relationship(Model('c','RNA_ORF','c','g','chr1','+',cds=[(10,28)],coding=[(10,28)]),[ref])=='overlapping_alternative_frame'
    assert non_nested(['PEPTIDE','EPTIDE','ACDEFGH','PEPTIDE'])==['ACDEFGH','PEPTIDE']


class MockEvidence:
    def __init__(self,hits,counts):self.hits=hits;self.multiplicity=counts
    def model_hits(self):return defaultdict(list,self.hits)
    def counts(self):return self.multiplicity


def hit(peptide,blocks,aa_start=0):
    return dict(peptide=peptide,canonical=peptide.replace('L','I'),chrom='chr1',strand='+',blocks=blocks,aa_start=aa_start,aa_end=aa_start+len(peptide),mapping_id=stable_id('L',peptide,blocks))


def config(tmp_path):
    for directory in ['mapping','orfs']:(tmp_path/directory).mkdir(exist_ok=True)
    return SimpleNamespace(output=tmp_path,gtf=['reference'],bam='reads',min_junction_support=2)


def test_boundary_events_need_diagnostic_peptides_and_not_absence(tmp_path):
    args=config(tmp_path)
    ref=Model('ref','A01','ref','g','chr1','+',cds=[(30,90)],coding=[(30,90)],protein='M'+'Q'*19)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,120)],coding=[(0,120)],protein='MPEPTIDEK'+'Q'*23+'ACDEFGHK',complete_start=True,complete_stop=True)
    hs=[hit('PEPTIDE',[(3,24)],1),hit('ACDEFGH',[(96,117)],32)]
    metadata={p:dict(quality='unknown') for p in ['PEPTIDE','ACDEFGH']}
    assess(args,[ref,candidate],MockEvidence({'c':hs},{'PEPTIDE':1,'ACDEFGH':1}),metadata,{})
    starts=table(tmp_path/'refinement/alternative_start_candidates.tsv');stops=table(tmp_path/'refinement/alternative_stop_candidates.tsv')
    assert len(starts)==1 and len(stops)==1
    assert all(r['confidence']!='high_confidence_candidate' for r in starts+stops)
    # Ordinary tryptic cleavage far from the N terminus is not alternative-start evidence.
    hs[0]['aa_start']=8;hs[0]['aa_end']=15
    assess(args,[ref,candidate],MockEvidence({'c':hs[:1]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    assert table(tmp_path/'refinement/alternative_start_candidates.tsv')==[]
    assert table(tmp_path/'refinement/alternative_stop_candidates.tsv')==[]


def test_novel_splice_requires_peptide_crossing_and_records_rna(tmp_path):
    args=config(tmp_path)
    ref=Model('r','A01','r','g','chr1','+',cds=[(0,30),(90,120)],coding=[(0,30),(90,120)],protein='M'+'Q'*19)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,33),(90,120)],coding=[(0,33),(90,120)],protein='M'+'Q'*20,complete_start=True,complete_stop=True)
    crossing=hit('PEPTIDE',[(21,33),(90,99)],7)
    assess(args,[ref,candidate],MockEvidence({'c':[crossing]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='qvalue_pass')},{('chr1',33,90,'+'):3})
    rows=table(tmp_path/'refinement/alternative_splice_candidates.tsv')
    assert len(rows)==1 and rows[0]['minimum_rna_support']=='3'
    assert rows[0]['confidence']=='rna_peptide_supported'
    assert 'full_length' in rows[0]['limitations']
    elsewhere=hit('PEPTIDE',[(0,21)],0)
    assess(args,[ref,candidate],MockEvidence({'c':[elsewhere]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{('chr1',33,90,'+'):3})
    assert table(tmp_path/'refinement/alternative_splice_candidates.tsv')==[]


def test_single_peptide_and_multimapping_never_high_confidence(tmp_path):
    args=config(tmp_path)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,60)],coding=[(0,60)],protein='M'+'Q'*19,complete_start=True,complete_stop=True)
    placements=[hit('PEPTIDE',[(0,21)])]
    assess(args,[candidate],MockEvidence({'c':placements},{'PEPTIDE':2}),{'PEPTIDE':dict(quality='qvalue_pass')},{})
    row=table(tmp_path/'orfs/novel_orf_candidates.tsv')[0]
    assert row['confidence']=='ambiguous'
    assert 'non_discriminating' in row['limitations']
