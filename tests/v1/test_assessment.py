from collections import defaultdict
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
import pysam

from gapms_v1.cli import main as cli_main
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
        path=args.output/'work/mock_rna_orfs.gtf';path.write_text(f'chr1\tx\tCDS\t131\t{len(sequence)}\t.\t+\t0\tParent=candidate;gene_id=n\n')
        return path
    monkeypatch.setattr('gapms_v1.workflow.reconstruct',reconstruction)
    monkeypatch.setattr('gapms_v1.workflow.tool_versions',lambda bam:{'synthetic_external_reconstruction':'mocked; BAM extraction and all mapping real'})
    output=tmp_path/'out';argv=['-a',str(assembly),'-p',str(peptides),'-b',str(bam),'-o',str(output),'-t','1']
    if with_annotation:argv+=['-g',str(reference)]
    assert main(argv)==0
    novel=table(output/'diagnostics/orfs/novel_orf_candidates.tsv')
    if with_annotation:
        assert len(novel)==1
        assert novel[0]['confidence']=='high_confidence_candidate'
        assert novel[0]['relationship']=='intergenic'
        assert novel[0]['distinct_diagnostic_peptides']=='2'
    else:
        assert novel==[]
        assert not (output/'review/novel_orfs.gff3').exists()
        assert len(table(output/'diagnostics/orfs/orf_peptide_evidence.tsv'))==1
        report=json.loads((output/'diagnostics/reports/summary.json').read_text())
        assert not report['assessment']['novelty_assessed']
    assert Path(output/'diagnostics/transcripts/rnaseq_junctions.tsv').exists()
    if with_annotation:
        assert len(table(output/'diagnostics/refinement/annotation_change_proposals.tsv'))==1
        assert len(table(output/'diagnostics/annotation/no_peptide_support_gene_models.tsv'))==1
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
    def model_hits(self, **kwargs):return defaultdict(list,self.hits)
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
    assert table(tmp_path/'orfs/novel_orf_candidates.tsv')==[]
    primary=table(tmp_path/'orfs/candidate_classification.tsv')
    assert len(primary)==1 and primary[0]['category']=='combined_refinement'
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
    assert table(tmp_path/'orfs/novel_orf_candidates.tsv')==[]
    assert table(tmp_path/'orfs/candidate_classification.tsv')[0]['category']=='alternative_splice'
    assert rows[0]['confidence']=='rna_peptide_supported'
    assert 'full_length' in rows[0]['limitations']
    elsewhere=hit('PEPTIDE',[(0,21)],0)
    assess(args,[ref,candidate],MockEvidence({'c':[elsewhere]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{('chr1',33,90,'+'):3})
    assert table(tmp_path/'refinement/alternative_splice_candidates.tsv')==[]


def test_single_peptide_and_multimapping_never_high_confidence(tmp_path):
    args=config(tmp_path)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,60)],coding=[(0,60)],protein='M'+'Q'*19,complete_start=True,complete_stop=True)
    placements=[hit('PEPTIDE',[(0,21)])]
    assess(args,[candidate,Model('unrelated','A01','unrelated','unrelated','chr1','+',cds=[(300,360)],coding=[(300,360)],protein='M'+'Q'*19)],MockEvidence({'c':placements},{'PEPTIDE':2}),{'PEPTIDE':dict(quality='qvalue_pass')},{})
    row=table(tmp_path/'orfs/novel_orf_candidates.tsv')[0]
    assert row['confidence']=='ambiguous'
    assert 'non_discriminating' in row['limitations']


def test_competing_orf_structures_prevent_high_confidence(tmp_path):
    args=config(tmp_path)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,90)],coding=[(0,90)],protein='M'+'Q'*29,complete_start=True,complete_stop=True)
    other=Model('other','RNA_ORF','o','g','chr1','+',cds=[(0,60),(90,120)],coding=[(0,60),(90,120)],protein='M'+'Q'*29,complete_start=True,complete_stop=True)
    placements=[hit('PEPTIDE',[(0,21)]),hit('ACDEFGH',[(30,51)],10)]
    evidence=MockEvidence({'c':placements,'other':placements},{'PEPTIDE':1,'ACDEFGH':1})
    assess(args,[candidate,other,Model('unrelated','A01','unrelated','unrelated','chr1','+',cds=[(300,360)],coding=[(300,360)],protein='M'+'Q'*19)],evidence,{p:dict(quality='qvalue_pass') for p in ['PEPTIDE','ACDEFGH']},{})
    rows=table(tmp_path/'orfs/novel_orf_candidates.tsv')
    assert rows and all(r['confidence']=='ambiguous' for r in rows)
    assert all('equally_compatible_coding_models' in r['limitations'] for r in rows)


def test_adjacent_blocks_are_one_genomic_placement(tmp_path):
    store=EvidenceStore(tmp_path/'evidence.sqlite')
    try:
        store.add('PEPTIDE','chr1','+',[(0,9),(9,21)],'PEPTIDE','annotation')
        store.add('PEPTIDE','chr1','+',[(0,21)],'PEPTIDE','contiguous')
        assert store.counts()=={'PEPTIDE':1}
        assert store.db.execute('SELECT blocks FROM hits').fetchone()[0]=='[[0,21]]'
    finally:store.db.close()


def test_unresolved_orf_residues_preclude_high_confidence(tmp_path):
    args=config(tmp_path)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,90)],coding=[(0,90)],protein='MPEPTIDEKACDEFGHKX',complete_start=True,complete_stop=True)
    placements=[hit('PEPTIDE',[(3,24)],1),hit('ACDEFGH',[(27,48)],9)]
    assess(args,[candidate,Model('unrelated','A01','unrelated','unrelated','chr1','+',cds=[(300,360)],coding=[(300,360)],protein='M'+'Q'*19)],MockEvidence({'c':placements},{'PEPTIDE':1,'ACDEFGH':1}),{p:dict(quality='qvalue_pass') for p in ['PEPTIDE','ACDEFGH']},{})
    row=table(tmp_path/'orfs/novel_orf_candidates.tsv')[0]
    assert row['confidence']!='high_confidence_candidate'
    assert 'unresolved_amino_acids_in_candidate_ORF' in row['limitations']


def test_alternative_frame_orf_is_not_boundary_refinement(tmp_path):
    args=config(tmp_path)
    ref=Model('ref','A01','ref','g','chr1','+',cds=[(31,91)],coding=[(31,91)],protein='M'+'Q'*19)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,120)],coding=[(0,120)],protein='MPEPTIDEK'+'Q'*31,complete_start=True,complete_stop=True)
    hs=[hit('PEPTIDE',[(3,24)],1)]
    assess(args,[ref,candidate],MockEvidence({'c':hs},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    assert len(table(tmp_path/'orfs/novel_orf_candidates.tsv'))==1
    assert table(tmp_path/'refinement/alternative_start_candidates.tsv')==[]
    assert table(tmp_path/'orfs/candidate_classification.tsv')[0]['category']=='novel_orf'


def test_incomplete_c_terminus_is_extension_not_stop(tmp_path):
    args=config(tmp_path)
    ref=Model('r','A01','r','g','chr1','+',cds=[(0,60)],coding=[(0,60)],protein='M'+'Q'*19)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,90)],coding=[(0,90)],protein='M'+'Q'*29,complete_start=True,complete_stop=False)
    placements=[hit('PEPTIDE',[(63,84)],21)]
    assess(args,[ref,candidate],MockEvidence({'c':placements},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    assert not table(tmp_path/'refinement/alternative_stop_candidates.tsv')
    row=table(tmp_path/'refinement/c_terminal_extension_candidates.tsv')[0]
    assert row['whole_model_status']=='local_evidence_only'
    assert 'incomplete_stop' in row['whole_model_review_reasons']
    assert 'termination_site_unresolved' in row['limitations']
    assert json.loads(row['local_peptide_blocks'])==[[63,84]]
    assert len(table(tmp_path/'refinement/local_evidence_only.tsv'))==1


def test_local_splice_evidence_does_not_validate_partial_whole_model(tmp_path):
    args=config(tmp_path)
    ref=Model('r','A01','r','g','chr1','+',cds=[(0,30),(90,150)],coding=[(0,30),(90,150)],protein='M'+'Q'*29)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(12,33),(90,150)],coding=[(12,33),(90,150)],protein='Q'*27,complete_stop=True)
    crossing=hit('PEPTIDE',[(21,33),(90,99)],3)
    assess(args,[ref,candidate],MockEvidence({'c':[crossing]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{('chr1',33,90,'+'):10})
    row=table(tmp_path/'refinement/alternative_splice_candidates.tsv')[0]
    assert row['confidence']=='rna_peptide_supported'
    assert row['whole_model_status']=='local_evidence_only'
    assert json.loads(row['diagnostic_junctions'])==[[33,90]]
    assert not table(tmp_path/'refinement/complete_models_for_manual_review.tsv')


@pytest.mark.parametrize('strand',['+','-'])
def test_boundaries_known_in_other_isoforms_are_not_new(tmp_path,strand):
    args=config(tmp_path)
    first=Model('a','A01','a','gene','chr1',strand,cds=[(0,60)],coding=[(0,60)],protein='M'+'Q'*19)
    second=Model('b','A01','b','gene','chr1',strand,cds=[(30,120)],coding=[(30,120)],protein='M'+'Q'*29)
    candidate=Model('c','RNA_ORF','c','rna','chr1',strand,cds=[(0,120)],coding=[(0,120)],protein='MPEPTIDEK'+'Q'*23+'ACDEFGHK',complete_start=True,complete_stop=True)
    placements=[hit('PEPTIDE',[(3,24)],1),hit('ACDEFGH',[(96,117)],32)]
    assess(args,[first,second,candidate],MockEvidence({'c':placements},{'PEPTIDE':1,'ACDEFGH':1}),{p:dict(quality='unknown') for p in ['PEPTIDE','ACDEFGH']},{})
    assert not table(tmp_path/'refinement/annotation_change_proposals.tsv')
    row=table(tmp_path/'orfs/annotation_isoform_checks.tsv')[0]
    assert json.loads(row['known_start_isoform_ids'])==(['a'] if strand=='+' else ['b'])
    assert json.loads(row['known_end_isoform_ids'])==(['b'] if strand=='+' else ['a'])


def test_exact_alternative_isoform_is_logged_without_any_proposal(tmp_path):
    args=config(tmp_path)
    a=Model('a','A01','a','g','chr1','+',cds=[(0,60)],coding=[(0,60)],protein='M'+'Q'*19)
    b=Model('b','A01','b','g','chr1','+',cds=[(0,30),(90,120)],coding=[(0,30),(90,120)],protein='M'+'Q'*19)
    candidate=Model('c','RNA_ORF','c','rna','chr1','+',cds=b.cds,coding=b.coding,protein=b.protein)
    placement=hit('PEPTIDE',[(21,30),(90,102)],7)
    assess(args,[a,b,candidate],MockEvidence({'b':[placement],'c':[placement]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    assert not table(tmp_path/'refinement/annotation_change_proposals.tsv')
    row=table(tmp_path/'orfs/annotation_isoform_checks.tsv')[0]
    assert row['interpretation']=='exact_annotated_isoform'
    assert json.loads(row['exact_annotated_isoform_ids'])==['b']


def test_recombining_known_boundaries_and_junctions_does_not_invent_features(tmp_path):
    args=config(tmp_path)
    def model(mid,source,blocks):
        return Model(mid,source,mid,'gene','chr1','+',cds=blocks,coding=blocks,protein='M'+'Q'*39,complete_start=True,complete_stop=True)
    a=model('a','A01',[(0,30),(60,90),(120,150)])
    b=model('b','A01',[(0,60),(90,120),(150,180)])
    c=model('c','RNA_ORF',[(0,30),(60,120),(150,180)])
    placement=hit('PEPTIDE',[(81,102)],17)
    assess(args,[a,b,c],MockEvidence({'c':[placement]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{('chr1',30,60,'+'):5,('chr1',120,150,'+'):5})
    assert not table(tmp_path/'refinement/annotation_change_proposals.tsv')
    row=table(tmp_path/'orfs/annotation_isoform_checks.tsv')[0]
    assert json.loads(row['new_junctions'])==[]
    assert row['exact_annotated_isoform_ids']=='[]'


@pytest.mark.parametrize('with_support',[False,True])
def test_candidate_exports_require_peptides_and_keep_ambiguous_support(tmp_path,with_support):
    from Bio import SeqIO
    args=config(tmp_path)
    models=[Model(mid,'RNA_ORF',mid,mid,'chr1','+',cds=[(i*100,i*100+60)],coding=[(i*100,i*100+60)],protein='M'+'Q'*19,complete_start=True,complete_stop=True) for i,mid in enumerate(['unique','shared','unsupported'])]
    hs={'unique':[hit('PEPTIDE',[(0,21)])],'shared':[hit('ACDEFGH',[(100,121)])]} if with_support else {}
    stats=assess(args,models,MockEvidence(hs,{'PEPTIDE':1,'ACDEFGH':2}),{p:dict(quality='unknown') for p in ['PEPTIDE','ACDEFGH']},{})
    expected={'unique','shared'} if with_support else set()
    assert {r['model_id'] for r in table(tmp_path/'orfs/orf_peptide_evidence.tsv')}==expected
    assert {r.id for r in SeqIO.parse(tmp_path/'orfs/candidate_orfs.faa','fasta')}==expected
    def gtf_ids(path):
        from gapms_v1.models import attributes
        return {attributes(line.rstrip().split('\t')[8])['transcript_id'] for line in path.read_text().splitlines()}
    assert gtf_ids(tmp_path/'orfs/candidate_orfs.gtf')==expected
    assert gtf_ids(tmp_path/'orfs/all_predicted_orfs.gtf')=={'unique','shared','unsupported'}
    assert len(table(tmp_path/'orfs/all_predicted_orf_evidence.tsv'))==3
    assert len(table(tmp_path/'orfs/orfs_without_peptide_support.tsv'))==3-len(expected)
    assert stats['predicted_orfs']==3 and stats['candidate_orfs']==len(expected)
    assert stats['orfs_without_peptide_support']==3-len(expected)
    assert stats['distinct_candidate_cds_paths']==len(expected)
    assert all(int(r['distinct_peptides'])>0 for r in table(tmp_path/'orfs/partial_or_ambiguous_orfs.tsv'))
    if with_support:
        assert next(r for r in table(tmp_path/'orfs/orf_peptide_evidence.tsv') if r['model_id']=='shared')['category']=='ambiguous'


def main(argv):
    """Scientific fixtures retain intermediates they independently inspect."""
    return cli_main([*argv, "--keep-work", "--non-tryptic"])


@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('coding',[True,False])
def test_known_junction_in_other_frame_or_gene_is_not_a_new_splice(tmp_path,strand,coding):
    args=config(tmp_path)
    def model(mid,source,blocks):
        return Model(mid,source,mid,mid,'chr1',strand,cds=blocks,
            coding=blocks if strand=='+' else list(reversed(blocks)),protein='M'+'Q'*21)
    ref=model('ref','A01',[(0,30),(90,120)])
    # Shares the candidate's junction but differs in phase throughout its CDS.
    other=model('other','A01',[(1,33),(90,121)])
    if not coding:
        other.exons=list(other.cds);other.cds=[];other.coding=[];other.protein=''
    candidate=model('candidate','RNA_ORF',[(0,33),(90,120)])
    h=hit('PEPTIDE',[(21,33),(90,99)],7);h['strand']=strand
    store=MockEvidence({'candidate':[h]},{'PEPTIDE':1})
    metadata={'PEPTIDE':dict(quality='unknown')}
    assert relationship(candidate,[other])==('overlapping_alternative_frame' if coding else 'UTR_associated')
    assess(args,[ref,candidate],store,metadata,{('chr1',33,90,strand):3})
    assert len(table(tmp_path/'refinement/alternative_splice_candidates.tsv'))==1
    assess(args,[ref,other,candidate],store,metadata,{('chr1',33,90,strand):3})
    assert not table(tmp_path/'refinement/alternative_splice_candidates.tsv')
    row=table(tmp_path/'orfs/annotation_isoform_checks.tsv')[0]
    assert json.loads(row['new_junctions'])==[]


def test_nontryptic_annotation_placement_still_excludes_false_novelty(tmp_path):
    args=config(tmp_path)
    ref=Model('ref','A01','ref','g','chr1','+',cds=[(0,90)],coding=[(0,90)],protein='M'+'Q'*29)
    candidate=Model('c','RNA_ORF','c','g','chr1','+',cds=[(0,60),(90,120)],coding=[(0,60),(90,120)],protein='M'+'Q'*29)
    placement=hit('PEPTIDE',[(3,24)],1)
    class DigestionEvidence(MockEvidence):
        def model_hits(self,all_mappings=False,model_ids=None):
            return defaultdict(list,{'c':[placement],**({'ref':[placement]} if all_mappings else {})})
    stats=assess(args,[ref,candidate],DigestionEvidence({}, {'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    assert table(tmp_path/'orfs/novel_orf_candidates.tsv')==[]
    checks=table(tmp_path/'orfs/annotation_isoform_checks.tsv')
    assert checks[0]['interpretation']=='existing_isoforms_explain_all_observed_peptides'
    assert stats['candidate_orfs']==1
