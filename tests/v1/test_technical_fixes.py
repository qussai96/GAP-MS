"""Regression tests for tool discovery, BAM indexing and output-fidelity fixes."""
import json
import os
from pathlib import Path

import pysam
import pytest

from gapms_v1.cli import main as cli_main
from gapms_v1.consensus_engine import feedback_model
from gapms_v1.rnaseq import find_tool, indexed_bam
from test_mapping import dna, table
import test_refinement
from test_refinement import run_case


def executable(path, text='#!/bin/sh\nexit 0\n'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text);path.chmod(0o755)
    return path


def test_transdecoder_helpers_found_next_to_wrapper(tmp_path, monkeypatch):
    # TransDecoder >=7 conda layout: only the wrapper is on PATH.
    install=tmp_path/'opt/transdecoder'
    wrapper=executable(install/'TransDecoder')
    longorfs=executable(install/'util/TransDecoder.LongOrfs')
    executable(install/'util/TransDecoder.Predict')
    bin_dir=tmp_path/'bin';bin_dir.mkdir();(bin_dir/'TransDecoder').symlink_to(wrapper)
    monkeypatch.setenv('PATH',str(bin_dir))
    assert find_tool('TransDecoder.LongOrfs')==str(longorfs)
    assert find_tool('TransDecoder.Predict')==str(install/'util/TransDecoder.Predict')
    assert find_tool('stringtie') is None
    # An executable on PATH always takes precedence.
    direct=executable(bin_dir/'TransDecoder.LongOrfs')
    assert find_tool('TransDecoder.LongOrfs')==str(direct)


def test_transdecoder_helpers_missing_are_reported(tmp_path, monkeypatch):
    executable(tmp_path/'opt/TransDecoder');monkeypatch.setenv('PATH',str(tmp_path/'opt'))
    assert find_tool('TransDecoder.LongOrfs') is None


def write_bam(path, length=320, left='', right=''):
    with pysam.AlignmentFile(path,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':length}]}) as f:
        for i in range(3):
            read=pysam.AlignedSegment();read.query_name=f'r{i}';read.reference_id=0;read.reference_start=0;read.mapping_quality=60
            read.cigartuples=[(0,len(left)),(3,20),(0,len(right))] if left else [(0,30)]
            read.query_sequence=(left+right) if left else 'A'*30
            read.query_qualities=pysam.qualitystring_to_array('I'*len(read.query_sequence));read.set_tag('NH',1);read.set_tag('XS','+');f.write(read)


def test_bam_index_built_in_work_directory_only(tmp_path):
    bam=tmp_path/'input/rna.bam';bam.parent.mkdir();write_bam(bam)
    work=tmp_path/'work/rna_index'
    linked=indexed_bam(bam,work)
    assert linked.parent==work and os.path.samefile(linked,bam)
    assert Path(str(linked)+'.bai').is_file()
    assert not list(bam.parent.glob('*.bai')) and not list(bam.parent.glob('*.csi'))
    with pysam.AlignmentFile(str(linked),'rb',index_filename=str(linked)+'.bai') as reader:
        assert reader.count('chr1')==3
    assert indexed_bam(bam,work)==linked  # resume reuses the link and index


def test_existing_bam_index_is_used_unchanged(tmp_path):
    bam=tmp_path/'rna.bam';write_bam(bam);pysam.index(str(bam))
    assert indexed_bam(bam,tmp_path/'work')==bam
    assert not (tmp_path/'work').exists()


def test_unlinkable_bam_fails_with_instruction(tmp_path, monkeypatch):
    bam=tmp_path/'rna.bam';write_bam(bam)
    def refuse(*args):raise OSError(18,'Invalid cross-device link')
    monkeypatch.setattr('gapms_v1.rnaseq.os.link',refuse)
    with pytest.raises(ValueError,match='samtools index'):
        indexed_bam(bam,tmp_path/'work')


def test_consensus_with_unindexed_bam_completes(tmp_path, monkeypatch):
    left=dna('M'+'A'*49);right=dna('A'*49+'*');sequence=left+'GT'+'A'*16+'AG'+right
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+sequence+'\n')
    annotations=[]
    for name in ['a','b']:
        g=tmp_path/(name+'.gtf');g.write_text('chr1\tx\tCDS\t1\t150\t.\t+\t0\tgene_id "g"; transcript_id "t";\n'+'chr1\tx\tCDS\t171\t320\t.\t+\t0\tgene_id "g"; transcript_id "t";\n');annotations.append(g)
    bam=tmp_path/'rna.bam';write_bam(bam,left=left,right=right)
    peptides=tmp_path/'p.txt';peptides.write_text('MAAAAAAAA\n')
    def reconstruction(args,reference,lengths,commands):
        (args.output/'orfs').mkdir(exist_ok=True);(args.output/'transcripts').mkdir(exist_ok=True)
        path=args.output/'work/mock_rna_orfs.gtf';path.write_text('chr1\tx\tCDS\t1\t150\t.\t+\t0\tParent=candidate;gene_id=n\n')
        return path
    monkeypatch.setattr('gapms_v1.workflow.reconstruct',reconstruction)
    monkeypatch.setattr('gapms_v1.workflow.tool_versions',lambda bam:{'reconstruction':'mocked'})
    output=tmp_path/'out'
    argv=['-a',str(genome),'-p',str(peptides),'-g',str(annotations[0]),'-g',str(annotations[1]),'-b',str(bam),'-o',str(output),'-t','1']
    assert cli_main(argv)==0
    manifest=json.loads((output/'diagnostics/run_manifest.json').read_text())
    assert manifest['status']=='completed'
    assert manifest['consensus_bam_index']=='built_in_work_directory_from_hard_link_of_input_BAM'
    assert json.loads((output/'diagnostics/consensus/report.json').read_text())['rna']['available']
    assert not Path(str(bam)+'.bai').exists()
    assert not (output/'diagnostics/work').exists()  # temporary link and index removed


@pytest.mark.parametrize('policy',['tryptic','non_tryptic'])
def test_rejected_local_hypothesis_is_not_a_placement(tmp_path, monkeypatch, policy):
    # run_case's wrapper always adds --non-tryptic; select the policy explicitly.
    flags=['--keep-work']+(['--non-tryptic'] if policy=='non_tryptic' else [])
    monkeypatch.setattr(test_refinement,'main',lambda argv:cli_main([*argv,*flags]))
    # The local donor shift explains PEPTVIDEK, but the residue before it is A,
    # so it is not trypsin-compatible in the proposed context.
    left=dna('MAPEPT');right=dna('IDEKACDEFGHK*');intron='GTGGT'+'C'*25+'AG'
    sequence='C'*90+left+intron+right+'C'*90;b=90+len(left);c=b+len(intron);end=c+len(right)
    out,rows=run_case(tmp_path,sequence,[[(90,b),(c,end)]],['PEPTVIDEK'])
    mapped=[r for r in table(out/'mapping/mapped_peptides.tsv') if r['peptide']=='PEPTVIDEK']
    guided=[r for r in mapped if 'annotation_guided' in r['evidence']]
    if policy=='tryptic':
        assert guided==[] and not [r for r in rows if r['event_type']=='alternative_splice']
        assert json.loads((out/'reports/summary.json').read_text())['annotation_guided_search']['nontryptic_local_matches_excluded']>=1
    else:
        assert len(guided)==1 and [r for r in rows if r['event_type']=='alternative_splice']


def test_events_and_models_report_identified_peptide_spelling(tmp_path):
    protein='MQQQQQQQQMAPEPTIDEKACDEFGHK*';sequence='C'*90+dna(protein)+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,90+3*len(protein))]],['MAPEPTLDEK','MQQQQQQQQ'])
    start=[r for r in rows if r['event_type']=='alternative_start'][0]
    assert json.loads(start['peptides'])==['MAPEPTLDEK']
    assert json.loads(start['mapping_multiplicity'])=={'MAPEPTLDEK':1}
    assert [c['peptide'] for c in json.loads(start['cleavage_context'])]==['MAPEPTLDEK']
    assert table(out/'refinement/annotation_change_proposals.tsv')[0]['peptides']==start['peptides']
    model=table(out/'annotation/gene_model_evidence.tsv')[0]
    assert set(json.loads(model['peptides']))=={'MAPEPTLDEK','MQQQQQQQQ'}
    gene=table(out/'annotation/gene_evidence.tsv')[0]
    assert gene['distinct_peptides']=='2' and gene['locus_discriminating_peptides']=='2'


class Reference:
    def __init__(self, sequence):self.sequence=sequence
    def fetch(self, chrom, a, b):return self.sequence[a:b]


@pytest.mark.parametrize('protein,start,stop,issues',[
    ('MPEPTIDEK*',True,True,[]),
    ('QPEPTIDEK*',False,True,[]),
    ('MPEPTIDEKA',True,False,[]),
    ('MPEP*TIDEK*',True,True,['internal_stop']),
])
def test_feedback_graph_model_completeness_from_sequence(protein, start, stop, issues):
    sequence='CC'+dna(protein)+'CC'
    row=dict(CDS_0based_halfopen=json.dumps([[2,2+3*len(protein)]]),contig='chr1',strand='+')
    model=feedback_model('GS.gx.t1',row,Reference(sequence))
    assert (model.complete_start,model.complete_stop,model.issues)==(start,stop,issues)
    assert not model.protein.endswith('*')


def test_feedback_graph_model_partial_codon_is_flagged():
    sequence=dna('MPEPTIDEK')+'A'
    row=dict(CDS_0based_halfopen=json.dumps([[0,len(sequence)]]),contig='chr1',strand='+')
    model=feedback_model('GS.gx.t1',row,Reference(sequence))
    assert not model.complete_stop and 'partial_terminal_codon' in model.issues
