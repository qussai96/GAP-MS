"""Frozen-engine integration, coordinate fidelity and evidence-universe tests."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest
from Bio.Seq import Seq

from gapms_v1.cli import main as cli_main, parser
from gapms_v1.consensus_engine import prepare,finish,write_evidence,normalize_sources
from gapms_v1.mapping import EvidenceStore
from gapms_v1.models import load_models
from gapms_v1._genesensus.peptides import read_bed
from gapms_v1._genesensus.pipeline import run as native_run
from gapms_v1._genesensus.annotation import read_annotation
from test_mapping import dna,table,independently_verify


def fixture(tmp_path,strand='+',compressed=False):
    sequence='C'*30+dna('MPEPTIDEK'+'ACDEFGHIK'*6+'*')+'C'*30
    if strand=='-':sequence=str(Seq(sequence).reverse_complement())
    genome=tmp_path/'genome.fa';genome.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(genome))
    annotations=[]
    for name in ['helixer','annevo']:
        text=f'chr1\tx\tCDS\t31\t{len(sequence)-30}\t.\t{strand}\t0\tgene_id "g"; transcript_id "t";\n'
        path=tmp_path/(name+'.gtf')
        if compressed:
            import gzip
            path=path.with_suffix('.gtf.gz')
            with gzip.open(path,'wt') as f:f.write(text)
        else:path.write_text(text)
        annotations.append(path)
    peptides=tmp_path/'peps.txt';peptides.write_text('PEPTIDEK\nACDEFGHIK\n')
    return genome,annotations,peptides


@pytest.mark.parametrize('strand',['+','-'])
def test_automatic_consensus_order_resume_and_native_equivalence(tmp_path,strand):
    genome,annotations,peptides=fixture(tmp_path,strand,compressed=True)
    out=tmp_path/'out'
    cli=['-a',str(genome),'-p',str(peptides),'-g',str(annotations[0]),'-g',str(annotations[1]),'-o',str(out),'-t','1']
    assert main(cli)==0
    summary=json.loads((out/'diagnostics/reports/summary.json').read_text())['consensus']
    assert summary['engine']=='GeneSensus' and summary['engine_version']=='0.6.0a3'
    assert summary['primary_models']==1 and summary['additional_isoforms']==0 and not summary['reference_used']
    independently_verify(genome,table(out/'diagnostics/mapping/mapped_peptides.tsv'))
    mm=read_annotation(out/'diagnostics/consensus/annotation.gff3');assert len(mm)==1
    native=tmp_path/'native'
    native_run(out/'diagnostics/work/consensus_inputs/config.json',native,threads=1,cache=out/'diagnostics/work/genesensus_cache')
    for name in ['annotation.gff3','proteins.faa','primary_annotation.gff3','additional_isoforms.gff3','graph_report.json']:
        assert (native/name).read_bytes()==(out/'diagnostics/consensus'/name).read_bytes()
    second=tmp_path/'second'
    assert main(['-a',str(genome),'-p',str(peptides),'-g',str(annotations[1]),'-g',str(annotations[0]),'-o',str(second),'-t','1'])==0
    assert (second/'diagnostics/consensus/annotation.gff3').read_bytes()==(out/'diagnostics/consensus/annotation.gff3').read_bytes()
    assert main(cli+['--resume'])==0
    assert '--consensus ' not in parser().format_help()


def test_full_ambiguity_universe_and_failed_identifications(tmp_path):
    store=EvidenceStore(tmp_path/'e.sqlite')
    for a in [0,100]:store.add('MPEPTIDEK','chr1','+',[(a,a+27)],'MPEPTIDEK','source')
    store.add('ACDEFGHIK','chr1','+',[(200,227)],'ACDEFGHIK','source')
    store.db.commit();bed=tmp_path/'e.bed'
    result=write_evidence(store,{'MPEPTIDEK':{'quality':'unknown'},'ACDEFGHIK':{'quality':'qvalue_fail'}},bed)
    assert result==dict(placements=2,excluded_failed_identification_classes=1)
    evidence,report=read_bed(bed,keep_ambiguous=True)
    assert len(bed.read_text().splitlines())==2
    assert 'ACDEFGHLK' not in bed.read_text()
    assert report['ambiguous_sequences_retained']==1
    store.db.close()


def test_private_engine_matches_frozen_release_bytes():
    import gapms_v1._genesensus as gs
    root=Path(gs.__file__).parent;record=json.loads((root/'VENDORED.json').read_text())
    assert record['version']=='0.6.0a3'
    assert set(p.name for p in root.glob('*.py'))==set(record['original_sha256'])
    for name,digest in record['original_sha256'].items():
        content=(root/name).read_bytes()
        if name=='rna_geometry.py':
            content=content.replace(b'from .rna import',b'from genesensus.rna import').replace(b'from .genome import',b'from genesensus.genome import')
        assert hashlib.sha256(content).hexdigest()==digest,name
    # Importing RNA geometry must work without any external genesensus installation.
    from gapms_v1._genesensus.rna_geometry import augment_geometry
    assert callable(augment_geometry)


def test_normalization_preserves_initial_phase_coordinates_and_gene_family(tmp_path):
    genome,ann,peps=fixture(tmp_path)
    with pysam.FastaFile(genome) as reference:
        models=[]
        for i,path in enumerate(ann):models+=load_models(path,str(i),reference,{'chr1':reference.lengths[0]})
        # A nonzero initial phase is not silently trimmed and reset before selection.
        m=models[0];m.phases[m.cds[0]]=1;m.issues=['inconsistent_internal_CDS_phase']
        predictions,audit=normalize_sources(SimpleNamespace(gtf=ann),models,tmp_path)
        normalized=read_annotation(predictions['prediction_01']['annotation'])
        assert normalized[0].cds==tuple(m.cds) and normalized[0].phase==1
        assert audit[0]['original_issues']==m.issues
        assert normalized[0].identifier==m.id


def test_single_annotation_keeps_previous_workflow(tmp_path):
    genome,ann,peps=fixture(tmp_path)
    out=tmp_path/'out'
    assert main(['-a',str(genome),'-p',str(peps),'-g',str(ann[0]),'-o',str(out),'-t','1'])==0
    assert not (out/'diagnostics/consensus').exists()


def test_consensus_rejects_unsupported_genetic_code_before_inference(tmp_path):
    genome,ann,peps=fixture(tmp_path)
    out=tmp_path/'out'
    assert main(['-a',str(genome),'-p',str(peps),'-g',str(ann[0]),'-g',str(ann[1]),'-o',str(out),'--genetic-code','2','-t','1'])==1


def test_adapter_uses_real_bam_geometry_without_extra_predictor_vote(tmp_path):
    from dataclasses import replace
    from gapms.protein_mapper import ProteinMatcher
    left=dna('M'+'A'*49);right=dna('A'*49+'*')
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+left+'GT'+'A'*16+'AG'+right+'\n');pysam.faidx(str(genome))
    annotations=[]
    for name in ['a','b']:
        g=tmp_path/(name+'.gtf');g.write_text('chr1\tx\tCDS\t1\t150\t.\t+\t0\tgene_id "g"; transcript_id "t";\n'+'chr1\tx\tCDS\t171\t320\t.\t+\t0\tgene_id "g"; transcript_id "t";\n');annotations.append(g)
    bam=tmp_path/'rna.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':320}]}) as f:
        for i in range(3):
            read=pysam.AlignedSegment();read.query_name=f'r{i}';read.reference_id=0;read.reference_start=0;read.mapping_quality=60;read.cigartuples=[(0,150),(3,20),(0,150)];read.query_sequence=left+right;read.query_qualities=pysam.qualitystring_to_array('I'*300);read.set_tag('NH',1);read.set_tag('XS','+');f.write(read)
    pysam.index(str(bam));peptides=tmp_path/'p.txt';peptides.write_text('MAAAAAAAA\n')
    output=tmp_path/'out';(output/'work').mkdir(parents=True)
    args=parser().parse_args(['-a',str(genome),'-p',str(peptides),'-g',str(annotations[0]),'-g',str(annotations[1]),'-b',str(bam),'-o',str(output),'-t','2'])
    store=EvidenceStore(output/'work/evidence.sqlite');store.add('MAAAAAAAA','chr1','+',[(0,27)],'MAAAAAAAA','genome');store.db.commit()
    metadata={'MAAAAAAAA':{'quality':'unknown'}}
    with pysam.FastaFile(genome) as reference:
        models=[m for i,g in enumerate(annotations) for m in load_models(g,str(i),reference,{'chr1':320})]
        models.append(replace(models[0],source='RNA_ORF',id='rna_orf'))
        ctx,graphs,reports,seconds=prepare(args,models,store,metadata,reference)
        result=finish(args,ctx,store,metadata,reference,{},reports,seconds,matcher=ProteinMatcher(['MAAAAAAAA']))
    config=json.loads(ctx['config_path'].read_text());assert config['rna']==dict(bam=str(bam),library='unstranded')
    assert len(config['predictions'])==2
    report=json.loads((output/'consensus/report.json').read_text());assert report['rna']['available']
    assert result['primary_models']==1
    assert (output/'consensus/annotation.gff3').read_bytes()==(output/'consensus/primary_annotation.gff3').read_bytes()
    store.db.close()


def main(argv):
    """Scientific fixtures retain intermediates they independently inspect."""
    return cli_main([*argv, "--keep-work"])
