"""Policy tests use the real default CLI, without the legacy-fixture wrapper."""
import json

import pytest
from Bio.Seq import Seq

from gapms_v1.cli import main, parser
from gapms_v1.consensus_engine import write_evidence
from gapms_v1.mapping import EvidenceStore
from peplocus.sequence import canonical
from test_mapping import dna, table


def test_policy_visible_and_default():
    assert '--non-tryptic' in parser().format_help()
    assert not parser().parse_args(['-a','a','-p','p','-o','o']).non_tryptic


@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('annotations',[0,1,2])
def test_selected_evidence_and_dual_exports(tmp_path,strand,annotations):
    peptide='ACDEFGHK'
    proteins=['MK'+peptide+'A*','MA'+peptide+'A*']
    blocks=[];sequence=''
    for protein in proteins:
        sequence+='N'*30
        coding=dna(protein)
        if strand=='-':coding=str(Seq(coding).reverse_complement())
        blocks.append((len(sequence),len(sequence)+len(coding)))
        sequence+=coding
    genome=tmp_path/'genome.fa';genome.write_text('>chr1\n'+sequence+'\n')
    peps=tmp_path/'p.txt';peps.write_text(peptide+'\n')
    gtfs=[]
    for i in range(annotations):
        path=tmp_path/f'{i}.gtf'
        path.write_text(''.join(f'chr1\ttest\tCDS\t{a+1}\t{b}\t.\t{strand}\t0\tgene_id "g{j}"; transcript_id "t{j}";\n' for j,(a,b) in enumerate(blocks)))
        gtfs+=['-g',str(path)]
    tracks=[]
    for nontryptic in (False,True):
        out=tmp_path/('all' if nontryptic else 'tryptic')
        argv=['-a',str(genome),'-p',str(peps),'-o',str(out),'-t','1',*gtfs,'--keep-work']
        if nontryptic:argv+=['--non-tryptic']
        assert main(argv)==0
        mappings=table(out/'diagnostics/mapping/mapped_peptides.tsv')
        assert len(mappings)==2
        assert sum(r['tryptic_compatible']=='True' for r in mappings)==1
        summary=json.loads((out/'diagnostics/mapping/mapping_summary.json').read_text())['mapping']
        assert summary['analysis_evidence']==('all_exact' if nontryptic else 'tryptic')
        assert summary['analysis_unique_peptides']==(0 if nontryptic else 1)
        assert summary['analysis_multiple_peptides']==(1 if nontryptic else 0)
        if annotations:
            models=table(out/'diagnostics/annotation/gene_model_evidence.tsv')
            assert sum(int(m['distinct_peptides']) for m in models)==annotations*(2 if nontryptic else 1)
            scores=table(out/'diagnostics/classification/scores.tsv')
            assert any(int(float(row['mapped_peptides']))==1 for row in scores)
            if not nontryptic:assert any(int(float(row['mapped_peptides']))==0 for row in scores)
        tracks.append([(out/name).read_text() for name in ('peptides.bed','peptides.tryptic.bed')])
    assert tracks[0]==tracks[1]


def test_consensus_filters_off_model_placements_and_il_aliases(tmp_path):
    metadata={p:dict(quality='unknown') for p in ('ACDEFGIK','ACDEFGLK')}
    for strict in (True,False):
        store=EvidenceStore(tmp_path/f'{strict}.db',require_tryptic=strict)
        for peptide in metadata:
            store.add(peptide,'chr1','+',[[0,24]],peptide,'genome',tryptic=True)
            store.add(peptide,'chr1','+',[[100,124]],peptide,'genome')
        path=tmp_path/f'{strict}.bed'
        result=write_evidence(store,metadata,path)
        assert result['placements']==(1 if strict else 2)
        assert store.counts()=={canonical('ACDEFGIK'):1 if strict else 2}
        store.db.close()


@pytest.mark.parametrize('with_annotation',[False,True])
def test_bam_branch_uses_same_default(tmp_path,monkeypatch,with_annotation):
    from test_assessment import make_bam
    sequence=dna('MKACDEFGHKA*')+'N'*30+dna('MAACDEFGHKA*')
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+sequence+'\n')
    peptides=tmp_path/'p';peptides.write_text('ACDEFGHK\n')
    bam=tmp_path/'r.bam';make_bam(bam,sequence)
    annotation=tmp_path/'a.gtf'
    annotation.write_text('chr1\tx\tCDS\t1\t36\t.\t+\t0\tgene_id "g"; transcript_id "t";\n')
    def reconstruction(args,reference,lengths,commands):
        (args.output/'orfs').mkdir(exist_ok=True)
        (args.output/'transcripts').mkdir(exist_ok=True)
        path=args.output/'work/mock.gtf'
        path.write_text(f'chr1\tx\tCDS\t67\t{len(sequence)}\t.\t+\t0\tgene_id "rna"; transcript_id "rna";\n')
        return path
    monkeypatch.setattr('gapms_v1.workflow.reconstruct',reconstruction)
    monkeypatch.setattr('gapms_v1.workflow.tool_versions',lambda bam:{'reconstruction':'mocked'})
    for unrestricted in (False,True):
        out=tmp_path/f'out{unrestricted}'
        argv=['-a',str(genome),'-p',str(peptides),'-b',str(bam),'-o',str(out),'-t','1']
        if with_annotation:argv+=['-g',str(annotation)]
        if unrestricted:argv+=['--non-tryptic']
        assert main(argv)==0
        rna=table(out/'diagnostics/orfs/all_predicted_orf_evidence.tsv')
        assert len(rna)==1
        assert int(rna[0]['distinct_peptides'])==int(unrestricted)


def test_possible_internal_methionine_start_retained_by_default(tmp_path):
    protein='MQQQQQQQQMAPEPTIDEKACDEFGHK*'
    genome=tmp_path/'g';genome.write_text('>chr1\n'+dna(protein)+'\n')
    gtf=tmp_path/'a.gtf';gtf.write_text(f'chr1\tx\tCDS\t1\t{3*len(protein)}\t.\t+\t0\tgene_id "g"; transcript_id "t";\n')
    peptides=tmp_path/'p';peptides.write_text('MAPEPTIDEK\n')
    out=tmp_path/'out'
    assert main(['-a',str(genome),'-p',str(peptides),'-g',str(gtf),'-o',str(out),'-t','1'])==0
    events=table(out/'diagnostics/refinement/annotation_guided_events.tsv')
    assert any(e['event_type']=='alternative_start' for e in events)
    assert 'trypsin-compatible exact placements' in (out/'report.html').read_text()


def test_shared_locus_does_not_transfer_cleavage_between_isoforms(tmp_path):
    peptide='ACDEFGHK'
    sequence=dna('MK')+'N'*30+dna('MA')+'N'*30+dna(peptide+'A*')
    genome=tmp_path/'g';genome.write_text('>chr1\n'+sequence+'\n')
    gtf=tmp_path/'a.gtf';lines=[]
    for model,first in [('tryptic',0),('nontryptic',36)]:
        for a,b in [(first,first+6),(72,len(sequence))]:
            lines.append(f'chr1\tx\tCDS\t{a+1}\t{b}\t.\t+\t0\tgene_id "g"; transcript_id "{model}";\n')
    gtf.write_text(''.join(lines));peps=tmp_path/'p';peps.write_text(peptide+'\n')
    for unrestricted in (False,True):
        out=tmp_path/f'out{unrestricted}'
        argv=['-a',str(genome),'-p',str(peps),'-g',str(gtf),'-o',str(out),'-t','1']
        if unrestricted:argv+=['--non-tryptic']
        assert main(argv)==0
        rows=table(out/'all_proteins_scores.tsv')
        assert sorted(int(r['mapped_peptides']) for r in rows)==([1,1] if unrestricted else [0,1])
        assert sorted(int(r['protein_specific_peptides']) for r in rows)==([0,0] if unrestricted else [0,1])


@pytest.mark.parametrize('protein,peptide,compatible',[
    ('MAPEPTIDEKA*','APEPTIDEK',True),  # initiator methionine removal
    ('MKPEPTIDEKA*','PEPTIDEK',False),  # K-P cannot be cleaved
    ('MKACDEFGHKP*','ACDEFGHK',False),  # terminal K-P cannot be cleaved
    ('MKACDEFGH*','ACDEFGH',True),      # translated protein stop
])
def test_terminal_exceptions_and_proline(tmp_path,protein,peptide,compatible):
    genome=tmp_path/'g';genome.write_text('>chr1\n'+dna(protein)+'\n')
    gtf=tmp_path/'a';gtf.write_text(f'chr1\tx\tCDS\t1\t{len(protein)*3}\t.\t+\t0\tgene_id "g"; transcript_id "t";\n')
    peptides=tmp_path/'p';peptides.write_text(peptide+'\n');out=tmp_path/'out'
    assert main(['-a',str(genome),'-p',str(peptides),'-g',str(gtf),'-o',str(out),'-t','1'])==0
    row=table(out/'all_proteins_scores.tsv')[0]
    assert int(row['mapped_peptides'])==int(compatible)
    assert bool((out/'peptides.tryptic.bed').read_text())==compatible
