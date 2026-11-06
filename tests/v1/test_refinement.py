import json
from pathlib import Path
import pytest
from Bio.Seq import Seq

from gapms_v1.cli import main as cli_main
from test_mapping import dna, table, independently_verify


def run_case(tmp_path,sequence,models,peptides,strand='+',extra=()):
    if strand=='-':
        length=len(sequence)
        sequence=str(Seq(sequence).reverse_complement())
        models=[[(length-b,length-a) for a,b in reversed(blocks)] for blocks in models]
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+sequence+'\n')
    annotation=tmp_path/'a.gtf';lines=[]
    for i,blocks in enumerate(models):
        cumulative=0
        for a,b in blocks:
            lines.append(f'chr1\tx\tCDS\t{a+1}\t{b}\t.\t{strand}\t{(-cumulative)%3}\tgene_id "gene"; transcript_id "t{i}";\n');cumulative+=b-a
    annotation.write_text(''.join(lines))
    queries=tmp_path/'p.txt';queries.write_text('\n'.join(peptides)+'\n');out=tmp_path/'out'
    assert main(['-a',str(genome),'-g',str(annotation),'-p',str(queries),'-o',str(out),'-t','1',*extra])==0
    independently_verify(genome,table(out/'diagnostics/mapping/mapped_peptides.tsv'))
    assert not (out/'diagnostics/orfs').exists()
    return out/'diagnostics',table(out/'diagnostics/refinement/annotation_guided_events.tsv')


@pytest.mark.parametrize('strand',['+','-'])
def test_internal_nontryptic_methionine_start_without_bam(tmp_path,strand):
    protein='MQQQQQQQQMAPEPTIDEKACDEFGHK*';sequence='C'*90+dna(protein)+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,90+3*len(protein))]],['MAPEPTIDEK','MQQQQQQQQ'],strand)
    starts=[r for r in rows if r['event_type']=='alternative_start']
    assert len(starts)==1
    row=starts[0];assert row['boundary_status']=='internal_ATG_nontryptic_N_terminus'
    assert row['upstream_peptide_support']=='1' # coexistence is allowed; absence is not required
    assert 'proteolytic_processing_not_excluded' in row['limitations']
    assert row['candidate_cds_blocks']=='[]' and row['candidate_model_id']==''
    assert row['evidence_origin']=='annotation_guided'
    assert len(table(out/'refinement/annotation_change_proposals.tsv'))==1
    assert 'start_codon' in (out/'refinement/annotation_guided_features.gtf').read_text()


@pytest.mark.parametrize('reason',['known_isoform','tryptic','multiple_loci','ambiguous_previous'])
def test_internal_start_filters(tmp_path,reason):
    prefix='MQQQQQQQK' if reason=='tryptic' else 'MQQQQQQQQ'
    protein=prefix+'MAPEPTIDEKACDEFGHK*';coding=dna(protein)
    if reason=='ambiguous_previous':coding=coding[:24]+'NNN'+coding[27:]
    sequence='C'*90+coding+'C'*90;blocks=[[(90,90+len(coding))]]
    if reason=='known_isoform':blocks.append([(117,90+len(coding))])
    if reason=='multiple_loci':sequence+='N'*90+coding
    out,rows=run_case(tmp_path,sequence,blocks,['MAPEPTIDEK'])
    assert not [r for r in rows if r['event_type']=='alternative_start']


@pytest.mark.parametrize('strand',['+','-'])
def test_n_terminal_extension_leaves_start_unresolved(tmp_path,strand):
    protein='MPEPTIDEKACDEFGHK*';sequence='C'*90+dna(protein)+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(117,90+3*len(protein))]],['PEPTIDEK'],strand)
    extensions=[r for r in rows if r['event_type']=='n_terminal_extension']
    assert len(extensions)==1 and extensions[0]['boundary_status']=='boundary_unresolved'
    assert extensions[0]['proposed_feature_blocks']=='[]'
    assert not [r for r in rows if r['event_type']=='alternative_start']


@pytest.mark.parametrize('strand',['+','-'])
def test_c_terminal_extension_infers_actual_stop(tmp_path,strand):
    protein='MQQQQQQQKACDEFGHK*';sequence='C'*90+dna(protein)+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,117)]],['ACDEFGHK'],strand)
    stops=[r for r in rows if r['event_type']=='alternative_stop']
    assert len(stops)==1 and stops[0]['boundary_status']=='genomic_stop_inferred'
    assert json.loads(stops[0]['proposed_feature_blocks'])
    assert 'not_directly_observed' in stops[0]['limitations']


def test_complete_existing_stop_cannot_be_crossed(tmp_path):
    sequence='C'*90+dna('MQQQQQQQK*ACDEFGHK*')+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,120)]],['ACDEFGHK'])
    assert not [r for r in rows if r['event_type'] in {'alternative_stop','c_terminal_extension'}]


def test_missing_stop_outside_window_is_unresolved(tmp_path):
    sequence='C'*90+dna('MQQQQQQQKACDEFGHKQQQQQ*')+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,117)]],['ACDEFGHK'],extra=['--refinement-window','30'])
    assert len(rows)==1 and rows[0]['event_type']=='c_terminal_extension'


@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('known',[False,True])
def test_local_junction_discovery_and_isoform_exclusion(tmp_path,strand,known):
    left=dna('MPEPT');right=dna('IDEKACDEFGHK*');intron='GTGGT'+'C'*25+'AG'
    sequence='C'*90+left+intron+right+'C'*90;b=90+len(left);c=b+len(intron);end=c+len(right)
    models=[[(90,b),(c,end)]]
    if known:models.append([(90,b+3),(c,end)])
    out,rows=run_case(tmp_path,sequence,models,['PEPTVIDEK'],strand)
    splice=[r for r in rows if r['event_type']=='alternative_splice']
    assert len(splice)==(0 if known else 1)
    if not known:
        assert splice[0]['boundary_status']=='junction_sequence_verified'
        assert json.loads(splice[0]['diagnostic_junctions'])
        assert json.loads((out/'reports/summary.json').read_text())['mapping']['mapped_peptides']==1


@pytest.mark.parametrize('strand',['+','-'])
def test_single_exon_skip_preserves_frame(tmp_path,strand):
    parts=[dna('MPEPT'),dna('AAAAAA'),dna('IDEKACDEFGHK*')];intron='GT'+'C'*26+'AG';sequence='C'*90;blocks=[]
    for i,part in enumerate(parts):
        a=len(sequence);sequence+=part;blocks.append((a,len(sequence)))
        if i<2:sequence+=intron
    sequence+='C'*90
    out,rows=run_case(tmp_path,sequence,[blocks],['PEPTIDEK'],strand)
    assert any(r['event_type']=='alternative_splice' for r in rows)


@pytest.mark.parametrize('strand',['+','-'])
def test_nearby_peptide_in_wrong_frame_is_not_an_extension(tmp_path,strand):
    prefix=dna('PEPTIDEK')+'A';ref=dna('MQQQQQQQK*');start=90+len(prefix)
    out,rows=run_case(tmp_path,'C'*90+prefix+ref+'C'*90,[[(start,start+len(ref))]],['PEPTIDEK'],strand)
    assert not rows


def test_unknown_sequence_does_not_establish_coding_connection(tmp_path):
    prefix=dna('PEPTIDEK')+'N'*6;ref=dna('MQQQQQQQK*');start=90+len(prefix)
    out,rows=run_case(tmp_path,'C'*90+prefix+ref+'C'*90,[[(start,start+len(ref))]],['PEPTIDEK'])
    assert not rows


def test_new_splice_placement_revokes_internal_start_uniqueness(tmp_path):
    protein='MQQQQQQQQMAPEPTIDEK*';first='C'*90+dna(protein)+'C'*90
    # The same peptide is encoded by a shifted junction at a second locus.
    left=dna('MMA');right=dna('EPTIDEK*');intron='GTGGT'+'C'*25+'AG'
    # +3 donor adds V, so use an internal MV... peptide with matching sequence.
    first='C'*90+dna('MQQQQQQQQMAPTVIDEK*')+'C'*90
    left=dna('MMAPT');right=dna('IDEK*')
    second_start=len(first)+90;b=second_start+len(left);c=b+len(intron);end=c+len(right)
    sequence=first+'C'*90+left+intron+right+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,90+len(dna('MQQQQQQQQMAPTVIDEK*')))],[(second_start,b),(c,end)]],['MAPTVIDEK'])
    mapped=table(out/'mapping/peptide_locus_summary.tsv')
    assert int(mapped[0]['loci'])>=2
    assert not [r for r in rows if r['event_type']=='alternative_start']


def test_bam_junction_evidence_strengthens_local_splice(tmp_path):
    import sqlite3,pysam
    from gapms_v1.cli import parser
    from gapms_v1.models import load_models
    from gapms_v1.mapping import EvidenceStore
    from gapms_v1.refinement import refine_annotations
    left=dna('MPEPT');right=dna('IDEKACDEFGHK*');intron='GTGGT'+'C'*25+'AG'
    b=90+len(left);c=b+len(intron);end=c+len(right)
    out,rows=run_case(tmp_path,'C'*90+left+intron+right+'C'*90,[[(90,b),(c,end)]],['PEPTVIDEK'])
    assert rows[0]['confidence']=='peptide_supported'
    args=parser().parse_args(['-a',str(tmp_path/'g.fa'),'-g',str(tmp_path/'a.gtf'),'-p',str(tmp_path/'p.txt'),'-o',str(out)])
    index=next((out/'work/reference').glob('*.fai'))
    store=object.__new__(EvidenceStore);store.db=sqlite3.connect(out/'work/evidence.sqlite')
    with pysam.FastaFile(str(tmp_path/'g.fa'),filepath_index=str(index)) as ref:
        models=load_models(tmp_path/'a.gtf','A01_a.gtf',ref,dict(zip(ref.references,ref.lengths)))
        strengthened,_=refine_annotations(args,models,store,{'PEPTVIDEK':dict(quality='unknown')},ref,{('chr1',b+3,c,'+'):4})
    store.db.close()
    assert strengthened[0]['confidence']=='rna_peptide_supported'
    assert strengthened[0]['minimum_rna_support']==4



def test_tryptic_explanation_in_another_isoform_blocks_start(tmp_path):
    a=dna('MQQQQQQQQ');b=dna('MKKKKKKKK');common=dna('MAPEPTIDEK*')
    start_a=90;start_b=start_a+len(a)+60;start_common=start_b+len(b)+60
    sequence='C'*90+a+'N'*60+b+'N'*60+common+'C'*90
    models=[[(start_a,start_a+len(a)),(start_common,start_common+len(common))],[(start_b,start_b+len(b)),(start_common,start_common+len(common))]]
    out,rows=run_case(tmp_path,sequence,models,['MAPEPTIDEK'])
    assert not [r for r in rows if r['event_type']=='alternative_start']
    stats=json.loads((out/'reports/summary.json').read_text())['annotation_guided_search']
    assert stats['start_matches_explained_by_another_isoform_cleavage']>=1



def test_known_stop_at_neighboring_gene_is_not_a_new_stop(tmp_path):
    first=dna('MQQQQQQQK');extension=dna('ACDEFGHK');neighbor=dna('MWWWWWWWK*')
    a=90;b=a+len(first);c=b+len(extension);end=c+len(neighbor)
    out,rows=run_case(tmp_path,'C'*90+first+extension+neighbor+'C'*90,[[(a,b)],[(c,end)]],['ACDEFGHK'])
    assert not [r for r in rows if r['event_type']=='alternative_stop']
    c_events=[r for r in rows if r['event_type']=='c_terminal_extension']
    assert c_events and c_events[0]['boundary_status']=='extension_with_already_annotated_end'


def main(argv):
    """Scientific fixtures retain intermediates they independently inspect."""
    return cli_main([*argv, "--keep-work"])
