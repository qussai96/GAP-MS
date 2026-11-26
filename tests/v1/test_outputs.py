import csv,json
from pathlib import Path
from urllib.parse import unquote
import pytest
from gapms_v1.cli import main
from gapms_v1.outputs import review_features,publish
from gapms_v1.assessment import EVENT_FIELDS
from test_mapping import dna,table


@pytest.mark.parametrize('multiple',[False,True])
def test_public_layout_no_redundant_annotation_and_default_cleanup(tmp_path,multiple):
    genome=tmp_path/'g.fa';genome.write_text('>c\n'+dna('MPEPTIDEK'+'ACDEFGHIK'*6+'*')+'\n')
    peptides=tmp_path/'p.tsv';peptides.write_text('sequence\nPEPTIDEK\nVVVVVVV\nAAA\n')
    ann=[]
    for name in ['a','b']:
        p=tmp_path/(name+'.gtf');p.write_text(f'c\tx\tCDS\t1\t{3*(9+9*6+1)}\t.\t+\t0\tgene_id "g"; transcript_id "t";\n');ann.append(p)
    out=tmp_path/'out';argv=['-a',str(genome),'-p',str(peptides),'-g',str(ann[0]),'-o',str(out),'-t','1']
    if multiple:argv+=['-g',str(ann[1])]
    assert main(argv)==0
    expected={'report.html','report.log','supported_genes.gtf','all_proteins_scores.tsv','peptides.bed','peptides.igv.bed','peptides.tryptic.bed','unmapped_peptides.tsv','diagnostics'}
    if multiple:expected|={'consensus.gff3','consensus_proteins.faa'}
    assert {p.name for p in out.iterdir()}==expected
    assert not (out/'diagnostics/work').exists()
    assert not (out/'annotation.gff3').exists()
    rows={r['peptide']:r for r in table(out/'unmapped_peptides.tsv')}
    assert rows['VVVVVVV']['status']=='unmapped_in_searched_space'
    assert rows['AAA']['status']=='not_searched_below_min_length'
    from gapms_v1.outputs import LEGACY_SCORE_FIELDS
    rows=table(out/'all_proteins_scores.tsv')
    with (out/'all_proteins_scores.tsv').open() as f:assert f.readline().strip().split('\t')==LEGACY_SCORE_FIELDS
    assert len(LEGACY_SCORE_FIELDS)==29
    assert all(r['supported'] in {'+','-'} for r in rows)
    assert all(r['confidence'] in {'high','moderate','no-support'} for r in rows)
    assert 'classifier_score' not in rows[0] if rows else True
    assert (out/'diagnostics/reports/overview.svg').exists()
    assert (out/'diagnostics/reports/overview.png').exists()
    figure=json.loads((out/'diagnostics/reports/overview_counts.json').read_text())
    assert figure['splice']+figure['internal']==figure['peptide_protein_pairs']
    assert 'data:image/svg+xml;base64,' in (out/'report.html').read_text()

    assert (out/'supported_genes.gtf').read_bytes()==(out/'diagnostics/classification/supported_genes.gtf').read_bytes()
    assert 'Results ready:' in (out/'report.log').read_text()
    assert (out/'peptides.bed').read_bytes()==(out/'diagnostics/mapping/mapped_peptides.bed').read_bytes()
    if multiple:assert (out/'consensus.gff3').read_bytes()==(out/'diagnostics/consensus/annotation.gff3').read_bytes()
    html=(out/'report.html').read_text();assert ('Multi-annotation consensus' if multiple else 'Single-annotation assessment') in html
    before=(out/'peptides.bed').read_bytes();assert main(argv+['--resume'])==0;assert (out/'peptides.bed').read_bytes()==before


def test_failed_run_retains_work_and_parent_log(tmp_path,monkeypatch):
    import gapms_v1.workflow as w
    def failure(args):
        (args.output/'work').mkdir(parents=True);(args.output/'work/proof').write_text('retain')
        raise RuntimeError('deliberate stage failure')
    monkeypatch.setattr(w,'_run',failure)
    out=tmp_path/'out'
    assert main(['-a','unused','-p','unused','-o',str(out)])==1
    assert (out/'diagnostics/work/proof').read_text()=='retain'
    assert 'deliberate stage failure' in (out/'report.log').read_text()
    assert not (out/'report.html').exists()


def test_local_proposal_gff_does_not_invent_transcript_boundaries(tmp_path):
    row={k:'' for k in EVENT_FIELDS};row.update(event_id='event1',event_type='alternative_stop',chrom='c',strand='-',confidence='peptide_supported',boundary_status='termination_unresolved',whole_model_status='local_evidence_only',local_peptide_blocks='[[12,30],[80,92]]',diagnostic_junctions='[[30,80]]',proposed_feature_blocks='[]')
    p=tmp_path/'proposals.gff3';review_features([row],p);lines=p.read_text().splitlines()[1:]
    assert {l.split('\t')[2] for l in lines}=={'peptide_supported_region','splice_junction'}
    assert not any('\tCDS\t' in l or '\tstop_codon\t' in l for l in lines)
    assert all('proposal_requires_review=true' in l for l in lines)
    row.update(event_type='alternative_start',proposed_feature_blocks='[[12,15]]')
    review_features([row],p);assert '\tstart_codon\t13\t15\t' in p.read_text()


def test_rescue_only_is_separate_from_exact_unmapped(tmp_path):
    g=tmp_path/'g.fa';g.write_text('>c\n'+dna('MPEPTIDEKACDEFGHIK*')+'\n')
    p=tmp_path/'p.txt';p.write_text('MAPEPTIDEKACDEFGHIK\n');out=tmp_path/'out'
    assert main(['-a',str(g),'-p',str(p),'-o',str(out),'--rescue-one-edit','-t','1'])==0
    rows=table(out/'unmapped_peptides.tsv');assert len(rows)==1
    assert rows[0]['status']=='rescue_only' and int(rows[0]['rescue_loci'])>=1
    assert not (out/'peptides.bed').read_text()
    assert (out/'diagnostics/mapping/rescue/rescue_mappings.tsv').exists()


def test_refresh_preserves_scores_and_deduplicates_novel_proposals(tmp_path):
    from gapms_v1.outputs import refresh_report
    genome=tmp_path/'g.fa';genome.write_text('>c\n'+dna('MPEPTIDEKACDEFGHIK*')+'\n')
    peptides=tmp_path/'p';peptides.write_text('PEPTIDEK\n')
    out=tmp_path/'out'
    assert main(['-a',str(genome),'-p',str(peptides),'-o',str(out),'-t','1'])==0
    diag=out/'diagnostics'
    row={k:'' for k in EVENT_FIELDS};row.update(event_id='same_event',event_type='novel_orf',chrom='c',start=0,end=27,strand='+',confidence='ambiguous',local_peptide_blocks='[[3,27]]',diagnostic_junctions='[]',candidate_cds_blocks='[]',peptides='["PEPTIDEK"]')
    import csv
    for relative in ['refinement/annotation_change_proposals.tsv','orfs/novel_orf_candidates.tsv']:
        target=diag/relative;target.parent.mkdir(exist_ok=True)
        with target.open('w') as f:
            w=csv.DictWriter(f,delimiter='\t',fieldnames=EVENT_FIELDS);w.writeheader();w.writerow(row)
    bed=(out/'peptides.bed').read_bytes()
    provenance=refresh_report(out)
    assert len(table(out/'review/candidates.tsv'))==1
    assert (out/'peptides.bed').read_bytes()==bed
    assert '1 proposal rows' in (out/'report.html').read_text()
    counts=json.loads((diag/'reports/overview_counts.json').read_text())
    assert counts['events']['novel_orf']==1
    assert counts['events_ambiguous']['novel_orf']==1
    assert counts['events_specific'].get('novel_orf',0)==0
    assert '0 supported by specific peptides; 1 ambiguous or not supported' in (out/'report.html').read_text()
    assert provenance['preserved_core_implementation_sha256']
