from types import SimpleNamespace
import pandas as pd
from gapms_v1.classification import classify,features
from gapms.extract import engineer_features


def frame():
    rows=[]
    for i in range(12):
        high=i<4;low=4<=i<8
        rows.append(dict(Protein=str(i),model_id=str(i),protein_length=100 if low else 300,
            mapped_peptides=5 if high else 0 if low else 2,protein_specific_peptides=2 if high else 0,
            gene_specific_peptides=3 if high else 0,sequence_coverage=.9 if high else .1,
            N_terminal_peptides=0,C_terminal_peptides=0,splice_peptides=0,internal_peptides=2,
            splice_sites=0,supported_splice_sites=0))
    return engineer_features(pd.DataFrame(rows))


def test_iterative_training_saves_model_and_preserves_seeds(tmp_path):
    result,audit=classify(frame(),tmp_path,1)
    assert audit['status']=='trained'
    assert result.loc[:3,'supported'].all()
    assert not result.loc[4:7,'supported'].any()
    assert result.classifier_score.between(0,1).all()
    assert (tmp_path/'model.json').exists()
    assert len(audit['iterations'])<=5


def test_insufficient_seeds_remain_unresolved_and_high_wins_overlap(tmp_path):
    df=frame().iloc[:1].copy();df.loc[0,['protein_length','mapped_peptides','gene_specific_peptides']]=[10,1,0]
    result,audit=classify(df,tmp_path,1)
    assert audit['status']=='insufficient_seeds'
    assert result.loc[0,'supported']
    assert result.loc[0,'seed_label']=='high'
    assert result.loc[0,'classifier_score']=='NA'


def test_il_dedup_repeated_motif_and_quality():
    model=SimpleNamespace(id='m',source='a',original_id='t',gene='g',chrom='c',strand='+',
        protein='MABCDEFGHIJK',cds=[(0,12),(30,54)],coding=[(0,12),(30,54)])
    hits={'m':[dict(peptide=p,canonical=c,aa_start=a,aa_end=b) for p,c,a,b in
        [('IL','LL',0,2),('LL','LL',0,2),('REP','REP',2,5),('REP','REP',6,9),('BAD','BAD',9,12)]]}
    meta={p:dict(quality='qvalue_fail' if p=='BAD' else 'unknown') for p in ['IL','LL','REP','BAD']}
    row=features([model],hits,meta).iloc[0]
    assert row.mapped_peptides==2
    assert row.splice_peptides==0
    assert row.supported_splice_sites==0
    assert row.N_terminal_peptides==1
    assert row.C_terminal_peptides==0


def test_features_match_legacy_extractor(tmp_path,monkeypatch):
    import pysam
    import gapms.extract as legacy
    from gapms_v1.models import load_models,write_models,write_proteins
    from test_mapping import dna
    genome=tmp_path/'g.fa';genome.write_text('>c\n'+dna('MPEPTIDEKACDEFGHIK*')+'\n')
    pysam.faidx(str(genome))
    ann=tmp_path/'a.gtf'
    ann.write_text('c\tx\tCDS\t1\t54\t.\t+\t0\tgene_id "g"; transcript_id "t";\n')
    with pysam.FastaFile(str(genome)) as ref:models=load_models(ann,'a',ref,{'c':54})
    m=models[0];write_models(models,ann);faa=tmp_path/'p.fa';write_proteins(models,faa)
    peptides=[(p,m.protein.index(p),m.protein.index(p)+len(p)) for p in ['PEPTIDEK','ACDEFGHIK']]
    mapping=tmp_path/'mapping.tsv';mapping.write_text('peptide\tprotein\tlocation\n'+''.join(f'{p}\t{m.id}\t{a+1}\n' for p,a,b in peptides))
    monkeypatch.setattr(legacy,'plot_mapped_percentage_bars',lambda *a,**k:None)
    legacy.extract_features(ann,faa,mapping,tmp_path)
    expected=pd.read_csv(tmp_path/'all_proteins_scores.tsv',sep='\t').iloc[0]
    hits={m.id:[dict(peptide=p,canonical=p.replace('I','L'),aa_start=a,aa_end=b) for p,a,b in peptides]}
    actual=features(models,hits,{p:dict(quality='unknown') for p,a,b in peptides}).iloc[0]
    from gapms.model_utils import ITERATIVE_FEATURE_COLS
    for column in ITERATIVE_FEATURE_COLS+['splice_sites','supported_splice_sites']:
        assert actual[column]==expected[column],column


def test_minus_strand_coding_is_already_in_translation_order():
    model=SimpleNamespace(id='m',source='a',original_id='t',gene='g',chrom='c',strand='-',
        protein='MABCDEFGHIJK',cds=[(0,24),(40,52)],coding=[(40,52),(0,24)])
    hits={'m':[dict(peptide='ABC',canonical='ABC',aa_start=2,aa_end=5)]}
    row=features([model],hits,{'ABC':dict(quality='unknown')}).iloc[0]
    assert row.splice_sites==1
    assert row.supported_splice_sites==1
    assert row.splice_peptides==1
