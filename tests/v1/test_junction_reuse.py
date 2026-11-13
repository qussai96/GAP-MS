import json
import pytest
from Bio.Seq import Seq
from gapms_v1.cli import main
from test_mapping import dna, table

@pytest.mark.parametrize('strand',['+','-'])
def test_precomputed_annotation_junctions_keep_zero_read_edges(tmp_path,monkeypatch,strand):
    coding=dna('MKACDEFGHKA*');cut=13;intron='GT'+'C'*28+'AG'
    seq=coding[:cut]+intron+coding[cut:];blocks=[(0,cut),(cut+len(intron),len(seq))]
    if strand=='-':
        seq=str(Seq(seq).reverse_complement());blocks=sorted((len(seq)-b,len(seq)-a) for a,b in blocks)
    genome=tmp_path/'g';genome.write_text('>chr1\n'+seq+'\n')
    gtf=tmp_path/'a';lines=[];used=0
    for a,b in blocks if strand=='+' else reversed(blocks):
        lines.append(f'chr1\tx\tCDS\t{a+1}\t{b}\t.\t{strand}\t{(-used)%3}\tgene_id "g"; transcript_id "t";\n');used+=b-a
    gtf.write_text(''.join(lines));peptides=tmp_path/'p';peptides.write_text('ACDEFGHK\n');out=tmp_path/'out'
    def duplicate_parse(*a,**kw):raise AssertionError('PepLocus reparsed already extracted GTF junctions')
    monkeypatch.setattr('peplocus.pipeline.read_annotation',duplicate_parse)
    assert main(['-a',str(genome),'-p',str(peptides),'-g',str(gtf),'-o',str(out),'-t','1','--keep-work'])==0
    assert (out/'peptides.tryptic.bed').read_text()
    manifest=json.loads((out/'diagnostics/run_manifest.json').read_text())
    options=next(c['options'] for c in manifest['commands'] if c.get('component')=='peplocus')
    assert 'gtf' not in options and options['junctions'].endswith('search_junctions.tsv')
    mapping=table(out/'diagnostics/mapping/mapped_peptides.tsv')
    assert any(len(json.loads(r['blocks']))==2 for r in mapping)
