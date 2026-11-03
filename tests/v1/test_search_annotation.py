from pathlib import Path
import pysam
from gapms_v1.cli import main
from gapms_v1.models import load_models
from test_mapping import dna, table


def test_search_annotation_adds_spliced_locus_without_becoming_target(tmp_path):
    peptide='PEPTIDE'
    coding=dna(peptide)
    sequence=coding+'N'*60+coding[:9]+'N'*60+coding[9:]
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+sequence+'\n')
    queries=tmp_path/'p.txt';queries.write_text(peptide+'\n')
    guide=tmp_path/'s.gff';guide.write_text('chr1\tx\tCDS\t82\t90\t.\t+\t0\tParent=known\nchr1\tx\tCDS\t151\t162\t.\t+\t0\tParent=known\n')
    out=tmp_path/'out'
    assert main(['-a',str(genome),'-p',str(queries),'--search-gtf',str(guide),'-o',str(out),'-t','1'])==0
    rows=table(out/'mapping/mapped_peptides.tsv')
    assert len(rows)==2 and {r['loci_for_peptide'] for r in rows}=={'2'}
    assert not (out/'annotation').exists()
    import json
    assert json.loads((out/'run_manifest.json').read_text())['mapping_only_models']==1


def test_mapping_only_reused_transcript_ids_are_locus_namespaced(tmp_path):
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+dna('MPEPTIDEK*')+'\n>chr2\n'+dna('MPEPTIDEK*')+'\n')
    pysam.faidx(str(genome))
    gff=tmp_path/'a.gff';gff.write_text('chr1\tx\tCDS\t1\t30\t.\t+\t0\tParent=reused\nchr2\tx\tCDS\t1\t30\t.\t+\t0\tParent=reused\n')
    with pysam.FastaFile(str(genome)) as ref:
        models=load_models(gff,'SEARCH',ref,dict(zip(ref.references,ref.lengths)),namespace_loci=True)
    assert len(models)==2 and len({m.id for m in models})==2
