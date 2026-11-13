import os
import sys
from gapms_v1.common import run_command


def test_child_environment_is_recorded_without_leaking_secrets(tmp_path):
    commands=[];env=dict(os.environ,PERL_HASH_SEED='0',PERL_PERTURB_KEYS='0',PRIVATE_TEST_VALUE='not-for-provenance')
    log=tmp_path/'run.log'
    run_command([sys.executable,'-c','import os; print(os.environ["PERL_HASH_SEED"], os.environ["PERL_PERTURB_KEYS"])'],tmp_path,log,commands,env=env)
    assert log.read_text().strip()=='0 0'
    assert commands[0]['environment']==dict(PERL_HASH_SEED='0',PERL_PERTURB_KEYS='0')
    assert 'PRIVATE_TEST_VALUE' not in str(commands)
    assert os.environ.get('PRIVATE_TEST_VALUE')!='not-for-provenance'


def test_transcript_ids_and_order_ignore_stringtie_scheduling(tmp_path):
    from gapms_v1.models import Model
    from gapms_v1.rnaseq import canonicalize_transcripts
    import csv
    first=[Model('x','RNA_transcripts','STRG.1.1','STRG.1','c','+',exons=[(0,9),(30,60)]),
           Model('y','RNA_transcripts','STRG.1.2','STRG.1','c','+',exons=[(0,12),(30,60)]),
           Model('z','RNA_transcripts','STRG.2.1','STRG.2','c','.',exons=[(90,120)])]
    second=[Model('b','RNA_transcripts','STRG.77.1','STRG.77','c','.',exons=[(90,120)]),
            Model('a','RNA_transcripts','STRG.88.2','STRG.88','c','+',exons=[(0,12),(30,60)]),
            Model('d','RNA_transcripts','STRG.88.3','STRG.88','c','+',exons=[(0,9),(30,60)]),
            Model('c','RNA_transcripts','STRG.88.1','STRG.88','c','+',exons=[(0,9),(30,60)])]
    outputs=[]
    for i,models in enumerate([first,second]):
        out=tmp_path/f'{i}.gtf';mapping=tmp_path/f'{i}.tsv'
        result=canonicalize_transcripts(models,out,mapping)
        assert len(result)==3 and len({m.gene for m in result})==2
        assert any(m.strand=='.' for m in result)
        with mapping.open() as f:records=list(csv.DictReader(f,delimiter='\t'))
        assert len(records)==len(models)
        assert all(r['original_transcript_id'].startswith('STRG.') and r['transcript_id'].startswith('RT') for r in records)
        outputs.append(out.read_bytes())
    assert outputs[0]==outputs[1]
