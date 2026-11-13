"""Release-boundary checks: minimal CLI and independent tryptic multiplicity."""
import csv
from gapms_v1.cli import parser, main
from gapms_v1.mapping import EvidenceStore


def test_default_help_is_small_but_advanced_controls_remain_available(capsys):
    basic=parser().format_help()
    assert '--assembly' in basic and '--gtf' in basic and '--bam' in basic
    assert '--consensus-beam' not in basic and '--min-orf-aa' not in basic
    assert main(['--help-all'])==0
    assert '--consensus-beam' in capsys.readouterr().out
    assert parser().parse_args(['-a','a','-p','p','-o','o']).threads==8


def test_dual_tracks_preserve_all_locations_and_score_subset_independently(tmp_path):
    p='ACDEFGK';alias='ACDEFGK'
    store=EvidenceStore(tmp_path/'hits.sqlite')
    store.add(p,'c','+',[(0,21)],p,'genome')
    store.add(p,'c','+',[(30,51)],p,'annotation',tryptic=True)
    result=store.export(tmp_path,[p],{p:dict(peptide_id='p',quality='unknown')})
    all_rows=[line.split('\t') for line in (tmp_path/'mapped_peptides.bed').read_text().splitlines()]
    subset=[line.split('\t') for line in (tmp_path/'tryptic_peptides.bed').read_text().splitlines()]
    assert len(all_rows)==2 and len(subset)==1
    assert all(r[4]=='0' for r in all_rows) and subset[0][4]=='1000'
    assert subset[0][1]=='30'
    assert result['multiple_peptides']==1 and result['tryptic_unique_peptides']==1
    store.db.close()
