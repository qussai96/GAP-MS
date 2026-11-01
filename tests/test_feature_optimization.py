"""Independent equivalence checks for vectorized feature calculations."""
from collections import Counter
import numpy as np
import pandas as pd
import pytest

from gapms.protein_projection import annotate_cds_coordinates, _annotate_cds_coordinates_reference
from gapms.extract import build_splice_junctions


def reference_junctions(cds):
    records=[];counts={}
    for protein,group in cds.groupby('Protein',sort=False):
        rows=list(group.sort_values('cds_nt_start').itertuples());count=0
        for a,b in zip(rows,rows[1:]):
            if a.Seqid!=b.Seqid:continue
            has_intron=(int(b.Start)>int(a.End)+1 if a.Strand=='+' else
                        int(a.Start)>int(b.End)+1 if a.Strand=='-' else False)
            if has_intron:records.append({'Protein':protein,'junction_nt':int(a.cds_nt_end)});count+=1
        counts[protein]=count
    return pd.DataFrame(records,columns=['Protein','junction_nt']),counts


@pytest.mark.parametrize('seed',range(8))
def test_cds_and_junctions_match_original(seed):
    rng=np.random.default_rng(seed);rows=[]
    for protein in range(40):
        start=1
        for exon in range(int(rng.integers(1,15))):
            length=int(rng.integers(1,90));end=start+length-1
            rows.append(dict(Protein=f'p{protein}',Seqid='chr'+str(protein%3),
                             Start=start,End=end,Strand='-' if protein%2 else '+',
                             Frame=rng.choice(['0','1','2','.','bad','1.0']),Gene=f'g{protein}'))
            start=end+int(rng.integers(1,50))
    df=pd.DataFrame(rows).sample(frac=1,random_state=seed).reset_index(drop=True)
    expected=_annotate_cds_coordinates_reference(df)
    actual=annotate_cds_coordinates(df)
    pd.testing.assert_frame_equal(actual,expected)
    observed,counts=build_splice_junctions(actual)
    wanted,wanted_counts=reference_junctions(expected)
    pd.testing.assert_frame_equal(observed,wanted,check_dtype=False)
    assert {k:counts.get(k,0) for k in wanted_counts}==wanted_counts


def test_duplicate_start_and_index_preserve_reference_sort():
    df=pd.DataFrame([dict(Protein='p',Seqid='c',Start=1,End=end,Strand='-',Frame=phase)
                     for end,phase in [(1,'2'),(5,'0'),(9,'1')]],index=[0,0,1])
    pd.testing.assert_frame_equal(annotate_cds_coordinates(df),_annotate_cds_coordinates_reference(df))


@pytest.mark.parametrize('columns', [True,False])
def test_missing_phase_and_empty_cds(columns):
    df=pd.DataFrame([dict(Protein='p',Seqid='c',Start=1,End=5,Strand='+',Frame='0')])
    if not columns:df=df.drop(columns='Frame')
    pd.testing.assert_frame_equal(annotate_cds_coordinates(df),_annotate_cds_coordinates_reference(df))
    pd.testing.assert_frame_equal(annotate_cds_coordinates(df.iloc[:0]),_annotate_cds_coordinates_reference(df.iloc[:0]))


def test_null_protein_and_no_junctions():
    df=pd.DataFrame([dict(Protein=None,Seqid='c',Start=1,End=1,Strand='+',Frame='2'),
                     dict(Protein=None,Seqid='c',Start=2,End=6,Strand='+',Frame='0')])
    actual=annotate_cds_coordinates(df)
    pd.testing.assert_frame_equal(actual,_annotate_cds_coordinates_reference(df))
    junctions,counts=build_splice_junctions(actual)
    assert junctions.empty and counts=={}
