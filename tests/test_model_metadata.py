"""Equivalence to the original projection and model-aware report evidence."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gapms.model_metadata import CodingRuns, attach_model_metadata
from gapms.protein_projection import map_codon_to_genome, map_peptides_to_genome
from gapms.compare_supp_ref import load_genomic_model_evidence
from tests.helpers import create_mini_gapms_dataset


@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('phase',[0,1,2])
def test_cached_runs_match_nucleotide_oracle(strand,phase):
    rng=np.random.default_rng(24)
    for case in range(20):
        rows=[];start=1
        for i in range(int(rng.integers(1,10))):
            end=start+int(rng.integers(1,30))
            rows.append(dict(Protein='p',Seqid='chr1',Start=start,End=end,Strand=strand,Frame=str(phase if (strand=='+' and i==0) else 0)))
            start=end+int(rng.integers(1,15))
        cds=pd.DataFrame(rows);index=CodingRuns(cds)
        for _ in range(20):
            a=int(rng.integers(1,max(2,index.length//3+2)))
            b=a+int(rng.integers(0,6))
            assert index.project(a,b)==map_codon_to_genome(a,b,cds,strand)


@pytest.mark.parametrize('strand',['+','-'])
def test_pathological_overlap_and_adjacent_segments(strand):
    cds=pd.DataFrame([dict(Protein='p',Seqid='c',Start=a,End=b,Strand=strand,Frame='0')
                      for a,b in [(1,9),(7,12),(13,21),(13,21)]])
    index=CodingRuns(cds)
    for a in range(1,12):
        for b in range(a,13):assert index.project(a,b)==map_codon_to_genome(a,b,cds,strand)


def test_complete_metadata_bed_and_evidence_are_identical(tmp_path):
    inputs=create_mini_gapms_dataset(tmp_path)
    old=tmp_path/'old';new=tmp_path/'new'
    map_peptides_to_genome(inputs['prediction_gtf'],inputs['protein_fasta'],inputs['mapping_file'],old)
    attach_model_metadata(inputs['prediction_gtf'],inputs['protein_fasta'],inputs['mapping_file'],new)
    assert (old/'peptides_mapped.bed').read_bytes()==(new/'peptides_mapped.bed').read_bytes()
    assert load_genomic_model_evidence(inputs['peptides_bed'],old/'peptides_mapped.bed')==load_genomic_model_evidence(inputs['peptides_bed'],new/'peptides_mapped.bed')


def test_reference_reports_match_with_old_and_fast_metadata(tmp_path):
    from gapms.compare_supp_ref import generate_annotation_report
    inputs=create_mini_gapms_dataset(tmp_path/'inputs')
    results=[]
    for label,method in [('old',map_peptides_to_genome),('fast',attach_model_metadata)]:
        output=tmp_path/label;output.mkdir()
        # Standard genomic coordinates are held fixed in both report runs.
        (output/'peptides_mapped.bed').write_bytes(inputs['peptides_bed'].read_bytes())
        method(inputs['prediction_gtf'],inputs['protein_fasta'],inputs['mapping_file'],output/'protein_projection')
        results.append(generate_annotation_report(
            output,inputs['prediction_gtf'],inputs['reference_gtf'],
            tmap_file=inputs['tmap_file'],all_proteins_scores=inputs['all_scores_tsv'],
            protein_fasta=inputs['protein_fasta']))
    for old,fast in zip(*results):pd.testing.assert_frame_equal(old,fast)
    for path in (tmp_path/'old/Compare_to_Reference').rglob('*'):
        if path.is_file():
            other=tmp_path/'fast/Compare_to_Reference'/path.relative_to(tmp_path/'old/Compare_to_Reference')
            assert path.read_bytes()==other.read_bytes()
