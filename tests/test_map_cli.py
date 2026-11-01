"""End-to-end tests of the bundled mapper and both BED conventions."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pysam
import pytest
from Bio.Seq import Seq


def make_inputs(tmp_path):
    codons = {'M':'ATG','P':'CCT','E':'GAA','T':'ACT','I':'ATT','D':'GAT','K':'AAA'}
    coding = ''.join(codons[a] for a in 'MPEPTIDEK')
    intron = 'GT' + 'C'*24 + 'AG'
    forward = 'N'*30 + coding[:13] + intron + coding[13:] + 'N'*30
    dna = forward + 'N'*30 + str(Seq(forward).reverse_complement())
    fasta = tmp_path/'assembly.fa'; fasta.write_text('>chr1\n'+dna+'\n')
    peptides = tmp_path/'peptides.tsv'
    peptides.write_text('Peptide Sequence\nMPEPTIDEK\nMPEPTIDER\n')
    gtf = tmp_path/'models.gff3'
    intervals = [(30,43,'+'),(71,85,'+')]
    intervals += [(len(dna)-b,len(dna)-a,'-') for a,b,_ in intervals.copy()]
    gtf.write_text(''.join(f'chr1\ttest\texon\t{a+1}\t{b}\t.\t{s}\t.\tParent=t{s}\n' for a,b,s in intervals))
    bam = tmp_path/'reads.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':len(dna)}]}) as out:
        for i in range(3):
            read=pysam.AlignedSegment(); read.query_name=f'r{i}';read.query_sequence=coding
            read.reference_id=0;read.reference_start=30;read.mapping_quality=60
            read.cigartuples=[(0,13),(3,28),(0,14)];read.set_tag('XS','+');out.write(read)
    return fasta,peptides,gtf,bam


@pytest.mark.parametrize('evidence',['gtf','bam','both'])
def test_map_cli_all_evidence_and_beds(tmp_path,evidence):
    genome,peptides,gtf,bam=make_inputs(tmp_path)
    output=tmp_path/'out'
    command=[sys.executable,'-m','gapms.cli','map','-a',str(genome),'-p',str(peptides),
             '-o',str(output),'--rescue-one-edit','--rescue-min-length','9']
    if evidence in {'gtf','both'}:command+=['-g',str(gtf)]
    if evidence in {'bam','both'}:command+=['-b',str(bam)]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    result=subprocess.run(command,capture_output=True,text=True,env=env)
    assert result.returncode==0,result.stdout+result.stderr
    assert (output/'peptides_mapped.bed').read_bytes()==(output/'peplocus/mappings.bed').read_bytes()
    assert (output/'peptides_mapped.igv.bed').read_bytes()==(output/'peplocus/mappings.igv_display.bed').read_bytes()
    rows=list(csv.DictReader((output/'peplocus/mappings.tsv').open(),delimiter='\t'))
    assert any(r['peptide']=='MPEPTIDEK' and r['strand']=='+' and len(json.loads(r['blocks']))==2 for r in rows)
    rescue=list(csv.DictReader((output/'peplocus/rescue_mappings.tsv').open(),delimiter='\t'))
    assert any(r['peptide']=='MPEPTIDER' and r['edit_type']=='substitution' for r in rescue)
    if evidence!='bam':
        assert any(r['strand']=='-' for r in rows)
        display=[l.split('\t') for l in (output/'peptides_mapped.igv.bed').read_text().splitlines() if not l.startswith('track')]
        assert any(r[3]=='KEDITPEPM' and r[5]=='-' for r in display)
    repeat=subprocess.run(command,capture_output=True,text=True,env=env)
    assert repeat.returncode!=0 and 'output already exists' in repeat.stderr


def test_workflow_adapter_uses_assembly_not_protein_projection(tmp_path):
    from gapms.peptides_to_genome import map_peptides_to_genome
    genome,peptides,gtf,bam=make_inputs(tmp_path)
    output=tmp_path/'out'
    # No protein or protein-mapping file is required for assembly mapping.
    map_peptides_to_genome(gtf,None,None,output,assembly=genome,peptides=peptides)
    assert (output/'peptides_mapped.bed').exists()
    assert (output/'peptides_mapped.igv.bed').exists()
    assert not (output/'protein_projection').exists()


def test_map_cli_does_not_import_classifier():
    code = "import sys; from gapms.map_cli import main\ntry: main(['--help'])\nexcept SystemExit as e: assert e.code == 0\nassert 'gapms.filter' not in sys.modules\nassert 'xgboost' not in sys.modules"
    command=[sys.executable,'-c',code]
    result=subprocess.run(command,capture_output=True,text=True)
    assert result.returncode==0 and '--assembly' in result.stdout


def test_reference_report_keeps_genomic_ambiguity_and_model_metadata(tmp_path):
    from gapms.compare_supp_ref import load_genomic_model_evidence, get_terminal_supporting_peptides
    genomic=tmp_path/'genomic.bed';projection=tmp_path/'model.bed'
    genomic.write_text('chr1\t0\t12\tMILI\t0\t+\t0\t12\t0\t1\t12\t0\n'
                       'chr1\t30\t42\tMILI\t0\t+\t30\t42\t0\t1\t12\t0\n')
    projection.write_text('chr1\t0\t12\tMLLL\t0\t+\t0\t12\t0\t1\t12\t0\tprotein=p1\tpep_start=1\tpep_end=4\tprot_len=10\n'
                          'chr1\t60\t72\tMLLL\t0\t+\t60\t72\t0\t1\t12\t0\tprotein=p2\n')
    evidence=load_genomic_model_evidence(genomic,projection)
    assert len(evidence['MLLL'])==2
    assert evidence['MLLL'][0]['protein']=='p1'
    assert evidence['MLLL'][1]['protein'] is None
    pred={'seqid':'chr1','strand':'+','cds_regions':[(1,30)]}
    assert get_terminal_supporting_peptides('p1',pred,evidence,'start')==set()
    genomic.write_text(genomic.read_text().splitlines()[0]+'\n')
    evidence=load_genomic_model_evidence(genomic,projection)
    assert get_terminal_supporting_peptides('p1',pred,evidence,'start')=={'MLLL'}


def test_model_compatibility_file_stays_separate(tmp_path):
    from gapms.peptides_to_genome import map_peptides_to_genome
    from gapms.compare_supp_ref import load_genomic_model_evidence
    genome,peptides,gtf,bam=make_inputs(tmp_path)
    gtf.write_text(gtf.read_text().replace('\texon\t','\tCDS\t'))
    proteins=tmp_path/'proteins.fa';proteins.write_text('>t+\nMPEPTIDEK\n')
    mapping=tmp_path/'proteomapper.tsv'
    mapping.write_text('Peptide\tProtein\tlocation\nMPEPTIDEK\tt+\t1\n')
    output=tmp_path/'out'
    map_peptides_to_genome(gtf,proteins,mapping,output,assembly=genome,
                           peptides=peptides,model_evidence=True)
    projected=output/'protein_projection/peptides_mapped.bed'
    assert projected.exists()
    assert all(len(line.split('\t'))==12 for line in (output/'peptides_mapped.bed').read_text().splitlines())
    evidence=load_genomic_model_evidence(output/'peptides_mapped.bed',projected)
    assert any(r['protein']=='t+' for r in evidence['MPEPTLDEK'])
    assert any(r['protein'] is None for r in evidence['MPEPTLDEK'])
