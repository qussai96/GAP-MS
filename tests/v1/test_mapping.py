import csv
import json
from pathlib import Path

import pytest
from Bio.Seq import Seq

from gapms_v1.cli import main
from gapms_v1.models import load_models, write_models
from gapms._peplocus.io import prepare_reference
import pysam


def table(path):
    with path.open() as stream:return list(csv.DictReader(stream,delimiter='\t'))


CODONS={'M':'ATG','A':'GCT','C':'TGT','D':'GAT','E':'GAA','F':'TTT','G':'GGT','H':'CAT','I':'ATT','K':'AAA','L':'CTG','N':'AAT','P':'CCT','Q':'CAA','R':'CGT','S':'TCT','T':'ACT','V':'GTT','W':'TGG','Y':'TAT','*':'TAA'}
def dna(protein):return ''.join(CODONS[a] for a in protein)


def test_genome_end_to_end_strands_frames_aliases_multiplicity(tmp_path):
    seq='CCC'+dna('MPEPTIDEK*')+'N'*37+str(Seq(dna('MPEPTIDEK*')).reverse_complement())+'N'*31+'A'+dna('MABCDEFG*'.replace('B','N'))
    assembly=tmp_path/'genome.fa';assembly.write_text('>chr1\n'+seq+'\n')
    peptides=tmp_path/'peptides.tsv';peptides.write_text('sequence\nPEPTIDE\nPEPTLDE\nVVVVVVV\nAAA\n')
    out=tmp_path/'out'
    assert main(['-a',str(assembly),'-p',str(peptides),'-o',str(out),'-t','1'])==0
    rows=table(out/'mapping/mapped_peptides.tsv')
    independently_verify(assembly,rows)
    assert {r['strand'] for r in rows}=={'+','-'}
    assert {r['peptide'] for r in rows}=={'PEPTIDE','PEPTLDE'}
    assert all(int(r['loci_for_peptide'])==2 for r in rows)
    assert [r['peptide'] for r in table(out/'mapping/excluded_peptides.tsv')]==['AAA']
    assert [r['peptide'] for r in table(out/'mapping/unmapped_peptides.tsv')]==['VVVVVVV']
    manifest=json.loads((out/'run_manifest.json').read_text());assert manifest['status']=='completed'
    assert main(['-a',str(assembly),'-p',str(peptides),'-o',str(out),'-t','1','--resume'])==0
    assert not Path(str(assembly)+'.fai').exists()
    peptides.write_text('sequence\nPEPTIDE\n')
    assert main(['-a',str(assembly),'-p',str(peptides),'-o',str(out),'-t','1','--resume'])==1


@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('cut',[10,11,12])
def test_annotation_split_codons_and_negative_strand(tmp_path,strand,cut):
    coding=dna('MPEPTIDEK*');intron='GT'+'C'*27+'AG'
    seq=coding[:cut]+intron+coding[cut:]
    blocks=[(0,cut),(cut+len(intron),len(seq))]
    if strand=='-':
        seq=str(Seq(seq).reverse_complement());blocks=sorted((len(seq)-b,len(seq)-a) for a,b in blocks)
    assembly=tmp_path/'genome.fa';assembly.write_text('>chr1\n'+seq+'\n')
    gtf=tmp_path/'model.gtf'
    ordered=blocks if strand=='+' else list(reversed(blocks));cumulative=0;lines=[]
    for a,b in ordered:
        phase=(-cumulative)%3
        lines.append(f'chr1\ttest\tCDS\t{a+1}\t{b}\t.\t{strand}\t{phase}\tgene_id "g"; transcript_id "t";\n');cumulative+=b-a
    gtf.write_text(''.join(lines))
    peptides=tmp_path/'peps.txt';peptides.write_text('PEPTIDE\n')
    output=tmp_path/'out'
    assert main(['-a',str(assembly),'-p',str(peptides),'-g',str(gtf),'-o',str(output),'-t','1'])==0
    rows=table(output/'mapping/mapped_peptides.tsv')
    independently_verify(assembly,rows)
    assert any(len(json.loads(r['blocks']))==2 and r['strand']==strand for r in rows)
    model=table(output/'annotation/gene_model_evidence.tsv')[0]
    assert model['category']=='peptide_supported' and model['distinct_peptides']=='1'
    assert Path(str(assembly)+'.fai').exists()==False


def test_multiple_annotation_sources_conflicts_and_no_support(tmp_path):
    assembly=tmp_path/'g.fa';assembly.write_text('>chr1\n'+dna('MPEPTIDEK*')+'N'*30+dna('MQQQQQQQK*')+'\n')
    p=tmp_path/'p.txt';p.write_text('PEPTIDE\n')
    a=tmp_path/'a.gtf';b=tmp_path/'b.gff3'
    a.write_text('chr1\tx\tCDS\t1\t30\t.\t+\t0\ttranscript_id "same"; gene_id "a";\n')
    b.write_text('chr1\tx\tCDS\t61\t90\t.\t+\t0\tParent=same;gene_id=b\n')
    out=tmp_path/'out';assert main(['-a',str(assembly),'-p',str(p),'-g',str(a),'-g',str(b),'-o',str(out),'-t','1'])==0
    rows=table(out/'annotation/gene_model_evidence.tsv')
    assert len({r['model_id'] for r in rows})==2
    assert {r['category'] for r in rows}=={'peptide_supported','no_observed_peptide_support'}
    assert len(table(out/'annotation/no_peptide_support_gene_models.tsv'))==1


def test_bad_inputs_fail_without_touching_inputs(tmp_path):
    p=tmp_path/'peps';p.write_text('PEPTIDE\n')
    assert main(['-a',str(tmp_path/'missing'),'-p',str(p),'-o',str(tmp_path/'out')])==1
    g=tmp_path/'g.fa';g.write_text('>chr1\nATG\n')
    ann=tmp_path/'a.gtf';ann.write_text('chr1\tx\tCDS\t1\t30\t.\t+\t0\tParent=t\n')
    assert main(['-a',str(g),'-p',str(p),'-g',str(ann),'-o',str(tmp_path/'out')])==1
    assert not Path(str(g)+'.fai').exists()


def independently_verify(assembly,rows,code=1):
    from Bio import SeqIO
    sequences={r.id:str(r.seq) for r in SeqIO.parse(assembly,'fasta')}
    for row in rows:
        blocks=json.loads(row['blocks'])
        coding=''.join(sequences[row['chrom']][a:b] for a,b in blocks)
        if row['strand']=='-':coding=str(Seq(coding).reverse_complement())
        translated=str(Seq(coding).translate(table=code))
        assert translated==row['matched_sequence']
        assert translated.replace('I','L')==row['peptide'].replace('I','L')


@pytest.mark.parametrize('strand',['+','-'])
def test_microexon_and_full_cds_projection(tmp_path,strand):
    coding=dna('MPEPTIDEK*');parts=[coding[:10],coding[10:12],coding[12:]]
    intron='GT'+'C'*28+'AG';sequence=intron.join(parts)
    blocks=[];position=0
    for part in parts:
        blocks.append((position,position+len(part)));position+=len(part)+len(intron)
    if strand=='-':
        sequence=str(Seq(sequence).reverse_complement());blocks=sorted((len(sequence)-b,len(sequence)-a) for a,b in blocks)
    assembly=tmp_path/'g.fa';assembly.write_text('>chr1\n'+sequence+'\n')
    annotation=tmp_path/'a.gff3';lines=[];cumulative=0
    for a,b in blocks if strand=='+' else reversed(blocks):
        lines.append(f'chr1\tx\tCDS\t{a+1}\t{b}\t.\t{strand}\t{(-cumulative)%3}\tParent=t\n');cumulative+=b-a
    annotation.write_text(''.join(lines));peptides=tmp_path/'p';peptides.write_text('PEPTIDE\n')
    out=tmp_path/'out'
    assert main(['-a',str(assembly),'-p',str(peptides),'-g',str(annotation),'-o',str(out),'-t','1'])==0
    rows=table(out/'mapping/mapped_peptides.tsv');assert any(len(json.loads(r['blocks']))==3 for r in rows)
    independently_verify(assembly,rows)


def test_nonstandard_genetic_code_and_multiple_frames(tmp_path):
    assembly=tmp_path/'g.fa';sequence='C'+'ATGCCGGAACCTACTTGGGATGAATAA'+'N'*21+'CC'+'ATGCCGGAACCTACTTGGGATGAATAA'
    assembly.write_text('>chr1\n'+sequence.replace('TGG','TGA')+'\n')
    peptide=tmp_path/'p';peptide.write_text('PEPTWDE\n')
    out=tmp_path/'out'
    assert main(['-a',str(assembly),'-p',str(peptide),'-o',str(out),'-t','1','--genetic-code','4'])==0
    rows=table(out/'mapping/mapped_peptides.tsv');assert len(rows)>=2
    independently_verify(assembly,rows,4)


def test_identification_metadata_and_bad_quality(tmp_path):
    assembly=tmp_path/'g.fa';assembly.write_text('>chr1\n'+dna('MPEPTIDEK*')+'\n')
    p=tmp_path/'p.tsv';p.write_text('sequence\tq-value\tPEP\tspectra\nPEPTIDE\t0.1\t0.2\t3\n')
    out=tmp_path/'out';assert main(['-a',str(assembly),'-p',str(p),'-o',str(out),'-t','1'])==0
    row=table(out/'mapping/mapped_peptides.tsv')[0];assert row['identification_quality']=='qvalue_fail'
    original=json.loads(table(out/'mapping/peptide_identifications.tsv')[0]['original_record']);assert original['PEP']=='0.2' and original['spectra']=='3'


def test_unstranded_rna_transcripts_are_extracted_plus_without_inventing_strand(tmp_path):
    assembly=tmp_path/'g.fa';assembly.write_text('>chr1\n'+dna('MPEPTIDEK*')+'\n')
    annotation=tmp_path/'rna.gtf';annotation.write_text('chr1\tStringTie\texon\t1\t30\t.\t.\t.\tgene_id "g"; transcript_id "t";\n')
    genome,index,lengths=prepare_reference(assembly,tmp_path/'cache')
    with pysam.FastaFile(genome,filepath_index=index) as reference:
        with pytest.raises(ValueError,match='strand'):
            load_models(annotation,'annotation',reference,lengths)
        models=load_models(annotation,'RNA_transcripts',reference,lengths,allow_unstranded=True)
        assert models[0].strand=='.' and models[0].protein==''
