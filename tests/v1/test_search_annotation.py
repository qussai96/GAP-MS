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


def test_search_annotation_does_not_suppress_candidate_relative_to_helixer(tmp_path,monkeypatch):
    from test_assessment import make_bam
    sequence=dna('MQQQQQQQK*')+'N'*100+dna('MPEPTIDEKACDEFGHK*')
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+sequence+'\n')
    queries=tmp_path/'p.tsv';queries.write_text('sequence\tq-value\nPEPTIDE\t0.001\nACDEFGH\t0.001\n')
    annotation=tmp_path/'helixer.gtf';annotation.write_text('chr1\tx\tCDS\t1\t30\t.\t+\t0\tParent=helixer\n')
    guide=tmp_path/'refseq.gtf';guide.write_text(f'chr1\tx\tCDS\t131\t{len(sequence)}\t.\t+\t0\tParent=known\n')
    bam=tmp_path/'r.bam';make_bam(bam,sequence)
    def reconstruction(args,reference,lengths,commands):
        (args.output/'orfs').mkdir(exist_ok=True);(args.output/'transcripts').mkdir(exist_ok=True)
        path=args.output/'orfs/candidate_orfs.gtf';path.write_text(guide.read_text());return path
    monkeypatch.setattr('gapms_v1.workflow.reconstruct',reconstruction)
    monkeypatch.setattr('gapms_v1.workflow.tool_versions',lambda bam:{'mock_reconstruction':True})
    out=tmp_path/'out'
    assert main(['-a',str(genome),'-p',str(queries),'-g',str(annotation),'--search-gtf',str(guide),'-b',str(bam),'-o',str(out),'-t','1'])==0
    candidates=table(out/'orfs/novel_orf_candidates.tsv')
    assert len(candidates)==1 and candidates[0]['relationship']=='intergenic'
    assert candidates[0]['confidence']=='high_confidence_candidate'
    assert len(table(out/'annotation/gene_model_evidence.tsv'))==1


def test_mapping_only_orphan_exon_is_reported_without_inventing_transcript(tmp_path):
    genome=tmp_path/'g.fa';genome.write_text('>chr1\n'+dna('MPEPTIDEK*')+'\n')
    pysam.faidx(str(genome))
    gff=tmp_path/'a.gff';gff.write_text('chr1\tx\tCDS\t1\t30\t.\t+\t0\tParent=known\nchr1\tx\texon\t1\t12\t.\t+\t.\tID=orphan\n')
    skipped=[]
    with pysam.FastaFile(str(genome)) as ref:
        models=load_models(gff,'SEARCH',ref,dict(zip(ref.references,ref.lengths)),namespace_loci=True,skipped_features=skipped)
    assert len(models)==1 and models[0].original_id=='known'
    assert len(skipped)==1 and skipped[0]['line']==2
    assert 'orphan_exon' in skipped[0]['reason']
