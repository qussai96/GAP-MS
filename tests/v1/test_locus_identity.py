import pysam
from gapms_v1.models import load_models


def test_reused_transcript_names_keep_locus_and_gene_parent_separate(tmp_path):
    genome=tmp_path/'g.fa';genome.write_text('>a\nATGAAATAA\n>b\nATGAAATAA\n');pysam.faidx(str(genome))
    gff=tmp_path/'g.gff';gff.write_text('a\tx\tmRNA\t1\t9\t.\t+\t.\tID=reused;Parent=geneA\na\tx\tCDS\t1\t9\t.\t+\t0\tParent=reused\nb\tx\tmRNA\t1\t9\t.\t+\t.\tID=reused;Parent=geneB\nb\tx\tCDS\t1\t9\t.\t+\t0\tParent=reused\n')
    with pysam.FastaFile(str(genome)) as ref:models=load_models(gff,'annotation',ref,dict(zip(ref.references,ref.lengths)))
    assert len(models)==2 and len({m.id for m in models})==2
    assert {(m.chrom,m.gene,m.original_id) for m in models}=={('a','geneA','reused'),('b','geneB','reused')}
