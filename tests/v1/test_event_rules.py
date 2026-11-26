"""Scientific event rules: locus identity, diagnostic evidence, starts, splicing, ORFs."""
from collections import defaultdict
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from gapms_v1.assessment import assess
from gapms_v1.cli import main as cli_main
from gapms_v1.mapping import EvidenceStore
from gapms_v1.models import Model, flag_noncanonical_introns
from gapms_v1.orf_projection import convert_transdecoder_to_genome_gtf
import test_refinement
from test_assessment import MockEvidence, config, hit
from test_mapping import dna, table
from test_refinement import run_case


def tryptic(monkeypatch):
    """run_case adds --non-tryptic by default; use the default tryptic policy."""
    monkeypatch.setattr(test_refinement,'main',lambda argv:cli_main([*argv,'--keep-work']))


def starts(rows):return [r for r in rows if r['event_type']=='alternative_start']


# ---------------------------------------------------------------- locus identity

def test_same_locus_paths_are_one_locus_distinct_sites_are_not(tmp_path):
    store=EvidenceStore(tmp_path/'e.sqlite')
    try:
        store.add('PEPTIDEK','chr1','+',[(0,24)],'PEPTIDEK','genome')
        # Last codon split by a junction 1 nt into the next exon: same origin.
        store.add('PEPTIDEK','chr1','+',[(0,23),(100,101)],'PEPTIDEK','annotation')
        assert store.counts()=={'PEPTIDEK':1} and store.placement_counts()=={'PEPTIDEK':2}
        store.add('PEPTIDEK','chr1','+',[(500,524)],'PEPTIDEK','genome')
        store.add('PEPTIDEK','chr1','-',[(0,24)],'PEPTIDEK','genome')
        assert store.counts()=={'PEPTIDEK':3}
    finally:store.db.close()


def test_identical_translation_junction_twin_keeps_annotated_peptide_unique(tmp_path):
    # Exon 2 starts with the same 6 nt as the intron: a donor/acceptor shift by
    # +6 translates identically. The annotated junction explains the peptide.
    exon1=dna('MSDEWHKAFNDYW');intron='GTGCAG'+'GT'+'C'*60+'CCCAG';exon2='GTGCAG'+dna('HEYFWKTNPDHWR*')
    sequence='C'*90+exon1+intron+exon2+'C'*90
    a=90;b=a+len(exon1);c=b+len(intron);d=c+len(exon2)
    out,rows=run_case(tmp_path,sequence,[[(a,b),(c,d)]],['AFNDYWVQHEYFWK'])
    assert not [r for r in rows if r['event_type']=='alternative_splice']
    model=table(out/'annotation/gene_model_evidence.tsv')[0]
    assert model['category']=='peptide_supported' and model['locus_discriminating_peptides']=='1'
    mapped=table(out/'mapping/mapped_peptides.tsv')
    assert {r['loci_for_peptide'] for r in mapped}=={'1'} and len(mapped)==2


# ---------------------------------------------------------------- alternative starts

@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('upstream',[False,True])
def test_internal_start_blocked_by_unique_upstream_peptide(tmp_path,strand,upstream):
    protein='MQQQQQQQQMAPEPTIDEKACDEFGHK*';sequence='C'*90+dna(protein)+'C'*90
    peptides=['MAPEPTIDEK']+(['MQQQQQQQQ'] if upstream else [])
    out,rows=run_case(tmp_path,sequence,[[(90,90+3*len(protein))]],peptides,strand)
    row,=starts(rows)
    if upstream:
        assert row['confidence']=='not_supported'
        assert 'conflicting_upstream_evidence' in json.loads(row['confidence_reasons'])
    else:
        assert row['confidence']=='peptide_supported'


def test_second_codon_start_is_initiator_methionine_ambiguity(tmp_path):
    protein='MMAPEPTIDEKACDEFGHK*';sequence='C'*90+dna(protein)+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,90+3*len(protein))]],['MAPEPTIDEK'])
    row,=starts(rows)
    assert row['boundary_status']=='initiator_methionine_ambiguity' and row['confidence']=='not_supported'


@pytest.mark.parametrize('strand',['+','-'])
def test_internal_start_with_initiator_methionine_removed(tmp_path,strand):
    protein='MQQQQQQQQMASDEFWYNKACDEFGHK*';sequence='C'*90+dna(protein)+'C'*90
    out,rows=run_case(tmp_path,sequence,[[(90,90+3*len(protein))]],['ASDEFWYNK'],strand)
    row,=starts(rows)
    assert row['boundary_status']=='internal_ATG_initiator_methionine_removed'
    assert row['confidence']=='peptide_supported'
    codon=json.loads(row['proposed_feature_blocks'])
    expected=[[117,120]] if strand=='+' else [[len(sequence)-120,len(sequence)-117]]
    assert codon==expected


def test_extension_peptide_directly_after_stop_without_methionine_is_rejected(tmp_path):
    upstream=dna('*VEPTIDEAK');model=dna('MQQQQQQQK*');start=90+len(upstream)
    out,rows=run_case(tmp_path,'C'*90+upstream+model+'C'*90,[[(start,start+len(model))]],['VEPTIDEAK'])
    row,=[r for r in rows if r['event_type']=='n_terminal_extension']
    assert row['boundary_status']=='upstream_stop_without_start' and row['confidence']=='not_supported'


def test_extension_start_blocked_by_peptide_further_upstream(tmp_path):
    # Q before M: not a trypsin site, so M is a start candidate.
    upstream=dna('WEPTIDEFKQMGSDEHWYNK');model=dna('MTRWQDKVNHGYFPDK*');start=90+len(upstream)
    out,rows=run_case(tmp_path,'C'*90+upstream+model+'C'*90,[[(start,start+len(model))]],['MGSDEHWYNK','WEPTIDEFK'])
    row,=starts(rows)
    assert row['boundary_status']=='ATG_supported_start_candidate'
    assert 'conflicting_upstream_evidence' in json.loads(row['confidence_reasons'])


# ---------------------------------------------------------------- strict uniqueness and decoys

@pytest.mark.parametrize('second_copy',[False,True])
def test_nontryptic_copy_elsewhere_makes_event_ambiguous(tmp_path,monkeypatch,second_copy):
    tryptic(monkeypatch)
    protein='MQQQQQQQKACDEFGHWK*';sequence='C'*90+dna(protein)+'C'*90
    if second_copy:sequence+=dna('QACDEFGHWKQ')+'C'*90  # preceded by Q: not trypsin-compatible
    out,rows=run_case(tmp_path,sequence,[[(90,117)]],['ACDEFGHWK'])
    row,=[r for r in rows if r['event_type']=='alternative_stop']
    assert row['confidence']==('ambiguous' if second_copy else 'peptide_supported')
    loci=json.loads(row['peptide_loci'])['ACDEFGHWK']
    assert loci['loci_all_exact']==(2 if second_copy else 1)
    search=json.loads((out/'reports/summary.json').read_text())['annotation_guided_search']
    assert set(search['decoy'])=={'alternative_start','alternative_stop','alternative_splice','n_terminal_extension','c_terminal_extension'}
    assert search['decoy']['alternative_stop']['supported_target_features']==(0 if second_copy else 1)
    assert 'm_initial_enrichment' in search['alternative_start_background']
    bins=search['decoy']['alternative_stop']['by_longest_peptide_length']
    assert set(bins)=={'<=7','8','9-10','11-14','15+'}
    assert bins['9-10']['supported_target_features']==(0 if second_copy else 1)
    if not second_copy:assert row['decoy_estimated_fdr']=='0.0'


# ---------------------------------------------------------------- RNA-derived features

def candidate(mid,blocks,**kw):
    return Model(mid,'RNA_ORF',mid,'rna','chr1','+',cds=blocks,coding=blocks,protein='M'+'Q'*39,complete_start=True,complete_stop=True,**kw)


def reference(mid,blocks):
    return Model(mid,'A01',mid,'g','chr1','+',cds=blocks,coding=blocks,protein='M'+'Q'*39)


def splice_rows(tmp_path):return table(tmp_path/'refinement/alternative_splice_candidates.tsv')


def test_peptide_in_retained_intron_is_reported(tmp_path):
    args=config(tmp_path)
    ref=reference('r',[(0,30),(90,120)]);cand=candidate('c',[(0,120)])
    assess(args,[ref,cand],MockEvidence({'c':[hit('PEPTIDE',[(36,57)],12)]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    row,=splice_rows(tmp_path)
    assert row['event_subtype']=='retained_intron' and row['confidence']=='peptide_supported'


def test_peptide_in_novel_exon_is_reported(tmp_path):
    args=config(tmp_path)
    ref=reference('r',[(0,30),(150,180)]);cand=candidate('c',[(0,30),(60,90),(150,180)])
    assess(args,[ref,cand],MockEvidence({'c':[hit('PEPTIDE',[(63,84)],11)]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    row,=splice_rows(tmp_path)
    assert row['event_subtype']=='novel_exon'


def test_peptide_in_shifted_reading_frame_is_reported(tmp_path):
    args=config(tmp_path)
    ref=reference('r',[(0,100)]);cand=candidate('c',[(0,30),(40,100)])
    assess(args,[ref,cand],MockEvidence({'c':[hit('PEPTIDE',[(49,70)],13)]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    row,=splice_rows(tmp_path)
    assert row['event_subtype']=='frame_change'


def test_isoforms_sharing_a_new_junction_support_one_event(tmp_path):
    args=config(tmp_path)
    ref=reference('r',[(0,30),(90,120)])
    c1=candidate('c1',[(0,33),(90,120)]);c2=candidate('c2',[(0,33),(90,150)])
    crossing=hit('PEPTIDE',[(21,33),(90,99)],7)
    assess(args,[ref,c1,c2],MockEvidence({'c1':[crossing],'c2':[crossing]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{('chr1',33,90,'+'):5})
    row,=splice_rows(tmp_path)
    assert set(json.loads(row['supporting_model_ids']))=={'c1','c2'}
    assert row['confidence']=='rna_peptide_supported'


def test_same_locus_annotated_placement_is_not_diagnostic(tmp_path):
    args=config(tmp_path)
    ref=reference('r',[(0,30),(90,120)]);cand=candidate('c',[(0,20),(90,120)])
    # Annotation explains PEPTIDE contiguously; the RNA path splits its last codon.
    hits={'r':[hit('PEPTIDE',[(0,21)])],'c':[hit('PEPTIDE',[(0,20),(90,91)])]}
    assess(args,[ref,cand],MockEvidence(hits,{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{})
    assert not table(tmp_path/'refinement/annotation_change_proposals.tsv')


def test_short_unique_peptide_supports_and_reports_its_length(tmp_path):
    # No length filter: the user judges short-peptide chance matches (decoy FDR by length).
    args=config(tmp_path)
    ref=reference('r',[(0,30),(90,120)]);cand=candidate('c',[(0,120)])
    assess(args,[ref,cand],MockEvidence({'c':[hit('PEPTID',[(36,54)],12)]},{'PEPTID':1}),{'PEPTID':dict(quality='unknown')},{})
    row,=splice_rows(tmp_path)
    assert row['confidence']=='peptide_supported'
    assert json.loads(row['peptide_loci'])['PEPTID']['length']==6


def test_gene_with_only_multilocus_peptides_is_ambiguous(tmp_path):
    args=config(tmp_path);args.bam=None
    ref=reference('r',[(0,60)])
    assess(args,[ref],MockEvidence({'r':[hit('PEPTIDE',[(0,21)])]},{'PEPTIDE':2}),{'PEPTIDE':dict(quality='unknown')},{})
    assert table(tmp_path/'annotation/gene_evidence.tsv')[0]['category']=='ambiguous'


# ---------------------------------------------------------------- ORF reconstruction hygiene

class Reference:
    def __init__(self,sequence):self.sequence=sequence
    def fetch(self,chrom,a,b):return self.sequence[a:b]


def test_noncanonical_introns_are_flagged_on_the_model_strand():
    sequence='A'*10+'GT'+'A'*16+'AG'+'A'*10+'CT'+'A'*16+'AC'+'A'*10
    plus=Model('p','RNA_ORF','p','g','chr1','+',cds=[(0,10),(30,40)])
    reverse=Model('r','RNA_ORF','r','g','chr1','+',cds=[(30,40),(60,70)])
    minus=Model('m','RNA_ORF','m','g','chr1','-',cds=[(30,40),(60,70)])
    flag_noncanonical_introns([plus,reverse,minus],Reference(sequence))
    assert plus.issues==[] and reverse.issues==['noncanonical_intron'] and minus.issues==[]


def test_reverse_orfs_on_stranded_transcripts_are_not_projected(tmp_path):
    transcripts=tmp_path/'t.gtf'
    transcripts.write_text(''.join([
        'chr1\tS\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "spliced";\n',
        'chr1\tS\texon\t201\t300\t.\t+\t.\tgene_id "g1"; transcript_id "spliced";\n',
        'chr1\tS\texon\t401\t600\t.\t.\t.\tgene_id "g2"; transcript_id "unstranded";\n']))
    gff=tmp_path/'td.gff3'
    rows=[]
    for tid,orf,strand in [('spliced','sense','+'),('spliced','reverse','-'),('unstranded','unstranded_reverse','-')]:
        rows.append(f'{tid}\tT\tmRNA\t10\t190\t.\t{strand}\t.\tID={orf};Parent={tid}.g\n')
        rows.append(f'{tid}\tT\tCDS\t10\t190\t.\t{strand}\t0\tID=cds.{orf};Parent={orf}\n')
    gff.write_text(''.join(rows))
    out=tmp_path/'orfs.gtf';convert_transdecoder_to_genome_gtf(gff,transcripts,out)
    ids={line.split('transcript_id "')[1].split('"')[0] for line in out.read_text().splitlines()}
    assert ids=={'sense','unstranded_reverse'}
    assert json.loads(out.with_suffix('.summary.json').read_text())['reverse_strand_orfs_on_stranded_transcripts_skipped']==1


@pytest.mark.parametrize('library,flag',[('unstranded',None),('firststrand','--rf'),('secondstrand','--fr')])
def test_stringtie_receives_library_strandedness(tmp_path,monkeypatch,library,flag):
    from gapms_v1 import rnaseq
    captured=[]
    class Stop(Exception):pass
    def run(command,*args,**kwargs):captured.append(command);raise Stop
    monkeypatch.setattr(rnaseq,'find_tool',lambda name:name)
    monkeypatch.setattr(rnaseq,'run_command',run)
    args=SimpleNamespace(output=tmp_path,bam=Path('r.bam'),threads=1,resume=False,rna_library=library,genetic_code=1,min_orf_aa=100)
    (tmp_path/'work').mkdir()
    with pytest.raises(Stop):rnaseq.reconstruct(args,None,{},[])
    command=[str(c) for c in captured[0]]
    assert ('--rf' in command)==(flag=='--rf') and ('--fr' in command)==(flag=='--fr')


def test_duplicated_exon_combinations_are_different_origins(tmp_path):
    # Tandem duplicated exons: the annotated pair (A1,B1) and an RNA path (A2,B2)
    # both encode PEPTIDE; a linking path (A2,B1) shares only one exon with each.
    # Sharing one exon is not the same origin: the peptide has two loci.
    args=config(tmp_path)
    ref=reference('r',[(100,129),(200,230)])
    c=candidate('c',[(0,29),(300,330)]);link=candidate('link',[(0,29),(200,230)])
    annotated=hit('PEPTIDE',[(115,129),(200,207)],5)
    new=hit('PEPTIDE',[(15,29),(300,307)],5);linking=hit('PEPTIDE',[(15,29),(200,207)],5)
    assess(args,[ref,c,link],MockEvidence({'r':[annotated],'c':[new],'link':[linking]},{'PEPTIDE':2}),{'PEPTIDE':dict(quality='unknown')},{})
    rows=table(tmp_path/'refinement/annotation_change_proposals.tsv')
    assert rows and all(r['confidence']=='ambiguous' for r in rows)


def test_gene_copies_bridged_by_readthrough_are_two_loci(tmp_path):
    store=EvidenceStore(tmp_path/'e.sqlite')
    try:
        # Copy 1, read-through (copy-1 short exon + copy-2 long exon), copy 2.
        for blocks in ([(1000,1009),(1400,1442)],[(1000,1009),(8200,8242)],[(7800,7809),(8200,8242)]):
            store.add('ISQVVDAIDIAQATMTK','chr1','+',blocks,'ISQVVDAIDIAQATMTK','genome')
        assert store.counts()=={'ISQVVDAIDIAQATMTK':2}
    finally:store.db.close()


class Genome:
    def __init__(self,sequence):self.sequence=sequence
    def fetch(self,chrom,a,b):return self.sequence[a:b]
    def get_reference_length(self,chrom):return len(self.sequence)


@pytest.mark.parametrize('upstream_codon,blocked',[('AAA',True),('TAA',False)])
def test_rna_start_preceded_by_trypsin_site_is_not_supported(tmp_path,upstream_codon,blocked):
    args=config(tmp_path)
    ref=reference('r',[(33,93)]);cand=candidate('c',[(3,123)])
    genome=Genome(upstream_codon+'ATG'+'C'*200)
    assess(args,[ref,cand],MockEvidence({'c':[hit('MPEPTIDEK',[(3,30)],0)]},{'MPEPTIDEK':1}),{'MPEPTIDEK':dict(quality='unknown')},{},reference=genome)
    row,=table(tmp_path/'refinement/alternative_start_candidates.tsv')
    if blocked:
        assert row['confidence']=='not_supported' and 'trypsin_site_before_proposed_start' in json.loads(row['confidence_reasons'])
    else:
        assert row['confidence']=='peptide_supported'


def test_peptide_explained_by_two_alternative_junctions_supports_neither(tmp_path):
    # Two donors translate identically into the same acceptor: the peptide proves
    # that one of the junctions exists, not which one.
    args=config(tmp_path)
    ref=reference('r',[(0,30),(90,120)])
    c1=candidate('c1',[(0,33),(90,120)]);c2=candidate('c2',[(0,45),(90,120)])
    # Short donor-side segments, shared acceptor side (majority): one origin, two paths.
    h1=hit('PEPTIDE',[(27,33),(90,105)],9);h2=hit('PEPTIDE',[(39,45),(90,105)],13)
    assess(args,[ref,c1,c2],MockEvidence({'c1':[h1],'c2':[h2]},{'PEPTIDE':1}),{'PEPTIDE':dict(quality='unknown')},{('chr1',33,90,'+'):5,('chr1',45,90,'+'):5})
    rows=splice_rows(tmp_path)
    assert len(rows)==2 and all(r['confidence']=='not_supported' for r in rows)
    assert all('peptide_path_ambiguous_within_locus' in json.loads(r['confidence_reasons']) for r in rows)


def test_min_length_default_and_no_discovery_length_option():
    from gapms_v1.cli import parser
    args=parser().parse_args(['-a','g.fa','-p','p.txt'])
    assert args.min_length==6 and not hasattr(args,'min_discovery_length')
    with pytest.raises(SystemExit):parser().parse_args(['-a','g.fa','-p','p.txt','--min-discovery-length','9'])
