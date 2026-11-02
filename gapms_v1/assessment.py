"""Conservative model assessment and event-specific annotation proposals."""
from collections import Counter, defaultdict
import json

from gapms._peplocus.sequence import canonical
from .common import stable_id, write_table
from .models import write_models, write_proteins

MODEL_FIELDS=['model_id','source','original_transcript_id','gene_id','chrom','start','end','strand','category','distinct_peptides','locus_discriminating_peptides','non_nested_discriminating_peptides','peptides','mapping_ids','complete_start','complete_stop','issues','quality','reasons']
EVENT_FIELDS=['event_id','event_type','candidate_model_id','reference_model_ids','annotation_sources','chrom','start','end','strand','candidate_cds_blocks','reference_cds_blocks','relationship','peptides','mapping_ids','distinct_diagnostic_peptides','mapping_multiplicity','rna_junctions','minimum_rna_support','confidence','reasons','limitations']


def non_nested(sequences):
    result=[]
    for sequence in sorted(set(sequences),key=lambda s:(-len(s),s)):
        if not any(sequence in longer for longer in result):result.append(sequence)
    return result


def interval_overlap(left,right):
    return any(max(a,c)<min(b,d) for a,b in left for c,d in right)


class ModelIndex:
    """Bin transcript extents; avoid candidate-by-entire-annotation comparisons."""
    def __init__(self,models):
        self.bins=defaultdict(list)
        for m in models:
            a=min(x for x,y in m.intervals);b=max(y for x,y in m.intervals)
            for bin_id in range(a//100000,(b-1)//100000+1):self.bins[(m.chrom,bin_id)].append(m)

    def overlap(self,model):
        a=min(x for x,y in model.intervals);b=max(y for x,y in model.intervals);found={}
        for bin_id in range(a//100000,(b-1)//100000+1):
            for m in self.bins.get((model.chrom,bin_id),[]):
                x=min(x for x,y in m.intervals);y=max(y for x,y in m.intervals)
                if max(a,x)<min(b,y):found[m.id]=m
        return sorted(found.values(),key=lambda m:m.id)


def relationship(candidate,references):
    if not references:return 'intergenic'
    same=[m for m in references if m.strand==candidate.strand]
    if not same:return 'antisense'
    overlapping=[m for m in same if interval_overlap(candidate.cds,m.cds)]
    if not overlapping:
        if any(interval_overlap(candidate.cds,m.exons) for m in same):return 'UTR_associated'
        return 'intronic'
    def frame_at(m,position):
        offset=0
        for a,b in m.coding:
            if a<=position<b:return (offset+(position-a if m.strand=='+' else b-1-position))%3
            offset+=b-a
    for m in overlapping:
        for a,b in candidate.cds:
            for c,d in m.cds:
                position=max(a,c)
                if position<min(b,d) and frame_at(candidate,position)==frame_at(m,position):return 'new_transcript_isoform'
    return 'overlapping_alternative_frame'


def boundary(model,which):
    segment=model.coding[0] if which=='start' else model.coding[-1]
    if which=='start':return segment[0] if model.strand=='+' else segment[1]
    return segment[1] if model.strand=='+' else segment[0]


def beyond(hit,reference,which):
    point=boundary(reference,which)
    left=min(a for a,b in hit['blocks']);right=max(b for a,b in hit['blocks'])
    if which=='start':return left<point if reference.strand=='+' else right>point
    return right>point if reference.strand=='+' else left<point


def assess(args,models,store,metadata,rna_junctions):
    hits=store.model_hits();counts=store.counts();by_id={m.id:m for m in models}
    qualities=defaultdict(list)
    for p,item in metadata.items():qualities[canonical(p)].append(item['quality'])
    def quality(sequences):
        if sequences and all('qvalue_pass' in qualities[p] for p in sequences):return 'qvalue_pass'
        if any('qvalue_fail' in qualities[p] and 'qvalue_pass' not in qualities[p] for p in sequences):return 'qvalue_fail'
        return 'unknown'
    records=[];categories={}
    for model in models:
        placements=hits[model.id];sequences={h['canonical'] for h in placements}
        unique={p for p in sequences if counts[p]==1}
        if not model.cds:category='no_CDS_annotation'
        elif not sequences:category='no_observed_peptide_support'
        elif model.issues:category='ambiguous'
        elif unique:category='peptide_supported'
        else:category='ambiguous'
        categories[model.id]=category
        records.append(dict(model_id=model.id,source=model.source,original_transcript_id=model.original_id,gene_id=model.gene,
            chrom=model.chrom,start=min(a for a,b in model.intervals),end=max(b for a,b in model.intervals),strand=model.strand,
            category=category,distinct_peptides=len(sequences),locus_discriminating_peptides=len(unique),
            non_nested_discriminating_peptides=len(non_nested(unique)),peptides=sorted(sequences),mapping_ids=sorted({h['mapping_id'] for h in placements}),
            complete_start=model.complete_start,complete_stop=model.complete_stop,issues=sorted(set(model.issues)),quality=quality(sequences),
            reasons='no_observed_evidence_is_not_evidence_of_incorrect_annotation' if not sequences else 'global_genomic_placement_multiplicity; I/L aliases counted once'))
    annotation=[m for m in models if m.source!='RNA_ORF'];candidates=[m for m in models if m.source=='RNA_ORF']
    if args.gtf:
        directory=args.output/'annotation';directory.mkdir(exist_ok=True)
        original_rows=[r for r in records if r['source']!='RNA_ORF']
        write_table(directory/'gene_model_evidence.tsv',MODEL_FIELDS,original_rows)
        for category,filename in [('no_observed_peptide_support','no_peptide_support_gene_models.tsv'),('ambiguous','ambiguous_gene_models.tsv')]:
            write_table(directory/filename,MODEL_FIELDS,(r for r in original_rows if r['category']==category))
        write_models([m for m in annotation if categories[m.id]=='peptide_supported'],directory/'supported_gene_models.gtf')
        genes=defaultdict(list)
        for r in original_rows:genes[(r['source'],r['gene_id'])].append(r)
        gene_rows=[]
        for (source,gene),rows in sorted(genes.items()):
            sequences={p for row in rows for p in row['peptides']}
            gene_rows.append(dict(source=source,gene_id=gene,models=len(rows),distinct_peptides=len(sequences),locus_discriminating_peptides=sum(counts[p]==1 for p in sequences),category='peptide_supported' if sequences else 'no_observed_peptide_support'))
        write_table(directory/'gene_evidence.tsv',['source','gene_id','models','distinct_peptides','locus_discriminating_peptides','category'],gene_rows)
    if args.bam:
        directory=args.output/'orfs'
        write_table(directory/'orf_peptide_evidence.tsv',MODEL_FIELDS,(r for r in records if r['source']=='RNA_ORF'))
        write_proteins(candidates,directory/'candidate_orfs.faa')
        write_table(directory/'partial_or_ambiguous_orfs.tsv',MODEL_FIELDS,(r for r in records if r['source']=='RNA_ORF' and (not r['complete_start'] or not r['complete_stop'] or r['category']=='ambiguous')))
    # Locus specificity alone does not distinguish alternative coding structures.
    sequence_models=defaultdict(set)
    for model in models:
        for hit in hits[model.id]:sequence_models[hit['canonical']].add(model.id)
    reference_index=ModelIndex(annotation)
    reference_loci=defaultdict(set)
    for m in annotation:
        for hit in hits[m.id]:reference_loci[hit['canonical']].add((m.chrom,m.strand,tuple(tuple(b) for b in hit['blocks'])))
    signature={m.signature() for m in annotation if m.coding}
    events=[]
    for candidate in candidates:
        if not hits[candidate.id] or candidate.signature() in signature:continue
        refs=reference_index.overlap(candidate)
        relation=relationship(candidate,refs)
        all_hits=hits[candidate.id]
        diagnostic=[h for h in all_hits if (candidate.chrom,candidate.strand,tuple(tuple(b) for b in h['blocks'])) not in reference_loci[h['canonical']]]
        rna_keys=[(candidate.chrom,a,b,candidate.strand) for a,b in candidate.junctions]
        minimum=min((rna_junctions.get(k,0) for k in rna_keys),default=0)
        def event(kind,placements,related,extra='',rna_required=False):
            if not placements:return
            sequences={h['canonical'] for h in placements}
            unique={p for p in sequences if counts[p]==1}
            coherent=not candidate.issues and bool(candidate.protein)
            compatible=set.intersection(*(sequence_models[p] for p in sequences)) - {candidate.id}
            competitors=sorted(mid for mid in compatible if by_id[mid].signature()!=candidate.signature())
            high=(len(non_nested(unique))>=2 and coherent and candidate.complete_start and candidate.complete_stop and quality(unique)=='qvalue_pass'
                  and (not rna_keys or minimum>=args.min_junction_support))
            if rna_required and minimum<args.min_junction_support:high=False
            if getattr(args,'splice_search_capped',False) or competitors:high=False
            confidence='high_confidence_candidate' if high else ('rna_peptide_supported' if minimum>=args.min_junction_support and unique and coherent else 'peptide_supported' if unique and coherent else 'ambiguous')
            if competitors:confidence='ambiguous'
            limitations=['candidate_annotation_change_requires_review','peptide_identification_quality_'+quality(sequences),'RNA_junctions_do_not_prove_full_length_connectivity']
            if any(counts[p]>1 for p in sequences):limitations.append('non_discriminating_peptides_present')
            if len(non_nested(unique))<2:limitations.append('fewer_than_two_non_nested_locus_discriminating_peptides')
            if competitors:limitations.append('equally_compatible_coding_models:'+','.join(competitors))
            if candidate.issues:limitations.extend(candidate.issues)
            if getattr(args,'splice_search_capped',False):limitations.append('splice_path_search_capped_uniqueness_incomplete')
            if kind=='alternative_start':limitations.extend(['proteolytic_processing_not_excluded','initiation_site_not_experimentally_confirmed'])
            if kind=='alternative_stop':limitations.append('stop_position_inferred_from_ORF_sequence')
            event_id=stable_id('E',kind,candidate.id,sorted(m.id for m in related))
            events.append(dict(event_id=event_id,event_type=kind,candidate_model_id=candidate.id,
                reference_model_ids=sorted(m.id for m in related),annotation_sources=sorted({m.source for m in related}),chrom=candidate.chrom,
                start=min(a for a,b in candidate.cds),end=max(b for a,b in candidate.cds),strand=candidate.strand,
                candidate_cds_blocks=sorted(candidate.coding),reference_cds_blocks={m.id:sorted(m.coding) for m in related},relationship=relation,
                peptides=sorted(sequences),mapping_ids=sorted({h['mapping_id'] for h in placements}),distinct_diagnostic_peptides=len(sequences),
                mapping_multiplicity={p:counts[p] for p in sorted(sequences)},rna_junctions=[dict(chrom=k[0],start=k[1],end=k[2],strand=k[3],reads=rna_junctions.get(k,0)) for k in rna_keys],
                minimum_rna_support=minimum,confidence=confidence,reasons=extra or 'direct_peptide_compatibility_with_candidate_coding_path',limitations=limitations))
        if diagnostic:event('novel_orf',diagnostic,refs)
        elif relation=='intergenic':event('novel_orf',all_hits,refs,'no_comparable_annotation_provided' if not annotation else 'intergenic_peptide_supported_ORF')
        same=[m for m in refs if m.strand==candidate.strand and m.coding and interval_overlap(m.cds,candidate.cds)]
        if same:
            # A peptide elsewhere in the protein never supports a changed boundary/junction.
            for ref in same:
                if boundary(candidate,'start')!=boundary(ref,'start') and candidate.complete_start:
                    terminal=[h for h in diagnostic if beyond(h,ref,'start') and (h['aa_start']==0 or (h['aa_start']==1 and candidate.protein.startswith('M') and candidate.protein[1:2] in 'ACGPSTV'))]
                    event('alternative_start',terminal,[ref],'diagnostic_N_terminal_peptide; ATG_start; tryptic_and_initiator_processing_context_recorded')
                if boundary(candidate,'end')!=boundary(ref,'end'):
                    terminal=[h for h in diagnostic if beyond(h,ref,'end')]
                    event('alternative_stop',terminal,[ref],'diagnostic_peptide_extends_beyond_annotated_coding_end')
            annotated_junctions={j for m in same for j in m.junctions}
            new_junctions=set(candidate.junctions)-annotated_junctions
            crossing=[h for h in diagnostic if any((a[1],b[0]) in new_junctions for a,b in zip(h['blocks'],h['blocks'][1:]))]
            if new_junctions:event('alternative_splice',crossing,same,'peptide_directly_crosses_junction_absent_from_overlapping_reference_CDS',True)
    events.sort(key=lambda e:(e['chrom'],e['start'],e['strand'],e['event_type'],e['event_id']))
    if args.bam:
        directory=args.output/'orfs'
        novel=[e for e in events if e['event_type']=='novel_orf']
        write_table(directory/'novel_orf_candidates.tsv',EVENT_FIELDS,novel)
        novel_models=[by_id[i] for i in sorted({e['candidate_model_id'] for e in novel})]
        write_models(novel_models,directory/'novel_orf_candidates.gtf');write_proteins(novel_models,directory/'novel_orf_candidates.faa')
    if args.gtf and args.bam:
        directory=args.output/'refinement';directory.mkdir(exist_ok=True)
        for kind in ['alternative_start','alternative_stop','alternative_splice']:
            selected=[e for e in events if e['event_type']==kind]
            write_table(directory/(kind+'_candidates.tsv'),EVENT_FIELDS,selected)
            write_models([by_id[i] for i in sorted({e['candidate_model_id'] for e in selected})],directory/(kind+'_candidates.gtf'))
        write_table(directory/'annotation_change_proposals.tsv',EVENT_FIELDS,events)
    # Traceable individual protein occurrences, including N-terminal cleavage context.
    occurrence_rows=[]
    for model in models:
        for h in hits[model.id]:
            previous=model.protein[h['aa_start']-1] if h['aa_start'] else ''
            first=model.protein[h['aa_start']:h['aa_start']+1]
            occurrence_rows.append(dict(model_id=model.id,mapping_id=h['mapping_id'],peptide=h['peptide'],aa_start=h['aa_start'],aa_end=h['aa_end'],previous_aa=previous,n_terminal_context='protein_N_terminus' if not previous else 'tryptic' if previous in 'KR' and first!='P' else 'nontryptic',initiator_methionine_processing=h['aa_start']==1 and model.protein.startswith('M') and first in 'ACGPSTV'))
    write_table(args.output/'mapping/model_peptide_occurrences.tsv',['model_id','mapping_id','peptide','aa_start','aa_end','previous_aa','n_terminal_context','initiator_methionine_processing'],occurrence_rows)
    return dict(models=len(models),annotation_models=len(annotation),candidate_orfs=len(candidates),model_categories=dict(Counter(categories.values())),annotation_events=dict(Counter(e['event_type'] for e in events)),event_confidences=dict(Counter(e['confidence'] for e in events)))
