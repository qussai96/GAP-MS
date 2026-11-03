"""Conservative model assessment and event-specific annotation proposals."""
from collections import Counter, defaultdict
import json

from gapms._peplocus.sequence import canonical
from .common import stable_id, write_table
from .models import write_models, write_proteins

MODEL_FIELDS=['model_id','source','original_transcript_id','gene_id','chrom','start','end','strand','category','distinct_peptides','locus_discriminating_peptides','non_nested_discriminating_peptides','peptides','mapping_ids','complete_start','complete_stop','issues','quality','reasons']
EVENT_FIELDS=['event_id','event_type','candidate_model_id','reference_model_ids','annotation_sources','chrom','start','end','strand','candidate_cds_blocks','reference_cds_blocks','relationship','peptides','mapping_ids','distinct_diagnostic_peptides','mapping_multiplicity','rna_junctions','minimum_rna_support','confidence','reasons','limitations','coding_structure_id','local_peptide_blocks','diagnostic_junctions','splice_chain_changed','whole_model_status','whole_model_review_reasons','uniqueness_scope']


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
            reasons='no_observed_evidence_is_not_evidence_of_incorrect_annotation' if not sequences else 'placement_multiplicity_within_searched_genome_and_splice_paths; I/L aliases counted once'))
    annotation=[m for m in models if m.source!='RNA_ORF'];candidates=[m for m in models if m.source=='RNA_ORF']
    if args.gtf:
        directory=args.output/'annotation';directory.mkdir(exist_ok=True)
        original_rows=[r for r in records if r['source']!='RNA_ORF']
        write_table(directory/'gene_model_evidence.tsv',MODEL_FIELDS,original_rows)
        for category,filename in [('no_observed_peptide_support','no_peptide_support_gene_models.tsv'),('ambiguous','ambiguous_gene_models.tsv')]:
            write_table(directory/filename,MODEL_FIELDS,(r for r in original_rows if r['category']==category))
        write_models([m for m in annotation if categories[m.id]=='peptide_supported'],directory/'supported_gene_models.gtf')
        genes=defaultdict(list)
        for r in original_rows:genes[(r['source'],r['chrom'],r['strand'],r['gene_id'])].append(r)
        gene_rows=[]
        for (source,chrom,strand,gene),rows in sorted(genes.items()):
            sequences={p for row in rows for p in row['peptides']}
            gene_rows.append(dict(source=source,chrom=chrom,strand=strand,gene_id=gene,models=len(rows),distinct_peptides=len(sequences),locus_discriminating_peptides=sum(counts[p]==1 for p in sequences),category='peptide_supported' if sequences else 'no_observed_peptide_support'))
        write_table(directory/'gene_evidence.tsv',['source','chrom','strand','gene_id','models','distinct_peptides','locus_discriminating_peptides','category'],gene_rows)
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
    events=[];isoform_checks=[]
    annotated_paths=defaultdict(list)
    for model in annotation:
        if model.coding:annotated_paths[model.signature()].append(model.id)
    for candidate in candidates:
        if not hits[candidate.id]:continue
        refs=reference_index.overlap(candidate)
        relation=relationship(candidate,refs)
        all_hits=hits[candidate.id]
        diagnostic=[h for h in all_hits if (candidate.chrom,candidate.strand,tuple(tuple(b) for b in h['blocks'])) not in reference_loci[h['canonical']]]
        same=[m for m in refs if m.coding and relationship(candidate,[m])=='new_transcript_isoform']
        known_starts=[m.id for m in same if boundary(m,'start')==boundary(candidate,'start')]
        known_ends=[m.id for m in same if boundary(m,'end')==boundary(candidate,'end')]
        annotated_junctions={j for m in same for j in m.junctions}
        new_junctions=set(candidate.junctions)-annotated_junctions
        exact=sorted(annotated_paths.get(candidate.signature(),[]))
        isoform_checks.append(dict(candidate_model_id=candidate.id,chrom=candidate.chrom,strand=candidate.strand,coding_structure_id=stable_id('CDS',candidate.signature()),exact_annotated_isoform_ids=exact,comparable_annotated_isoform_ids=sorted(m.id for m in same),known_start_isoform_ids=sorted(known_starts),known_end_isoform_ids=sorted(known_ends),new_junctions=sorted(new_junctions),observed_peptide_sequences=len({h['canonical'] for h in all_hits}),diagnostic_peptide_sequences=len({h['canonical'] for h in diagnostic}),interpretation='exact_annotated_isoform' if exact else 'existing_isoforms_explain_all_observed_peptides' if not diagnostic else 'diagnostic_evidence_requires_feature_specific_review'))
        if exact:continue
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
                  and set(candidate.protein)<=set('ACDEFGHIKLMNPQRSTVWY')
                  and (not rna_keys or minimum>=args.min_junction_support))
            if rna_required and minimum<args.min_junction_support:high=False
            if getattr(args,'splice_search_capped',False) or competitors:high=False
            confidence='high_confidence_candidate' if high else ('rna_peptide_supported' if minimum>=args.min_junction_support and unique and coherent else 'peptide_supported' if unique and coherent else 'ambiguous')
            if competitors:confidence='ambiguous'
            limitations=['candidate_annotation_change_requires_review','peptide_identification_quality_'+quality(sequences),'RNA_junctions_do_not_prove_full_length_connectivity']
            if any(counts[p]>1 for p in sequences):limitations.append('non_discriminating_peptides_present')
            if len(non_nested(unique))<2:limitations.append('fewer_than_two_non_nested_locus_discriminating_peptides')
            if competitors:limitations.append('equally_compatible_coding_models:'+','.join(competitors))
            if set(candidate.protein)-set('ACDEFGHIKLMNPQRSTVWY*'):limitations.append('unresolved_amino_acids_in_candidate_ORF')
            if candidate.issues:limitations.extend(candidate.issues)
            if getattr(args,'splice_search_capped',False):limitations.append('splice_path_search_capped_uniqueness_incomplete')
            if kind=='alternative_start':limitations.extend(['proteolytic_processing_not_excluded','initiation_site_not_experimentally_confirmed'])
            if kind=='alternative_stop':limitations.append('stop_position_inferred_from_ORF_sequence')
            if kind=='c_terminal_extension':limitations.append('termination_site_unresolved_no_complete_stop_codon')
            limitations.append('uniqueness_is_relative_to_searched_splice_paths_not_all_possible_splicing')
            review_reasons=[]
            if not candidate.complete_start:review_reasons.append('incomplete_start')
            if not candidate.complete_stop:review_reasons.append('incomplete_stop')
            if candidate.issues:review_reasons.append('candidate_sequence_or_phase_issues')
            if confidence=='ambiguous':review_reasons.append('ambiguous_diagnostic_evidence')
            if any(rna_junctions.get(k,0)<args.min_junction_support for k in rna_keys):review_reasons.append('coding_junction_without_required_RNA_support')
            if len(non_nested(unique))<2:review_reasons.append('fewer_than_two_independent_diagnostic_peptides')
            genes={(m.source,m.gene) for m in related if m.strand==candidate.strand and interval_overlap(candidate.cds,m.cds)}
            if len(genes)>1:review_reasons.append('overlaps_multiple_annotation_gene_ids_review_fusion_or_isoforms')
            if any(any(sum(max(0,min(b,d)-max(a,c)) for c,d in candidate.coding)<b-a for a,b in m.coding) for m in related if m.coding):review_reasons.append('compared_annotation_CDS_not_fully_preserved_local_support_does_not_validate_deletions')
            changed=any(set(candidate.junctions)!=set(m.junctions) for m in related if m.coding)
            if kind in {'alternative_stop','c_terminal_extension'} and changed:limitations.append('terminal_change_accompanied_by_splice_chain_change')
            # This is a context model, never a whole-gene replacement inferred from local evidence.
            status='local_evidence_only' if review_reasons else 'complete_model_for_manual_review'
            diagnostic_junctions=sorted({(a[1],b[0]) for h in placements for a,b in zip(h['blocks'],h['blocks'][1:]) if (a[1],b[0]) in candidate.junctions and not any((a[1],b[0]) in m.junctions for m in related)})
            event_id=stable_id('E',kind,candidate.id,sorted(m.id for m in related))
            events.append(dict(event_id=event_id,event_type=kind,candidate_model_id=candidate.id,
                reference_model_ids=sorted(m.id for m in related),annotation_sources=sorted({m.source for m in related}),chrom=candidate.chrom,
                start=min(a for a,b in candidate.cds),end=max(b for a,b in candidate.cds),strand=candidate.strand,
                candidate_cds_blocks=sorted(candidate.coding),reference_cds_blocks={m.id:sorted(m.coding) for m in related},relationship=relation,
                peptides=sorted(sequences),mapping_ids=sorted({h['mapping_id'] for h in placements}),distinct_diagnostic_peptides=len(sequences),
                mapping_multiplicity={p:counts[p] for p in sorted(sequences)},rna_junctions=[dict(chrom=k[0],start=k[1],end=k[2],strand=k[3],reads=rna_junctions.get(k,0)) for k in rna_keys],
                minimum_rna_support=minimum,confidence=confidence,coding_structure_id=stable_id('CDS',candidate.signature()),local_peptide_blocks=sorted({tuple(b) for h in placements for b in h['blocks']}),diagnostic_junctions=diagnostic_junctions,splice_chain_changed=changed,whole_model_status=status,whole_model_review_reasons=review_reasons,uniqueness_scope='searched_assembly_and_supplied_splice_paths',reasons=extra or 'direct_peptide_compatibility_with_candidate_coding_path',limitations=limitations))
        if diagnostic and relation!='new_transcript_isoform':event('novel_orf',diagnostic,refs)
        elif relation=='intergenic':event('novel_orf',all_hits,refs,'no_comparable_annotation_provided' if not annotation else 'intergenic_peptide_supported_ORF')
        if same:
            # A peptide elsewhere in the protein never supports a changed boundary/junction.
            for ref in same:
                if not known_starts and boundary(candidate,'start')!=boundary(ref,'start') and candidate.complete_start:
                    terminal=[h for h in diagnostic if beyond(h,ref,'start') and (h['aa_start']==0 or (h['aa_start']==1 and candidate.protein.startswith('M') and candidate.protein[1:2] in 'ACGPSTV'))]
                    event('alternative_start',terminal,[ref],'diagnostic_N_terminal_peptide; ATG_start; tryptic_and_initiator_processing_context_recorded')
                if not known_ends and boundary(candidate,'end')!=boundary(ref,'end'):
                    terminal=[h for h in diagnostic if beyond(h,ref,'end')]
                    event('alternative_stop' if candidate.complete_stop else 'c_terminal_extension',terminal,[ref],'diagnostic_peptide_extends_beyond_annotated_coding_end')
            crossing=[h for h in diagnostic if any((a[1],b[0]) in new_junctions for a,b in zip(h['blocks'],h['blocks'][1:]))]
            if new_junctions:event('alternative_splice',crossing,same,'peptide_directly_crosses_junction_absent_from_overlapping_reference_CDS',True)
    if args.bam:
        write_table(args.output/'orfs/annotation_isoform_checks.tsv',['candidate_model_id','chrom','strand','coding_structure_id','exact_annotated_isoform_ids','comparable_annotated_isoform_ids','known_start_isoform_ids','known_end_isoform_ids','new_junctions','observed_peptide_sequences','diagnostic_peptide_sequences','interpretation'],isoform_checks)
    # A coding model has one primary classification; component events remain auditable.
    grouped=defaultdict(list)
    for row in events:grouped[row['candidate_model_id']].append(row)
    primary=[]
    for model_id,rows in sorted(grouped.items()):
        kinds=sorted({row['event_type'] for row in rows})
        primary.append(dict(candidate_model_id=model_id,category=kinds[0] if len(kinds)==1 else 'combined_refinement',component_changes=kinds,event_ids=[row['event_id'] for row in rows]))
    if args.bam:
        write_table(args.output/'orfs/candidate_classification.tsv',['candidate_model_id','category','component_changes','event_ids'],primary)
    events.sort(key=lambda e:(e['chrom'],e['start'],e['strand'],e['event_type'],e['event_id']))
    if args.bam:
        directory=args.output/'orfs'
        novel=[e for e in events if e['event_type']=='novel_orf']
        write_table(directory/'novel_orf_candidates.tsv',EVENT_FIELDS,novel)
        novel_models=[by_id[i] for i in sorted({e['candidate_model_id'] for e in novel})]
        write_models(novel_models,directory/'novel_orf_candidates.gtf');write_proteins(novel_models,directory/'novel_orf_candidates.faa')
    if args.gtf and args.bam:
        directory=args.output/'refinement';directory.mkdir(exist_ok=True)
        for kind in ['alternative_start','alternative_stop','alternative_splice','c_terminal_extension']:
            selected=[e for e in events if e['event_type']==kind]
            write_table(directory/(kind+'_candidates.tsv'),EVENT_FIELDS,selected)
            write_models([by_id[i] for i in sorted({e['candidate_model_id'] for e in selected})],directory/(kind+'_candidates.gtf'))
        write_table(directory/'annotation_change_proposals.tsv',EVENT_FIELDS,events)
        write_table(directory/'local_evidence_only.tsv',EVENT_FIELDS,(e for e in events if e['whole_model_status']=='local_evidence_only'))
        write_table(directory/'complete_models_for_manual_review.tsv',EVENT_FIELDS,(e for e in events if e['whole_model_status']=='complete_model_for_manual_review'))
        with (directory/'local_peptide_evidence.bed').open('w') as stream:
            for e in events:
                for i,(a,b) in enumerate(e['local_peptide_blocks'],1):
                    stream.write(f"{e['chrom']}\t{a}\t{b}\t{e['event_id']}.{i}\t0\t{e['strand']}\n")
        (directory/'README.txt').write_text('Candidate GTF files contain full RNA models for context and manual review only. They are not corrected annotations or automatic replacements. Local peptide blocks and diagnostic junctions in the event tables delimit the direct evidence. local_peptide_evidence.bed contains separate supported blocks (BED6), not reconstructed CDSs. Complete models still require independent full-structure validation. IDs are run-scoped; join with the run/crop. Coding structure IDs collapse duplicate CDS paths within an assembly; retain run identity across assemblies.\n')
    # Traceable individual protein occurrences, including N-terminal cleavage context.
    occurrence_rows=[]
    for model in models:
        for h in hits[model.id]:
            previous=model.protein[h['aa_start']-1] if h['aa_start'] else ''
            first=model.protein[h['aa_start']:h['aa_start']+1]
            occurrence_rows.append(dict(model_id=model.id,mapping_id=h['mapping_id'],peptide=h['peptide'],aa_start=h['aa_start'],aa_end=h['aa_end'],previous_aa=previous,n_terminal_context='protein_N_terminus' if not previous else 'tryptic' if previous in 'KR' and first!='P' else 'nontryptic',initiator_methionine_processing=h['aa_start']==1 and model.protein.startswith('M') and first in 'ACGPSTV'))
    write_table(args.output/'mapping/model_peptide_occurrences.tsv',['model_id','mapping_id','peptide','aa_start','aa_end','previous_aa','n_terminal_context','initiator_methionine_processing'],occurrence_rows)
    return dict(isoform_checks=dict(Counter(row['interpretation'] for row in isoform_checks)),distinct_candidate_cds_paths=len({m.signature() for m in candidates if m.coding}),whole_model_statuses=dict(Counter(e['whole_model_status'] for e in events)),models=len(models),annotation_models=len(annotation),candidate_orfs=len(candidates),model_categories=dict(Counter(categories.values())),candidate_categories=dict(Counter(r['category'] for r in primary)),annotation_events=dict(Counter(e['event_type'] for e in events)),event_confidences=dict(Counter(e['confidence'] for e in events)))
