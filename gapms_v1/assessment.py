"""Conservative model assessment and event-specific annotation proposals."""
from collections import Counter, defaultdict
import json

from peplocus.sequence import canonical
from .common import stable_id, write_table
from .models import write_models, write_proteins
from .mapping import ClusterIndex
from .events import SUPPORTED, feature_consistent, feature_id, judge, merge_events, quality_lookup

MODEL_FIELDS=['model_id','source','original_transcript_id','gene_id','chrom','start','end','strand','category','distinct_peptides','locus_discriminating_peptides','non_nested_discriminating_peptides','peptides','mapping_ids','complete_start','complete_stop','issues','quality','reasons']
EVENT_FIELDS=['event_id', 'event_type', 'event_subtype', 'confidence', 'confidence_reasons', 'candidate_model_id', 'supporting_model_ids', 'reference_model_ids', 'annotation_sources', 'chrom', 'start', 'end', 'strand', 'candidate_cds_blocks', 'reference_cds_blocks', 'relationship', 'peptides', 'mapping_ids', 'distinct_diagnostic_peptides', 'mapping_multiplicity', 'rna_junctions', 'minimum_rna_support', 'reasons', 'limitations', 'coding_structure_id', 'local_peptide_blocks', 'diagnostic_junctions', 'splice_chain_changed', 'whole_model_status', 'whole_model_review_reasons', 'uniqueness_scope', 'evidence_origin', 'boundary_status', 'cleavage_context', 'upstream_peptide_support', 'proposed_feature_blocks', 'peptide_loci', 'decoy_estimated_fdr']


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


def frame_at(model,position):
    """Codon phase of a genomic position within a model's CDS, or None."""
    offset=0
    for a,b in model.coding:
        if a<=position<b:return (offset+(position-a if model.strand=='+' else b-1-position))%3
        offset+=b-a
    return None


def first_coding_base(blocks,strand):
    """The peptide's first nucleotide in translation direction."""
    return blocks[0][0] if strand=='+' else blocks[-1][1]-1


def intron_containing(blocks,models):
    """Smallest annotated CDS intron (gap) of any model that contains all blocks."""
    left=blocks[0][0];right=blocks[-1][1];best=None
    for m in models:
        cds=sorted(m.cds)
        for (a,b),(c,d) in zip(cds,cds[1:]):
            if b<=left and right<=c and (best is None or c-b<best[1]-best[0]):best=(b,c)
    return best


def supporting(event):
    ids=event.get('supporting_model_ids') or []
    if isinstance(ids,str):ids=json.loads(ids) if ids else []
    return [i for i in ids if i]+([event['candidate_model_id']] if event.get('candidate_model_id') else [])


def relationship(candidate,references):
    if not references:return 'intergenic'
    same=[m for m in references if m.strand==candidate.strand]
    if not same:return 'antisense'
    overlapping=[m for m in same if interval_overlap(candidate.cds,m.cds)]
    if not overlapping:
        if any(interval_overlap(candidate.cds,m.exons) for m in same):return 'UTR_associated'
        return 'intronic'
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


def upstream_residue(reference, model, table=1):
    """In-frame residue encoded immediately before the model's first codon (contiguous genome)."""
    if reference is None or not model.coding:return None
    from peplocus.sequence import reverse_complement, translate
    a,b=model.coding[0]
    if model.strand=='+':
        if a<3:return None
        codon=reference.fetch(model.chrom,a-3,a)
    else:
        if b+3>reference.get_reference_length(model.chrom):return None
        codon=reverse_complement(reference.fetch(model.chrom,b,b+3))
    return translate(codon.upper(),table) or None


def assess(args,models,store,metadata,rna_junctions,guided_events=None,reference=None):
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
            non_nested_discriminating_peptides=len(non_nested(unique)),peptides=sorted({h['peptide'] for h in placements}),mapping_ids=sorted({h['mapping_id'] for h in placements}),
            complete_start=model.complete_start,complete_stop=model.complete_stop,issues=sorted(set(model.issues)),quality=quality(sequences),
            reasons='no_observed_evidence_is_not_evidence_of_incorrect_annotation' if not sequences else 'placement_multiplicity_within_searched_genome_and_splice_paths; I/L aliases counted once'))
    annotation=[m for m in models if m.source!='RNA_ORF'];candidates=[m for m in models if m.source=='RNA_ORF']
    supported_candidates=[m for m in candidates if hits[m.id]]
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
            sequences={canonical(p) for row in rows for p in row['peptides']}
            discriminating=sum(counts.get(p,0)==1 for p in sequences)
            # Gene support follows the model rule: at least one single-locus peptide.
            gene_rows.append(dict(source=source,chrom=chrom,strand=strand,gene_id=gene,models=len(rows),distinct_peptides=len(sequences),locus_discriminating_peptides=discriminating,category='peptide_supported' if discriminating else 'ambiguous' if sequences else 'no_observed_peptide_support'))
        write_table(directory/'gene_evidence.tsv',['source','chrom','strand','gene_id','models','distinct_peptides','locus_discriminating_peptides','category'],gene_rows)
    if args.bam:
        directory=args.output/'orfs'
        rna_records=[r for r in records if r['source']=='RNA_ORF']
        write_table(directory/'all_predicted_orf_evidence.tsv',MODEL_FIELDS,rna_records)
        write_table(directory/'orf_peptide_evidence.tsv',MODEL_FIELDS,(r for r in rna_records if r['distinct_peptides']>0))
        write_table(directory/'orfs_without_peptide_support.tsv',MODEL_FIELDS,(r for r in rna_records if r['distinct_peptides']==0))
        write_models(supported_candidates,directory/'candidate_orfs.gtf')
        write_proteins(supported_candidates,directory/'candidate_orfs.faa')
        write_models(candidates,directory/'all_predicted_orfs.gtf')
        write_proteins(candidates,directory/'all_predicted_orfs.faa')
        write_table(directory/'partial_or_ambiguous_orfs.tsv',MODEL_FIELDS,(r for r in rna_records if r['distinct_peptides']>0 and (not r['complete_start'] or not r['complete_stop'] or r['category']=='ambiguous')))
        (directory/'README.txt').write_text('candidate_orfs.gtf, candidate_orfs.faa and orf_peptide_evidence.tsv contain only RNA-derived ORFs with at least one exact/I-L-equivalent peptide placement. Shared or ambiguous support is retained and labeled; peptide support does not establish identification quality, isoform specificity, novelty or full-model correctness. all_predicted_orfs.gtf, all_predicted_orfs.faa and all_predicted_orf_evidence.tsv retain the complete reconstruction for audit. orfs_without_peptide_support.tsv lists predictions without observed peptide support; absence of support does not establish that an ORF is incorrect. partial_or_ambiguous_orfs.tsv contains only peptide-supported candidates.\n')
    # Diagnostic = not explained by any supplied annotation at the same locus
    # (any placement sharing nucleotides), whatever the enzyme policy.
    strict=store.counts(all_mappings=True)
    qual=quality_lookup(metadata)
    capped=getattr(args,'splice_search_capped',False)
    reference_index=ModelIndex(annotation)
    annotated=ClusterIndex()
    annotation_hits=store.model_hits(all_mappings=True,model_ids={m.id for m in annotation}) if annotation and candidates else {}
    for m in annotation:
        for hit in annotation_hits.get(m.id,[]):annotated.add(hit['canonical'],m.chrom,m.strand,hit['blocks'],explaining=True)
    if annotation and candidates:
        # Every exact placement links locus clusters (e.g. duplicated exons).
        candidate_sequences={h['canonical'] for c in candidates for h in hits[c.id]}
        for c in candidates:
            for h in hits[c.id]:annotated.add(h['canonical'],c.chrom,c.strand,h['blocks'])
        if hasattr(store,'placements'):
            for sequence,items in store.placements(candidate_sequences).items():
                for chrom,strand,blocks in items:annotated.add(sequence,chrom,strand,blocks)
    events=[];isoform_checks=[];unassigned=Counter()
    annotated_paths=defaultdict(list)
    annotated_junctions_by_locus=defaultdict(set)
    for model in annotation:
        annotated_junctions_by_locus[model.chrom,model.strand].update(model.junctions)
        annotated_junctions_by_locus[model.chrom,model.strand].update(
            (a[1],b[0]) for a,b in zip(model.exons,model.exons[1:]) if a[1]<b[0])
    for model in annotation:
        if model.coding:annotated_paths[model.signature()].append(model.id)
    for candidate in candidates:
        # Without a supplied annotation, sequence support cannot establish novelty.
        # Keep peptide-supported RNA ORFs in candidate_orfs and their evidence table.
        if not hits[candidate.id] or not annotation:continue
        refs=reference_index.overlap(candidate)
        relation=relationship(candidate,refs)
        all_hits=hits[candidate.id]
        diagnostic=[h for h in all_hits if not annotated.explains(h['canonical'],candidate.chrom,candidate.strand,h['blocks'])]
        same=[m for m in refs if m.coding and relationship(candidate,[m])=='new_transcript_isoform']
        known_starts=[m.id for m in same if boundary(m,'start')==boundary(candidate,'start')]
        known_ends=[m.id for m in same if boundary(m,'end')==boundary(candidate,'end')]
        # A known splice boundary is not new merely because its annotation
        # has a different reading frame, gene ID or comparison-model subset.
        annotated_junctions=annotated_junctions_by_locus[candidate.chrom,candidate.strand]
        new_junctions=set(candidate.junctions)-annotated_junctions
        exact=sorted(annotated_paths.get(candidate.signature(),[]))
        isoform_checks.append(dict(candidate_model_id=candidate.id,chrom=candidate.chrom,strand=candidate.strand,coding_structure_id=stable_id('CDS',candidate.signature()),exact_annotated_isoform_ids=exact,comparable_annotated_isoform_ids=sorted(m.id for m in same),known_start_isoform_ids=sorted(known_starts),known_end_isoform_ids=sorted(known_ends),new_junctions=sorted(new_junctions),observed_peptide_sequences=len({h['canonical'] for h in all_hits}),diagnostic_peptide_sequences=len({h['canonical'] for h in diagnostic}),interpretation='exact_annotated_isoform' if exact else 'existing_isoforms_explain_all_observed_peptides' if not diagnostic else 'diagnostic_evidence_requires_feature_specific_review'))
        if exact or not diagnostic:continue
        rna_keys=[(candidate.chrom,a,b,candidate.strand) for a,b in candidate.junctions]
        minimum=min((rna_junctions.get(k,0) for k in rna_keys),default=0)
        total=sum(b-a for a,b in candidate.coding)
        # Every diagnostic placement is assigned to exactly one feature.
        features=defaultdict(list)
        if relation!='new_transcript_isoform':
            reason='intergenic_peptide_supported_ORF' if relation=='intergenic' else 'diagnostic_peptides_outside_annotated_coding_frame'
            features[('novel_orf','',('orf',),reason)]=list(diagnostic)
        else:
            for h in diagnostic:
                blocks=[tuple(b) for b in h['blocks']]
                crossed=sorted({(a[1],b[0]) for a,b in zip(blocks,blocks[1:]) if (a[1],b[0]) in new_junctions})
                first=first_coding_base(blocks,candidate.strand)
                ref_frames=[f for f in (frame_at(m,first) for m in same) if f is not None]
                n_terminal=h['aa_start']==0 or (h['aa_start']==1 and candidate.protein.startswith('M') and candidate.protein[1:2] in 'ACGPSTV')
                if crossed:
                    features[('alternative_splice','new_junction',('junctions',tuple(crossed)),'peptide_crosses_unannotated_junction')].append(h)
                elif all(beyond(h,m,'start') for m in same):
                    if candidate.complete_start and not known_starts and n_terminal:
                        codon=tuple(tuple(b) for b in candidate.project(0,3))
                        features[('alternative_start','',('codon',codon),'ORF_N_terminal_peptide')].append(h)
                    else:
                        features[('n_terminal_extension','',('start',boundary(candidate,'start'),'ORF_start_not_directly_observed'),'ORF_start_not_directly_observed')].append(h)
                elif all(beyond(h,m,'end') for m in same):
                    if candidate.complete_stop and not known_ends:
                        codon=tuple(tuple(b) for b in candidate.project(total-3,total))
                        features[('alternative_stop','',('codon',codon),'diagnostic_peptide_extends_beyond_annotated_coding_end')].append(h)
                    else:
                        features[('c_terminal_extension','',('end',boundary(candidate,'end'),'boundary_unresolved'),'diagnostic_peptide_extends_beyond_annotated_coding_end')].append(h)
                elif not any(interval_overlap(blocks,m.cds) for m in same):
                    host=intron_containing(blocks,same)
                    if host is None:
                        unassigned['outside_comparable_CDS_span']+=1;continue
                    exon=[b for b in candidate.coding if interval_overlap(blocks,[b])]
                    subtype='retained_intron' if any(interval_overlap(exon,m.cds) for m in same) else 'novel_exon'
                    features[('alternative_splice',subtype,('intron',host,subtype),'diagnostic_peptide_in_annotated_intron')].append(h)
                elif ref_frames and all(f!=frame_at(candidate,first) for f in ref_frames):
                    # Comparable CDSs cover the peptide's first base, all in another frame.
                    features[('alternative_splice','frame_change',('region',tuple(blocks)),'diagnostic_peptide_in_alternative_reading_frame')].append(h)
                else:
                    unassigned['same_frame_within_annotated_CDS']+=1
        for (kind,subtype,key,reason),placements in features.items():
            sequences={h['canonical'] for h in placements}
            blockers=[]
            if candidate.issues:
                blockers.append('noncanonical_intron_in_candidate' if 'noncanonical_intron' in candidate.issues else 'candidate_sequence_or_phase_issues')
            if not candidate.protein:blockers.append('candidate_sequence_or_phase_issues')
            if kind=='alternative_start' and any(h['aa_start']==0 for h in placements):
                # An ORF-initial M peptide preceded in frame by K/R may be a tryptic
                # product of a longer protein, not a protein N-terminus.
                previous=upstream_residue(reference,candidate,getattr(args,'genetic_code',1))
                if previous in {'K','R'}:blockers.append('trypsin_site_before_proposed_start')
            junction_key=key[1] if key[0]=='junctions' else ()
            host=key[1] if key[0]=='intron' else None
            consistent={h['canonical'] for h in placements if all(
                feature_consistent(kind,subtype,candidate.strand,x['blocks'],annotated.members(x['canonical'],candidate.chrom,candidate.strand,x['blocks']),junction_key,host)
                for x in placements if x['canonical']==h['canonical'])}
            confidence,specific,reasons=judge(sequences,strict,qual,blockers,consistent)
            crossed=sorted({(a[1],b[0]) for h in placements for a,b in zip(h['blocks'],h['blocks'][1:]) if (a[1],b[0]) in new_junctions})
            event_rna=min((rna_junctions.get((candidate.chrom,a,b,candidate.strand),0) for a,b in crossed),default=minimum) if crossed else minimum
            if confidence=='peptide_supported':
                high=(len(non_nested(specific))>=2 and candidate.complete_start and candidate.complete_stop and quality(specific)=='qvalue_pass'
                      and set(candidate.protein)<=set('ACDEFGHIKLMNPQRSTVWY') and (not rna_keys or minimum>=args.min_junction_support) and not capped)
                if kind=='alternative_splice' and event_rna<args.min_junction_support:high=False
                if high:confidence='high_confidence_candidate'
                elif (crossed or rna_keys) and event_rna>=args.min_junction_support:confidence='rna_peptide_supported'
            limitations=['candidate_annotation_change_requires_review','peptide_identification_quality_'+quality(sequences),'RNA_junctions_do_not_prove_full_length_connectivity']
            if any(strict.get(p,0)>1 for p in sequences):limitations.append('non_discriminating_peptides_present')
            if len(non_nested(specific))<2:limitations.append('fewer_than_two_non_nested_locus_discriminating_peptides')
            if set(candidate.protein)-set('ACDEFGHIKLMNPQRSTVWY*'):limitations.append('unresolved_amino_acids_in_candidate_ORF')
            if candidate.issues:limitations.extend(candidate.issues)
            if capped:limitations.append('splice_path_search_capped_uniqueness_incomplete')
            if kind=='alternative_start':
                limitations.extend(['proteolytic_processing_not_excluded','initiation_site_not_experimentally_confirmed'])
                if reference is None:limitations.append('upstream_cleavage_context_not_checked')
                else:limitations.append('upstream_cleavage_context_from_contiguous_genome')
            if kind=='alternative_stop':limitations.append('stop_position_inferred_from_ORF_sequence')
            if kind=='c_terminal_extension':limitations.append('termination_site_unresolved_no_complete_stop_codon')
            if kind=='n_terminal_extension':limitations.append('initiation_site_unresolved_ORF_start_not_observed')
            limitations.append('uniqueness_is_relative_to_searched_splice_paths_not_all_possible_splicing')
            related=same if kind!='novel_orf' else refs
            review_reasons=[]
            if not candidate.complete_start:review_reasons.append('incomplete_start')
            if not candidate.complete_stop:review_reasons.append('incomplete_stop')
            if candidate.issues:review_reasons.append('candidate_sequence_or_phase_issues')
            if confidence not in SUPPORTED:review_reasons.append('diagnostic_evidence_not_specific')
            if any(rna_junctions.get(k,0)<args.min_junction_support for k in rna_keys):review_reasons.append('coding_junction_without_required_RNA_support')
            if len(non_nested(specific))<2:review_reasons.append('fewer_than_two_independent_diagnostic_peptides')
            genes={(m.source,m.gene) for m in related if m.strand==candidate.strand and interval_overlap(candidate.cds,m.cds)}
            if len(genes)>1:review_reasons.append('overlaps_multiple_annotation_gene_ids_review_fusion_or_isoforms')
            if any(any(sum(max(0,min(b,d)-max(a,c)) for c,d in candidate.coding)<b-a for a,b in m.coding) for m in related if m.coding):review_reasons.append('compared_annotation_CDS_not_fully_preserved_local_support_does_not_validate_deletions')
            changed=any(set(candidate.junctions)!=set(m.junctions) for m in related if m.coding)
            if kind in {'alternative_stop','c_terminal_extension'} and changed:limitations.append('terminal_change_accompanied_by_splice_chain_change')
            # This is a context model, never a whole-gene replacement inferred from local evidence.
            status='local_evidence_only' if review_reasons else 'complete_model_for_manual_review'
            feature=[list(b) for b in key[1]] if key[0]=='codon' else []
            events.append(dict(event_id=feature_id(kind,candidate.chrom,candidate.strand,key if kind!='novel_orf' else ('orf',candidate.id)),event_type=kind,event_subtype=subtype,
                candidate_model_id=candidate.id,supporting_model_ids=[candidate.id],
                reference_model_ids=sorted(m.id for m in related),annotation_sources=sorted({m.source for m in related}),chrom=candidate.chrom,
                start=min(a for a,b in candidate.cds),end=max(b for a,b in candidate.cds),strand=candidate.strand,
                candidate_cds_blocks=sorted(candidate.coding),reference_cds_blocks={m.id:sorted(m.coding) for m in related},relationship=relation,
                peptides=sorted({h['peptide'] for h in placements}),mapping_ids=sorted({h['mapping_id'] for h in placements}),distinct_diagnostic_peptides=len(sequences),
                mapping_multiplicity={h['peptide']:strict.get(h['canonical'],0) for h in sorted(placements,key=lambda h:h['peptide'])},
                peptide_loci={h['peptide']:dict(loci_all_exact=strict.get(h['canonical'],0),loci_analysis=counts.get(h['canonical'],0),specific=h['canonical'] in specific,length=len(h['peptide'])) for h in placements},
                rna_junctions=[dict(chrom=k[0],start=k[1],end=k[2],strand=k[3],reads=rna_junctions.get(k,0)) for k in rna_keys],
                minimum_rna_support=event_rna,confidence=confidence,confidence_reasons=reasons,coding_structure_id=stable_id('CDS',candidate.signature()),
                local_peptide_blocks=sorted({tuple(b) for h in placements for b in h['blocks']}),diagnostic_junctions=crossed,splice_chain_changed=changed,
                whole_model_status=status,whole_model_review_reasons=review_reasons,uniqueness_scope='genomic_loci_all_exact_placements_and_searched_paths',
                reasons=reason,limitations=limitations,proposed_feature_blocks=feature,boundary_status=reason,evidence_origin='RNA_ORF',decoy_estimated_fdr=''))
    # Novel ORFs from different RNA isoforms that share any diagnostic placement are one novel locus.
    novel=[e for e in events if e['event_type']=='novel_orf']
    parent={e['event_id']:e['event_id'] for e in novel}
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    owner={}
    for e in novel:
        for mapping_id in e['mapping_ids']:
            if mapping_id in owner:parent[find(e['event_id'])]=find(owner[mapping_id])
            else:owner[mapping_id]=e['event_id']
    roots=defaultdict(list)
    for e in novel:roots[find(e['event_id'])].append(e)
    for members in roots.values():
        anchor=min(m for e in members for m in e['mapping_ids'])
        for e in members:e['event_id']=feature_id('novel_orf',e['chrom'],e['strand'],('placement',anchor))
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
    rna_event_count=len(events)
    # One row per biological feature: isoforms, sources and evidence origins
    # that support the same start, stop, junction or locus are merged.
    events=merge_events(events+list(guided_events or []))
    events.sort(key=lambda e:(e['chrom'],e['start'],e['strand'],e['event_type'],e['event_id']))
    if args.bam:
        directory=args.output/'orfs'
        novel=[e for e in events if e['event_type']=='novel_orf']
        write_table(directory/'novel_orf_candidates.tsv',EVENT_FIELDS,novel)
        novel_models=[by_id[i] for i in sorted({i for e in novel for i in supporting(e)}) if i in by_id]
        write_models(novel_models,directory/'novel_orf_candidates.gtf');write_proteins(novel_models,directory/'novel_orf_candidates.faa')
    if args.gtf:
        directory=args.output/'refinement';directory.mkdir(exist_ok=True)
        for kind in ['alternative_start','alternative_stop','alternative_splice','n_terminal_extension','c_terminal_extension']:
            selected=[e for e in events if e['event_type']==kind]
            write_table(directory/(kind+'_candidates.tsv'),EVENT_FIELDS,selected)
            context_models=[by_id[i] for i in sorted({i for e in selected for i in supporting(e)}) if i in by_id]
            context_path=directory/(kind+'_rna_model_context.gtf')
            if context_models:write_models(context_models,context_path)
            else:context_path.unlink(missing_ok=True)
        write_table(directory/'annotation_change_proposals.tsv',EVENT_FIELDS,events)
        write_table(directory/'local_evidence_only.tsv',EVENT_FIELDS,(e for e in events if e['whole_model_status']=='local_evidence_only'))
        write_table(directory/'complete_models_for_manual_review.tsv',EVENT_FIELDS,(e for e in events if e['whole_model_status']=='complete_model_for_manual_review'))
        with (directory/'local_peptide_evidence.bed').open('w') as stream:
            for e in events:
                for i,(a,b) in enumerate(e['local_peptide_blocks'],1):
                    stream.write(f"{e['chrom']}\t{a}\t{b}\t{e['event_id']}.{i}\t0\t{e['strand']}\n")
        (directory/'README.txt').write_text('Annotation-guided events have local features in annotation_guided_features.gtf; their candidate_model_id and full-CDS blocks are empty because no full transcript was inferred. RNA and annotation-guided rows can support the same local feature and are not counts of distinct genes. Files ending _rna_model_context.gtf contain full RNA models for context and manual review only. They are not corrected annotations or automatic replacements. Local peptide blocks and diagnostic junctions in the event tables delimit the direct evidence. local_peptide_evidence.bed contains separate supported blocks (BED6), not reconstructed CDSs. Complete models still require independent full-structure validation. IDs are run-scoped; join with the run/crop. Coding structure IDs collapse duplicate CDS paths within an assembly; retain run identity across assemblies.\n')
    # Traceable individual protein occurrences, including N-terminal cleavage context.
    occurrence_rows=[]
    for model in models:
        for h in hits[model.id]:
            previous=model.protein[h['aa_start']-1] if h['aa_start'] else ''
            first=model.protein[h['aa_start']:h['aa_start']+1]
            occurrence_rows.append(dict(model_id=model.id,mapping_id=h['mapping_id'],peptide=h['peptide'],aa_start=h['aa_start'],aa_end=h['aa_end'],previous_aa=previous,n_terminal_context='protein_N_terminus' if not previous else 'tryptic' if previous in 'KR' and first!='P' else 'nontryptic',initiator_methionine_processing=h['aa_start']==1 and model.protein.startswith('M') and first in 'ACGPSTV'))
    write_table(args.output/'mapping/model_peptide_occurrences.tsv',['model_id','mapping_id','peptide','aa_start','aa_end','previous_aa','n_terminal_context','initiator_methionine_processing'],occurrence_rows)
    return dict(novelty_assessed=bool(annotation),rna_event_rows=rna_event_count,annotation_guided_event_rows=len(guided_events or []),isoform_checks=dict(Counter(row['interpretation'] for row in isoform_checks)),distinct_candidate_cds_paths=len({m.signature() for m in supported_candidates if m.coding}),distinct_predicted_cds_paths=len({m.signature() for m in candidates if m.coding}),whole_model_statuses=dict(Counter(e['whole_model_status'] for e in events)),models=len(models),annotation_models=len(annotation),predicted_orfs=len(candidates),candidate_orfs=len(supported_candidates),orfs_without_peptide_support=len(candidates)-len(supported_candidates),model_categories=dict(Counter(categories.values())),candidate_categories=dict(Counter(r['category'] for r in primary)),annotation_events=dict(Counter(e['event_type'] for e in events)),event_confidences=dict(Counter(e['confidence'] for e in events)),
                supported_events=dict(Counter(e['event_type'] for e in events if e['confidence'] in SUPPORTED)),unassigned_diagnostic_placements=dict(unassigned))
