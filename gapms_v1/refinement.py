"""Bounded annotation-guided peptide evidence; no RNA reconstruction required.

Hypotheses are local coding paths, not inferred full-length transcripts. All
local hypotheses are recorded as possible placements before multiplicity is
evaluated, so a second possible origin always revokes uniqueness. Events are
only proposed for peptides that no supplied annotation explains anywhere, and
are judged by the shared rules in :mod:`gapms_v1.events`. A reversed-peptide
decoy set is searched through the identical local hypothesis space to estimate
the false-event rate of each event type.
"""
from collections import Counter, defaultdict
from itertools import product
import json
import time

import ahocorasick

from peplocus.digestion import METHIONINE_REMOVAL, cleavage_context
from peplocus.sequence import canonical, join_blocks, translate, project
from .assessment import EVENT_FIELDS, boundary
from .common import LOG, stable_id, write_table
from .events import SUPPORTED, feature_consistent, feature_id, judge, quality_lookup
from .mapping import ClusterIndex, LocusIndex, locus_key

STANDARD = set('ACDEFGHIKLMNPQRSTVWY')


def ordered(blocks, strand):
    blocks=json.loads(locus_key('',strand,blocks)[2])
    return [tuple(b) for b in (blocks if strand=='+' else reversed(blocks))]


def sliced(model, start, end):
    return ordered(project(start,end,model.coding,model.strand),model.strand)


def cleavage(protein, start, end, native_start=False, native_end=False):
    previous=protein[start-1] if start else ''
    following=protein[end:end+1]
    n_ok=(start==0 and native_start) or bool(previous in {'K','R'} and protein[start:start+1]!='P')
    c_ok=(native_end and (end==len(protein) or following=='*')) or bool(protein[end-1:end] in {'K','R'} and following in STANDARD-{'P'})
    return dict(previous_aa=previous,following_aa=following,n_terminus='protein_boundary' if start==0 and native_start else 'tryptic' if n_ok else 'nontryptic' if previous in STANDARD else 'upstream_stop_boundary' if previous=='*' else 'unknown',c_terminus='protein_boundary' if native_end and (end==len(protein) or following=='*') else 'tryptic' if c_ok else 'nontryptic' if following else 'unknown',tryptic_n=n_ok,tryptic_c=c_ok)


LENGTH_BINS={'<=7':(0,7),'8':(8,8),'9-10':(9,10),'11-14':(11,14),'15+':(15,10**6)}


def length_bin(length):
    return next(label for label,(low,high) in LENGTH_BINS.items() if low<=length<=high)


def decoy(sequence):
    """Reverse all but the C-terminal residue (keeps tryptic C-termini)."""
    return sequence[-2::-1]+sequence[-1]


def codon_point(codon, strand, which):
    """Genomic boundary of a start (5') or stop (3') codon."""
    if which=='start':return codon[0][0] if strand=='+' else codon[-1][1]
    return codon[-1][1] if strand=='+' else codon[0][0]


def start_background(annotations, associations):
    """Run-level control: is an M-initial non-tryptic peptide more frequent than chance?

    Among internal peptides whose N-terminus is not trypsin-compatible, the
    fraction beginning with Met (or following Met with an iMet-removal residue)
    is compared with the proteome Met frequency. Without enrichment, internal
    M-initial peptides are indistinguishable from random non-tryptic fragments.
    """
    proteins={m.id:m.protein for m in annotations if m.protein}
    residues=Counter()
    for protein in set(proteins.values()):residues.update(protein)
    total=sum(residues.values())
    met=residues['M']/total if total else 0.0
    seen=set();internal=0;m_initial=0;after_m=0;removal_residue=0
    for mid,peptide,start in associations:
        protein=proteins.get(mid)
        if not protein or start<2 or peptide in seen:continue
        previous=protein[start-1]
        if previous in 'KR*' or previous not in STANDARD:continue
        seen.add(peptide);internal+=1
        if peptide[0]=='M':m_initial+=1
        if peptide[0] in METHIONINE_REMOVAL:
            removal_residue+=1
            if previous=='M':after_m+=1
    m_fraction=m_initial/internal if internal else None
    removal_fraction=after_m/removal_residue if removal_residue else None
    return dict(internal_nontryptic_peptides=internal,proteome_met_frequency=round(met,4),
                m_initial_fraction=None if m_fraction is None else round(m_fraction,4),
                m_initial_enrichment=None if m_fraction is None or not met else round(m_fraction/met,2),
                methionine_before_removal_residue_fraction=None if removal_fraction is None else round(removal_fraction,4),
                methionine_before_removal_residue_enrichment=None if removal_fraction is None or not met else round(removal_fraction/met,2),
                interpretation='Enrichment near 1 means internal start candidates are not distinguishable from random non-tryptic fragments in this dataset.')


def refine_annotations(args, annotations, store, metadata, reference, rna_junctions):
    started=time.monotonic()
    annotations=[m for m in annotations if m.source not in {'RNA_ORF','consensus_graph'} and m.coding]
    stats=Counter();pending=[];decoy_pending=[]
    radius=getattr(args,'refinement_window',3000)
    splice_radius=getattr(args,'refinement_splice_window',30)
    cap=getattr(args,'max_refinement_paths',64)
    anchor=getattr(args,'min_refinement_anchor',9)
    max_peptide=max(map(len,metadata))
    max_nt=3*(max_peptide+1)
    lookup={m.id:m for m in annotations}
    known_paths=set();annotated=LocusIndex();ordinary_n=LocusIndex();internal=[];occurrences=defaultdict(list);all_assoc=[]
    aliases=defaultdict(list)
    for p in metadata:
        if len(p)>=args.min_length:aliases[canonical(p)].append(p)
    # One streaming pass; RNA/consensus associations do not explain annotation novelty.
    for mid,peptide,sequence,chrom,strand,encoded,start,end in store.db.execute('SELECT * FROM associations'):
        if mid not in lookup:continue
        blocks=json.loads(encoded)
        known_paths.add((sequence,chrom,strand,encoded));annotated.add(sequence,chrom,strand,blocks)
        m=lookup[mid];protein=m.protein
        occurrences[mid].append((start,sequence));all_assoc.append((mid,peptide,start))
        context=cleavage(protein,start,end,True,m.complete_stop)
        n_terminal=start==0 or (start==1 and protein.startswith('M') and protein[1:2] in METHIONINE_REMOVAL)
        if n_terminal or context['tryptic_n']:ordinary_n.add(sequence,chrom,strand,blocks)
        if m.issues or start<1:continue
        # A start is either the M that begins the peptide, or the M removed before it.
        if protein[start]=='M' and protein[start-1] not in 'KR':
            initiator,mode=start,'internal_ATG_nontryptic_N_terminus'
        elif start>=2 and protein[start-1]=='M' and protein[start] in METHIONINE_REMOVAL:
            initiator,mode=start-1,'internal_ATG_initiator_methionine_removed'
        else:continue
        if protein[initiator-1:initiator] not in STANDARD and initiator>0:continue  # unknown residue before the start
        if store.require_tryptic and not context['tryptic_c']:continue
        codon=m.project(initiator*3,initiator*3+3)
        if join_blocks(reference,m.chrom,m.strand,codon).upper()!='ATG':continue
        internal.append(dict(model=m,peptide=peptide,sequence=sequence,blocks=blocks,start=start,initiator=initiator,mode=mode,context=context,codon=codon))
    explained=set(annotated.items)
    eligible=set(aliases)-explained
    stats['unexplained_query_sequences']=len(eligible)
    automaton=None
    if aliases:
        automaton=ahocorasick.Automaton()
        for sequence in sorted(aliases):automaton.add_word(sequence,(sequence,False))
        decoys=0
        for sequence in sorted(eligible):
            d=decoy(sequence)
            if d!=sequence and d not in aliases and not automaton.exists(d):
                automaton.add_word(d,(d,True));decoys+=1
        stats['decoy_sequences']=decoys
        automaton.make_automaton()
    known_boundaries={which:{(m.chrom,m.strand,boundary(m,which)) for m in annotations} for which in ('start','end')}
    def is_known_boundary(model, point, which):
        return (model.chrom,model.strand,point) in known_boundaries[which]

    def add(kind,model,peptide,blocks,context,key,status,feature=(),edges=(),subtype='',upstream=(),is_decoy=False):
        chrom,strand,encoded=locus_key(model.chrom,model.strand,blocks)
        sequence=canonical(peptide)
        compatible=cleavage_context(peptide,context.get('previous_aa',''),context.get('following_aa',''))['enzyme_compatible']
        if store.require_tryptic and not compatible:
            # A rejected hypothesis is not a genomic placement: keep it out of
            # every mapping track and multiplicity count.
            if not is_decoy:stats['nontryptic_local_matches_excluded']+=1
            return
        item=dict(kind=kind,model=model,peptide=peptide,sequence=sequence,blocks=json.loads(encoded),encoded=encoded,context=context,
                  key=key,status=status,feature=[list(b) for b in feature],edges=[list(e) for e in edges],subtype=subtype,upstream=list(upstream))
        if is_decoy:
            decoy_pending.append(item);return
        if kind!='alternative_start' and (sequence,chrom,strand,encoded) in known_paths:return
        matched=translate(join_blocks(reference,chrom,strand,json.loads(encoded)),args.genetic_code)
        if canonical(matched)!=sequence:raise RuntimeError('annotation-guided projection verification failed')
        # Every hypothesis is a possible origin and enters multiplicity.
        store.add(peptide,chrom,strand,json.loads(encoded),matched,'annotation_guided',tryptic=compatible)
        # Only peptides that no annotation explains anywhere become events.
        if sequence in eligible:pending.append(item)

    def emit(kind,model,sequence,is_decoy,blocks,context,key,status,**extra):
        if is_decoy:
            add(kind,model,sequence,blocks,context,key,status,is_decoy=True,**extra)
        else:
            for peptide in aliases[sequence]:add(kind,model,peptide,blocks,context,key,status,**extra)

    def scan(model, segments, kind, old_offset):
        if automaton is None:return
        protein=translate(join_blocks(reference,model.chrom,model.strand,segments),args.genetic_code)
        if kind=='n_terminal_extension':
            # Only the stop-free path connected to the existing coding start.
            stop=protein.rfind('*',0,old_offset//3)
            valid_start=stop+1
            stop_after=protein.find('*',old_offset//3)
            valid_end=stop_after if stop_after>=0 else len(protein)
        else:
            valid_start=0;stop_after=protein.find('*');valid_end=stop_after if stop_after>=0 else len(protein)
        matches=[]
        for end,(sequence,is_decoy) in automaton.iter(canonical(protein)):
            end+=1;start=end-len(sequence)
            if start<valid_start or end>valid_end:continue
            if kind=='n_terminal_extension' and (start*3>=old_offset or 'X' in protein[start:old_offset//3]):continue
            if kind=='c_terminal_extension' and (end*3<=old_offset or 'X' in protein[old_offset//3:end]):continue
            matches.append((start,end,sequence,is_decoy))
        for start,end,sequence,is_decoy in matches:
            blocks=project(start*3,end*3,segments,model.strand)
            context=cleavage(protein,start,end,False,stop_after==end)
            if kind=='n_terminal_extension':
                previous=protein[start-1] if start else ''
                initiator=None
                if protein[start:start+1]=='M' and previous not in {'K','R'} and start>0:initiator,status=start,'ATG_supported_start_candidate'
                elif previous=='M' and protein[start:start+1] in METHIONINE_REMOVAL:initiator,status=start-1,'ATG_supported_start_initiator_methionine_removed'
                if initiator is not None:
                    codon=project(initiator*3,initiator*3+3,segments,model.strand)
                    if join_blocks(reference,model.chrom,model.strand,codon).upper()=='ATG' and not is_known_boundary(model,codon_point(codon,model.strand,'start'),'start'):
                        # Any other observed peptide further upstream in this
                        # stop-free window contradicts initiation at this M.
                        upstream=sorted({s for a,_,s,d in matches if not d and a<initiator})
                        emit('alternative_start',model,sequence,is_decoy,blocks,context,('codon',tuple(map(tuple,codon))),status,feature=codon,upstream=upstream)
                        continue
                if previous=='*':
                    emit('n_terminal_extension',model,sequence,is_decoy,blocks,context,('start',boundary(model,'start'),'upstream_stop_without_start'),'upstream_stop_without_start')
                else:
                    emit('n_terminal_extension',model,sequence,is_decoy,blocks,context,('start',boundary(model,'start'),'boundary_unresolved'),'boundary_unresolved')
            else:
                if stop_after>=0:
                    codon=project(stop_after*3,(stop_after+1)*3,segments,model.strand)
                    if not is_known_boundary(model,codon_point(codon,model.strand,'end'),'end'):
                        emit('alternative_stop',model,sequence,is_decoy,blocks,context,('codon',tuple(map(tuple,codon))),'genomic_stop_inferred',feature=codon)
                        continue
                    if not is_decoy:stats['extensions_with_already_annotated_stop']+=1
                    status='extension_with_already_annotated_end'
                else:status='boundary_unresolved'
                emit('c_terminal_extension',model,sequence,is_decoy,blocks,context,('end',boundary(model,'end'),status),status)

    # Downstream starts are a context test of existing placements, not a new path.
    for item in internal:
        m=item['model'];point=codon_point(item['codon'],m.strand,'start')
        if is_known_boundary(m,point,'start'):
            stats['internal_starts_already_annotated']+=1;continue
        if ordinary_n.explains(item['sequence'],m.chrom,m.strand,item['blocks']):
            stats['start_matches_explained_by_another_isoform_cleavage']+=1;continue
        status=item['mode'] if item['initiator']>1 else 'initiator_methionine_ambiguity'
        upstream=sorted({s for offset,s in occurrences[m.id] if offset<item['initiator']})
        record=dict(kind='alternative_start',model=m,peptide=item['peptide'],sequence=item['sequence'],blocks=item['blocks'],
                    encoded=locus_key(m.chrom,m.strand,item['blocks'])[2],context=item['context'],key=('codon',tuple(map(tuple,item['codon']))),
                    status=status,feature=[list(b) for b in item['codon']],edges=[],subtype='internal_start',upstream=upstream,internal=True)
        pending.append(record)
    representatives={}
    for m in annotations:
        if not m.issues and m.protein:representatives.setdefault(m.signature(),m)
    LOG.info('Annotation-guided refinement: %s coding paths, %s unexplained peptide sequences',len(representatives),len(eligible))
    for m in representatives.values():
        if automaton is None:break
        total=sum(b-a for a,b in m.coding)
        first=m.coding[0];last=m.coding[-1]
        extent=first[0] if m.strand=='+' else reference.get_reference_length(m.chrom)-first[1]
        extension=min(radius,extent)//3*3
        if extension:
            prefix=(first[0]-extension,first[0]) if m.strand=='+' else (first[1],first[1]+extension)
            segments=ordered([prefix]+sliced(m,0,min(total,max_nt)),m.strand)
            scan(m,segments,'n_terminal_extension',extension)
        if not m.complete_stop:
            extent=reference.get_reference_length(m.chrom)-last[1] if m.strand=='+' else last[0]
            extension=min(radius,extent)
            if extension:
                suffix=(last[1],last[1]+extension) if m.strand=='+' else (last[0]-extension,last[0])
                start=max(0,(total-max_nt)//3*3)
                segments=ordered(sliced(m,start,total)+[suffix],m.strand)
                scan(m,segments,'c_terminal_extension',total-start)
    known_edges={(m.chrom,m.strand,a,b) for m in annotations for a,b in m.junctions}
    seen=set()
    if automaton is not None and not getattr(args,'no_refinement_splice',False):
        for m in representatives.values():
            genomic=sorted(m.coding);offsets={};offset=0
            for block in m.coding:offsets[block]=offset;offset+=block[1]-block[0]
            for i,left in enumerate(genomic[:-1]):
                for skip in (0,1):
                    k=i+1+skip
                    if k>=len(genomic):continue
                    right=genomic[k];a,b=left;c,d=right
                    skipped=sum(v-u for u,v in genomic[i+1:k])
                    context_key=(m.chrom,m.strand,a,b,c,d,offsets[left if m.strand=='+' else right]%3,skipped%3)
                    if context_key in seen:continue
                    seen.add(context_key)
                    low_motif,high_motif=('GT','AG') if m.strand=='+' else ('CT','AC')
                    low_start=max(a+1,b-splice_radius);low_end=min(d-1,b+splice_radius)
                    high_start=max(a+1,c-splice_radius);high_end=min(d-1,c+splice_radius)
                    low_seq=reference.fetch(m.chrom,low_start,min(reference.get_reference_length(m.chrom),low_end+2)).upper()
                    high_base=max(0,high_start-2);high_seq=reference.fetch(m.chrom,high_base,high_end).upper()
                    donors=[x for x in range(low_start,low_end+1) if low_seq[x-low_start:x-low_start+2]==low_motif]
                    acceptors=[x for x in range(high_start,high_end+1) if high_seq[x-2-high_base:x-high_base]==high_motif]
                    variants=[(x,y) for x,y in product(donors,acceptors) if y-x>=20 and (x-b-(y-c)-skipped)%3==0 and (m.chrom,m.strand,x,y) not in known_edges]
                    variants.sort(key=lambda xy:(abs(xy[0]-b)+abs(xy[1]-c),xy))
                    if len(variants)>cap:stats['splice_contexts_capped']+=1
                    subtype='exon_skip' if skip else 'junction_shift'
                    for x,y in variants[:cap]:
                        stats['local_splice_paths_searched']+=1
                        if m.strand=='+':
                            start=max(a,x-max_nt);start+=(-(offsets[left]+start-a))%3
                            segments=[(start,x),(y,min(d,y+max_nt))]
                        else:
                            end=min(d,y+max_nt);end-=(-(offsets[right]+d-end))%3
                            segments=[(y,end),(max(a,x-max_nt),x)]
                        if any(u>=v for u,v in segments):continue
                        protein=translate(join_blocks(reference,m.chrom,m.strand,segments),args.genetic_code)
                        for end,(sequence,is_decoy) in automaton.iter(canonical(protein)):
                            end+=1;start=end-len(sequence)
                            if '*' in protein[:end]:continue
                            blocks=project(start*3,end*3,segments,m.strand)
                            if len(blocks)!=2 or min(v-u for u,v in blocks)<anchor:continue
                            context=cleavage(protein,start,end,False,False)
                            emit('alternative_splice',m,sequence,is_decoy,blocks,context,('junctions',((x,y),)),'junction_sequence_verified',edges=[(x,y)],subtype=subtype)
    store.db.commit()
    strict=store.counts(all_mappings=True)
    qualities=quality_lookup(metadata)
    grouped={}
    for p in pending:
        m=p['model']
        if p.get('internal') and strict.get(p['sequence'],0)!=1:
            # Multi-locus internal starts are explained by the annotation at one
            # of their loci; they are counted, not proposed.
            stats['nonunique_internal_start_matches']+=1;continue
        event_id=feature_id(p['kind'],m.chrom,m.strand,p['key'])
        grouped.setdefault(event_id,[]).append(p)
    # Decoy events: identical local filters; genome-wide uniqueness is not
    # evaluated for decoys, so the estimate is conservative.
    # Each decoy feature is binned by its longest peptide, like target features.
    decoy_features=defaultdict(dict)
    stats['decoy_local_matches']=len(decoy_pending)
    for p in decoy_pending:
        if p['status'] in {'upstream_stop_without_start'}:continue
        fid=feature_id(p['kind'],p['model'].chrom,p['model'].strand,p['key'])
        decoy_features[p['kind']][fid]=max(decoy_features[p['kind']].get(fid,0),len(p['sequence']))
    # Locus clusters of every placement (searched, annotated and hypothetical).
    clusters=ClusterIndex()
    for sequence,items in store.placements({p['sequence'] for items in grouped.values() for p in items}).items():
        for chrom,strand,blocks in items:clusters.add(sequence,chrom,strand,blocks)
    background=start_background(annotations,all_assoc)
    def weak(key):
        value=background[key]
        return value is not None and value<2
    weak_background={'internal_ATG_nontryptic_N_terminus':weak('m_initial_enrichment'),
                     'internal_ATG_initiator_methionine_removed':weak('methionine_before_removal_residue_enrichment')}
    rows=[]
    for event_id,items in sorted(grouped.items()):
        first=items[0];kind=first['kind'];m=first['model'];chrom,strand=m.chrom,m.strand
        seqs=sorted({p['sequence'] for p in items});refs={p['model'].id:p['model'] for p in items}
        blocks=sorted({tuple(b) for p in items for b in p['blocks']})
        statuses=sorted({p['status'] for p in items})
        blockers=[s for s in statuses if s in {'initiator_methionine_ambiguity','upstream_stop_without_start'}]
        upstream=sorted({u for p in items for u in p['upstream'] if strict.get(u,0)==1 and qualities.get(u)!={'qvalue_fail'} and u not in seqs})
        if kind=='alternative_start' and upstream:blockers.append('conflicting_upstream_evidence')
        consistent={p['sequence'] for p in items if all(
            feature_consistent(kind,first['subtype'] if kind=='alternative_splice' else '',strand,x['blocks'],clusters.members(x['sequence'],chrom,strand,x['blocks']),[tuple(e) for e in x['edges']])
            for x in items if x['sequence']==p['sequence'])}
        confidence,specific,reasons=judge(seqs,strict,qualities,blockers,consistent)
        support_edges=sorted({(a[1],b[0]) for p in items for a,b in zip(p['blocks'],p['blocks'][1:]) if a[1]<b[0]})
        rna_minimum=min((rna_junctions.get((chrom,a,b,strand),0) for a,b in support_edges),default=0)
        # RNA support is evidence for the junction itself, so it promotes only splice features.
        if kind=='alternative_splice' and confidence=='peptide_supported' and support_edges and rna_minimum>=args.min_junction_support:confidence='rna_peptide_supported'
        q={p for s in seqs for p in qualities.get(s,{'unknown'})}
        q='qvalue_pass' if 'qvalue_pass' in q else 'qvalue_fail' if q=={'qvalue_fail'} else 'unknown'
        limits=['local_evidence_does_not_validate_full_transcript','bounded_annotation_guided_search_not_exhaustive','peptide_identification_quality_'+q,'uniqueness_is_genomic_loci_over_all_exact_placements_and_searched_paths','coding_connection_is_hypothetical_until_independently_supported']
        if not all(p['context']['tryptic_n'] and p['context']['tryptic_c'] for p in items):limits.append('not_fully_tryptic_in_proposed_local_context')
        if kind=='alternative_start':
            limits+=['proteolytic_processing_not_excluded','translation_initiation_not_confirmed']
            if first.get('internal') and any(weak_background.get(s) for s in statuses):limits.append('internal_start_rate_not_above_random_fragment_background_in_this_run')
        if kind=='alternative_stop':limits.append('stop_inferred_from_assembly_not_directly_observed')
        if stats['splice_contexts_capped']:limits.append('annotation_guided_splice_search_capped')
        status=statuses[0] if len(statuses)==1 else ','.join(statuses)
        feature=[list(b) for b in first['feature']]
        edges=sorted({tuple(e) for p in items for e in p['edges']})
        upstream_count=len({s for p in items for s in p['upstream']}) if kind=='alternative_start' else None
        rows.append(dict(event_id=event_id,event_type=kind,event_subtype=first['subtype'] or ('initiator_methionine_removed' if 'initiator_methionine_removed' in status else ''),
            candidate_model_id='',supporting_model_ids=[],reference_model_ids=sorted(refs),annotation_sources=sorted({x.source for x in refs.values()}),
            chrom=chrom,start=blocks[0][0],end=blocks[-1][1],strand=strand,candidate_cds_blocks=[],reference_cds_blocks={x.id:sorted(x.coding) for x in refs.values()},
            relationship='annotation_guided_local_change',peptides=sorted({p['peptide'] for p in items}),
            mapping_ids=sorted({stable_id('L',p['sequence'],chrom,strand,p['encoded']) for p in items}),distinct_diagnostic_peptides=len(seqs),
            mapping_multiplicity={p['peptide']:strict.get(p['sequence'],0) for p in sorted(items,key=lambda x:x['peptide'])},
            peptide_loci={p['peptide']:dict(loci_all_exact=strict.get(p['sequence'],0),specific=p['sequence'] in specific,length=len(p['peptide'])) for p in items},
            rna_junctions=[dict(chrom=chrom,start=a,end=b,strand=strand,reads=rna_junctions.get((chrom,a,b,strand),0)) for a,b in support_edges],
            minimum_rna_support=rna_minimum,confidence=confidence,confidence_reasons=reasons,reasons=status,limitations=limits,coding_structure_id='',
            local_peptide_blocks=blocks,diagnostic_junctions=[list(e) for e in edges],splice_chain_changed=bool(edges),whole_model_status='local_evidence_only',
            whole_model_review_reasons=['no_full_transcript_inferred'],uniqueness_scope='genomic_loci_all_exact_placements_and_searched_paths',evidence_origin='annotation_guided',
            boundary_status=status,cleavage_context=[dict(peptide=p['peptide'],reference_model_id=p['model'].id,**p['context']) for p in items],
            upstream_peptide_support=upstream_count,proposed_feature_blocks=feature,
            longest_specific_peptide=max((len(s) for s in specific),default=0)))
    def fdr(decoys,targets):return None if not targets else round(min(1.0,decoys/targets),4)
    decoy_summary={}
    for kind in ['alternative_start','alternative_stop','alternative_splice','n_terminal_extension','c_terminal_extension']:
        # Internal starts test existing placements, not the local search space.
        searched=[r for r in rows if r['event_type']==kind and r['confidence'] in SUPPORTED and r['event_subtype']!='internal_start']
        decoys=decoy_features.get(kind,{})
        by_length={}
        for label,(low,high) in LENGTH_BINS.items():
            t=sum(low<=r['longest_specific_peptide']<=high for r in searched)
            d=sum(low<=n<=high for n in decoys.values())
            by_length[label]=dict(supported_target_features=t,decoy_features=d,estimated_fdr=fdr(d,t))
        decoy_summary[kind]=dict(supported_target_features=len(searched),decoy_features_passing_local_filters=len(decoys),
                                 estimated_fdr=fdr(len(decoys),len(searched)),by_longest_peptide_length=by_length)
    for row in rows:
        # Each row carries the FDR of its event type and peptide-length class.
        if row['event_subtype']=='internal_start' or row['confidence'] not in SUPPORTED:row['decoy_estimated_fdr']=''
        else:row['decoy_estimated_fdr']=decoy_summary[row['event_type']]['by_longest_peptide_length'][length_bin(row['longest_specific_peptide'])]['estimated_fdr']
        del row['longest_specific_peptide']
    directory=args.output/'refinement';directory.mkdir(exist_ok=True)
    write_table(directory/'annotation_guided_events.tsv',EVENT_FIELDS,rows)
    with (directory/'annotation_guided_features.gtf').open('w') as stream:
        for row in rows:
            feature='start_codon' if row['event_type']=='alternative_start' else 'stop_codon' if row['event_type']=='alternative_stop' else 'peptide_supported_region'
            for a,b in row['proposed_feature_blocks'] or row['local_peptide_blocks']:
                attrs='; '.join(k+' '+json.dumps(str(v)) for k,v in [('gene_id',row['annotation_sources'][0]+':'+lookup[row['reference_model_ids'][0]].gene),('transcript_id',row['event_id']),('event_id',row['event_id']),('event_type',row['event_type']),('confidence',row['confidence']),('evidence_origin','annotation_guided'),('interpretation','local_candidate_requires_review')])+';'
                stream.write(f'{row["chrom"]}\tGAPMS_v1\t{feature}\t{a+1}\t{b}\t.\t{row["strand"]}\t.\t{attrs}\n')
    stats.update(events=len(rows),internal_start_matches_examined=len(internal),
                 supported_events=sum(r['confidence'] in SUPPORTED for r in rows))
    result=dict(stats,seconds=round(time.monotonic()-started,3),refinement_window_nt=radius,splice_window_nt=splice_radius,max_paths_per_context=cap,
                min_peptide_anchor_nt=anchor,decoy=decoy_summary,alternative_start_background=background,
                scope='contiguous terminal extensions, internal ATG starts (incl. initiator-Met removal), nearby canonical GT-AG junction shifts and single-exon skips; no arbitrary new-exon or full-transcript inference')
    LOG.info('Annotation-guided refinement completed: %s local events (%s supported), %s splice paths',len(rows),result['supported_events'],stats['local_splice_paths_searched'])
    return rows,result
