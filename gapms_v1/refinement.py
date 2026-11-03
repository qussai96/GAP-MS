"""Bounded annotation-guided peptide evidence; no RNA reconstruction required.

Hypotheses are local coding paths, not inferred full-length transcripts. New
paths are inserted before multiplicity and confidence are evaluated.
"""
from collections import Counter, defaultdict
from itertools import product
import json
import time

import ahocorasick

from gapms._peplocus.sequence import canonical, join_blocks, translate, project
from .assessment import EVENT_FIELDS, ModelIndex, boundary, relationship
from .common import LOG, stable_id, write_table
from .mapping import locus_key


def ordered(blocks, strand):
    blocks=json.loads(locus_key('',strand,blocks)[2])
    return [tuple(b) for b in (blocks if strand=='+' else reversed(blocks))]


def sliced(model, start, end):
    return ordered(project(start,end,model.coding,model.strand),model.strand)


def cleavage(protein, start, end, native_start=False, native_end=False):
    previous=protein[start-1] if start else ''
    following=protein[end:end+1]
    standard=set('ACDEFGHIKLMNPQRSTVWY')
    n_ok=(start==0 and native_start) or bool(previous in {'K','R'} and protein[start:start+1]!='P')
    c_ok=(native_end and (end==len(protein) or following=='*')) or bool(protein[end-1:end] in {'K','R'} and following in standard-{'P'})
    return dict(previous_aa=previous,following_aa=following,n_terminus='protein_boundary' if start==0 and native_start else 'tryptic' if n_ok else 'nontryptic' if previous in standard else 'upstream_stop_boundary' if previous=='*' else 'unknown',c_terminus='protein_boundary' if native_end and (end==len(protein) or following=='*') else 'tryptic' if c_ok else 'nontryptic' if following else 'unknown',tryptic_n=n_ok,tryptic_c=c_ok)


def refine_annotations(args, annotations, store, metadata, reference, rna_junctions):
    started=time.monotonic()
    annotations=[m for m in annotations if m.source!='RNA_ORF' and m.coding]
    stats=Counter();pending=[]
    radius=getattr(args,'refinement_window',3000)
    splice_radius=getattr(args,'refinement_splice_window',30)
    cap=getattr(args,'max_refinement_paths',64)
    anchor=getattr(args,'min_refinement_anchor',9)
    max_peptide=max(map(len,metadata))
    max_nt=3*(max_peptide+1)
    lookup={m.id:m for m in annotations}
    index=ModelIndex(annotations)
    known_paths=set();internal=[];model_occurrences=defaultdict(list)
    aliases=defaultdict(list)
    for p in metadata:
        if len(p)>=args.min_length:aliases[canonical(p)].append(p)
    # One streaming pass; RNA associations do not explain annotation novelty.
    for mid,peptide,sequence,chrom,strand,encoded,start,end in store.db.execute('SELECT * FROM associations'):
        if mid not in lookup:continue
        known_paths.add((sequence,chrom,strand,encoded))
        model_occurrences[mid].append((start,sequence))
        m=lookup[mid]
        if m.issues or start<=0 or not m.protein[start:].startswith('M'):continue
        context=cleavage(m.protein,start,end,True,m.complete_stop)
        if context['n_terminus']!='nontryptic' or not context['tryptic_c']:continue
        codon=m.project(start*3,start*3+3)
        if join_blocks(reference,m.chrom,m.strand,codon).upper()!='ATG':continue
        internal.append((m,peptide,sequence,json.loads(encoded),start,context,codon))
    eligible=set(aliases)
    explained_sequences={key[0] for key in known_paths}
    eligible-=explained_sequences
    for sequence,chrom,strand,encoded in store.db.execute('SELECT canonical,chrom,strand,blocks FROM hits'):
        if (sequence,chrom,strand,encoded) not in known_paths:eligible.add(sequence)
    stats['unexplained_query_sequences']=len(eligible)
    automaton=None
    if aliases:
        automaton=ahocorasick.Automaton()
        for sequence in sorted(aliases):automaton.add_word(sequence,sequence)
        automaton.make_automaton()
    comparable_cache={}
    def comparable(model):
        if model.id not in comparable_cache:
            comparable_cache[model.id]=[m for m in index.overlap(model) if m.strand==model.strand and relationship(model,[m])=='new_transcript_isoform']
        return comparable_cache[model.id]
    def is_known_boundary(model, point, which):
        return any(boundary(m,which)==point for m in comparable(model))
    def add(kind,model,peptide,blocks,context,diagnostic_edges=(),status='unresolved',feature=(),upstream=None):
        chrom,strand,encoded=locus_key(model.chrom,model.strand,blocks)
        sequence=canonical(peptide)
        if kind!='alternative_start' and (sequence,chrom,strand,encoded) in known_paths:return
        matched=translate(join_blocks(reference,chrom,strand,json.loads(encoded)),args.genetic_code)
        if canonical(matched)!=sequence:raise RuntimeError('annotation-guided projection verification failed')
        store.add(peptide,chrom,strand,json.loads(encoded),matched,'annotation_guided')
        pending.append(dict(kind=kind,model=model,peptide=peptide,sequence=sequence,blocks=json.loads(encoded),encoded=encoded,context=context,edges=tuple(diagnostic_edges),status=status,feature=list(feature),upstream=upstream))
    def scan(model, segments, kind, old_offset, terminal_stop=False):
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
        for end,sequence in automaton.iter(canonical(protein)):
            end+=1;start=end-len(sequence)
            if start<valid_start or end>valid_end:continue
            if kind=='n_terminal_extension' and start*3>=old_offset:continue
            if kind=='n_terminal_extension' and 'X' in protein[start:old_offset//3]:continue
            if kind=='c_terminal_extension' and 'X' in protein[old_offset//3:end]:continue
            if kind=='c_terminal_extension' and end*3<=old_offset:continue
            blocks=project(start*3,end*3,segments,model.strand)
            context=cleavage(protein,start,end,False,stop_after==end)
            event_kind=kind;status='boundary_unresolved';feature=[]
            if kind=='n_terminal_extension' and protein[start:start+1]=='M' and context['n_terminus'] in {'nontryptic','upstream_stop_boundary'} and context['tryptic_c']:
                codon=project(start*3,start*3+3,segments,model.strand)
                if join_blocks(reference,model.chrom,model.strand,codon).upper()=='ATG':
                    point=codon[0][0] if model.strand=='+' else codon[-1][1]
                    if not is_known_boundary(model,point,'start'):
                        event_kind='alternative_start';status='ATG_supported_start_candidate';feature=list(codon)
            if kind=='c_terminal_extension' and stop_after>=0:
                codon=project(stop_after*3,(stop_after+1)*3,segments,model.strand)
                point=codon[-1][1] if model.strand=='+' else codon[0][0]
                if not is_known_boundary(model,point,'end'):
                    event_kind='alternative_stop';status='genomic_stop_inferred';feature=list(codon)
                else:
                    stats['extensions_with_already_annotated_stop']+=1
            for peptide in aliases[sequence]:add(event_kind,model,peptide,blocks,context,status=status,feature=feature)
    # Downstream starts are a special context test, not a new sequence/path.
    for m,peptide,sequence,blocks,start,context,codon in internal:
        point=codon[0][0] if m.strand=='+' else codon[-1][1]
        if is_known_boundary(m,point,'start'):
            stats['internal_starts_already_annotated']+=1;continue
        upstream=len({p for offset,p in model_occurrences[m.id] if offset<start})
        add('alternative_start',m,peptide,blocks,context,status='internal_ATG_nontryptic_N_terminus',feature=codon,upstream=upstream)
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
                    skipped=sum(v-u for u,v in genomic[i+1:k])
                    variants=[(x,y) for x,y in product(donors,acceptors) if y-x>=20 and (x-b-(y-c)-skipped)%3==0 and (m.chrom,m.strand,x,y) not in known_edges]
                    variants.sort(key=lambda xy:(abs(xy[0]-b)+abs(xy[1]-c),xy))
                    if len(variants)>cap:stats['splice_contexts_capped']+=1
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
                        for end,sequence in automaton.iter(canonical(protein)):
                            end+=1;start=end-len(sequence)
                            if '*' in protein[:end]:continue
                            blocks=project(start*3,end*3,segments,m.strand)
                            if len(blocks)!=2 or min(v-u for u,v in blocks)<anchor:continue
                            context=cleavage(protein,start,end,False,False)
                            for peptide in aliases[sequence]:add('alternative_splice',m,peptide,blocks,context,[(x,y)],'junction_sequence_verified')
    store.db.commit()
    counts=store.counts()
    grouped={}
    for p in pending:
        # Internal starts require locus specificity after all new paths were added.
        if p['kind']=='alternative_start' and counts[p['sequence']]!=1:
            stats['nonunique_internal_start_matches']+=1;continue
        m=p['model']
        locus_anchor=boundary(m,'start' if p['kind']=='n_terminal_extension' else 'end') if not p['edges'] and not p['feature'] else None
        key=(p['kind'],m.chrom,m.strand,m.source,m.gene,tuple(p['edges']),p['status'],tuple(map(tuple,p['feature'])),locus_anchor)
        group=grouped.setdefault(key,[]);group.append(p)
    rows=[]
    for key,items in sorted(grouped.items(),key=lambda x:str(x[0])):
        kind,chrom,strand,source,gene,edges,status,feature,locus_anchor=key
        seqs=sorted({p['sequence'] for p in items});refs={p['model'].id:p['model'] for p in items}
        blocks=sorted({tuple(b) for p in items for b in p['blocks']});unique=any(counts[p]==1 for p in seqs)
        qualities={metadata[p['peptide']]['quality'] for p in items}
        q='qvalue_pass' if 'qvalue_pass' in qualities else 'qvalue_fail' if 'qvalue_fail' in qualities else 'unknown'
        support_edges=sorted({(a[1],b[0]) for p in items for a,b in zip(p['blocks'],p['blocks'][1:]) if a[1]<b[0]})
        minimum=min((rna_junctions.get((chrom,a,b,strand),0) for a,b in support_edges),default=0)
        confidence='ambiguous' if not unique or q=='qvalue_fail' else 'rna_peptide_supported' if support_edges and minimum>=args.min_junction_support else 'peptide_supported'
        limits=['local_evidence_does_not_validate_full_transcript','bounded_annotation_guided_search_not_exhaustive','peptide_identification_quality_'+q,'uniqueness_relative_to_searched_paths','coding_connection_is_hypothetical_until_independently_supported']
        if not all(p['context']['tryptic_n'] and p['context']['tryptic_c'] for p in items):limits.append('not_fully_tryptic_in_proposed_local_context')
        if kind=='alternative_start':limits+=['proteolytic_processing_not_excluded','translation_initiation_not_confirmed','upstream_peptide_absence_does_not_establish_absence_of_translation']
        if kind=='alternative_stop':limits.append('stop_inferred_from_assembly_not_directly_observed')
        if stats['splice_contexts_capped']:limits.append('annotation_guided_splice_search_capped')
        ident=stable_id('AG',key)
        rows.append(dict(event_id=ident,event_type=kind,candidate_model_id='',reference_model_ids=sorted(refs),annotation_sources=sorted({m.source for m in refs.values()}),chrom=chrom,start=blocks[0][0],end=blocks[-1][1],strand=strand,candidate_cds_blocks=[],reference_cds_blocks={m.id:sorted(m.coding) for m in refs.values()},relationship='annotation_guided_local_change',peptides=seqs,mapping_ids=sorted({stable_id('L',p['sequence'],chrom,strand,p['encoded']) for p in items}),distinct_diagnostic_peptides=len(seqs),mapping_multiplicity={p:counts[p] for p in seqs},rna_junctions=[dict(chrom=chrom,start=a,end=b,strand=strand,reads=rna_junctions.get((chrom,a,b,strand),0)) for a,b in support_edges],minimum_rna_support=minimum,confidence=confidence,reasons=status,limitations=limits,coding_structure_id='',local_peptide_blocks=blocks,diagnostic_junctions=list(edges),splice_chain_changed=bool(edges),whole_model_status='local_evidence_only',whole_model_review_reasons=['no_full_transcript_inferred'],uniqueness_scope='searched_assembly_and_supplied_splice_paths',evidence_origin='annotation_guided',boundary_status=status,cleavage_context=[dict(peptide=p['sequence'],reference_model_id=p['model'].id,**p['context']) for p in items],upstream_peptide_support=max((p['upstream'] for p in items if p['upstream'] is not None),default=None),proposed_feature_blocks=list(feature)))
    directory=args.output/'refinement';directory.mkdir(exist_ok=True)
    write_table(directory/'annotation_guided_events.tsv',EVENT_FIELDS,rows)
    with (directory/'annotation_guided_features.gtf').open('w') as stream:
        for row in rows:
            feature='start_codon' if row['event_type']=='alternative_start' else 'stop_codon' if row['event_type']=='alternative_stop' else 'peptide_supported_region'
            for a,b in row['proposed_feature_blocks'] or row['local_peptide_blocks']:
                attrs='; '.join(k+' '+json.dumps(str(v)) for k,v in [('gene_id',row['annotation_sources'][0]+':'+lookup[row['reference_model_ids'][0]].gene),('transcript_id',row['event_id']),('event_id',row['event_id']),('event_type',row['event_type']),('evidence_origin','annotation_guided'),('interpretation','local_candidate_requires_review')])+';'
                stream.write(f'{row["chrom"]}\tGAPMS_v1\t{feature}\t{a+1}\t{b}\t.\t{row["strand"]}\t.\t{attrs}\n')
    stats.update(events=len(rows),internal_start_matches_examined=len(internal))
    result=dict(stats,seconds=round(time.monotonic()-started,3),refinement_window_nt=radius,splice_window_nt=splice_radius,max_paths_per_context=cap,min_peptide_anchor_nt=anchor,scope='contiguous terminal extensions, internal ATG starts, nearby canonical GT-AG junction shifts and single-exon skips; no arbitrary new-exon or full-transcript inference')
    LOG.info('Annotation-guided refinement completed: %s local events, %s splice paths',len(rows),stats['local_splice_paths_searched'])
    return rows,result
