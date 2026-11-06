"""Native observed RNA-chain geometry; references here are genomic DNA only.

Caller owns input fingerprints and cache publication. This adapter preserves the
junction template counts from extract_junctions and adds unique-chain geometry.
"""
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor
from .rna import alignment_index,transcript_strand
from .genome import Genome

_GENOME=None

def _initialize(path,index):
    global _GENOME
    _GENOME=Genome.from_index(path,index)

def fragments(start,cigar,settings):
    pos=left=start;blocks=[];valid=[]
    for i,(op,n) in enumerate(cigar):
        if op==3:
            blocks.append((left,pos))
            before=cigar[i-1][1] if i and cigar[i-1][0] in (0,7,8) else 0
            after=cigar[i+1][1] if i+1<len(cigar) and cigar[i+1][0] in (0,7,8) else 0
            valid.append(settings['min_intron']<=n<=settings['max_intron'] and min(before,after)>=settings['min_overhang'])
            pos+=n;left=pos
        elif op in (0,2,7,8):
            pos+=n
    blocks.append((left,pos));run=[blocks[0]];answer=[]
    for ok,block in zip(valid,blocks[1:]):
        if ok:
            run.append(block)
        else:
            if len(run)>1:answer.append(tuple(run))
            run=[block]
    if len(run)>1:answer.append(tuple(run))
    return answer

def _accepted(read,settings):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:return False
    if settings['skip_duplicates'] and read.is_duplicate:return False
    nh=read.get_tag('NH') if read.has_tag('NH') else None
    if nh is not None and nh!=1:return False
    if read.mapping_quality==255:return nh==1
    return read.mapping_quality>=settings['min_mapq']

def _orientations(read,chrom,blocks,settings,genome):
    strand=transcript_strand(read,settings['library'])
    if strand!='.':return [strand]
    xs=read.get_tag('XS') if read.has_tag('XS') else None
    if xs in ('+','-'):return [xs]
    votes=set()
    for left,right in zip(blocks,blocks[1:]):
        motif=(genome.fetch(chrom,left[1],left[1]+2).upper(),genome.fetch(chrom,right[0]-2,right[0]).upper())
        if motif in {('GT','AG'),('GC','AG'),('AT','AC')}:votes.add('+')
        elif motif in {('CT','AC'),('CT','GC'),('GT','AT')}:votes.add('-')
    return sorted(votes) if len(votes)==1 else ['+','-']

def _work(args,genome=None):
    import pysam
    bam,chrom,settings,index,prior,minimum=args
    genome=genome if genome is not None else _GENOME
    junctions=Counter();exons=Counter();all_exons=Counter();observed=set();seen=set();stats=Counter()
    def templates(a,b,strand):return max(prior.get((a,b,strand),0),prior.get((a,b,'.'),0))
    with pysam.AlignmentFile(str(bam),'rb',index_filename=index) as reader:
        for read in reader.fetch(chrom):
            stats['records']+=1
            if not _accepted(read,settings):continue
            for chain in fragments(read.reference_start,read.cigartuples or [],settings):
                # Unstranded compatibility preserves the same collapse convention
                # as the matched control. Known-strand libraries retain strand.
                strand=transcript_strand(read,settings['library'])
                key=chain if strand=='.' else (strand,chain)
                if key in seen:continue
                seen.add(key);stats['unique_geometry_chains']+=1
                introns=[(left[1],right[0]) for left,right in zip(chain,chain[1:])]
                observed.update(introns)
                for orient in _orientations(read,chrom,chain,settings,genome):
                    stats['emitted_alignment_chains']+=1
                    for a,b in introns:junctions[a,b,orient]+=1
                    for i,(a,b) in enumerate(chain[1:-1],1):
                        all_exons[a,b,orient]+=1
                        if templates(*introns[i-1],orient)>=minimum and templates(*introns[i],orient)>=minimum:
                            exons[a,b,orient]+=1
    assert observed=={(a,b) for a,b,s in prior},('RNA geometry differs from junction extraction',chrom)
    return chrom,junctions,exons,all_exons,dict(stats)

def augment_geometry(bam,genome,prior,threads=1,minimum_templates=3):
    """Add structure counts in a second indexed pass; never assemble transcripts.

    `genome` is an opened genesensus.genome.Genome. `prior` is the immutable
    output of extract_junctions. Counts represent distinct genomic chain
    geometries, not independent templates; original template counts are kept.
    """
    if isinstance(threads,bool) or not isinstance(threads,int) or threads<1:raise ValueError('threads must be a positive integer')
    if isinstance(minimum_templates,bool) or not isinstance(minimum_templates,int) or minimum_templates<1:raise ValueError('minimum_templates must be a positive integer')
    if prior.get('format_version')!=1 or prior.get('coordinate_system')!='0-based-half-open':raise ValueError('Invalid RNA junction cache format')
    lengths=prior['references']
    if not all(genome.lengths.get(c)==n for c,n in lengths.items()):raise ValueError('RNA/genome reference mismatch')
    settings=dict(prior['settings']);index=alignment_index(bam,lengths);local=defaultdict(dict)
    for chrom,a,b,strand,n in prior['junctions']:
        key=(a,b,strand);local[chrom][key]=max(local[chrom].get(key,0),n)
    args=[(str(bam),c,settings,index,local[c],minimum_templates) for c in sorted(lengths)]
    if threads==1:
        results=[_work(arg,genome) for arg in args]
    else:
        with ProcessPoolExecutor(max_workers=threads,initializer=_initialize,initargs=(str(genome.path),genome.index)) as executor:
            results=list(executor.map(_work,args))
    data=dict(prior)
    data['rna_chain_junctions']=[[chrom,a,b,strand,n] for chrom,junctions,exons,all_exons,stats in results for (a,b,strand),n in sorted(junctions.items())]
    data['observed_internal_exon_counts']=[[chrom,a,b,strand,n] for chrom,junctions,exons,all_exons,stats in results for (a,b,strand),n in sorted(exons.items())]
    data['observed_internal_exon_counts_all']=[[chrom,a,b,strand,n] for chrom,junctions,exons,all_exons,stats in results for (a,b,strand),n in sorted(all_exons.items())]
    data['observed_internal_exons']=[row[:-1] for row in data['observed_internal_exon_counts']]
    data['alignment_geometry_provenance']=dict(method='native indexed observed-chain extraction',reference_annotations_used=False,transcript_assembly=False,junction_geometry_unchanged=True,internal_exon_min_templates=minimum_templates,unfiltered_exon_count_policy='All retained observed internal exon geometries, scoring source-exact exons only; not discovery',depth_policy='One vote per unique genomic chain geometry; template counts unchanged',source_bam=str(bam),stats=dict(sum((Counter(stats) for c,j,e,all_e,stats in results),Counter())))
    return data
