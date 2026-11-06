"""Qualified local RNA observations for consensus event corroboration.

A union of matched bases supports coverage, template counts support junctions
and fully spanned retained introns. These are local, not full-transcript phases.
"""
from collections import defaultdict
import pysam


def union(intervals):
    out=[]
    for a,b in sorted(intervals):
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],b))
        else:out.append((a,b))
    return out


def subtract(left,right):
    out=[];right=union(right)
    for a,b in union(left):
        cursor=a
        for x,y in right:
            if y<=cursor:continue
            if x>=b:break
            if x>cursor:out.append((cursor,min(x,b)))
            cursor=max(cursor,y)
            if cursor>=b:break
        if cursor<b:out.append((cursor,b))
    return out


class RNAEvidence:
    def __init__(self):
        self.coverage=defaultdict(list);self.junction_templates=defaultdict(set);self.blocks=defaultdict(list);self.reads=0

    @classmethod
    def from_bam(cls,path,loci,min_mapq=20,anchor=8,library='unstranded',index_path=None):
        evidence=cls()
        if path is None:return evidence
        regions=defaultdict(list)
        for locus in loci:
            c=locus[0];regions[c.model.chrom].append((min(d.blocks[0][0] for d in locus),max(d.blocks[-1][1] for d in locus)))
        index_filename=None
        with pysam.AlignmentFile(path,'rb') as probe:
            indexed=probe.has_index()
        if not indexed:
            if index_path is None:raise ValueError('consensus RNA evidence needs a BAM index or a writable run-local index path')
            index_filename=str(index_path)
            pysam.index('-c','-o',index_filename,str(path))
        with pysam.AlignmentFile(path,'rb',index_filename=index_filename) as bam:
            for chrom,intervals in sorted(regions.items()):
                for left,right in union(intervals):
                    for read in bam.fetch(chrom,left,right):
                        if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate
                            or read.mapping_quality<max(20,min_mapq) or (read.has_tag('NH') and read.get_tag('NH')!=1)
                            or (read.mapping_quality==255 and not(read.has_tag('NH') and read.get_tag('NH')==1))):continue
                        strand='.'
                        if library!='unstranded':
                            reverse=read.is_reverse
                            if read.is_paired and read.is_read2:reverse=not reverse
                            if library=='firststrand':reverse=not reverse
                            strand='-' if reverse else '+'
                        token=(read.get_tag('RG') if read.has_tag('RG') else '',read.query_name)
                        pos=read.reference_start;blocks=[];splice=[];current=None;ops=read.cigartuples or []
                        for i,(op,n) in enumerate(ops):
                            if op in (0,7,8):
                                if current is None:current=pos
                                pos+=n
                            else:
                                if current is not None:blocks.append((current,pos));current=None
                                if op==3:
                                    left_anchor=ops[i-1][1] if i and ops[i-1][0] in (0,7,8) else 0
                                    right_anchor=ops[i+1][1] if i+1<len(ops) and ops[i+1][0] in (0,7,8) else 0
                                    if min(left_anchor,right_anchor)>=anchor and 20<=n<=500000:splice.append((pos,pos+n))
                                if op in (2,3):pos+=n
                        if current is not None:blocks.append((current,pos))
                        evidence.coverage[chrom,strand].extend(blocks)
                        for a,b in splice:evidence.junction_templates[chrom,a,b,strand].add(token)
                        for a,b in blocks:evidence.blocks[chrom,strand].append((a,b,token))
                        evidence.reads+=1
        evidence.coverage={k:union(v) for k,v in evidence.coverage.items()}
        return evidence

    def count(self,chrom,strand,a,b):
        return len(self.junction_templates.get((chrom,a,b,strand),set())|self.junction_templates.get((chrom,a,b,'.'),set()))

    def supports_change(self,candidate,primary):
        from .consensus import junctions,overlap_nt
        c,p=candidate,primary;ch,st=c.model.chrom,c.model.strand
        new=set(junctions(c.blocks))-set(junctions(p.blocks))
        retained=[(a,b) for a,b in junctions(p.blocks) if any(x<=a and b<=y for x,y in c.blocks)]
        changed=subtract(c.blocks,p.blocks)
        if not(new or retained):return False
        if any(self.count(ch,st,a,b)<3 for a,b in new):return False
        coverage=union(self.coverage.get((ch,st),[])+self.coverage.get((ch,'.'),[]))
        if changed and overlap_nt(changed,coverage)<.95*sum(b-a for a,b in changed):return False
        for a,b in retained:
            templates={token for x,y,token in self.blocks.get((ch,st),[])+self.blocks.get((ch,'.'),[]) if x<=a-8 and y>=b+8}
            if len(templates)<3:return False
        return True
