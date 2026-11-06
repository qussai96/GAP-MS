"""Splice-context qualification of existing original-source quality; no novel-protein q."""
from .local_molecular import LocalMolecular,frame_at
from .intervals import merge

class IntrinsicSourceMolecular(LocalMolecular):
    def __init__(self,originals,parameters,prior,genome):
        super().__init__(originals,parameters);self.prior=prior;self.genome=genome
    def splice(self,model):
        total=0.;width=self.p['boundary_window_nt']
        for a,b in model.introns:
            factor,_=self.prior.score(self.genome,model,a,b);source={}
            for c in self.originals:
                m=c.model
                if ((m.contig,m.strand)!=(model.contig,model.strand) or (a,b) not in m.introns or frame_at(model,a-1)!=frame_at(m,a-1) or frame_at(model,b)!=frame_at(m,b)):continue
                flanks=merge((max(x,u),min(y,v)) for x,y in m.cds for u,v in ((max(0,a-width),a),(b,b+width)) if max(x,u)<min(y,v))
                value=self.correction(c,flanks)-c.q*(1-factor)
                for s in c.sources:source[s]=max(source.get(s,-float('inf')),value)
            if source:total+=(b-a)*sum(source.values())/len(source)
        return total
