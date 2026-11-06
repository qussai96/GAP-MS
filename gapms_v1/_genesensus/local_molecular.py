"""Source-context molecular corrections, never reconstructed-protein scores."""
from .molecular import local_risk
from .intervals import merge


def frame_at(model, position):
    used = 0
    for a, b in model.cds:
        if a <= position < b:
            offset = used + position-a if model.strand == '+' else model.length-used-(position-a)-1
            return (offset-model.phase) % 3
        used += b-a
    return None


class LocalMolecular:
    def __init__(self, originals, parameters):
        self.originals = [c for c in originals if c.sources and c.q is not None]
        self.p = parameters

    def correction(self, candidate, regions):
        risk, known = local_risk(candidate, regions, self.p['localizer_gate'])
        q = candidate.q
        return q - self.p['risk_strength']*(1-q)*(risk if known else 0.)

    def coding(self, model, regions=None):
        """Integrate source-capped, same-frame corrections over coding segments."""
        originals = [c for c in self.originals
                     if (c.model.contig,c.model.strand)==(model.contig,model.strand)]
        total = 0.
        for left,right in (model.cds if regions is None else regions):
            boundaries = {left,right}
            for c in originals:
                for a,b in c.model.cds:
                    if max(a,left)<min(b,right):
                        boundaries.update((max(a,left),min(b,right)))
            points = sorted(boundaries)
            for a,b in zip(points,points[1:]):
                source = {}
                for c in originals:
                    if (any(x<=a and b<=y for x,y in c.model.cds)
                            and frame_at(model,a)==frame_at(c.model,a)):
                        value = self.correction(c,((a,b),))
                        for s in c.sources:
                            source[s] = max(source.get(s,-float('inf')),value)
                if source:
                    total += (b-a)*sum(source.values())/len(source)
        return total

    def splice(self, model):
        total = 0.
        width = self.p['boundary_window_nt']
        for a,b in model.introns:
            source = {}
            for c in self.originals:
                m = c.model
                if ((m.contig,m.strand)!=(model.contig,model.strand)
                        or (a,b) not in m.introns
                        or frame_at(model,a-1)!=frame_at(m,a-1)
                        or frame_at(model,b)!=frame_at(m,b)):
                    continue
                flanks = merge((max(x,u),min(y,v)) for x,y in m.cds
                               for u,v in ((max(0,a-width),a),(b,b+width))
                               if max(x,u)<min(y,v))
                value = self.correction(c,flanks)
                for s in c.sources:
                    source[s] = max(source.get(s,-float('inf')),value)
            if source:
                total += (b-a)*sum(source.values())/len(source)
        return total

    def features(self, model):
        return self.coding(model),self.splice(model)

    def graph_corrections(self, models):
        """Edge priors on supplied geometry, without guessing a novel frame.

        Search state phase is checked during graph extension. Corrections are
        indexed by edge and source reading frame to avoid promoting an
        incompatible translation during pruning.
        """
        bounds = sorted({p for m in models for iv in m.cds for p in iv})
        result = {}
        for a,b in zip(bounds,bounds[1:]):
            frames = {}
            for c in self.originals:
                if any(x<=a and b<=y for x,y in c.model.cds):
                    frame = frame_at(c.model,a if c.model.strand=='+' else b-1)
                    source = frames.setdefault(frame,{})
                    value = self.correction(c,((a,b),))
                    for s in c.sources:
                        source[s] = max(source.get(s,-float('inf')),value)
            for frame,source in frames.items():
                result['coding',a,b,frame] = self.p.get('molecular_weight',1.)*sum(source.values())/len(source)
        for a,b in sorted({iv for m in models for iv in m.introns}):
            frames = {}
            for c in self.originals:
                if (a,b) not in c.model.introns:
                    continue
                m=c.model
                position=a-1 if m.strand=='+' else b
                frame=(frame_at(m,position)+1)%3
                width=self.p['boundary_window_nt']
                flanks=merge((max(x,u),min(y,v)) for x,y in m.cds
                             for u,v in ((max(0,a-width),a),(b,b+width))
                             if max(x,u)<min(y,v))
                source=frames.setdefault(frame,{})
                for s in c.sources:
                    source[s]=max(source.get(s,-float('inf')),self.correction(c,flanks))
            for frame,source in frames.items():
                result['splice',a,b,frame]=self.p.get('molecular_weight',1.)*sum(source.values())/len(source)
        return result
