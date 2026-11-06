"""Direct measured candidate quality, with conditional local error adjustments.

One additional bounded quality channel; no structural admission contribution.
The source-context channel remains the controlled comparator. This is a
protein plausibility measurement, not calibrated transcript correctness.
"""
from .molecular import local_risk
from .local_molecular import LocalMolecular
from .intervals import merge

class DirectMolecular(LocalMolecular):
    def __init__(self,candidates,parameters):
        super().__init__(candidates,parameters)
        self.bykey={c.model.key:c for c in candidates}

    def features(self,model):
        candidate=self.bykey[model.key]
        if candidate.q is None:
            # Missing is unobserved. Do not substitute an inherited protein q.
            return 0.,0.
        coding=model.length*self.correction(candidate,model.cds)
        splice=0.
        width=self.p['boundary_window_nt']
        for a,b in model.introns:
            flanks=merge((max(x,u),min(y,v)) for x,y in model.cds
                         for u,v in ((max(0,a-width),a),(b,b+width))
                         if max(x,u)<min(y,v))
            splice+=(b-a)*self.correction(candidate,flanks)
        return coding,splice
