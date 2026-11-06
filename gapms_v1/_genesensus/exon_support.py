"""Observed exon pairing and source-frame support for novel graph paths."""

class ExonSupport:
    def __init__(self, source_models, support, peptides=None):
        self.protein=support.protein_exon_frames
        self.peptides=peptides or {}
        self.observed=set()
        for models in source_models.values():
            for model in models:
                self.observed.update(self.tokens(model))
        self.rna={(chrom,strand,iv) for (chrom,strand),intervals in support.rna_internal.items() for iv in intervals}

    @staticmethod
    def tokens(model):
        ordered=model.cds if model.strand=='+' else tuple(reversed(model.cds))
        used=-model.phase
        for i,iv in enumerate(ordered):
            role='single' if len(ordered)==1 else 'first' if i==0 else 'last' if i==len(ordered)-1 else 'internal'
            yield model.contig,model.strand,role,iv,used%3
            used+=iv[1]-iv[0]

    def supports(self, model, matched_peptides=()):
        for token in self.tokens(model):
            if token in self.observed:
                continue
            chrom,strand,role,iv,frame=token
            if (chrom,strand,iv,frame) in self.protein:
                continue
            if any(iv in self.peptides[key].blocks[1:-1] for key in matched_peptides):
                continue
            # RNA does not assert coding frame; the complete ORF has been validated
            # on the genome. It must supply the observed whole-exon geometry.
            if role!='single' and (chrom,strand,iv) in self.rna:
                continue
            return False
        return True
