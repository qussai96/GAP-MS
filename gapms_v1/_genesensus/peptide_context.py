"""Positive regional peptide evidence with a conserved multi-placement budget.

Fractional weights are ambiguity allocations, not calibrated probabilities.
The dict interface deliberately exposes only unique genomic placements to old
hard splice / exon / ORF-veto guards. Ambiguous observations are soft evidence.
"""
from collections import defaultdict
from .peptides import Peptide, compatible


class PeptideEvidence(dict):
    def strong_translation(self,candidate):
        """Tiered positive support for the existing low-score ORF veto only.

        This is a conservative heuristic, not a protein-level FDR guarantee.
        Regional peptide observations remain available for reconstruction.
        """
        sequences={seq for seq in candidate.peptides if len(seq)>=9}
        if any(len(self[seq].blocks)>1 and min(b-a for a,b in self[seq].blocks)>=3
               for seq in sequences):return True
        nonnested=[seq for seq in sequences if not any(seq!=other and seq in other for other in sequences)]
        return len(nonnested)>=2

    def __init__(self, unique, placements):
        super().__init__(unique)
        self.placements = {seq: tuple(Peptide(seq,*loc) for loc in sorted(locations))
                           for seq,locations in sorted(placements.items())}
        self.matches = {}
        self.bins = defaultdict(set)
        for seq, locations in self.placements.items():
            for index, peptide in enumerate(locations):
                for a,b in peptide.blocks:
                    for cell in range(a//10000,(b-1)//10000+1):
                        self.bins[peptide.contig,peptide.strand,cell].add((seq,index))

    def attach(self, candidates, genome):
        for candidate in candidates:
            model=candidate.model
            possible=set()
            for a,b in model.cds:
                for cell in range(a//10000,(b-1)//10000+1):
                    possible.update(self.bins.get((model.contig,model.strand,cell),()))
            if not possible:
                self.matches[candidate.uid]=()
                continue
            protein=genome.assess(model)[2]
            matched=[]
            for seq,index in sorted(possible):
                peptide=self.placements[seq][index]
                if compatible(model,peptide,protein):
                    matched.append((seq,index))
                    if len(self.placements[seq])==1:
                        candidate.peptides.add(seq)
            self.matches[candidate.uid]=tuple(matched)

    def features(self, candidate, ambiguous=True):
        """At most one channel-vote/base; split credit capped by observed coding nt.

        Each sequence has a total placement budget of one over the full mapping
        universe, including locations outside the supplied annotation inventory.
        Overlaps and repeated placements never multiply per-base coverage above
        one. A split peptide divides its splice budget among its junctions.
        """
        events=defaultdict(float); junctions=defaultdict(lambda:defaultdict(float)); sequences=set()
        for seq,index in self.matches.get(candidate.uid,()):
            locations=self.placements[seq]
            if not ambiguous and len(locations)>1: continue
            peptide=locations[index]; weight=1./len(locations)
            sequences.add(seq)
            for a,b in peptide.blocks:
                events[a]+=weight; events[b]-=weight
            jj=[(a[1],b[0]) for a,b in zip(peptide.blocks,peptide.blocks[1:])]
            if jj:
                for iv in jj:
                    for a,b in peptide.blocks:
                        junctions[iv][a]+=weight/len(jj)
                        junctions[iv][b]-=weight/len(jj)
        coding=0.; active=0.; previous=None
        for position,delta in sorted(events.items()):
            if previous is not None:
                coding+=(position-previous)*min(1.,max(0.,active))
            active+=delta; previous=position
        splice=0.
        for a,b in candidate.model.introns:
            covered=0.; active=0.; previous=None
            for position,delta in sorted(junctions.get((a,b),{}).items()):
                if previous is not None:
                    covered+=(position-previous)*min(1.,max(0.,active))
                active+=delta; previous=position
            splice+=min(b-a,covered)
        return coding,splice,len(sequences)

    def graph_observations(self, models, support):
        """Unique placements only; never import excluded-predictor junctions.

        The mapper's complete HAT universe remains intact for uniqueness. A
        split path becomes a search hint only if each junction already occurs
        among the supplied predictions or independently supplied alignments.
        """
        if not models: return ()
        chrom,strand=models[0].contig,models[0].strand
        lo=min(m.cds[0][0] for m in models); hi=max(m.cds[-1][1] for m in models)
        jj={iv for m in models for iv in m.introns}
        possible=set()
        for m in models:
            for a,b in m.cds:
                for cell in range(a//10000,(b-1)//10000+1):
                    possible.update(self.bins.get((chrom,strand,cell),()))
        result=[]
        for seq,index in sorted(possible):
            if len(self.placements[seq]) != 1: continue
            p=self.placements[seq][index]
            if p.blocks[0][0]<lo or p.blocks[-1][1]>hi: continue
            introns=[(a[1],b[0]) for a,b in zip(p.blocks,p.blocks[1:])]
            if any(iv not in jj and not support.supported_junction(models[0],iv) for iv in introns): continue
            result.append(p)
        return tuple(result)

    def conservative_repair(self, candidates, chosen, valid, weights, support, p):
        """One-for-one local repair, retaining incumbent direct observations.

        No reward for missing peptides. An intervention needs newly explained,
        unique, locally contradictory observations: one >=9 aa split peptide
        with >=3 nt on each block, or two nonnested >=9 aa peptide sequences.
        These conservative fixed identification gates are not FDR estimates.
        A peptide outside the incumbent's span cannot on its own force extension.
        Baseline strong RNA states and all incumbent unique peptides are kept.
        """
        from .intervals import overlap
        chosen=set(chosen); events=[]
        while True:
            proposals=[]
            for i,c in enumerate(candidates):
                if i in chosen or not valid[i] or not c.complete or not c.canonical: continue
                removed=[j for j in chosen if overlap(c.model.cds,candidates[j].model.cds)]
                if len(removed)!=1: continue
                j=removed[0]; old=candidates[j]
                # Nonoverlapping retained models cannot be silently fused or lost.
                if not set(old.peptides)<=set(c.peptides): continue
                old_m=old.model; new_m=c.model
                if any(support.rna_count(old_m,iv)>=p['minimum_rna_templates']
                       and iv not in new_m.introns for iv in old_m.introns): continue
                if any(support.rna_exon_counts.get((old_m.contig,old_m.strand,*iv),0)>=p['minimum_rna_templates']
                       and iv not in new_m.cds for iv in old_m.cds): continue
                new=[]; junction=False
                for seq in sorted(set(c.peptides)-set(old.peptides)):
                    if len(seq)<9: continue
                    peptide=self[seq]
                    # Evidence must contradict a state inside the competing
                    # locus, not merely add an unsupported outer boundary.
                    inside=old_m.cds[0][0]<=peptide.blocks[0][0] and peptide.blocks[-1][1]<=old_m.cds[-1][1]
                    if not inside: continue
                    new.append(seq)
                    if len(peptide.blocks)>1 and min(b-a for a,b in peptide.blocks)>=3:
                        junction=True
                nonnested=[seq for seq in new if not any(seq!=other and seq in other for other in new)]
                if not junction and len(nonnested)<2: continue
                proposals.append((len(nonnested),weights[i],c.uid,i,j,new,junction))
            if not proposals: break
            _,_,_,i,j,new,junction=max(proposals,key=lambda x:(x[0],x[1],x[2]))
            chosen.remove(j); chosen.add(i)
            events.append(dict(event='peptide_contrast_repair',accepted=True,candidate=candidates[i].uid,
                incumbent=candidates[j].uid,new_unique_sequences=new,
                directly_splice_supported=junction,baseline_utility=weights[j],
                replacement_utility=weights[i],preserved_unique_peptides=True,
                preserved_strong_RNA_states=True))
        return chosen,events
