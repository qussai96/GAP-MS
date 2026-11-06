"""Reference-free, additive secondary-isoform admissions. Frozen pilot policy."""
from .intervals import intersection_length,subtract,disputed
from .exon_support import ExonSupport
from .molecular import local_risk
GATE=.437721461

def shared_frame_bases(a,b):
 def segments(m):
  used=-m.phase;out=[]
  for x,y in m.cds if m.strand=='+' else reversed(m.cds):
   out.append((x,y,(used-x)%3 if m.strand=='+' else (used+y-1)%3));used+=y-x
  return out
 total=0
 for x,y,p in segments(a):
  for u,v,q in segments(b):
   if p==q:total+=max(0,min(y,v)-max(x,u))
 return total

def relationship(alt,primary):
 if alt.contig!=primary.contig or alt.strand!=primary.strand:return False,'different_contig_or_strand',0.
 shared=intersection_length(alt.cds,primary.cds)
 fraction=min(shared/alt.length,shared/primary.length)
 if fraction<.5:return False,'insufficient_shared_coding_anchor',fraction
 if shared_frame_bases(alt,primary)!=shared:return False,'shared_CDS_frame_disagreement',fraction
 if alt.introns==primary.introns:return False,'same_intron_chain_terminal_variant',fraction
 if not(set(alt.introns)&set(primary.introns) or set(alt.cds)&set(primary.cds)):return False,'no_shared_exact_exon_or_junction',fraction
 return True,'compatible_coding_alternative',fraction

def routes(candidate,parent,support,use_scores,use_rna,diagnostics,gate=GATE):
 m=candidate.model;new=sorted(set(m.introns)-set(parent.introns));lost=sorted(set(parent.introns)-set(m.introns));extra=subtract(m.cds,parent.cds);nextra=sum(b-a for a,b in extra)
 q=min(candidate.gene_scores) if use_scores and candidate.gene_scores else None
 qpositive=q is not None and q>=gate
 rna=False;counts=[];coverage=None;retention_observations=[]
 if use_rna:
  counts=[support.rna_count(m,iv) for iv in new]
  coverage=support.rna_coding_overlap(m) if not nextra else intersection_length(extra,support.rna_blocks.get((m.contig,m.strand),()))/nextra
  # rna_coding_overlap returns a fraction; no extra coding bases requires no new-base assertion.
  if not nextra:coverage=1.
  if not new and lost:
   for a,b in lost:
    observed=[support.rna_exon_counts_all.get((m.contig,m.strand,x,y),0) for x,y in m.cds[1:-1] if x<=a and y>=b]
    retention_observations.append(max(observed,default=0))
  event_supported=(bool(new) and all(n>=3 for n in counts)) or (bool(lost) and not new and all(n>=3 for n in retention_observations))
  rna=bool(event_supported and coverage>=.95)
 risk,known=local_risk(candidate,disputed(m,parent,9),gate) if use_scores else (None,False)
 spurious=bool(use_scores and q is not None and q<gate and any(cl=='cryptic_orf' and posterior>=.5 for cl,posterior in diagnostics))
 risk_veto=bool(use_scores and q is not None and q<gate and known and risk>=.5)
 return dict(q=q,q_positive=qpositive,distinct_predictor_sources=len(candidate.sources),RNA_event_supported=rna,new_junctions=new,lost_junctions=lost,alternative_junction_templates=counts,alternative_coding_coverage=coverage,retention_exon_observations=retention_observations,residue_risk=risk,residue_risk_observed=known,spurious_ORF_veto=spurious,residue_risk_veto=risk_veto)

def qualifies(features,policy):
 molecular=features['q_positive'] and (policy=='score_positive' or features['distinct_predictor_sources']>=2)
 supported=molecular or features['RNA_event_supported']
 if policy=='risk_guarded' and (features['spurious_ORF_veto'] or features['residue_risk_veto']):supported=False
 return bool(supported)
