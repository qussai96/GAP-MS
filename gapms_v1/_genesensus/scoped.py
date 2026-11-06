"""Context-dependent local conflict decisions; compatibility admissions/fallback.

No empirical distinction: retain the existing molecular/consensus decision.
Local peptide, RNA or alignment observations: apply coverage-aware local guards.
No source is chosen by peptide counts or canonical-ORF fractions.
"""
from .adaptive import compare_adaptive
from .selection import compare

def compare_scoped(c,old,peptides,support,chains,votes,n_sources,p):
    local=compare_adaptive(c,old,peptides,support,chains,votes,n_sources,p)
    if any(k in local['channels'] for k in ('peptide','rna','protein')):
        local['decision_mode']='local_empirical'
        return local
    result=compare(c,old,peptides,support,chains,votes,n_sources,p)
    result.update(policy='scoped',decision_mode='compatibility_without_local_empirical')
    return result
