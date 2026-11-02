"""Summarize completed production runs and blinded holdout benchmarks."""
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve();jobs=json.loads((root/'jobs.json').read_text());production=[];holdout=[]
for j in jobs:
 out=Path(j['output']);m=json.loads((out/'run_manifest.json').read_text());assert m['status']=='completed'
 s=json.loads((out/'reports/summary.json').read_text());c=s['mapping'];n=c['mapped_peptides'];a=s['assessment']
 production.append(dict(crop=j['crop'],species=' '.join(j['crop'].split('_')[-2:]),job_id=j['job_id'],peptides=c['peptides'],mapped=n,mapped_percent=100*n/c['peptides'],unique=c['unique_peptides'],multiple=c['multiple_peptides'],ambiguous_percent_mapped=100*c['multiple_peptides']/n,mapping_rows=c['mapping_rows'],minutes=s['seconds']/60,implementation_sha256=m['implementation_sha256'],**{'primary_'+k:v for k,v in a['candidate_categories'].items()}))
 h=json.loads((root/'heldout'/j['crop']/'recovery_summary.json').read_text())
 with (root/'heldout'/j['crop']/'target_recovery.tsv').open() as f:targets=list(csv.DictReader(f,delimiter='\t'))
 selection=json.loads((root/'heldout'/j['crop']/'selection.json').read_text())
 h['targets_below_100aa']=sum(len(m['protein'])<100 for m in selection['target_models'])
 h['targets_without_complete_start_or_stop']=sum(not(m['complete_start'] and m['complete_stop']) for m in selection['target_models'])
 h['genes_with_peptide_supported_RNA_candidate']=sum(int(r['supported_candidate_models'])>0 for r in targets)
 h['genes_all_expected_peptide_loci_recovered']=sum(r['expected_peptide_loci']==r['recovered_peptide_loci'] for r in targets)
 h['single_exon_exact_cds']=sum(r['exact_cds_recovered']=='True' and int(r['truth_junctions'])==0 for r in targets)
 h['multi_exon_exact_cds']=sum(r['exact_cds_recovered']=='True' and int(r['truth_junctions'])>0 for r in targets)
 holdout.append(h)
def save(name,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (root/(name+'.tsv')).open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
 (root/(name+'.json')).write_text(json.dumps(rows,indent=2)+'\n')
save('production_results',production);save('holdout_results',holdout)
sums={k:sum(r[k] for r in holdout) for k in ['heldout_genes','multi_exon_targets','expected_peptide_loci','recovered_peptide_loci','exact_cds_recovered','exact_splice_chain_recovered','start_recovered','stop_recovered','junction_tp','junction_fn','junction_fp','genes_with_peptide_supported_RNA_candidate','genes_all_expected_peptide_loci_recovered','single_exon_exact_cds','multi_exon_exact_cds']}
sums.update(peptide_coordinate_recall=sums['recovered_peptide_loci']/sums['expected_peptide_loci'],full_cds_recall=sums['exact_cds_recovered']/sums['heldout_genes'],full_splice_chain_recall=sums['exact_splice_chain_recovered']/sums['multi_exon_targets'],splice_junction_recall=sums['junction_tp']/(sums['junction_tp']+sums['junction_fn']),splice_junction_precision=sums['junction_tp']/(sums['junction_tp']+sums['junction_fp']))
(root/'holdout_pooled.json').write_text(json.dumps(sums,indent=2)+'\n')
lines=['# RefSeq peptide runs and annotation-removal results','', 'All nine production runs completed using the previous assembly, Helixer annotation, BAM and mapping settings, with the supplied RefSeq FragPipe peptide files. The same-frame refinement classification fix affects proposal labels, not genomic mapping.','', '| Species | Searched | Mapped | Mapped % | Multiply mapped % of mapped | Runtime (min) |','|---|---:|---:|---:|---:|---:|']
for r in production:lines.append(f"| {r['species']} | {r['peptides']:,} | {r['mapped']:,} | {r['mapped_percent']:.2f} | {r['ambiguous_percent_mapped']:.2f} | {r['minutes']:.1f} |")
lines+=['','Production counts above refer to normalized query sequences; the ambiguity figures additionally collapse I/L-equivalent query aliases.','', '## Blinded annotation-removal benchmark','', '| Species | Genes | Exact peptide coordinates | Exact full CDS | Exact splice chain (multi-exon) | Junction recall | Junction precision |','|---|---:|---:|---:|---:|---:|---:|']
for r in holdout:lines.append(f"| {' '.join(r['crop'].split('_')[-2:])} | {r['heldout_genes']} | {r['recovered_peptide_loci']}/{r['expected_peptide_loci']} | {r['exact_cds_recovered']}/{r['heldout_genes']} | {r['exact_splice_chain_recovered']}/{r['multi_exon_targets']} | {r['junction_recall']:.2%} | {r['junction_precision']:.2%} |")
lines+=['',f"Pooled exact peptide-coordinate recall: **{sums['peptide_coordinate_recall']:.2%}**. Exact full-CDS recovery: **{sums['full_cds_recall']:.2%}**. Full coding splice-chain recovery among multi-exon targets: **{sums['full_splice_chain_recall']:.2%}**.",'','These are different tasks. Correct peptide placements do not establish complete CDS structures. Full-CDS recovery requires one peptide-supported RNA model to match the complete strand-aware coding path; merely reporting an ORF at the locus does not count.','', 'Read [the protocol](README.md) before interpreting these values. Targets are deterministic, stratified samples of peptide-supported annotated genes, not a random sample of all genes. Reference peptide identification predates masking. Candidate junctions unmatched to the selected target isoform can be legitimate alternative isoforms; this precision is concordance with withheld annotations, not a genome-wide biological false-discovery rate.','', '## Figures','', 'Each box represents the crop-level distribution at one peptide length, with individual crop dots and mapped-query sample sizes. All observed lengths are included; lengths below seven were not searched.','', '![Casanovo ambiguity](figures/casanovo/ambiguity_all_lengths.png)','', '![RefSeq ambiguity](figures/refseq/ambiguity_all_lengths.png)','', 'Vector PDF/SVG files and crop-by-length count/percentage TSV files accompany each PNG.','']
notes=[f"The weaker ORF-presence criterion would count **{sums['genes_with_peptide_supported_RNA_candidate']}/{sums['heldout_genes']} genes ({100*sums['genes_with_peptide_supported_RNA_candidate']/sums['heldout_genes']:.1f}%)** as rediscovered because they have at least one peptide-supported overlapping RNA candidate. Only **{sums['exact_cds_recovered']}/{sums['heldout_genes']} ({100*sums['full_cds_recall']:.1f}%)** match the complete target CDS. Individual splice-junction recall is **{sums['splice_junction_recall']:.2%}**, and concordance precision against the selected target structures is **{sums['splice_junction_precision']:.2%}**."]
if (root/'heldout_truth_audit.json').exists():
 audits=json.loads((root/'heldout_truth_audit.json').read_text());assert len(audits)==len(jobs) and all(r['status']=='passed' for r in audits)
 misses=json.loads((root/'heldout_missed_peptide_loci.json').read_text())
 absent=sum(r['reason']=='junction_absent_from_search_evidence' for r in misses)
 notes.append(f"All **{sums['expected_peptide_loci']:,} expected peptide-coordinate paths** passed independent nucleotide projection and BioPython translation checks. **{absent}/{len(misses)} unrecovered paths** require junctions absent from the search evidence. Retained-annotation overlap checks and cached RNA/input identity checks passed for all crops.")
if (root/'RNA_structure_comparison.json').exists():
 rna=json.loads((root/'RNA_structure_comparison.json').read_text());shares=[r['shared_old_percent'] for r in rna]
 notes.append(f"Fresh projected TransDecoder CDS structure multisets shared **{min(shares):.1f}–{max(shares):.1f}%** with the previous predictions. Assembled RNA exon structures were identical in all crops (see RNA_exon_structure_comparison.json). Fresh ORF counts are therefore not a perfectly isolated peptide-list comparison; the held-out benchmark uses fixed RNA predictions to avoid this confound.")
if (root/'validation.json').exists():
 checks=json.loads((root/'validation.json').read_text());assert all(r['status']=='validated' for r in checks)
 notes.append(f"All production validations passed, including independent translation of **{sum(r['independently_translated_sample'] for r in checks):,} sampled peptide mappings** and **{sum(r['independently_translated_orfs'] for r in checks):,} sampled whole candidate coding paths**.")
report='\n'.join(lines).replace('## Figures','\n\n'.join(notes)+'\n\n## Figures')
(root/'RESULTS.md').write_text(report)
print(json.dumps(sums,indent=2))
