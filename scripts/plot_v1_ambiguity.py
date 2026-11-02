"""Crop-level mapping ambiguity at every observed peptide length."""
import argparse,csv,json
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('jobs',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--label',default='RefSeq FragPipe');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
jobs=json.loads(a.jobs.read_text());rows=[];maximum=7
for j in jobs:
 counts=defaultdict(Counter);seen={}
 with (Path(j['output'])/'mapping/peptide_locus_summary.tsv').open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   key=r['peptide'].replace('I','L')
   if key in seen:
    assert seen[key]==(r['status'],r['loci']);continue
   seen[key]=(r['status'],r['loci']);counts[int(r['length'])][r['status']]+=1
 maximum=max(maximum,max(counts))
 for length,c in sorted(counts.items()):
  mapped=c['unique']+c['multiple'];total=sum(c.values())
  rows.append(dict(crop=j['crop'],species=' '.join(j['crop'].split('_')[-2:]),length=length,unique=c['unique'],multiple=c['multiple'],unmapped=c['unmapped'],mapped=mapped,searched=total,ambiguous_percent_mapped=100*c['multiple']/mapped if mapped else '',ambiguous_percent_searched=100*c['multiple']/total))
with (a.out/'ambiguity_by_length.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
fig,(ax,nax)=plt.subplots(2,1,figsize=(max(15,maximum*.29),7.4),gridspec_kw={'height_ratios':[5,1]},sharex=True)
positions=list(range(7,maximum+1));data=[[r['ambiguous_percent_mapped'] for r in rows if r['length']==x and r['mapped']] for x in positions]
valid=[(x,d) for x,d in zip(positions,data) if d]
ax.boxplot([d for x,d in valid],positions=[x for x,d in valid],widths=.68,whis=(0,100),manage_ticks=False,showfliers=False,patch_artist=True,boxprops={'facecolor':'#e6e6e6','edgecolor':'#555'},medianprops={'color':'black','linewidth':1.3},whiskerprops={'color':'#777'},capprops={'color':'#777'})
colors=plt.get_cmap('tab10').colors
for i,j in enumerate(jobs):
 points=[r for r in rows if r['crop']==j['crop'] and r['mapped']]
 ax.scatter([r['length']+(i-4)*.065 for r in points],[r['ambiguous_percent_mapped'] for r in points],s=11,color=colors[i],label=' '.join(j['crop'].split('_')[-2:]),zorder=3,alpha=.8)
 nax.plot([r['length'] for r in points],[r['mapped'] for r in points],color=colors[i],alpha=.65,lw=.9)
ax.set(ylabel='Peptides with >1 genomic location (%)',ylim=(-2,103),title=f'{a.label}: ambiguity by peptide length across nine crops')
ax.grid(axis='y',alpha=.2);ax.legend(ncol=3,fontsize=9,loc='upper right');ax.spines[['top','right']].set_visible(False)
nax.set(yscale='log',ylabel='Mapped\nqueries',xlabel='Peptide length (amino acids)',xlim=(6.4,maximum+.6));nax.grid(axis='y',alpha=.2);nax.spines[['top','right']].set_visible(False)
step=1 if maximum<=50 else 5;ticks=list(range(7,maximum+1,step));nax.set_xticks(ticks,labels=[str(x) for x in ticks]);nax.tick_params(labelsize=8)
fig.text(.065,.015,'Each dot = one crop. Boxes: median and interquartile range; whiskers: minimum–maximum across crops.\nDenominator = mapped query sequences at that length; I/L-equivalent query sequences and genomic placements are deduplicated. Lengths <7 were not searched. Empty lengths have no box.',fontsize=9)
fig.tight_layout(rect=(0,.07,1,1))
for ext in ['png','pdf','svg']:
 output=a.out/f'ambiguity_all_lengths.{ext}';fig.savefig(output,dpi=220)
 if ext=='svg':output.write_text('\n'.join(line.rstrip() for line in output.read_text().splitlines())+'\n')
(a.out/'figure_methods.json').write_text(json.dumps(dict(label=a.label,minimum_searched=7,maximum_observed=maximum,percentage='100 * multiple / (unique + multiple)',unit='one crop at one peptide length; I/L-canonical distinct query sequences',missing='no mapped queries: omitted, not imputed as zero',boxes='unweighted crop-level median/IQR; min/max whiskers',input_jobs=str(a.jobs.resolve())),indent=2)+'\n')
print(a.out, 'lengths 7–'+str(maximum))
