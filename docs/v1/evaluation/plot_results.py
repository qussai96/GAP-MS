from pathlib import Path
from collections import Counter
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root=Path(__file__).parent
rows=json.loads((root/'release_results.json').read_text());validation=json.loads((root/'ambiguity_by_length.json').read_text())
plt.rcParams['pdf.fonttype']=42
plt.rcParams['ps.fonttype']=42
fig,axes=plt.subplots(3,1,figsize=(11,12),constrained_layout=True)
labels=[r['species'] for r in rows];x=np.arange(len(rows))
unique=np.array([100*r['unique']/r['peptides'] for r in rows]);multiple=np.array([100*r['multiple']/r['peptides'] for r in rows]);unmapped=100-unique-multiple
ax=axes[0]
ax.bar(x,unique,label='One genomic path',color='#377eb8');ax.bar(x,multiple,bottom=unique,label='Multiple genomic paths',color='#e69f00');ax.bar(x,unmapped,bottom=unique+multiple,label='Unmapped in searched space',color='#cccccc')
ax.set_xticks(x,labels,rotation=20,ha='right',fontsize=8);ax.set_ylabel('Searched peptide queries (%)');ax.set_ylim(0,100);ax.set_title('A. Exact I/L-equivalent mapping; minimum length 7 aa',loc='left',fontsize=11);ax.legend(frameon=False,ncol=3,fontsize=8,loc='upper center')
ax=axes[1];colors=plt.cm.tab10.colors;pooled=Counter()
for i,record in enumerate(validation):
 counts=Counter({(r['length'],r['status']):r['peptides'] for r in record['ambiguity_by_length']});pooled.update(counts)
 lengths=[n for n in range(7,26) if counts[n,'unique']+counts[n,'multiple']>=100]
 fractions=[100*counts[n,'multiple']/(counts[n,'unique']+counts[n,'multiple']) for n in lengths]
 ax.plot(lengths,fractions,label=labels[i],color=colors[i],linewidth=1.3)
ax.set_xticks(range(7,26,2));ax.set_xlim(7,25);ax.set_ylim(0,100);ax.set_ylabel('Multiple paths among mapped queries (%)');ax.set_xlabel('Peptide length (aa)');ax.set_title('B. Length does not resolve all duplicated genomic loci',loc='left',fontsize=11);ax.legend(frameon=False,ncol=3,fontsize=7,loc='upper right')
ax=axes[2];bottom=np.zeros(len(rows))
parts=[('StringTie','stringtie_minutes','#999999'),('LongOrfs','longorfs_minutes','#e69f00'),('TransDecoder Predict','transdecoder_predict_minutes','#984ea3'),('Genome/path mapping','genome_mapping_minutes','#377eb8')]
for label,key,color in parts:
 values=np.array([r[key] for r in rows]);ax.bar(x,values,bottom=bottom,label=label,color=color);bottom+=values
other=np.maximum(0,np.array([r['initial_minutes'] for r in rows])-bottom);ax.bar(x,other,bottom=bottom,label='RNA extraction, association and other work',color='#4daf4a')
ax.set_xticks(x,labels,rotation=20,ha='right',fontsize=8);ax.set_ylabel('Initial full run (minutes)');ax.set_title('C. One observed run per crop; audited refreshes excluded',loc='left',fontsize=11);ax.legend(frameon=False,ncol=3,fontsize=7,loc='upper left')
for ax in axes:
 ax.spines[['top','right']].set_visible(False);ax.set_axisbelow(True);ax.grid(axis='y',alpha=.15)
fig.savefig(root/'evaluation.png',dpi=180);fig.savefig(root/'evaluation.pdf');plt.close(fig)
strata=[]
for low,high in [(7,8),(9,11),(12,14),(15,1000)]:
 unique=sum(v for (length,status),v in pooled.items() if low<=length<=high and status=='unique')
 multiple=sum(v for (length,status),v in pooled.items() if low<=length<=high and status=='multiple')
 strata.append(dict(length=f'{low}-{high}' if high<1000 else '15+',mapped=unique+multiple,multiple=multiple,multiple_percent=round(100*multiple/(unique+multiple),3)))
(root/'length_strata.json').write_text(json.dumps(strata,indent=2)+'\n');print(json.dumps(strata,indent=2))
