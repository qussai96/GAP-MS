"""Offline, exportable two-panel result figure with explicit counting units."""
import base64
from .events import SUPPORTED
from collections import Counter
from .common import write_json


def overview(root,diag,args,summary,events):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .outputs import assessed_scores
    rows=assessed_scores(diag,args)
    mapping=summary['mapping'];total=summary.get('input',{}).get('unique_peptides',mapping['peptides']);mapped=mapping['mapped_peptides']
    pairs=sum(int(r['mapped_peptides']) for r in rows)
    n=sum(int(r['N_terminal_peptides']) for r in rows)
    c=sum(int(r['C_terminal_peptides']) for r in rows)
    splice=sum(int(r['splice_peptides']) for r in rows)
    internal=sum(int(r['internal_peptides']) for r in rows)
    # Each review proposal is counted once; a proposal is not a distinct gene.
    unique_events=list({e['event_id']:e for e in events}.values())
    counts=Counter(e['event_type'] for e in unique_events)
    ambiguous=Counter(e['event_type'] for e in unique_events if e['confidence'] not in SUPPORTED)
    specific=counts-ambiguous
    categories=[('novel_orf','Novel ORF'),('alternative_start','Alternative start'),('n_terminal_extension','N-terminal extension'),
                ('alternative_stop','Alternative stop'),('c_terminal_extension','C-terminal extension'),('alternative_splice','Alternative splice')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none'})
    fig,(left,right)=plt.subplots(1,2,figsize=(15,5.8),gridspec_kw={'width_ratios':[1.15,1]})
    teal='#198678';gray='#e1e7eb';blue='#3e77ae';orange='#d9933a'
    data=[('Mapped',mapped,total-mapped,total,teal),('N-terminal',n,max(0,pairs-n),pairs,blue),
          ('C-terminal',c,max(0,pairs-c),pairs,blue),
          ('Splice',splice,max(0,pairs-splice),pairs,orange),
          ('Internal',internal,max(0,pairs-internal),pairs,'#8a74a6')]
    for y,(label,a,b,den,color) in enumerate(data):
        pct=100*a/den if den else 0
        left.barh(y,pct,color=color,height=.56)
        left.barh(y,100-pct,left=pct,color=gray,height=.56)
        left.text(102,y,f'{pct:.1f}%' if den else 'N/A',va='center',fontsize=11)
        if label=='Mapped':left.text(0,-.58,f'{mapped:,} mapped / {total:,} peptides',fontsize=11,color='#344451')
    left.set(yticks=range(len(data)),yticklabels=[d[0] for d in data],xlim=(0,117),xticks=[0,25,50,75,100],xlabel='Percent')
    left.set_ylim(len(data)-.35,-1);left.set_title('A  Peptide support',loc='left',fontweight='bold',pad=18)
    from matplotlib.patches import Patch
    left.legend(handles=[Patch(color=teal,label='Mapped'),Patch(color=blue,label='Terminal'),Patch(color=orange,label='Splice'),Patch(color='#8a74a6',label='Internal'),Patch(color=gray,label='Remaining')],loc='upper left',bbox_to_anchor=(-.03,-.18),ncol=3,frameon=False,fontsize=10)
    values=[counts[k] for k,label in categories]
    specific_values=[specific[k] for k,label in categories]
    ambiguous_values=[ambiguous[k] for k,label in categories]
    right.barh(range(len(values)),specific_values,color=teal,height=.58,label='Specific support')
    right.barh(range(len(values)),ambiguous_values,left=specific_values,color='#b8bfc5',height=.58,label='Ambiguous or not supported')
    right.legend(loc='upper left',bbox_to_anchor=(-.03,-.18),frameon=False,fontsize=10)
    right.set(yticks=range(len(values)),yticklabels=[label for k,label in categories],xlabel='Proposals')
    right.invert_yaxis();right.set_title('B  Changes to review',loc='left',fontweight='bold',pad=18)
    maximum=max(values,default=0);right.set_xlim(0,max(1,maximum)*1.25)
    for y,value in enumerate(values):right.text(value+max(1,maximum)*.025,y,f'{value:,}',va='center')
    from matplotlib.ticker import MaxNLocator
    right.xaxis.set_major_locator(MaxNLocator(integer=True,nbins=4))
    for ax in [left,right]:
        ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0);ax.set_axisbelow(True);ax.grid(axis='x',alpha=.15)
    fig.subplots_adjust(left=.15,right=.97,top=.85,bottom=.24,wspace=.85)
    directory=diag/'reports';svg=directory/'overview.svg';png=directory/'overview.png'
    fig.savefig(svg);fig.savefig(png,dpi=180);plt.close(fig)
    write_json(directory/'overview_counts.json',dict(peptides=total,mapped=mapped,peptide_protein_pairs=pairs,N_terminal=n,C_terminal=c,splice=splice,internal=internal,events=dict(counts),events_specific=dict(specific),events_ambiguous=dict(ambiguous),annotation='consensus' if len(args.gtf)>1 else 'supplied' if args.gtf else 'RNA_ORFs'))
    data=base64.b64encode(svg.read_bytes()).decode('ascii')
    return '<h2>Overview</h2><img style="width:100%;height:auto" alt="Peptide mapping, terminal and splice support, and review proposal counts" src="data:image/svg+xml;base64,'+data+'"><p><a href="diagnostics/reports/overview.svg">SVG</a> · <a href="diagnostics/reports/overview.png">PNG</a></p>'
