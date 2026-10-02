# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""A concise, evidence-focused replacement for Figure 3."""
import os, sys, csv
import numpy as np
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK, ASM06, ASM08, read_fasta
from fig_style import plt, savefig, C, gene_arrow

OUT = os.path.join(WORK, '09_figures', 'Fig3_shared_plasmid_revised')
ANN = os.path.join(WORK, '04_plasmidome', 'annotation')
P06 = read_fasta(ASM06['Plasmid1'])[0][1]
P08 = [x[1] for x in read_fasta(ASM08['Plasmid1']) if 'Plasmid1' in x[0]][0]
L = len(P06)
rot = (P06 + P06).find(P08)
P08n = P08[L-rot:] + P08[:L-rot] if rot else P08

genes=[]
with open(os.path.join(ANN, 'pl1_genes.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['region'] == 'pl1_smbu06': genes.append(r)
mges=[]
with open(os.path.join(ANN, 'mge_inventory.tsv'), encoding='utf-8') as f:
    next(f)
    for line in f:
        q=line.rstrip().split('\t')
        if len(q)>=5 and q[0]=='smbu06' and q[1]=='Plasmid1': mges.append((int(q[3]),int(q[4])))

def colour(r):
    p=(r['product'] or '').lower(); tag=r['locus_tag']
    if tag in ('smbu06GL002641','smbu06GL002679'): return C['bacteriocin']
    if tag == 'smbu06GL002642': return C['immunity']
    if 'transposase' in p or 'integrase' in p: return C['mge']
    return C['core'] if r['product'] != 'hypothetical protein' else C['other']

fig=plt.figure(figsize=(7.2,6.1))
gs=gridspec.GridSpec(3,1,height_ratios=[1.05,0.72,1.45],hspace=.55)

# a: readable whole-plasmid map
ax=fig.add_subplot(gs[0]); ax.set_xlim(0,L/1000); ax.set_ylim(-1.05,1.35); ax.axis('off')
ax.plot([0,L/1000],[0,0],color='#666666',lw=4,solid_capstyle='round')
for r in genes:
    a,b=int(r['start'])/1000,int(r['end'])/1000
    col=colour(r)
    if col != C['other']:
        ax.plot([a,b],[0,0],color=col,lw=6,solid_capstyle='butt')
for a,b in mges:
    ax.plot([(a+b)/2000],[-.22],marker='|',ms=10,color=C['mge'],mew=1.4)
labels=[(36.1,'enterocin P-like',C['bacteriocin'],.54),(36.4,'candidate\nimmunity',C['immunity'],-.74),(69.7,'hiracin JM79-like',C['bacteriocin'],.54)]
for x,t,c,y in labels:
    ax.annotate(t,(x,0),xytext=(x,y),ha='center',va='center',fontsize=7,color=c,
                arrowprops=dict(arrowstyle='-',lw=.55,color=c))
ax.text(0,1.18,'a',fontweight='bold',fontsize=10)
ax.text(4,1.18,'Shared circular plasmid backbone: 128,837 bp',fontsize=9,fontweight='bold')
ax.text(L/1000,1.18,'smbu06 = smbu08',fontsize=8,ha='right',color='#4D4D4D')
ax.text(L/2000,-.48,'Grey backbone: plasmid sequence; ticks: transposase/integrase annotations',ha='center',fontsize=6.8,color='#555555')

# b: one compact identity proof
ax=fig.add_subplot(gs[1]); ax.set_xlim(0,L/1000); ax.set_ylim(99.65,100.12)
x=np.arange(0,L,1000)/1000
ident=[100*sum(a==b for a,b in zip(P06[i:i+1000],P08n[i:i+1000]))/len(P06[i:i+1000]) for i in range(0,L,1000)]
ax.fill_between(x,99.65,ident,color='#D9D9D9'); ax.plot(x,ident,color='#4D4D4D',lw=1.2)
ax.set_yticks([99.7,100.0]); ax.set_ylabel('1-kb identity (%)',fontsize=7); ax.set_xlabel('smbu06 plasmid coordinate (kb)',fontsize=7); ax.tick_params(labelsize=6)
ax.text(.01,.12,'b',transform=ax.transAxes,fontweight='bold',fontsize=10)
ax.text(.99,.12,'Circular origin normalized; 128,837/128,837 bp; 0 SNPs; 0 indels',transform=ax.transAxes,ha='right',fontsize=6.8)

# c: locus neighborhood
ax=fig.add_subplot(gs[2]); lo,hi=30,76
for r in genes:
    a,b=int(r['start'])/1000,int(r['end'])/1000
    if b < lo or a > hi: continue
    gene_arrow(ax,a,b,.26 if r['strand']=='+' else -.26,r['strand'],colour(r))
for a,b in mges:
    a,b=a/1000,b/1000
    if lo<=a<=hi: ax.plot([a,b],[0,0],color=C['mge'],lw=1.2)
for x,t,c,y in labels:
    if lo<=x<=hi:
        ax.annotate(t,(x,.26 if y>0 else -.26),xytext=(x,y),ha='center',fontsize=7,color=c,
                    arrowprops=dict(arrowstyle='-',lw=.55,color=c))
ax.set_xlim(lo,hi); ax.set_ylim(-1.05,1.05); ax.set_yticks([]); ax.set_xlabel('position in shared plasmid (kb)',fontsize=7); ax.tick_params(labelsize=6)
ax.set_title('Bacteriocin-associated neighbourhood',loc='left',fontsize=9,fontweight='bold',pad=8)
ax.text(.01,.91,'c',transform=ax.transAxes,fontweight='bold',fontsize=10)
ax.text(.99,.91,'Candidate gene calls are sequence annotations; activity was not tested.',transform=ax.transAxes,ha='right',fontsize=6.5,color='#555555')

savefig(fig, os.path.join(OUT,'Fig3_revised'))
print('Fig3 revised complete')
