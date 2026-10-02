import os, sys
import numpy as np
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK, ASM06, ASM08, read_fasta
from fig_style import plt, savefig, C

BLUE, ORANGE = '#2B6CB0', C['prophage']
s06 = read_fasta(ASM06['Chromosome1'])[0][1]
s08 = read_fasta(ASM08['Chromosome1'])[0][1]
attL, attJ = s06[1255782:1255928], s08[1255782:1255928]
sites = [i for i in range(146) if attL[i] != attJ[i]]

fig = plt.figure(figsize=(7.2, 3.82))
gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.0], hspace=0.26, wspace=0.34)

# a: all resolved chromosome-scale information on one clean axis.
ax = fig.add_subplot(gs[0, :]); ax.axis('off'); ax.set_xlim(-0.22, 2.95); ax.set_ylim(-0.36, 1.42)
ax.text(-0.02, 1.05, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold')
ax.text(0, 1.24, 'Chromosome-scale comparison', fontsize=8.5, fontweight='bold')
for y, L, label in [(0.80, 2.676906, 'smbu06 chromosome'), (0.12, 2.638187, 'smbu08 chromosome')]:
    ax.plot([0, L], [y, y], color=BLUE, lw=4.2, solid_capstyle='round')
    ax.text(-0.05, y, label, ha='right', va='center', fontsize=6.8)
    ax.text(L+0.03, y, f'{L:.6f} Mb', color=BLUE, va='center', fontsize=6.0)
x0,x1=1.255782,1.294501
ax.add_patch(Rectangle((x0,0.60),x1-x0,0.40,color=ORANGE))
ax.text((x0+x1)/2,1.10,'38.7-kb att-flanked region',ha='center',fontsize=6.4,color=ORANGE)
ax.plot([x0,x0],[-0.08,0.32],color=ORANGE,lw=1.5); ax.text(x0,-0.23,'attJ',ha='center',fontsize=6,color=ORANGE)
ax.plot(2.404143,0.12,'o',color='#222222',ms=3.7); ax.text(2.404143,-0.23,'1 SNP',ha='center',fontsize=6)
ax.text(0,-0.34,'2,638,186 bp of smbu08 aligned colinearly to smbu06 (100.00% coverage).',fontsize=6.1)

# b: zoomed structural logic
ax = fig.add_subplot(gs[1,0]); ax.axis('off'); ax.set_xlim(0,100); ax.set_ylim(-0.50,2.05)
ax.text(-0.10,1.05,'b',transform=ax.transAxes,fontsize=10,fontweight='bold')
ax.set_title('Resolved structural difference',fontsize=7.3,loc='left',pad=3)
ax.plot([10,94],[1.45,1.45],color=BLUE,lw=4,solid_capstyle='round')
ax.add_patch(Rectangle((43,1.20),14,0.50,color=ORANGE)); ax.text(50,1.92,'38.7 kb',ha='center',fontsize=6.3,color=ORANGE)
ax.text(43,0.98,'attL',ha='center',fontsize=5.8); ax.text(57,0.98,'attR',ha='center',fontsize=5.8); ax.text(7,1.45,'smbu06',ha='right',va='center',fontsize=6.4)
ax.plot([10,94],[0.36,0.36],color=BLUE,lw=4,solid_capstyle='round'); ax.plot([43,43],[0.06,0.68],color=ORANGE,lw=1.5)
ax.text(43,-0.21,'attJ',ha='center',fontsize=6.1,color=ORANGE); ax.text(7,0.36,'smbu08',ha='right',va='center',fontsize=6.4)
ax.text(50,-0.52,'One resolved 38.7-kb difference plus one SNP',ha='center',fontsize=5.9)

# c: att evidence and compact facts card
ax = fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_xlim(-10,154); ax.set_ylim(-0.60,2.05)
ax.text(-0.13,1.05,'c',transform=ax.transAxes,fontsize=10,fontweight='bold')
ax.set_title('attL and attJ sequence evidence',fontsize=7.3,loc='left',pad=3)
ax.plot(range(146),[1.18]*146,color=BLUE,lw=4,alpha=.55,solid_capstyle='butt')
ax.plot(range(146),[0.30]*146,color=ORANGE,lw=4,alpha=.55,solid_capstyle='butt')
for i in sites: ax.plot([i,i],[.95,1.42],color='#333333',lw=.52)
ax.text(-6,1.18,'attL',ha='right',va='center',fontsize=6.2); ax.text(-6,.30,'attJ',ha='right',va='center',fontsize=6.2)
ax.text(73,-.08,f'{len(sites)}/146 positions differ',ha='center',fontsize=6.2)
ax.text(0,-.34,'Interpretation',fontsize=6.1,fontweight='bold'); ax.text(0,-.55,'Gene content and att architecture are consistent with\na PBSX-like phage-related element.',fontsize=5.7,va='top')

fig.text(.5,.008,'This figure reports resolved sequence structure only; it does not infer induction, replication mechanism, or excision timing.',ha='center',fontsize=6)
savefig(fig, os.path.join(WORK,'09_figures','Fig2_chromosome_revised','Fig2_revised'))
