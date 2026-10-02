# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""fig11b_pangenome_plot.py — 从 pangenome_summary.tsv 重绘 Figure 11（不重算聚类）"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p3lib import WORK
from fig_style import plt, savefig

OUT = os.path.join(WORK, '09_figures')
SUM = os.path.join(WORK, '12_functional_annotation', 'summary')

vals = {}
shared = {}
names = []
with open(os.path.join(SUM, 'pangenome_summary.tsv'), encoding='utf-8') as f:
    for row in csv.reader(f, delimiter='\t'):
        if len(row) < 2:
            continue
        k, v = row[0], row[1]
        if k == 'genomes':
            names = v.split(';')
        elif k.startswith('genes_in_'):
            g = k[len('genes_in_'):]
            parts = v.replace(' core=', '\tcore=').replace(' accessory=', '\taccessory=').replace(' unique=', '\tunique=').split('\t')
            total = int(parts[0])
            dd = {}
            for p in parts[1:]:
                kk, vv = p.split('=')
                dd[kk] = int(vv)
            vals[g] = (total, dd)
        elif k.startswith('shared_'):
            g1, g2 = k[len('shared_'):].split('_vs_')
            shared[(g1, g2)] = int(v)
        elif k == 'total_clusters':
            ntot = int(v)


DISPLAY = {'smbu06': 'smbu06', 'smbu08': 'smbu08', '194': '194', 'IDCC2105': 'IDCC 2105',
           'CX2-6_2': 'CX 2-6_2', 'E_faecium_64-3': 'E. faecium'}
SHORT = [DISPLAY.get(n, n) for n in names]

fig = plt.figure(figsize=(7.4, 6.2))
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.42)

ax = fig.add_subplot(gs[0, 0])
core_n = [vals[n][1]['core'] for n in names]
acc_n = [vals[n][1]['accessory'] for n in names]
uniq_n = [vals[n][1]['unique'] for n in names]
y = np.arange(len(names))
ax.barh(y, core_n, color='#4DAF4A', label='core (all %d genomes)' % len(names), edgecolor='black', lw=0.3)
ax.barh(y, acc_n, left=core_n, color='#FFD92F', label='accessory (2–%d)' % (len(names) - 1), edgecolor='black', lw=0.3)
ax.barh(y, uniq_n, left=[c + a for c, a in zip(core_n, acc_n)], color='#E41A1C', label='unique (1)', edgecolor='black', lw=0.3)
for i, (c, a, u) in enumerate(zip(core_n, acc_n, uniq_n)):
    ax.text(c / 2, i, str(c), ha='center', va='center', fontsize=5.0, color='white')
    ax.text(c + a / 2, i, str(a), ha='center', va='center', fontsize=5.0)
    ax.text(c + a + u + 40, i, str(u), ha='left', va='center', fontsize=5.0, color='#E41A1C')
ax.set_yticks(y)
ax.set_yticklabels(SHORT, fontsize=6)
ax.invert_yaxis()
ax.set_xlim(0, 3050)
ax.set_xlabel('number of genes', fontsize=6)
ax.tick_params(labelsize=5.5)
ax.legend(fontsize=5.0, loc='lower right')
ax.set_title('a  Pangenome composition per genome', fontsize=7.5)
ax.text(0.01, 0.97, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

ax2 = fig.add_subplot(gs[0, 1])
S = np.zeros((len(names), len(names)))
for i, a in enumerate(names):
    for j, b in enumerate(names):
        if i == j:
            S[i, j] = vals[a][0]
        else:
            S[i, j] = shared.get((a, b), shared.get((b, a), 0))
im = ax2.imshow(S, cmap='YlGnBu')
ax2.set_xticks(range(len(names)))
ax2.set_xticklabels(SHORT, rotation=45, ha='right', fontsize=5.5)
ax2.set_yticks(range(len(names)))
ax2.set_yticklabels(SHORT, fontsize=6)
for i in range(len(names)):
    for j in range(len(names)):
        ax2.text(j, i, '%d' % S[i, j], ha='center', va='center', fontsize=5.0,
                 color='white' if S[i, j] > S.max() * 0.6 else 'black')
ax2.set_title('b  Shared orthologue clusters (pairwise)', fontsize=7.5)
ax2.text(-0.16, 1.03, 'b', transform=ax2.transAxes, fontsize=9, fontweight='bold')
cb = fig.colorbar(im, ax=ax2, shrink=0.75)
cb.ax.tick_params(labelsize=5.5)

fig.text(0.5, 0.005, ('Proteins predicted uniformly with pyrodigal; orthologue clusters built by all-vs-all BLASTP (>=60% identity, >=60% query coverage) and connected components ({n} clusters total).\n'
         'Comparison panel includes the strains used in the reference-style comparison (IDCC 2105, CX 2-6_2).').format(n=ntot),
         fontsize=5.6, ha='center', va='bottom')
savefig(fig, os.path.join(OUT, 'Fig11_pangenome', 'Fig11'))
print('fig11b done')
