# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""figS2_nisin.py — Supplementary Fig. S2：nisin 双层检索与阳性/阴性对照"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p3lib import WORK
from fig_style import plt, savefig, C
import matplotlib.gridspec as gridspec

OUT = os.path.join(WORK, '09_figures')
SEA = os.path.join(WORK, '05_bacteriocins', 'search')

# 读取 tblastn 最佳命中
best = {}
with open(os.path.join(SEA, 'nisin_tblastn_hits.tsv'), encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if p[0] != 'tblastn':
            continue
        q = p[2].split('|')[0]
        t = p[3]
        key = (q, t)
        pid, qcov = float(p[4]), float(p[7])
        if key not in best or pid > best[key][0]:
            best[key] = (pid, qcov)

targets = ['F44_positive_control', 'IL1403_negative_control', 'smbu06_complete', 'smbu08_complete']
qnames = ['nisA', 'nisB', 'nisC', 'nisI', 'nisP', 'nisR', 'nisT', 'nisE', 'nisF', 'nisG', 'nisK']

fig = plt.figure(figsize=(7.2, 5.0))
gs = gridspec.GridSpec(1, 1, bottom=0.36, top=0.94)

ax = fig.add_subplot(gs[0, 0])
M = np.zeros((len(qnames), len(targets)))
for i, q in enumerate(qnames):
    for j, t in enumerate(targets):
        M[i, j] = best.get((q, t), (0, 0))[0]
im = ax.imshow(M, cmap='RdYlGn', vmin=0, vmax=100, aspect='auto')
ax.set_xticks(range(len(targets)))
ax.set_xticklabels(['F44\n(positive control)', 'IL1403\n(negative control)',
                    'smbu06', 'smbu08'], fontsize=6)
ax.set_yticks(range(len(qnames)))
ax.set_yticklabels(qnames, fontsize=6)
for i in range(len(qnames)):
    for j in range(len(targets)):
        v = M[i, j]
        ax.text(j, i, '%.1f' % v if v > 0 else '—', ha='center', va='center', fontsize=5.2,
                color='black' if v < 70 else 'white')
ax.set_title('Nisin reference proteins (tblastn, best identity %)\n'
             'cluster call: complete = nisA+nisB+nisC ≥80% identity / ≥80% qcov', fontsize=7)
cb = fig.colorbar(im, ax=ax, shrink=0.8)
cb.set_label('identity (%)', fontsize=6)
cb.ax.tick_params(labelsize=5.5)

fig.text(0.5, 0.02,
         'F44 positive control: complete nisin cluster (nisA 92.98%, nisB/C/I/P/R/T 99.6–100%). IL1403 negative control: not detected.\n'
         'smbu06/smbu08: not detected — only generic ABC-family homology (nisF 45.81%, nisT 34.45%) and a weak nisR-family response-regulator hit (34.96%);\n'
         'the F44 nisin-cluster nucleotide region (14,401 bp) has zero BLASTN hits in either genome. Conclusion: no complete nisin locus in smbu06/smbu08.',
         fontsize=5.8, ha='center', va='bottom')

savefig(fig, os.path.join(OUT, 'FigS2_nisin', 'FigS2'))
print('figS2 done')
