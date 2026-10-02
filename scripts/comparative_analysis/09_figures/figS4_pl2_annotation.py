# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""figS4_pl2_annotation.py — Supplementary Fig. S4：Plasmid2 与 38.7 kb 元件完整注释"""
import os
import sys
import csv

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK
from fig_style import plt, savefig, C, gene_arrow
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

OUT = os.path.join(WORK, '09_figures')
ANN = os.path.join(WORK, '04_plasmidome', 'annotation')

MODC = {'lysogeny': '#7570B3', 'DNA_packaging': '#1B9E77', 'head': '#D95F02',
        'tail': '#E7298A', 'lysis': '#E6AB02', 'recombination_DNA': '#66A2E5',
        'regulation': '#A6761D', 'hypothetical': C['other']}

rows = []
with open(os.path.join(ANN, 'prophage_annotation_full.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        r['start'] = int(r['start']); r['end'] = int(r['end'])
        rows.append(r)

fig = plt.figure(figsize=(7.2, 6.2))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 1.15], hspace=0.28)

# a) 元件基因图（smbu06 坐标系）
ax = fig.add_subplot(gs[0, 0])
lo, hi = 1254000, 1300000
drawn = [(r['start'], r['end'], r['strand'], r['locus_tag'], r['functional_module'],
          r['nr_description'] or r['product']) for r in rows if r['end'] > lo and r['start'] < hi]
for s0, s1, strand, lt, mod, desc in drawn:
    y = 0.45 if strand == '+' else -0.45
    gene_arrow(ax, s0 / 1000, s1 / 1000, y, strand, MODC.get(mod, C['other']), height=0.3)
# 标注关键基因
KEYG = {'smbu06GL001245': 'integrase', 'smbu06GL001272': 'terminase (PBSX)', 'smbu06GL001273': 'portal',
        'smbu06GL001284': 'tape measure (1,483 aa)', 'smbu06GL001292': 'holin',
        'smbu06GL001293': 'amidase (endolysin)'}
_k = 0
for s0, s1, strand, lt, mod, desc in drawn:
    if lt in KEYG:
        y = 0.45 if strand == '+' else -0.45
        yoff = (1.02, 1.20)[_k % 2] if y > 0 else (-1.02, -1.20)[_k % 2]
        _k += 1
        ax.annotate(KEYG[lt], xy=((s0 + s1) / 2000, y), xytext=((s0 + s1) / 2000, yoff),
                    ha='center', fontsize=4.8, arrowprops=dict(arrowstyle='-', lw=0.35))
# att 标记
for x, lab in [(1255.783, 'attL'), (1294.501, 'attR')]:
    ax.plot([x, x], [-0.85, 0.85], color='black', lw=0.6, ls='--')
    ax.text(x, 0.9, lab, fontsize=5.0, ha='center')
ax.set_xlim(lo / 1000, hi / 1000)
ax.set_ylim(-1.25, 1.25)
ax.set_yticks([])
ax.set_xlabel('smbu06 chromosome position (kb)', fontsize=6)
ax.tick_params(labelsize=5.5)
ax.set_title('a  The 38,719-bp PBSX-like element (51 CDS; modules colour-coded)', fontsize=7)
ax.legend(handles=[mpatches.Patch(facecolor=c, label=k) for k, c in MODC.items()],
          loc='lower right', ncol=4, fontsize=4.6)
ax.text(0.01, 0.94, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold')

# b) 基因表（按坐标；两列布局）
ax2 = fig.add_subplot(gs[1, 0])
ax2.axis('off')
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
items = [r for r in rows]
half = (len(items) + 1) // 2
for col, sub in enumerate([items[:half], items[half:]]):
    y = 1.0
    for r in sub:
        desc = (r['nr_description'] or r['product'])[:44]
        ax2.text(0.01 + col * 0.5, y, '%s  %s–%s(%s) %s' % (
            r['locus_tag'].replace('smbu06GL00', 'GL'), r['start'], r['end'],
            r['strand'], desc), fontsize=4.0, va='top', family='DejaVu Sans')
        y -= 0.043
ax2.set_title('b  Gene inventory of the element (NR-derived descriptions; full table in Supplementary Table S5)',
              fontsize=7)
ax2.text(0.0, 1.04, 'b', transform=ax2.transAxes, fontsize=9, fontweight='bold')

savefig(fig, os.path.join(OUT, 'FigS4_pl2_annotation', 'FigS4'))
print('figS4 done')
