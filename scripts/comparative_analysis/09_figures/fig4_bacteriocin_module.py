# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig4_bacteriocin_module.py — Figure 4：共享质粒细菌素模块（enterocin P-like / 免疫 / hiracin-like）
a) 641–642 位点放大（基因箭头）
b) 679 位点放大
c) 与已知细菌素的对齐（差异位点标注）
d) IS 分布与模块区（密度对照）
e) 候选性质表
"""
import os
import sys
import csv

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p3lib import WORK, read_fasta
from fig_style import plt, savefig, C, gene_arrow
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

OUT = os.path.join(WORK, '09_figures')
ANN = os.path.join(WORK, '04_plasmidome', 'annotation')
SEA = os.path.join(WORK, '05_bacteriocins', 'search')

genes = []
with open(os.path.join(ANN, 'pl1_genes.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['region'] == 'pl1_smbu06':
            r['start'] = int(r['start']); r['end'] = int(r['end'])
            genes.append(r)
mge = []
with open(os.path.join(ANN, 'mge_inventory.tsv'), encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 9 and p[0] == 'smbu06' and p[1] == 'Plasmid1':
            mge.append(dict(start=int(p[3]), end=int(p[4]), desc=(p[7] or p[6])))

KEY = {'smbu06GL002641': ('enterocin P-like', C['bacteriocin']),
       'smbu06GL002642': ('Sakacin-A immunity factor', C['immunity']),
       'smbu06GL002679': ('hiracin JM79-like', C['bacteriocin'])}
# 标注避让偏移（dx kb, dy, 对齐）
LABEL_OFF = {'smbu06GL002641': (-1.2, 1.30, 'right'),
             'smbu06GL002642': (1.0, 1.30, 'left'),
             'smbu06GL002679': (-1.0, -1.20, 'right')}

fig = plt.figure(figsize=(7.2, 6.4))
gs = gridspec.GridSpec(3, 2, height_ratios=[0.8, 0.75, 1.0], hspace=0.52, wspace=0.32)


def locus_panel(ax, lo, hi, title):
    drawn = [(g['start'], g['end'], g['strand'], g['locus_tag'], g['product']) for g in genes
             if g['end'] > lo and g['start'] < hi]
    # 行分配（简单避让：按 strand 两条轨道 + 溢出到扩展轨道）
    for s0, s1, strand, lt, prod in drawn:
        if lt in KEY:
            col = KEY[lt][1]
        elif 'transposase' in prod.lower() or 'integrase' in prod.lower():
            col = C['mge']
        elif prod != 'hypothetical protein':
            col = C['core']
        else:
            col = C['other']
        y = 0.30 if strand == '+' else -0.30
        gene_arrow(ax, s0 / 1000, s1 / 1000, y, strand, col)
        if lt in KEY:
            nm = KEY[lt][0]
            dx, dy, ha = LABEL_OFF.get(lt, (0, 1.05 if y > 0 else -1.05, 'center'))
            ax.annotate('%s\n%s (%s aa)' % (lt.replace('smbu06GL', 'GL'), nm, (s1 - s0 + 1) // 3 - 1),
                        xy=((s0 + s1) / 2000, y), xytext=((s0 + s1) / 2000 + dx, dy),
                        ha=ha, fontsize=5.2, color=col,
                        arrowprops=dict(arrowstyle='-', lw=0.4, color=col))
    ax.set_xlim(lo / 1000, hi / 1000)
    ax.set_ylim(-1.35, 1.45)
    ax.set_yticks([])
    ax.set_xlabel('position in shared plasmid (kb)', fontsize=6)
    ax.tick_params(labelsize=5.5)
    ax.set_title(title, fontsize=7)
    return drawn


ax = fig.add_subplot(gs[0, :])
locus_panel(ax, 34500, 37500, 'a  enterocin P-like locus (GL002641–GL002642, 10-bp intergenic gap)')
ax.text(0.01, 0.94, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

ax2 = fig.add_subplot(gs[1, :])
locus_panel(ax2, 68000, 70600, 'b  hiracin JM79-like locus (GL002679)')
ax2.text(0.01, 0.94, 'b', transform=ax2.transAxes, fontsize=9, fontweight='bold', va='top')

# c) 对齐
def read_align(path):
    with open(path, encoding='utf-8') as f:
        for line in f:
            p = line.rstrip('\n').lstrip('#').split('\t')
            if len(p) >= 9 and p[0].replace('.', '').isdigit():
                return float(p[0]), p[7], p[8]
    return None, None, None


ax3 = fig.add_subplot(gs[2, 0])
ax3.axis('off')
rows = [
    ('GL002641 vs enterocin P (O30434)', *read_align(os.path.join(SEA, 'alignment_GL002641_vs_enterocinP.txt'))),
    ('GL002679 vs hiracin JM79', *read_align(os.path.join(SEA, 'alignment_GL002679_vs_hiracinJM79.txt'))),
]
y = 0.95
for name, pid, qs, ss in rows:
    ax3.text(0.0, y, '%s — %.1f%% identity' % (name, pid), fontsize=6.0, fontweight='bold')
    y -= 0.115
    for label, seq in [('cand', qs), ('hmlg', ss)]:
        ax3.text(0.0, y, label, fontsize=5.2, family='DejaVu Sans Mono', style='italic')
        for j, ch in enumerate(seq):
            same = (j < len(ss)) and (ss[j] == qs[j])
            col = 'black' if same else 'red'
            ax3.text(0.065 + j * 0.0140, y, ch, fontsize=5.4, color=col, family='DejaVu Sans Mono')
        y -= 0.095
    # 标记 mismatch 位置
    mm = [j for j in range(min(len(qs), len(ss))) if qs[j] != ss[j]]
    ax3.text(0.0, y, 'diff at positions: %s' % (', '.join(str(m + 1) for m in mm) or 'none'),
             fontsize=5.2, color='red')
    y -= 0.155
ax3.text(0.0, 1.06, 'c  Alignment with characterised bacteriocins', transform=ax3.transAxes,
         fontsize=7, fontweight='bold', va='bottom')

# d) IS 分布
ax4 = fig.add_subplot(gs[2, 1])
for m in mge:
    ax4.plot([m['start'] / 1000, m['end'] / 1000], [0.0, 0.0], color=C['mge'], lw=1.2, alpha=0.85)
ax4.add_patch(Rectangle((31.073, -0.35), 74.688 - 31.073, 0.7, facecolor=C['bacteriocin'],
                        alpha=0.18, edgecolor=C['bacteriocin'], lw=0.6))
for lt, (nm, col) in [('smbu06GL002641', KEY['smbu06GL002641']), ('smbu06GL002642', KEY['smbu06GL002642']),
                      ('smbu06GL002679', KEY['smbu06GL002679'])]:
    g = [x for x in genes if x['locus_tag'] == lt][0]
    ax4.plot((g['start'] + g['end']) / 2000, 0.0, marker='v', ms=3.5, color=col)
ax4.set_xlim(0, 128.8)
ax4.set_ylim(-0.6, 0.6)
ax4.set_yticks([])
ax4.set_xlabel('shared plasmid position (kb)', fontsize=6)
ax4.tick_params(labelsize=5.5)
ax4.set_title('d  IS elements (grey) across the plasmid', fontsize=7)
ax4.text(0.02, 0.9, 'module region 31.1–74.7 kb:\n0.32 IS/kb (14/43.6 kb)\nplasmid average: 0.31 IS/kb (40/128.8 kb)',
         transform=ax4.transAxes, fontsize=5.4, va='top')

# e) 性质
ax5 = fig.add_subplot(gs[2, 1])
ax5.remove()
ax5 = fig.add_axes([0.55, 0.055, 0.42, 0.0])
ax5.axis('off')

fig.text(0.02, 0.035, 'e  Candidates: GL002641 (54 aa; pI 6.74; GRAVY −0.21) — enterocin P-like, YGNGV→YDNGI variant;  '
         'GL002642 (88 aa) — Sakacin-A immunity factor;  GL002679 (67 aa; pI 9.14; GRAVY +0.31) — hiracin JM79-like (59.3%).',
         fontsize=6.0, va='top')
fig.text(0.02, 0.012, 'Sequence-level homology only; no expression, processing, secretion or activity data are generated in this study.',
         fontsize=5.6, style='italic', color='#444444')

savefig(fig, os.path.join(OUT, 'Fig4_bacteriocin_module', 'Fig4'))
print('fig4 done')
