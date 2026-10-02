# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""fig8_cazy.py — Figure 8：CAZy 碳水化合物活性酶谱（参考论文 Fig5 风格）"""
import csv
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK
from fig_style import plt, savefig, C
import matplotlib.gridspec as gridspec

OUT = os.path.join(WORK, '09_figures')
SUM = os.path.join(WORK, '12_functional_annotation', 'summary')


def load(name):
    with open(os.path.join(SUM, name), encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


fig = plt.figure(figsize=(7.4, 6.4))
gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.15], hspace=0.42, wspace=0.3)

# a) 5-class bars (both strains)
ax = fig.add_subplot(gs[0, 0])
c5 = load('cazy_5class.tsv')
classes = ['GH', 'GT', 'CBM', 'CE', 'AA', 'PL']
x = range(len(classes))
w = 0.38
for i, s in enumerate(['smbu06', 'smbu08']):
    row = [r for r in c5 if r['strain'] == s][0]
    vals = [int(row[c]) for c in classes]
    ax.bar([xx + (i - 0.5) * w for xx in x], vals, width=w,
           color=C['smbu06'] if s == 'smbu06' else C['smbu08'], label=s, edgecolor='black', lw=0.3)
    for xx, v in zip(x, vals):
        ax.text(xx + (i - 0.5) * w, v + 0.8, str(v), ha='center', fontsize=4.6)
ax.set_xticks(list(x))
ax.set_xticklabels(['GH', 'GT', 'CBM', 'CE', 'AA', 'PL'], fontsize=6)
ax.set_ylabel('gene counts', fontsize=6)
ax.tick_params(labelsize=5.5)
ax.legend(fontsize=5)
ax.set_title('a  CAZy classes in both genomes', fontsize=7.5)
ax.text(0.01, 0.97, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

# b) donut of class distribution (smbu06)
ax2 = fig.add_subplot(gs[0, 1])
row = [r for r in c5 if r['strain'] == 'smbu06'][0]
vals = [int(row[c]) for c in classes]
from matplotlib.colors import LinearSegmentedColormap
cols = ['#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', '#B15928']
wedges, _, autot = ax2.pie(vals, colors=cols, autopct='%1.0f%%', pctdistance=0.78,
                           textprops=dict(fontsize=5.2), startangle=90,
                           wedgeprops=dict(width=0.42, edgecolor='white', lw=0.5))
ax2.legend(wedges, ['%s (%d)' % (c, v) for c, v in zip(classes, vals)], fontsize=5.0,
           loc='center left', bbox_to_anchor=(0.92, 0.5))
ax2.set_title('b  Distribution (smbu06)', fontsize=7.5)
ax2.text(-1.25, 1.15, 'b', fontsize=9, fontweight='bold')

# c) top families
ax3 = fig.add_subplot(gs[1, :])
fam = load('cazy_families.tsv')
f6 = {}
for r in fam:
    if r['strain'] == 'smbu06':
        f6[r['family']] = int(r['n_genes'])
top = sorted(f6.items(), key=lambda x: -x[1])[:18][::-1]


def cl(fam):
    m = re.match(r'([A-Z]+)', fam)
    return {'GH': '#1F78B4', 'GT': '#33A02C', 'CBM': '#E31A1C', 'CE': '#FF7F00',
            'AA': '#6A3D9A', 'PL': '#B15928'}.get(m.group(1) if m else '', '#999999')


ax3.bar(range(len(top)), [v for _, v in top], color=[cl(k) for k, _ in top], edgecolor='black', lw=0.3)
ax3.set_xticks(range(len(top)))
ax3.set_xticklabels([k for k, _ in top], fontsize=5.2, rotation=38, ha='right')
ax3.set_ylabel('gene counts', fontsize=6)
ax3.tick_params(labelsize=5.5)
ax3.set_title('c  Top CAZy families (smbu06): glycoside hydrolases (GH) and glycosyltransferases (GT) dominate', fontsize=7.5)
import matplotlib.patches as mpatches
ax3.legend(handles=[mpatches.Patch(facecolor=v, label=k) for k, v in
                    [('GH', '#1F78B4'), ('GT', '#33A02C'), ('CBM', '#E31A1C'), ('CE', '#FF7F00')]],
           fontsize=5, loc='upper right')

fig.text(0.5, 0.005, 'CAZy annotation from the sequencing provider pipeline (dbCAN-style; see Methods for thresholds). '
         'GH = glycoside hydrolases; GT = glycosyltransferases; CBM = carbohydrate-binding modules; CE = carbohydrate esterases; AA = auxiliary activities; PL = polysaccharide lyases.',
         fontsize=5.6, ha='center', va='bottom')
savefig(fig, os.path.join(OUT, 'Fig8_cazy', 'Fig8'))
print('fig8 done')
