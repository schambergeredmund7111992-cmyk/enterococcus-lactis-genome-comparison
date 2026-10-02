# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""fig7_functional.py — Figure 7：功能基因组景观（KEGG / COG / GO）"""
import csv
import os
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


L1COL = {'Metabolism': '#66C2A5', 'Genetic Information Processing': '#FC8D62',
         'Environmental Information Processing': '#8DA0CB', 'Cellular Processes': '#E78AC3',
         'Human Diseases': '#A6D854', 'Organismal Systems': '#FFD92F',
         'Unclassified': '#B3B3B3', 'Drug Development': '#B3B3B3'}

fig = plt.figure(figsize=(7.4, 9.4))
gs = gridspec.GridSpec(3, 1, height_ratios=[1.25, 1.25, 1.0], hspace=0.36)

# a) KEGG Level2 top-15
ax = fig.add_subplot(gs[0])
k2 = [r for r in load('kegg_level2.tsv') if r['strain'] == 'smbu06']
k2.sort(key=lambda r: -int(r['n_genes']))
top = k2[:15][::-1]
ax.barh(range(len(top)), [int(r['n_genes']) for r in top],
        color=[L1COL.get(r['level1'], '#999999') for r in top], edgecolor='black', lw=0.3)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([('%s — %s' % (r['level1'][:22], r['level2']))[:52] for r in top], fontsize=5.0)
ax.set_xlabel('number of genes (smbu06)', fontsize=6)
ax.tick_params(labelsize=5.5)
ax.set_title('a  KEGG pathway classification (top 15 of %d Level-2 categories)' % len(k2), fontsize=7.5)
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(facecolor=v, label=k) for k, v in L1COL.items() if k != 'Unclassified'],
          fontsize=4.4, loc='lower right', ncol=2)
ax.text(0.99, 0.62, 'total KEGG-annotated genes: %d' % sum(int(r['n_genes']) for r in load('kegg_level2.tsv') if r['strain'] == 'smbu06'),
        transform=ax.transAxes, ha='right', fontsize=5.4)

# b) COG
ax2 = fig.add_subplot(gs[1])
cog = [r for r in load('cog.tsv') if r['strain'] == 'smbu06']
cog.sort(key=lambda r: r['cog_code'])
FCOL = {'METABOLISM': '#66C2A5', 'CELLULAR': '#FC8D62', 'INFORMATION': '#8DA0CB'}
ax2.bar(range(len(cog)), [int(r['n_genes']) for r in cog],
        color=[FCOL.get(r['first_class'].upper()[:11], '#999999') for r in cog],
        edgecolor='black', lw=0.3)
ax2.set_xticks(range(len(cog)))
ax2.set_xticklabels([r['cog_code'] for r in cog], fontsize=5.2)
for i, r in enumerate(cog):
    ax2.text(i, int(r['n_genes']) + 1.5, r['cog_code'], fontsize=4.0, ha='center')
ax2.set_ylabel('number of genes', fontsize=6)
ax2.tick_params(labelsize=5.5)
ax2.set_title('b  COG functional categories (smbu06; %d annotated genes)' % sum(int(r['n_genes']) for r in cog), fontsize=7.5)
ax2.legend(handles=[mpatches.Patch(facecolor=v, label=k) for k, v in FCOL.items()], fontsize=4.6, loc='upper right')

# c) GO level2
ax3 = fig.add_subplot(gs[2])
go = load('go_level2.tsv')
go = [r for r in go if r['strain'] == 'smbu06']
onto_col = {'biological_process': '#66C2A5', 'cellular_component': '#FC8D62', 'molecular_function': '#8DA0CB'}
onto_lab = {'biological_process': 'BP', 'cellular_component': 'CC', 'molecular_function': 'MF'}
x = 0
xt = []
for onto in ['biological_process', 'cellular_component', 'molecular_function']:
    rows = sorted([r for r in go if r['ontology'] == onto], key=lambda r: -int(r['n_genes']))[:10][::-1]
    for r in rows:
        ax3.barh(x, int(r['n_genes']), color=onto_col[onto], edgecolor='black', lw=0.25)
        xt.append(('%s  %s' % (onto_lab[onto], r['go_class']))[:46])
        x += 1
ax3.set_yticks(range(len(xt)))
ax3.set_yticklabels(xt, fontsize=4.4)
ax3.set_xlabel('number of genes (smbu06, top 10 per ontology)', fontsize=6)
ax3.tick_params(labelsize=5.5)
ax3.set_title('c  GO classification (Level-2 terms)', fontsize=7.5)
ax3.legend(handles=[mpatches.Patch(facecolor=v, label=k) for k, v in onto_col.items()], fontsize=4.6, loc='lower right')

fig.text(0.5, 0.005, 'Functional categories from the sequencing provider annotation tables (KEGG/COG/GO); smbu08 values are identical except for ≤2 genes per category.\n'
         'COG codes: full descriptions in Supplementary Table S7.', fontsize=5.6, ha='center', va='bottom')
savefig(fig, os.path.join(OUT, 'Fig7_functional', 'Fig7'))
print('fig7 done')
