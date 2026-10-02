# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""fig9_safety.py — Figure 9：安全性相关注释（CARD / ARDB / VFDB，如实报告）"""
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


fig = plt.figure(figsize=(7.4, 6.6))
gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.0], hspace=0.42, wspace=0.34)

# a) CARD strict hits
ax = fig.add_subplot(gs[0, 0])
card = [r for r in load('card_hits.tsv') if r['strain'] == 'smbu06']
card.sort(key=lambda r: float(r['best_identities']))
names = ['%s (%s)' % (r['best_aro'], r['gene'].replace('smbu06GL', 'GL')) for r in card]
idents = [float(r['best_identities']) for r in card]
mec = [r['resistance_mechanism'] for r in card]
MECCOL = {'antibiotic efflux': '#FC8D62', 'antibiotic inactivation': '#E41A1C',
          'antibiotic target protection': '#377EB8'}
ax.barh(range(len(card)), idents, color=[MECCOL.get(m, '#999999') for m in mec], edgecolor='black', lw=0.3)
ax.set_yticks(range(len(card)))
ax.set_yticklabels(names, fontsize=5.6)
for i, (v, r) in enumerate(zip(idents, card)):
    ax.text(v + 1, i, '%.1f%%' % v, va='center', fontsize=5.0)
ax.set_xlim(0, 112)
ax.set_xlabel('identity to CARD reference (%)', fontsize=6)
ax.tick_params(labelsize=5.5)
ax.set_title('a  CARD strict hits (4 genes, both strains)', fontsize=7.5)
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(facecolor=v, label=k) for k, v in MECCOL.items()], fontsize=4.6,
          loc='upper left', bbox_to_anchor=(0.22, 1.0))
ax.text(0.01, 0.97, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

# b) ARDB: identity vs database threshold
ax2 = fig.add_subplot(gs[0, 1])
ardb = [r for r in load('ardb_hits.tsv') if r['strain'] == 'smbu06']
ids = [float(r['identity']) for r in ardb]
mins = [float(r['min_identity']) for r in ardb if r['min_identity'] not in ('', '--')]
thr = mins[0] if mins else 80
ax2.scatter(range(len(ids)), ids, color=C['mge'], s=14, zorder=3, label='ARDB best hits')
ax2.axhline(thr, color=C['smbu08'], lw=1.0, ls='--', label='database Min_Identity threshold (%g%%)' % thr)
ax2.set_xlabel('ARDB hit index (n=%d)' % len(ids), fontsize=6)
ax2.set_ylabel('identity (%)', fontsize=6)
ax2.set_ylim(0, 100)
ax2.tick_params(labelsize=5.5)
ax2.legend(fontsize=4.8, loc='upper left')
ax2.set_title('b  ARDB hits are below the database threshold', fontsize=7.5)
ax2.text(0.01, 0.97, 'b', transform=ax2.transAxes, fontsize=9, fontweight='bold', va='top')

# c) VFDB categories
ax3 = fig.add_subplot(gs[1, 0])
vcat = {}
with open(os.path.join(SUM, 'vfdb_hits.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['strain'] == 'smbu06':
            cat = r['description'].split(']')[0].lstrip('[').strip() if r['description'] else 'unknown'
            vcat[cat] = vcat.get(cat, 0) + 1
top = sorted(vcat.items(), key=lambda x: -x[1])[:8][::-1]
ax3.barh(range(len(top)), [v for _, v in top], color='#B3B3B3', edgecolor='black', lw=0.3)
ax3.set_yticks(range(len(top)))
ax3.set_yticklabels([k[:38] for k, _ in top], fontsize=5.0)
ax3.set_xlabel('VFDB hits (smbu06, n=%d; all "Predicted")' % sum(vcat.values()), fontsize=6)
ax3.tick_params(labelsize=5.5)
ax3.set_title('c  VFDB category counts', fontsize=7.5)
ax3.text(0.01, 0.97, 'c', transform=ax3.transAxes, fontsize=9, fontweight='bold', va='top')

# d) 说明框
ax4 = fig.add_subplot(gs[1, 1])
ax4.axis('off')
vf = [r for r in load('vfdb_hits.tsv') if r['strain'] == 'smbu06']
vmax = max(float(r['identity']) for r in vf)
txt = (
    'Safety-relevant annotations (delivered provider tables):\n\n'
    "• CARD strict: four resistance genes at high identity were\n"
    "  detected in BOTH strains (a): AAC(6')-Ii 99.5% (aminoglycoside\n"
    '  inactivation), msrC 97.2%, eatAv 96.0% (target protection),\n'
    '  efrA 82.6% (efflux). All four are located on the chromosome\n'
    '  (not on the shared plasmid or the prophage).\n'
    '• ARDB: {n_ar} hits at 41–{mx}% identity — all below the database\n'
    '  Min_Identity threshold ({thr}%), i.e. not counted as resistance\n'
    '  determinants (b).\n'
    '• VFDB: {n_vf} hits, all flagged "Predicted"; the majority match\n'
    '  housekeeping proteins (max identity {vmax}%) (c).\n\n'
    'Interpretation: the genomes carry several known resistance\n'
    'genes; this should be stated precisely in any safety discussion\n'
    'and differs from a "no resistance genes" claim.'
).format(n_ar=len(ardb), mx=int(max(ids)), thr=int(thr), n_vf=sum(vcat.values()), vmax=round(vmax, 1))
ax4.text(0.0, 1.0, txt, fontsize=5.8, va='top', linespacing=1.5)

savefig(fig, os.path.join(OUT, 'Fig9_safety', 'Fig9'))
print('fig9 done')
