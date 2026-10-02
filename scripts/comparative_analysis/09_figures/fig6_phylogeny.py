# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig6_phylogeny.py — Figure 6：宿主核心树 + 分布/基因/质粒热图；质粒内容相关树
数据：core_tree.treefile（IQ-TREE）；gene_presence_matrix；pl1_distribution；nisin_panel_status；
      plasmid_gene_content.nwk（Jaccard 平均连接）
"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p3lib import WORK
from fig_style import plt, savefig, C
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

# vendored Biopython（Newick 解析）
sys.path.insert(0, r'<PROJECT_ROOT>\分析与初稿\tools\python_packages')
from Bio import Phylo
import io as _io

OUT = os.path.join(WORK, '09_figures')
TRE = os.path.join(WORK, '07_phylogenomics', 'trees')
PAN = os.path.join(WORK, '06_public_panel', 'analysis')

tree_path = os.path.join(TRE, 'core_tree_alrt.treefile')
if not os.path.exists(tree_path):
    tree_path = os.path.join(TRE, 'core_tree.treefile')
if not os.path.exists(tree_path):
    print('NO treefile yet — run build_tree first')
    sys.exit(1)
tree = Phylo.read(tree_path, 'newick')

_dedup_meta_pre = {m['accession'] for m in csv.DictReader(
    open(os.path.join(WORK, '06_public_panel', 'metadata', 'panel_metadata.tsv'), encoding='utf-8'), delimiter='\t')}
KEEP = _dedup_meta_pre | {'smbu06_project', 'smbu08_project', 'E_lactis_194',
                          'L_lactis_14B4', 'L_lactis_MG1363', 'L_lactis_IL1403'}
_drop = [t.name for t in tree.get_terminals() if t.name not in KEEP]
for name in _drop:
    try:
        tree.prune(name)
    except Exception:
        pass
tree.ladderize()
print('tree leaves after prune:', len(tree.get_terminals()))

# 数据
meta = {m['accession']: m for m in csv.DictReader(
    open(os.path.join(WORK, '06_public_panel', 'metadata', 'panel_metadata.tsv'), encoding='utf-8'), delimiter='\t')}
pres = {}
for r in csv.DictReader(open(os.path.join(PAN, 'gene_presence_matrix.tsv'), encoding='utf-8'), delimiter='\t'):
    pres.setdefault(r['accession'], {})[r['gene']] = r
dist = {r['accession']: r for r in csv.DictReader(open(os.path.join(PAN, 'pl1_distribution.tsv'), encoding='utf-8'), delimiter='\t')}
nimap = {r['accession']: r['cluster_call'] for r in csv.DictReader(open(os.path.join(PAN, 'nisin_panel_status.tsv'), encoding='utf-8'), delimiter='\t')}

SPECIES_COL = {'Enterococcus faecium': '#E41A1C', 'Enterococcus faecalis': '#377EB8',
               'Enterococcus lactis': '#4DAF4A', 'Enterococcus hirae': '#984EA3',
               'Enterococcus mundtii': '#FF7F00', 'Enterococcus durans': '#A65628',
               'Lactococcus lactis subsp. lactis': '#999999'}
DEFAULT_SPECIES = '#DDDDDD'
CLS_COL = {'full_length_near_identical': '#2166AC', 'full_length_divergent': '#4393C3',
           'partial_backbone': '#92C5DE', 'fragment': '#D1E5F0', 'trace': '#F7F7F7',
           'no_plasmid_hit': '#FFFFFF'}


def leaf_name(t):
    if t.name:
        return t.name
    return ''


fig = plt.figure(figsize=(7.4, 10.0))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[3.4, 1.2], hspace=0.18, wspace=0.32)

# a) 树 + 热图
ax = fig.add_subplot(gs[0, 0])
leaves = tree.get_terminals()
n = len(leaves)
maxx = max((c.branch_length or 0) for c in tree.find_clades())
def draw(clade, x, y):
    """返回该 clade 的 y 中心；对 ≥15 叶的节点标注 UFBoot 支持度"""
    if not clade.clades:
        return y
    x2 = x + (clade.branch_length or 0)
    ys_ = []
    yy = y
    for ch in clade.clades:
        yc = draw(ch, x2, yy)
        ys_.append(yc)
        yy = yc + 1
    ymid = (ys_[0] + ys_[-1]) / 2
    ax.plot([x, x2], [ys_[0], ys_[0]], color='black', lw=0.25)
    ax.plot([x, x2], [ys_[-1], ys_[-1]], color='black', lw=0.25)
    ax.plot([x2, x2], [ys_[0], ys_[-1]], color='black', lw=0.25)
    try:
        nleaf = len(clade.get_terminals())
    except Exception:
        nleaf = 0
    # 仅标注大分支（≥30 叶）且非近末端节点，避免梳状区标签互叠
    if clade.confidence is not None and nleaf >= 30 and x2 < 0.80 * maxx:
        ax.text(x2, ymid + 0.8, '%d' % round(clade.confidence), fontsize=3.0,
                va='bottom', ha='center', color='#333333')
    return (ys_[0] + ys_[-1]) / 2


draw(tree.root, 0.0, 0)
# 叶标签与热图
x0 = maxx * 1.02
colw = maxx * 0.085
cols = ['pl1', 'ent', 'imm', 'hir', 'nisin', 'sp']
head = ['plasmid', '641', '642', '679', 'nisin', 'species']
for ci, h in enumerate(head):
    ax.text(x0 + colw * (ci + 0.5), -1.6, h, fontsize=5.0, ha='center', rotation=0)
for i, t in enumerate(reversed(leaves)):
    name = leaf_name(t)
    acc = name
    y = i
    lab = name + ('*' if name in ('smbu06_project', 'smbu08_project') else '')
    ax.text(x0 + colw * 6 + 0.02, y, lab, fontsize=3.6, va='center')
    # pl1 class
    d = dist.get(acc, {})
    cls = d.get('class', 'no_plasmid_hit')
    ax.add_patch(plt.Rectangle((x0, y - 0.45), colw, 0.9, facecolor=CLS_COL.get(cls, 'white'), edgecolor='none'))
    pg = pres.get(acc, {})
    for ci, g in enumerate(['enterocinP_like', 'sakacinA_immunity', 'hiracinJM79_like'], start=1):
        v = pg.get(g, {})
        on = v.get('presence_70_70') == 'yes'
        soft = v.get('presence_60_70') == 'yes' and not on
        col = C['bacteriocin'] if on else ('#F5B7B1' if soft else 'white')
        ax.add_patch(plt.Rectangle((x0 + colw * ci, y - 0.45), colw, 0.9, facecolor=col,
                                   edgecolor='#BBBBBB', lw=0.1))
    # nisin
    nc = nimap.get(acc, 'not_detected')
    ax.add_patch(plt.Rectangle((x0 + colw * 4, y - 0.45), colw, 0.9,
                               facecolor=C['core'] if nc == 'complete_cluster' else 'white',
                               edgecolor='#BBBBBB', lw=0.1))
    # species
    sp = meta.get(acc, {}).get('organism', '')
    ax.add_patch(plt.Rectangle((x0 + colw * 5, y - 0.45), colw, 0.9,
                               facecolor=SPECIES_COL.get(sp, DEFAULT_SPECIES), edgecolor='none'))
    # 高亮两株
    if name in ('smbu06_project', 'smbu08_project'):
        ax.add_patch(plt.Rectangle((0, y - 0.5), x0 + colw * 6, 1.0, facecolor='#FFFF00', alpha=0.25, zorder=0))
ax.set_xlim(-maxx * 0.02, x0 + colw * 6 + maxx * 0.32)
ax.set_ylim(-2.6, n + 1.2)
ax.axis('off')
ax.text(0.0, 1.004, 'a  Core-genome ML tree (72 markers, LG+G, UFBoot 1000) with distribution heatmap',
        transform=ax.transAxes, fontsize=7, fontweight='bold', va='bottom')
legend_items = [Patch(facecolor=v, label=k.replace('_', ' ')) for k, v in CLS_COL.items() if k != 'no_plasmid_hit']
legend_items += [Patch(facecolor=v, label=k) for k, v in list(SPECIES_COL.items())[:4]]
ax.legend(handles=legend_items, loc='lower right', fontsize=4.2, ncol=2, bbox_to_anchor=(1.02, -0.02))

# b) 质粒内容树
ax2 = fig.add_subplot(gs[0, 1])
ptree_path = os.path.join(TRE, 'plasmid_gene_content.nwk')
if os.path.exists(ptree_path):
    pt = Phylo.read(ptree_path, 'newick')
    pt.ladderize()
    pleaves = pt.get_terminals()
    dists = {}
    mat = {}
    with open(os.path.join(TRE, 'plasmid_gene_content_matrix.tsv'), encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            mat[r['accession']] = float(r['distance_to_smbu06'])
            dists[r['accession']] = int(r['n_genes'])
    # 只画 top 45 最相关 + 两项目株
    order = sorted([l.name for l in pleaves if l.name not in ('smbu06_project', 'smbu08_project')],
                   key=lambda a: mat.get(a, 1.0))[:45]
    keep = set(order) | {'smbu06_project', 'smbu08_project'}
    pt2 = pt
    for l in pleaves:
        if l.name not in keep:
            l.name = None

    def draw2(clade, x, y):
        chs = [c for c in clade.clades]
        if not chs:
            return y
        x2 = x + (clade.branch_length or 0)
        ys_ = []
        yy = y
        for ch in chs:
            yc = draw2(ch, x2, yy)
            ys_.append(yc)
            yy = yc + 1
        ax2.plot([x, x2], [ys_[0], ys_[0]], color='black', lw=0.25)
        ax2.plot([x, x2], [ys_[-1], ys_[-1]], color='black', lw=0.25)
        ax2.plot([x2, x2], [ys_[0], ys_[-1]], color='black', lw=0.25)
        return (ys_[0] + ys_[-1]) / 2

    draw2(pt.root, 0.0, 0)
    nameds = [t for t in pt.get_terminals() if t.name]
    maxx2 = max((c.branch_length or 0) for c in pt.find_clades())
    for i, t in enumerate(reversed(nameds)):
        col = C['smbu08'] if t.name.startswith('smbu') else '#666666'
        ax2.text(maxx2 * 1.03, i, '%s (%d genes)' % (t.name, dists.get(t.name, 0)),
                 fontsize=3.4, va='center', color=col)
    ax2.set_xlim(-maxx2 * 0.02, maxx2 * 1.75)
    ax2.set_ylim(-1, len(nameds) + 1)
    ax2.axis('off')
    ax2.text(0.0, 1.004, 'b  Shared-plasmid gene-content relatedness (45 closest genomes, Jaccard/average linkage)',
             transform=ax2.transAxes, fontsize=6.6, fontweight='bold', va='bottom')

# c) 分布摘要条
ax3 = fig.add_subplot(gs[1, :])
classes = ['full_length_near_identical', 'partial_backbone', 'fragment', 'trace', 'no_plasmid_hit']
counts = {c: sum(1 for a in meta if dist.get(a, {}).get('class', 'no_plasmid_hit') == c) for c in classes}
bar = [counts[c] for c in classes]
labels = ['exact\n(≥99% id)', 'partial backbone\n(≥30% cov)', 'fragment\n(≥5%)', 'trace', 'no match']
cols = [CLS_COL[c] for c in classes]
b = ax3.bar(range(len(bar)), bar, color=cols, edgecolor='black', lw=0.4)
for i, v in enumerate(bar):
    ax3.text(i, v + 1.5, str(v), ha='center', fontsize=6)
ax3.set_xticks(range(len(labels)))
ax3.set_xticklabels(labels, fontsize=5.6)
ax3.set_ylabel('genomes (of 159)', fontsize=6)
ax3.tick_params(labelsize=5.5)
ax3.set_title('c  Shared-plasmid distribution in the deduplicated public panel', fontsize=7)
for i, g in enumerate(['641', '642', '679']):
    pass
ax3.text(0.99, 0.85, 'cassette carriers: 17 (641) / 9 (642) / 11 (679)\nnisin clusters: 5 (the added positive controls)',
         transform=ax3.transAxes, ha='right', fontsize=5.4)

savefig(fig, os.path.join(OUT, 'Fig6_phylogeny', 'Fig6'))
print('fig6 done; host leaves:', n)
