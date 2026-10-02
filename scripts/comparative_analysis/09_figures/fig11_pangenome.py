# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig11_pangenome.py — Figure 11：泛基因组比较（Venn/UpSet 思想）
6 基因组：smbu06、smbu08、194、IDCC 2105、CX 2-6_2（参考论文比较株）、E. faecium 64/3
方法：pyrodigal 统一预测蛋白 → 全对全 BLASTP（id≥60%、cov≥60%）→ 连通分量聚类 →
      core（6/6）/ accessory（2–5）/ unique（1）；成对共享簇热图。
输出：Figure 11 与 pangenome_summary.tsv（簇级与基因级计数）
"""
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pyrodigal
from p3lib import WORK, IN06, IN08, BLASTP, MAKEDB, read_fasta, write_fasta, run, open_log
from fig_style import plt, savefig, C

OUT = os.path.join(WORK, '09_figures')
FA = os.path.join(WORK, '12_functional_annotation')
TMP = os.path.join(FA, 'tmp_pangenome')
os.makedirs(TMP, exist_ok=True)
LOG = open_log(os.path.join(FA, 'logs', 'fig11_pangenome.log'))


def find_fna(acc):
    for d in [WORK + r'\06_public_panel\downloads\\' + acc,
              WORK + r'\02_taxonomy\refs\\' + acc]:
        for f in glob.glob(d + r'\**\*.fna', recursive=True):
            return f
    return None


SOURCES = [
    ('smbu06', IN06 + r'\smbu06\2.Assembly\smbu06.Complete.genome.fasta'),
    ('smbu08', IN08 + r'\smbu08\2.Assembly\smbu08.Complete.genome.fasta'),
    ('194', WORK + r'\02_taxonomy\refs\GCF_056582645.1\ncbi_dataset\data\GCF_056582645.1\GCF_056582645.1_ASM5658264v1_genomic.fna'),
    ('IDCC2105', find_fna('GCF_023612275.1') or find_fna('GCA_023612275.1')),
    ('CX2-6_2', find_fna('GCF_019343125.1') or find_fna('GCA_019343125.1')),
    ('E_faecium_64-3', find_fna('GCF_001298485.1')),
]
genomes = {}
for name, path in SOURCES:
    if path and os.path.exists(path):
        genomes[name] = path
    else:
        LOG.write('MISSING: %s\n' % name)
names = list(genomes)

# 1) pyrodigal 统一预测
prot_files = {}
gene_index = {}   # (genome, gene_id) -> cluster 用
for n in names:
    fna = genomes[n]
    seqs = read_fasta(fna)
    gf = pyrodigal.GeneFinder(meta=False, min_gene=90)
    prot = []
    gi = 0
    for hdr, s in seqs:
        gf.train(s)
        for g in gf.find_genes(s):
            gi += 1
            aa = g.translate()
            prot.append(('%s_g%05d' % (n, gi), aa))
    pf = os.path.join(TMP, n + '.faa')
    write_fasta(pf, prot)
    prot_files[n] = pf
    LOG.write('%s: %d proteins predicted\n' % (n, len(prot)))

# 2) 全对全 blastp
parent = {}
def find(x):
    while parent.get(x, x) != x:
        parent[x] = parent.get(parent[x], parent[x])
        x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb

edges = 0
for i, a in enumerate(names):
    run([MAKEDB, '-in', prot_files[a], '-dbtype', 'prot', '-out', os.path.join(TMP, a + '_db')], log=None, check=True)
for a in names:
    for b in names:
        if a == b:
            continue
        r = run([BLASTP, '-query', prot_files[a], '-db', os.path.join(TMP, b + '_db'),
                 '-outfmt', '6 qseqid sseqid pident length qlen', '-evalue', '1e-5',
                 '-max_target_seqs', '3', '-num_threads', '6', '-seg', 'no'], log=None)
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            q, s, pid, alen, qlen = line.split('\t')
            if float(pid) >= 60 and int(alen) >= 0.6 * int(qlen):
                union(q, s)
                edges += 1
LOG.write('edges: %d\n' % edges)

# 3) 簇统计
clusters = defaultdict(set)
all_genes = []
for n in names:
    for gid, aa in read_fasta(prot_files[n]):
        root = find(gid)
        clusters[root].add(n)
        all_genes.append((n, gid, root))

# 保存簇→基因组集合
core = {k: v for k, v in clusters.items() if len(v) == len(names)}
acc = {k: v for k, v in clusters.items() if 2 <= len(v) < len(names)}
uniq = {k: v for k, v in clusters.items() if len(v) == 1}
LOG.write('clusters=%d core=%d accessory=%d unique=%d\n' % (len(clusters), len(core), len(acc), len(uniq)))

with open(os.path.join(FA, 'summary', 'pangenome_summary.tsv'), 'w', encoding='utf-8') as f:
    f.write('item\tvalue\n')
    f.write('genomes\t%s\n' % ';'.join(names))
    f.write('total_clusters\t%d\n' % len(clusters))
    f.write('core_clusters(in_all_%d)\t%d\n' % (len(names), len(core)))
    f.write('accessory_clusters(2_to_%d)\t%d\n' % (len(names) - 1, len(acc)))
    f.write('unique_clusters\t%d\n' % len(uniq))
    # 每基因组的基因级构成
    for n in names:
        genes = [(gid, root) for (nn, gid, root) in all_genes if nn == n]
        nc = sum(1 for _, r in genes if len(clusters[r]) == len(names))
        na = sum(1 for _, r in genes if 2 <= len(clusters[r]) < len(names))
        nu = sum(1 for _, r in genes if len(clusters[r]) == 1)
        f.write('genes_in_%s\t%d core=%d accessory=%d unique=%d\n' % (n, len(genes), nc, na, nu))
    # 成对共享簇
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            shared = sum(1 for k, v in clusters.items() if a in v and b in v)
            f.write('shared_%s_vs_%s\t%d\n' % (a, b, shared))

# 4) 绘图

DISPLAY = {'smbu06': 'smbu06', 'smbu08': 'smbu08', '194': '194', 'IDCC2105': 'IDCC 2105',
           'CX2-6_2': 'CX 2-6_2', 'E_faecium_64-3': 'E. faecium'}
SHORT = [DISPLAY.get(n, n) for n in names]

fig = plt.figure(figsize=(7.4, 6.2))
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.42)

# a) 每基因组构成堆叠条
ax = fig.add_subplot(gs[0, 0])
core_n = []
acc_n = []
uniq_n = []
for n in names:
    genes = [(gid, root) for (nn, gid, root) in all_genes if nn == n]
    core_n.append(sum(1 for _, r in genes if len(clusters[r]) == len(names)))
    acc_n.append(sum(1 for _, r in genes if 2 <= len(clusters[r]) < len(names)))
    uniq_n.append(sum(1 for _, r in genes if len(clusters[r]) == 1))
y = np.arange(len(names))
ax.barh(y, core_n, color='#4DAF4A', label='core (all %d genomes)' % len(names), edgecolor='black', lw=0.3)
ax.barh(y, acc_n, left=core_n, color='#FFD92F', label='accessory (2–%d)' % (len(names) - 1), edgecolor='black', lw=0.3)
ax.barh(y, uniq_n, left=[c + a for c, a in zip(core_n, acc_n)], color='#E41A1C', label='unique (1)', edgecolor='black', lw=0.3)
for i, (c, a, u) in enumerate(zip(core_n, acc_n, uniq_n)):
    ax.text(c / 2, i, str(c), ha='center', va='center', fontsize=5.0, color='white')
    ax.text(c + a / 2, i, str(a), ha='center', va='center', fontsize=5.0)
    ax.text(c + a + u + 30, i, str(u), ha='left', va='center', fontsize=5.0, color='#E41A1C')
ax.set_yticks(y)
ax.set_yticklabels(SHORT, fontsize=6)
ax.invert_yaxis()
ax.set_xlabel('number of genes', fontsize=6)
ax.tick_params(labelsize=5.5)
ax.legend(fontsize=5.0, loc='lower right')
ax.set_title('a  Pangenome composition per genome', fontsize=7.5)
ax.text(0.01, 0.97, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

# b) 成对共享热图
ax2 = fig.add_subplot(gs[0, 1])
S = np.zeros((len(names), len(names)))
for i in range(len(names)):
    for j in range(len(names)):
        if i == j:
            S[i, j] = sum(1 for k, v in clusters.items() if names[i] in v)
        else:
            S[i, j] = sum(1 for k, v in clusters.items() if names[i] in v and names[j] in v)
im = ax2.imshow(S, cmap='YlGnBu')
ax2.set_xticks(range(len(names)))
ax2.set_xticklabels(SHORT, rotation=30, ha='right', fontsize=6)
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
         'Comparison panel includes the strains used in the reference-style comparison (IDCC 2105, CX 2-6_2).').format(n=len(clusters)),
         fontsize=5.6, ha='center', va='bottom')
savefig(fig, os.path.join(OUT, 'Fig11_pangenome', 'Fig11'))
print('fig11 done; clusters=%d core=%d' % (len(clusters), len(core)))
