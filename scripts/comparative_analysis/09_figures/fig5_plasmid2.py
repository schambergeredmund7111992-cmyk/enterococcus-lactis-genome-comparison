# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig5_plasmid2.py — Figure 5：smbu08 特异 77.4 kb 序列（切离态前噬菌体二聚体）
a) 切离模型示意（整合态 → 环化 → 二聚体）
b) Plasmid2 基因图（102 CDS，按模块着色，标注两拷贝）
c) 深度剖面对照（染色体 vs 切离单元；illumina+nanopore）
d) 与 smbu06 染色体区段的对齐（revcomp）
"""
import os
import sys
import csv

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p3lib import WORK, ASM06, ASM08, read_fasta, revcomp
from fig_style import plt, savefig, C, gene_arrow
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyArrowPatch, Arc

OUT = os.path.join(WORK, '09_figures')
ANN = os.path.join(WORK, '04_plasmidome', 'annotation')
MP = os.path.join(WORK, '01_assembly_qc', 'mapping')

pl2_genes = []
with open(os.path.join(ANN, 'pl2_annotation_full.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        r['start'] = int(r['start']); r['end'] = int(r['end'])
        pl2_genes.append(r)

fig = plt.figure(figsize=(7.2, 6.6))
gs = gridspec.GridSpec(3, 2, height_ratios=[0.8, 1.05, 0.95], hspace=0.62, wspace=0.3)

# a) 切离模型
ax = fig.add_subplot(gs[0, :])
ax.axis('off')
ax.set_xlim(0, 10)
ax.set_ylim(0, 2.4)
# 染色体
ax.plot([0.2, 3.2], [1.6, 1.6], color=C['smbu06'], lw=2.6)
ax.add_patch(Rectangle((1.35, 1.45), 0.55, 0.3, facecolor=C['prophage'], edgecolor='none'))
ax.text(1.62, 2.0, 'smbu06 chromosome\n(integrated prophage, 38,719 bp)', fontsize=5.6, ha='center', color=C['prophage'])
# 箭头
ax.add_patch(FancyArrowPatch((3.5, 1.6), (4.6, 1.6), arrowstyle='-|>', mutation_scale=8, lw=0.8))
ax.text(4.05, 1.85, 'excision\n(int/att)', fontsize=5.4, ha='center')
# 环
th = np.linspace(0, 2 * np.pi, 100)
ax.plot(5.6 + 0.55 * np.cos(th), 1.6 + 0.55 * np.sin(th) * 0.5, color=C['plasmid2'], lw=2.0)
ax.text(5.6, 0.75, 'excised unit U\n38,719 bp circle', fontsize=5.4, ha='center', color=C['plasmid2'])
# 二聚体
ax.add_patch(FancyArrowPatch((6.6, 1.6), (7.6, 1.6), arrowstyle='-|>', mutation_scale=8, lw=0.8))
ax.text(7.1, 1.85, 'replication /\ntandem duplication', fontsize=5.4, ha='center')
for i in range(2):
    pts = [(8.2 + i * 0.9, 1.1), (8.2 + i * 0.9, 2.1), (9.0 + i * 0.9, 2.1), (9.0 + i * 0.9, 1.1)]
    ax.plot([p[0] for p in pts[:2]], [p[1] for p in pts[:2]], color=C['plasmid2'], lw=2.0)
    ax.plot([p[0] for p in pts[2:]], [p[1] for p in pts[2:]], color=C['plasmid2'], lw=2.0)
    ax.annotate('', xy=pts[2], xytext=pts[1], arrowprops=dict(arrowstyle='-|>', color=C['plasmid2'], lw=1.0))
    ax.plot([pts[0][0], pts[3][0]], [pts[0][1], pts[3][1]], color=C['plasmid2'], lw=1.0, ls=':')
ax.text(8.65, 0.75, 'smbu08 Plasmid2\n77,437 bp (2,600× depth\nin source culture)', fontsize=5.4,
        ha='center', color=C['plasmid2'])
ax.text(0.02, 0.95, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold')

# b) pl2 基因图
ax2 = fig.add_subplot(gs[1, :])
MODC = {'lysogeny': '#7570B3', 'DNA_packaging': '#1B9E77', 'head': '#D95F02',
        'tail': '#E7298A', 'lysis': '#E6AB02', 'recombination_DNA': '#66A2E5',
        'regulation': '#A6761D', 'hypothetical': C['other']}
for r in pl2_genes:
    col = MODC.get(r['functional_module'], C['other'])
    y = 0.30 if r['strand'] == '+' else -0.30
    gene_arrow(ax2, r['start'] / 1000, r['end'] / 1000, y, r['strand'], col, fontsize=4.2)
ax2.axvline(38.719, color='black', lw=0.7, ls='--')
ax2.text(38.719, 1.05, 'unit boundary\n(1-bp deletion in 2nd copy)', fontsize=5.0, ha='center')
ax2.text(19, 0.95, 'copy 1 (38,719 bp)', fontsize=5.6, ha='center')
ax2.text(58, 0.95, 'copy 2 (38,718 bp)', fontsize=5.6, ha='center')
ax2.set_xlim(-1, 78.5)
ax2.set_ylim(-1.0, 1.35)
ax2.set_yticks([])
ax2.set_xlabel('Plasmid2 position (kb)', fontsize=6)
ax2.tick_params(labelsize=5.5)
import matplotlib.patches as mpatches
ax2.legend(handles=[mpatches.Patch(facecolor=c, label=k) for k, c in MODC.items()],
           loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=8, fontsize=4.6,
           handlelength=1.0, columnspacing=0.9, handletextpad=0.3)
ax2.text(0.01, 0.93, 'b', transform=ax2.transAxes, fontsize=9, fontweight='bold')

# c) 深度剖面
ax3 = fig.add_subplot(gs[2, 0])
from p3lib import read_fasta as rf
ref = rf(os.path.join(MP, 'combined_ref.fa'))
offs = {}
pos = 0
for name, s in ref:
    offs[name] = pos
    pos += len(s)
L2 = {}
for name, s in ref:
    L2[name] = len(s)
d8i = np.load(os.path.join(MP, 'depth_smbu08_illumina.npz'))['depth']
d8n = np.load(os.path.join(MP, 'depth_smbu08_nanopore.npz'))['depth']


def get(depth, name):
    return depth[offs[name]:offs[name] + L2[name]]


w = 500
pl2_i = get(d8i, 'smbu08|Plasmid2')
chr_i = get(d8i, 'smbu08|Chromosome1')
x = np.arange(0, len(pl2_i), w) / 1000
ax3.plot(x, [pl2_i[i:i + w].mean() for i in range(0, len(pl2_i), w)], color=C['plasmid2'], lw=0.9,
         label='Plasmid2 (excised unit), Illumina')
xc = np.arange(0, len(chr_i), w) / 1000
ax3.plot(xc, [chr_i[i:i + w].mean() for i in range(0, len(chr_i), w)], color=C['smbu08'], lw=0.7,
         label='chromosome, Illumina')
ax3.set_xlabel('position (kb)', fontsize=6)
ax3.set_ylabel('depth (×)', fontsize=6)
ax3.tick_params(labelsize=5.5)
ax3.legend(fontsize=5.0, loc='upper right')
ax3.set_title('Read depth: excised unit ≫ chromosome', fontsize=7)
ax3.text(0.02, 0.92, 'c', transform=ax3.transAxes, fontsize=9, fontweight='bold')

# d) 与 smbu06 染色体比对（脚本内运行 blastn，结果亦写入模块 03）
ax4 = fig.add_subplot(gs[2, 1])
sq = read_fasta(ASM06['Chromosome1'])[0][1]
from p3lib import BLASTN as _B
hsp_file = os.path.join(WORK, '03_replicon_compare', 'pl2_vs_chr06.hsp.tsv')
if not os.path.exists(hsp_file):
    from p3lib import run as _run
    pl2fa = os.path.join(WORK, '03_replicon_compare', 'alignments', 'smbu08_pl2.fa')
    r = _run([_B, '-task', 'dc-megablast', '-query', pl2fa, '-db', os.path.join(WORK, '03_replicon_compare', 'alignments', 'smbu06_chr_db'),
              '-outfmt', '6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
              '-evalue', '1e-5', '-max_hsps', '50'])
    with open(hsp_file, 'w') as fh:
        fh.write('sseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\n')
        fh.write(r.stdout)
hsp = []
with open(hsp_file, encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 9:
            continue
        hsp.append((int(p[5]), int(p[6]), int(p[7]), int(p[8]), float(p[1])))
for qs, qe, ss, se, pid in hsp:
    if abs(qe - qs) < 1000:
        continue
    ax4.plot([qs / 1000, qe / 1000], [ss / 1e6, se / 1e6], color=C['plasmid2'], lw=1.2)
ax4.axvspan(1255.783, 1294.501, color=C['prophage'], alpha=0.25)
ax4.text(1275, 0.6, 'prophage region\nin smbu06 chromosome (Mb)', fontsize=5.2, color=C['prophage'])
ax4.set_xlabel('Plasmid2 position (kb)', fontsize=6)
ax4.set_ylabel('smbu06 chromosome (Mb)', fontsize=6)
ax4.tick_params(labelsize=5.5)
ax4.set_title('Plasmid2 maps to the prophage region (reverse strand)', fontsize=7)
ax4.text(0.02, 0.92, 'd', transform=ax4.transAxes, fontsize=9, fontweight='bold')

savefig(fig, os.path.join(OUT, 'Fig5_plasmid2', 'Fig5'))
print('fig5 done')
