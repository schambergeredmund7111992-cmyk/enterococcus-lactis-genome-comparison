# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig2_chromosome.py — Figure 2：染色体比较（近等基因背景 + 38.7 kb 前噬菌体缺失 + 1 SNP）
a) 全染色体 dotplot（megablast HSP）
b) 共线性示意图（两株线性图，缺失区/att/SNP 标注）
c) 连接区 att 结构（attL vs attJ(==attR) 差异位点图）
d) 统计摘要框
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import subprocess
from p3lib import WORK, ASM06, ASM08, BLASTN, MAKEDB, read_fasta, write_fasta, run
from fig_style import plt, savefig, C, gene_arrow
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

OUT = os.path.join(WORK, '09_figures')
ALN = os.path.join(WORK, '03_replicon_compare', 'alignments')

s = read_fasta(ASM06['Chromosome1'])[0][1]
q = read_fasta(ASM08['Chromosome1'])[0][1]
Ls, Lq = len(s), len(q)

# blastn HSPs（dotplot 用）
cmd = [BLASTN, '-task', 'megablast', '-query', ALN + r'\smbu08_chr.fa', '-db', ALN + r'\smbu06_chr_db',
       '-outfmt', '6 pident length qstart qend sstart send', '-evalue', '1e-5', '-max_hsps', '2000',
       '-num_threads', '8']
r = run(cmd)
hsps = []
for line in r.stdout.strip().split('\n'):
    if line:
        p = line.split('\t')
        hsps.append((int(p[2]), int(p[3]), int(p[4]), int(p[5]), float(p[0]), int(p[1])))
big = [h for h in hsps if h[5] >= 50000]

fig = plt.figure(figsize=(7.2, 6.8))
gs = gridspec.GridSpec(2, 2, height_ratios=[1.15, 0.85], hspace=0.38, wspace=0.3)

# a) dotplot
ax = fig.add_subplot(gs[0, 0])
for qs, qe, ss, se, pid, al in hsps:
    if al < 2000:
        continue
    lw = 1.4 if al >= 50000 else 0.35
    alpha = 1.0 if al >= 50000 else 0.5
    col = C['shared'] if al >= 50000 else C['mge']
    ax.plot([qs / 1e6, qe / 1e6], [ss / 1e6, se / 1e6], color=col, lw=lw, alpha=alpha)
ax.plot([0, 2.64], [0, 2.64], ls=':', color='gray', lw=0.5)
# 缺失指示
ax.axvline(1.2558, color=C['prophage'], lw=0.5, ls='--')
ax.plot([1.2558, 1.2558], [1.2558, 1.2945], color=C['prophage'], lw=2.2)
ax.annotate('38,719-bp prophage\nabsent in smbu08', xy=(1.26, 1.28), xytext=(0.35, 1.7),
            fontsize=5.8, color=C['prophage'],
            arrowprops=dict(arrowstyle='->', lw=0.5, color=C['prophage']))
ax.annotate('1 SNP\n(C→T)', xy=(2.404, 2.443), xytext=(1.9, 2.05), fontsize=5.5,
            arrowprops=dict(arrowstyle='->', lw=0.5))
ax.set_xlabel('smbu08 chromosome (Mb)', fontsize=6)
ax.set_ylabel('smbu06 chromosome (Mb)', fontsize=6)
ax.tick_params(labelsize=5.5)
ax.set_title('Whole-chromosome alignment (BLASTN)', fontsize=7)
ax.text(0.03, 0.95, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

# b) 共线性示意
ax2 = fig.add_subplot(gs[0, 1])
def chr_track(y, length, color, label):
    ax2.plot([0, length / 1e6], [y, y], color=color, lw=2.2)
    ax2.text(-0.05, y, label, ha='right', va='center', fontsize=6)
chr_track(1.0, Ls, C['smbu06'], 'smbu06')
chr_track(0.0, Lq, C['smbu08'], 'smbu08')
# prophage block
ax2.add_patch(Rectangle((1.2558, 0.88), 0.0387, 0.24, facecolor=C['prophage'], edgecolor='none'))
ax2.text(1.275, 1.28, 'prophage (38,719 bp)', ha='center', fontsize=5.6, color=C['prophage'])
# att markers（交错高度避免重叠）
ax2.plot([1.2558, 1.2558], [1.02, 1.14], color='black', lw=0.6)
ax2.plot([1.2945, 1.2945], [1.02, 1.20], color='black', lw=0.6)
ax2.text(1.2350, 1.15, 'attL', ha='right', fontsize=5.2)
ax2.text(1.3150, 1.21, 'attR', ha='left', fontsize=5.2)
ax2.plot([1.2558, 1.2558], [-0.14, -0.02], color='black', lw=0.6)
ax2.text(1.2558, -0.19, 'attJ (== attR)', ha='center', fontsize=5.2)
# SNP
ax2.plot([2.442862 / 1e6, 2.404143 / 1e6], [0.94, 0.06], color='black', lw=0.5, ls='--')
ax2.plot(2.404143 / 1e6, 0.0, marker='v', color='black', ms=3)
ax2.text(2.404143 / 1e6 - 0.03, -0.19, '1 SNP', ha='right', fontsize=5.6)
ax2.set_xlim(-0.25, 2.85)
ax2.set_ylim(-0.3, 1.45)
ax2.axis('off')
ax2.set_title('Synteny: one prophage excision + one SNP', fontsize=7)
ax2.text(0.03, 0.95, 'b', transform=ax2.transAxes, fontsize=9, fontweight='bold', va='top')

# c) att 结构差异
ax3 = fig.add_subplot(gs[1, 0])
attL = s[1255782:1255928]
attJ = q[1255782:1255928]
diff = [i for i in range(146) if attL[i] != attJ[i]]
ax3.plot(range(146), [1] * 146, color=C['smbu06'], lw=3, solid_capstyle='butt', alpha=0.5)
for i in diff:
    ax3.plot([i, i], [0.88, 1.12], color='black', lw=0.5)
ax3.text(-3, 1.0, 'attL\n(smbu06)', ha='right', va='center', fontsize=5.4)
ax3.plot(range(146), [0] * 146, color=C['smbu08'], lw=3, solid_capstyle='butt', alpha=0.5)
ax3.text(-3, 0.0, 'attJ\n(smbu08)', ha='right', va='center', fontsize=5.4)
ax3.set_xlim(-30, 150)
ax3.set_ylim(-0.5, 1.5)
ax3.axis('off')
ax3.set_title('Junction repeat (146 bp): attL vs attJ (== attR) — %d/146 differences' % len(diff), fontsize=7)
ax3.text(0.01, 0.95, 'c', transform=ax3.transAxes, fontsize=9, fontweight='bold', va='top')
ax3.text(70, -0.42, 'position in repeat (bp)', ha='center', fontsize=5.5)

# d) 摘要
ax4 = fig.add_subplot(gs[1, 1])
ax4.axis('off')
txt = (
    'smbu06 chromosome : 2,676,906 bp\n'
    'smbu08 chromosome : 2,638,187 bp\n'
    'Difference       : −38,719 bp (one excisable element)\n'
    '                   + 1 SNP (smbu08 2,404,143 C→T;\n'
    '                   smbu06 2,442,862) + junction-repeat\n'
    '                   variation confined to attL vs attJ\n'
    '\n'
    'Colinear backbone: 2,638,186/2,638,187 bp covered\n'
    '(100.00% of smbu08) in two BLASTN HSPs;\n'
    'identity 100.000% (back segment) / 99.998% (front)\n'
    '\n'
    'smbu06 is byte-identical to E. lactis 194\n'
    '(GCF_056582645.1) chromosome (0 rotation needed);\n'
    'ANIb vs 194 = 100.00%; vs L. lactis 14B4 = 83.81%'
)
ax4.text(0.0, 1.0, txt, fontsize=6.0, va='top', linespacing=1.55, family='DejaVu Sans')
ax4.text(0.0, 1.06, 'd', transform=ax4.transAxes, fontsize=9, fontweight='bold', va='bottom')

savefig(fig, os.path.join(OUT, 'Fig2_chromosome', 'Fig2'))
print('fig2 done; big HSPs:', len(big))
