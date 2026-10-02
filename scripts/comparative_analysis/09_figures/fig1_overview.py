# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig1_overview.py — Figure 1：两株闭环复制子概览与组装统计
a) smbu06：染色体环 + Plasmid1 环（标注前噬菌体区、细菌素模块）
b) smbu08：染色体环 + Plasmid1 环 + Plasmid2 环（标注 attJ、切离二聚体深度）
c) 组装与测序统计表
"""
import os
import sys
import csv

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p3lib import WORK, read_fasta
from fig_style import plt, savefig, C
from matplotlib.patches import Patch

OUT = os.path.join(WORK, '09_figures')

stats = {}
with open(os.path.join(WORK, '01_assembly_qc', 'replicon_stats.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        stats[(r['strain'], r['replicon'])] = r

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.32, wspace=0.45)


def ring_panel(ax, replicons, title, notes, tag):
    """replicons: list of (name, length, color, radius, features[(a0,a1,color)])"""
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    for name, L, color, r, feats, lw in replicons:
        th = np.linspace(0, 2 * np.pi, 400)
        ax.plot(th, [r] * 400, color=color, lw=lw, solid_capstyle='butt')
        for a0, a1, fcol in feats:
            t = np.linspace(a0 / L * 2 * np.pi, a1 / L * 2 * np.pi, 60)
            ax.plot(t, [r] * 60, color=fcol, lw=lw + 1.6, solid_capstyle='butt')
        ax.text(2 * np.pi * 0.02, r + 0.09, '%s\n%s bp' % (name, format(L, ',')), fontsize=5.2,
                color=color, va='center')
    ax.set_ylim(0.35, 1.75)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.set_title(title, fontsize=7.5, pad=6)
    for i, (txt, col) in enumerate(notes):
        ax.text(0.5, -0.02 - i * 0.07, txt, transform=ax.transAxes, fontsize=5.4, color=col, ha='center')
    ax.text(0.03, 1.02, tag, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')


# a) smbu06
ax = fig.add_subplot(gs[0, 0], projection='polar')
ring_panel(ax, [
    ('Chromosome1', 2676906, C['smbu06'], 1.20,
     [(1255782, 1294501, C['prophage']), (898687, 899844, C['mge'])], 1.6),
    ('Plasmid1', 128837, C['shared'], 0.70,
     [(31073, 74688, C['bacteriocin'])] + [(r['start'], r['end'], C['mge']) for r in []], 1.6),
], 'smbu06 (2 replicons)',
    [('prophage 38,719 bp (1,255,783–1,294,501)', C['prophage']),
     ('bacteriocin module (31.1–74.7 kb)', C['bacteriocin'])], 'a')

# b) smbu08
ax2 = fig.add_subplot(gs[0, 1], projection='polar')
ring_panel(ax2, [
    ('Chromosome1', 2638187, C['smbu08'], 1.20,
     [(1255782, 1255928, C['prophage']), (898687 - 38719, 899844 - 38719, C['mge'])], 1.6),
    ('Plasmid1', 128837, C['shared'], 0.78,
     [(31073, 74688, C['bacteriocin'])], 1.6),
    ('Plasmid2', 77437, C['plasmid2'], 0.48,
     [(0, 38719, C['plasmid2']), (38720, 77437, C['plasmid2'])], 1.6),
], 'smbu08 (3 replicons)',
    [('attJ (prophage excised)', C['prophage']),
     ('Plasmid2 = excised-prophage dimer (2,600× depth)', C['plasmid2'])], 'b')

# c) 统计表
ax3 = fig.add_subplot(gs[:, 2])
ax3.axis('off')
rows = [
    ['', 'smbu06', 'smbu08'],
    ['Chromosome (bp)', '2,676,906', '2,638,187'],
    ['Plasmid1 (bp)', '128,837', '128,837'],
    ['Plasmid2 (bp)', '—', '77,437'],
    ['Total (bp)', '2,805,743', '2,844,461'],
    ['GC (%)', '38.4', '38.3'],
    ['CDS (chr/pl1/pl2)', '2,600/145/—', '2,547/146/102'],
    ['tRNA (chr)', '69', '70'],
    ['rRNA operons', '6', '6'],
    ['Illumina reads (pairs)', '4,173,470', '4,271,009'],
    ['Read length (bp)', '150 × 2', '150 × 2'],
    ['Chromosome depth (×)', '460', '440'],
    ['Plasmid1 depth (×)', '570', '600'],
    ['Plasmid2 depth (×)', '—', '2,600'],
    ['Hybrid assembly', 'Illumina + Nanopore', 'Illumina + Nanopore'],
]
tbl = ax3.table(cellText=[r[1:] for r in rows[1:]], rowLabels=[r[0] for r in rows[1:]],
                colLabels=rows[0][1:], loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(5.8)
tbl.scale(1.0, 1.32)
for (ri, ci), cell in tbl.get_celld().items():
    cell.set_linewidth(0.4)
    if ri == 0:
        cell.set_text_props(fontweight='bold')
ax3.set_title('Assembly and sequencing statistics', fontsize=7.5)
ax3.text(0.0, 0.99, 'c', transform=ax3.transAxes, fontsize=9, fontweight='bold', va='top')

fig.text(0.5, 0.012,
         'Both strains share a byte-identical 128,837-bp plasmid (grey). smbu06 carries a 38.7-kb PBSX-like prophage that is excised in smbu08,\n'
         'where the excised unit persists as a high-copy tandem-dimer plasmid (Plasmid2).',
         fontsize=6.0, ha='center', va='bottom')

savefig(fig, os.path.join(OUT, 'Fig1_overview', 'Fig1'))
print('fig1 done')
