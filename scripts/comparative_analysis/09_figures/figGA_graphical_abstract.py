# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""figGA_graphical_abstract.py — Graphical abstract（水平转移仅作待检验假说）"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK
from fig_style import plt, savefig, C
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Patch

OUT = os.path.join(WORK, '09_figures')
fig = plt.figure(figsize=(7.2, 3.6))
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('off')
ax.set_xlim(0, 100)
ax.set_ylim(0, 50)

# 左：两株基因组
ax.text(1, 47.5, 'Graphical abstract', fontsize=10, fontweight='bold')

# smbu06
ax.plot([4, 24], [35, 35], color=C['smbu06'], lw=3)
ax.add_patch(Rectangle((12.6, 33.9), 3.6, 2.2, facecolor=C['prophage'], edgecolor='none'))
ax.text(14.4, 38.6, '38,719-bp PBSX-like\nelement (integrated)', fontsize=5.6, ha='center', color=C['prophage'])
ax.text(4, 32.4, 'smbu06 chromosome 2,676,906 bp', fontsize=6)
ax.text(4, 30.2, '= E. lactis 194 chromosome\n  (byte-identical; 0 SNPs)', fontsize=5.6, color='#555555')
# smbu08
ax.plot([4, 24], [18, 18], color=C['smbu08'], lw=3)
ax.plot(12.7, 17.2, marker='|', color=C['prophage'], ms=8, mew=1.6)
ax.text(14.6, 20.4, 'element excised\n(attJ junction)', fontsize=5.6, ha='left', color=C['prophage'])
ax.text(4, 15.4, 'smbu08 chromosome 2,638,187 bp', fontsize=6)
ax.text(4, 13.2, 'difference = 38,719-bp element + 1 SNP', fontsize=5.6, color='#555555')

# 中：共享质粒 + 盒
ax.add_patch(FancyArrowPatch((25, 35), (31, 35), arrowstyle='-', lw=2.2, color=C['shared']))
ax.add_patch(FancyArrowPatch((25, 18), (31, 18), arrowstyle='-', lw=2.2, color=C['shared']))
circ = Circle((36, 26.5), 7.5, fill=False, edgecolor=C['shared'], lw=2.2)
ax.add_patch(circ)
ax.text(36, 39.2, 'shared plasmid 128,837 bp\nBYTE-IDENTICAL in both genomes\n(canonical SHA-256 + full-length BLASTN; 0 SNPs)',
        fontsize=5.8, ha='center', color=C['shared'])
# 盒位置
ax.add_patch(Rectangle((37.8, 25.3), 2.6, 2.4, facecolor=C['bacteriocin'], edgecolor='none', zorder=3))
ax.text(36, 12.9, 'bacteriocin cassette (−20 kb apart in total):\nenterocin P-like (94.4% vs enterocin P, motif YDNGI)\nSakacin-A immunity factor (10 bp downstream)\nhiracin JM79-like (59.3% vs hiracin JM79)\nno conjugation module',
        fontsize=5.4, ha='center', color=C['bacteriocin'])

# 右：切离态 + 面板分布
ax.add_patch(FancyArrowPatch((43.8, 26.5), (49, 26.5), arrowstyle='-|>', mutation_scale=8, lw=0.9))
ax.plot(51.5, 26.5, marker='o', ms=26, mfc='none', mec=C['plasmid2'], mew=1.8)
ax.text(51.5, 35.2, 'excised state in smbu08:\ntandem-dimer replicon\n77,437 bp at 4.1–5.0×\nchromosome depth', fontsize=5.6, ha='center', color=C['plasmid2'])

# 面板分布条
ax.text(66, 40.5, 'Public panel (159 deduplicated genomes):', fontsize=6.2, fontweight='bold')
rows = [('exact plasmid (100%)', 1, 159, C['shared']),
        ('plasmid backbone ≥30% shared', 34, 159, '#8899AA'),
        ('enterocin P-like gene', 17, 159, C['bacteriocin']),
        ('immunity gene', 9, 159, C['immunity']),
        ('hiracin JM79-like gene', 11, 159, C['bacteriocin']),
        ('nisin cluster', 5, 159, C['core'])]
for i, (lab, v, tot, col) in enumerate(rows):
    y = 36.5 - i * 4.1
    ax.barh(y, 22, left=66, color='#EEEEEE', height=2.2)
    ax.barh(y, 22 * v / tot, left=66, color=col, height=2.2)
    ax.text(65.4, y, lab, fontsize=5.4, va='center', ha='right')
    ax.text(66 + 22 * v / tot + 0.6, y, str(v), fontsize=5.4, va='center', ha='left')
ax.text(89, 36.5 - 6 * 4.1 - 2, '5/159 nisin clusters = the\ndeliberately added controls only', fontsize=5.0, va='top')
# 假说标注
ax.add_patch(Rectangle((3, 3.5), 94, 5.6, fill=False, edgecolor='#888888', lw=0.7, ls='--'))
ax.text(4, 6.3, 'OPEN HYPOTHESIS (not tested here): inter-plasmid/backbone mobilisation of the cassette, and any direction or timing of transfer — require population data and cannot be inferred from these static genomes.',
        fontsize=5.2, va='center', color='#444444')

savefig(fig, os.path.join(OUT, 'GraphicalAbstract', 'GraphicalAbstract'))
print('graphical abstract done')
