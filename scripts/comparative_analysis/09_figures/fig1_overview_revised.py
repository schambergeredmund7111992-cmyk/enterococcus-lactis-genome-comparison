# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""Revised Figure 1: concise replicon overview and structural comparison.

This version intentionally moves detailed assembly statistics to Supplementary
Table S1 and avoids claiming that the smbu08 circular component is an
autonomously replicating plasmid before long-read junction validation.
"""
import os
import sys
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK
from fig_style import plt, savefig, C

OUT = os.path.join(WORK, '09_figures', 'Fig1_overview_revised')

# Colour encodes replicon/feature type, not strain identity.
CHROM = '#2B6CB0'
SHARED = C['shared']
DIFF = C['prophage']
CIRCULAR = C['plasmid2']
BACT = C['bacteriocin']

fig = plt.figure(figsize=(7.2, 6.25))
gs = fig.add_gridspec(2, 2, height_ratios=[1.10, 1.22], hspace=0.43, wspace=0.28)


def ring(ax, items, title, subtitle, panel):
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    for length, radius, base_color, highlights, lw in items:
        theta = np.linspace(0, 2 * np.pi, 720)
        ax.plot(theta, np.repeat(radius, len(theta)), color=base_color, lw=lw, solid_capstyle='round')
        for start, end, colour in highlights:
            segment = np.linspace(start / length * 2 * np.pi, end / length * 2 * np.pi, 120)
            ax.plot(segment, np.repeat(radius, len(segment)), color=colour, lw=lw + 2.0,
                    solid_capstyle='round', zorder=3)
    ax.set_ylim(0.25, 1.42)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.set_title(title, loc='left', fontsize=8.3, fontweight='bold', pad=6)
    ax.text(0.01, 0.94, subtitle, transform=ax.transAxes, fontsize=6.5, va='top')
    ax.text(-0.10, 1.08, panel, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')


ax_a = fig.add_subplot(gs[0, 0], projection='polar')
ring(ax_a, [
    (2676906, 1.10, CHROM, [(1255782, 1294501, DIFF)], 3.2),
    (128837, 0.65, SHARED, [(31073, 74688, BACT)], 3.2),
], 'smbu06', '2.677-Mb chromosome + 128.8-kb shared plasmid', 'a')

ax_b = fig.add_subplot(gs[0, 1], projection='polar')
ring(ax_b, [
    (2638187, 1.10, CHROM, [(1255782, 1255928, DIFF)], 3.2),
    (128837, 0.72, SHARED, [(31073, 74688, BACT)], 3.2),
    (77437, 0.39, CIRCULAR, [(0, 38719, CIRCULAR), (38720, 77437, CIRCULAR)], 3.2),
], 'smbu08', '2.638-Mb chromosome + shared plasmid + 77.4-kb component', 'b')

legend = [
    Line2D([0], [0], color=CHROM, lw=3.4, label='Chromosome'),
    Line2D([0], [0], color=SHARED, lw=3.4, label='Shared 128.8-kb plasmid'),
    Line2D([0], [0], color=DIFF, lw=3.4, label='38.7-kb differential region / attJ'),
    Line2D([0], [0], color=CIRCULAR, lw=3.4, label='smbu08 77.4-kb circular component'),
    Line2D([0], [0], color=BACT, lw=3.4, label='Bacteriocin-associated region'),
]
fig.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.5, 0.615), ncol=3,
           fontsize=6.1, handlelength=1.8, columnspacing=1.2, frameon=False)

# Panel c: the structural point of the figure, using schematic rather than a dense table.
ax_c = fig.add_subplot(gs[1, :])
ax_c.set_xlim(0, 108)
ax_c.set_ylim(-0.55, 3.55)
ax_c.axis('off')
ax_c.text(-0.01, 1.04, 'c', transform=ax_c.transAxes, fontsize=10, fontweight='bold', va='bottom')
ax_c.text(0, 3.28, 'Structural difference between the near-isogenic genomes', fontsize=8.3, fontweight='bold')

# smbu06 chromosomal state
ax_c.text(0, 2.45, 'smbu06 chromosome', fontsize=7, ha='left', va='center')
ax_c.plot([24, 86], [2.45, 2.45], color=CHROM, lw=4, solid_capstyle='round')
ax_c.add_patch(Rectangle((52, 2.25), 8.0, 0.40, facecolor=DIFF, edgecolor='none', zorder=3))
ax_c.text(56, 2.86, '38.7-kb phage-related region', fontsize=6.4, ha='center', color=DIFF)
ax_c.text(52, 2.02, 'attL', fontsize=5.8, ha='center')
ax_c.text(60, 2.02, 'attR', fontsize=5.8, ha='center')

# smbu08 chromosomal state
ax_c.text(0, 1.42, 'smbu08 chromosome', fontsize=7, ha='left', va='center')
ax_c.plot([24, 86], [1.42, 1.42], color=CHROM, lw=4, solid_capstyle='round')
ax_c.plot([56, 56], [1.12, 1.72], color=DIFF, lw=1.5)
ax_c.text(56, 0.90, 'attJ', fontsize=6.2, ha='center', color=DIFF)
ax_c.annotate('', xy=(56, 1.74), xytext=(56, 2.18), arrowprops=dict(arrowstyle='-|>', color=DIFF, lw=0.8))

# smbu08 circular component model
ax_c.text(0, 0.18, 'smbu08 circular component', fontsize=7, ha='left', va='center')
ax_c.add_patch(Rectangle((35, -0.06), 24, 0.48, facecolor=CIRCULAR, edgecolor='white', linewidth=0.8))
ax_c.add_patch(Rectangle((59, -0.06), 24, 0.48, facecolor=CIRCULAR, edgecolor='white', linewidth=0.8))
ax_c.text(47, 0.18, 'unit 1', fontsize=6.4, color='white', ha='center', va='center')
ax_c.text(71, 0.18, 'unit 2', fontsize=6.4, color='white', ha='center', va='center')
ax_c.add_patch(FancyArrowPatch((83, 0.18), (35, 0.18), connectionstyle='arc3,rad=0.52',
                                arrowstyle='-|>', mutation_scale=9, lw=0.9, color=CIRCULAR))
ax_c.text(59, -0.42, '77.4 kb; two-unit assembly model (topology unresolved)', fontsize=6.3, ha='center', color=CIRCULAR)

ax_c.text(89, 2.45, '2,676,906 bp', fontsize=6.4, ha='left', va='center', color=CHROM)
ax_c.text(89, 1.42, '2,638,187 bp', fontsize=6.4, ha='left', va='center', color=CHROM)

fig.text(0.5, 0.006,
         'Detailed assembly and sequencing statistics are provided in Supplementary Table S1. '
         'Long-read support excludes a junction-level assembly artefact; molecular topology remains unresolved.',
         fontsize=6.0, ha='center', va='bottom')

savefig(fig, os.path.join(OUT, 'Fig1_revised'))
print('revised Fig1 complete')
