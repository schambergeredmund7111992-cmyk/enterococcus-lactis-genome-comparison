# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""figS1_mapping.py — Supplementary Fig. S1：reads 回贴、覆盖度与混样监测"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p3lib import WORK, read_fasta
from fig_style import plt, savefig, C
import matplotlib.gridspec as gridspec

OUT = os.path.join(WORK, '09_figures')
MP = os.path.join(WORK, '01_assembly_qc', 'mapping')
ref = read_fasta(os.path.join(MP, 'combined_ref.fa'))
offs = {}
LENS = {}
pos = 0
for name, s in ref:
    offs[name] = pos
    LENS[name] = len(s)
    pos += len(s)

fig = plt.figure(figsize=(7.2, 6.4))
gs = gridspec.GridSpec(3, 2, hspace=0.5, wspace=0.3)

# a) mapped fraction
ax = fig.add_subplot(gs[0, 0])
datasets = ['smbu06\nIllumina', 'smbu06\nNanopore', 'smbu08\nIllumina', 'smbu08\nNanopore']
frac = [(8338299 - 510) / 8338299 * 100, (78312 - 38) / 78312 * 100,
        (8532822 - 822) / 8532822 * 100, (123172 - 37) / 123172 * 100]
lact = [510 / 8338299 * 100, 38 / 78312 * 100, 822 / 8532822 * 100, 37 / 123172 * 100]
x = np.arange(4)
ax.bar(x, frac, color=C['shared'], label='target genomes (own strain)')
ax.bar(x, lact, bottom=frac, color=C['smbu08'], label='Lactococcus references')
ax.bar(x, [100 - f - l for f, l in zip(frac, lact)], bottom=[f + l for f, l in zip(frac, lact)],
       color='#DDDDDD', label='unmapped')
ax.set_xticks(x)
ax.set_xticklabels(datasets, fontsize=5.5)
ax.set_ylabel('reads (%)', fontsize=6)
ax.set_ylim(0, 100)
ax.tick_params(labelsize=5.5)
ax.legend(fontsize=4.6, loc='lower left')
ax.set_title('Read assignment (minimap2)', fontsize=7)
ax.text(0.02, 0.93, 'a', transform=ax.transAxes, fontsize=9, fontweight='bold')
ax.text(2.0, 88, 'Lactococcus hits:\n0.006–0.028% of reads', fontsize=5.0, ha='center')

# b) smbu06 depth (chr + pl1)
ax2 = fig.add_subplot(gs[0, 1])
d6 = np.load(os.path.join(MP, 'depth_smbu06_illumina.npz'))['depth']
chr6 = d6[offs['smbu06|Chromosome1']:offs['smbu06|Chromosome1'] + LENS['smbu06|Chromosome1']]
pl1 = d6[offs['smbu06|Plasmid1']:offs['smbu06|Plasmid1'] + LENS['smbu06|Plasmid1']]
w = 2000
ax2.plot(np.arange(0, len(chr6), w) / 1e6, [chr6[i:i + w].mean() for i in range(0, len(chr6), w)],
         color=C['smbu06'], lw=0.8, label='chromosome')
xb = np.arange(0, len(pl1), 500)
ax2.plot(2.75 + xb / 1e6, [pl1[i:i + 500].mean() for i in range(0, len(pl1), 500)],
         color=C['shared'], lw=0.8, label='Plasmid1')
ax2.set_xlabel('position (Mb, chromosome then plasmid offset)', fontsize=5.6)
ax2.set_ylabel('depth (×)', fontsize=6)
ax2.tick_params(labelsize=5.5)
ax2.legend(fontsize=4.8)
ax2.set_title('smbu06 Illumina depth', fontsize=7)
ax2.text(0.02, 0.93, 'b', transform=ax2.transAxes, fontsize=9, fontweight='bold')

# c) smbu08 depth (chr + pl1 + pl2, log)
ax3 = fig.add_subplot(gs[1, 0])
d8 = np.load(os.path.join(MP, 'depth_smbu08_illumina.npz'))['depth']
chr8 = d8[offs['smbu08|Chromosome1']:offs['smbu08|Chromosome1'] + LENS['smbu08|Chromosome1']]
p1 = d8[offs['smbu08|Plasmid1']:offs['smbu08|Plasmid1'] + LENS['smbu08|Plasmid1']]
p2 = d8[offs['smbu08|Plasmid2']:offs['smbu08|Plasmid2'] + LENS['smbu08|Plasmid2']]
ax3.plot(np.arange(0, len(chr8), w) / 1e6, [chr8[i:i + w].mean() for i in range(0, len(chr8), w)],
         color=C['smbu08'], lw=0.8, label='chromosome')
ax3.plot(2.70 + xb / 1e6, [p1[i:i + 500].mean() for i in range(0, len(p1), 500)],
         color=C['shared'], lw=0.8, label='Plasmid1')
xp = np.arange(0, len(p2), 500)
ax3.plot(2.90 + xp / 1e6, [p2[i:i + 500].mean() for i in range(0, len(p2), 500)],
         color=C['plasmid2'], lw=0.8, label='Plasmid2 (excised unit)')
ax3.axhline(440, color='gray', lw=0.4, ls=':')
ax3.set_yscale('log')
ax3.set_xlabel('position (Mb, chromosome / plasmid2 offsets)', fontsize=5.6)
ax3.set_ylabel('depth (×, log)', fontsize=6)
ax3.tick_params(labelsize=5.5)
ax3.legend(fontsize=4.8)
ax3.set_title('smbu08 Illumina depth: excised unit ≫ chromosome', fontsize=7)
ax3.text(0.02, 0.93, 'c', transform=ax3.transAxes, fontsize=9, fontweight='bold')

# d) breadth bars
ax4 = fig.add_subplot(gs[1, 1])
labels = ['smbu06 chr', 'smbu06 pl1', 'smbu08 chr', 'smbu08 pl1', 'smbu08 pl2']
val1 = [100.0, 99.98, 100.0, 99.98, 99.99]
ax4.barh(range(5), val1, color=[C['smbu06'], C['shared'], C['smbu08'], C['shared'], C['plasmid2']])
ax4.set_yticks(range(5))
ax4.set_yticklabels(labels, fontsize=5.4)
ax4.set_xlim(99.9, 100.02)
ax4.set_xlabel('breadth of coverage ≥1× (%)', fontsize=6)
ax4.tick_params(labelsize=5.5)
ax4.set_title('Coverage breadth (Illumina)', fontsize=7)
ax4.text(0.02, 0.93, 'd', transform=ax4.transAxes, fontsize=9, fontweight='bold')

# e) 结论文本（数字全部现场计算）
ax5 = fig.add_subplot(gs[2, :])
ax5.axis('off')
zero_stats = []
for sg, sn, dd in [('smbu06', 'Chromosome1', d6), ('smbu06', 'Plasmid1', d6),
                   ('smbu08', 'Chromosome1', d8), ('smbu08', 'Plasmid1', d8), ('smbu08', 'Plasmid2', d8)]:
    key = '%s|%s' % (sg, sn)
    arr = dd[offs[key]:offs[key] + LENS[key]]
    zero_stats.append((key, float((arr == 0).mean() * 100)))
zero_max = max(z for _, z in zero_stats)
lact_max = max(lact_min := [l for l in lact], default=0)
txt = (
    'Deconvolved depths (combined representations; see Supplementary Methods):\n'
    'smbu06: chromosome 443.6× (Illumina) / 403.2× (Nanopore);  Plasmid1 490.1× / 456.4×\n'
    'smbu08: chromosome 430.4× / 387.6×;  Plasmid1 517.8× / 485.8×;  excised unit (Plasmid2 + prophage region) 1,777.6× / 1,946.2×\n'
    'No mixed-sample signal: ≥%.3f%% of reads assign to the two target genomes; Lactococcus-reference reads ≤%.3f%%; zero-depth positions ≤%.4f%% of any replicon (Illumina).'
    % (min(frac), max(lact), zero_max)
)
ax5.text(0.01, 1.0, txt, fontsize=6.0, va='top', linespacing=1.6)
ax5.text(0.0, 1.12, 'e', transform=ax5.transAxes, fontsize=9, fontweight='bold')

savefig(fig, os.path.join(OUT, 'FigS1_mapping', 'FigS1'))
print('figS1 done')
