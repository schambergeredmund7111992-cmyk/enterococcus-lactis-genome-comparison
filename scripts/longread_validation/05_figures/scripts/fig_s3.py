# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig_s3.py — Supplementary Figure S3（本会话实现，自包含；数据取自各分析输出）
实线 = 直接 read 证据；虚线 = 模型推断。
A 结构模型与样本混合状态（保留态多数 / 切离态少数 / 质粒）
B 关键 junction 的跨越 read 示意（read ID、长度、左右锚定、mapQ、identity）
C IGV 风格 read 排布（attJ 与质粒 junction）
D 区域深度与相对拷贝数（两种算法 + bootstrap CI95）
E 竞争模型支持对照（strong/ambiguous）
F 染色体状态计数：smbu08 vs smbu06 阳性对照（混合群体核心证据）
输出：05_figures/{png(600dpi),pdf,svg}/FigS3_junction_reads.*
"""
import os
import re
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
import matplotlib.gridspec as gridspec

ATTJ = r'<PROJECT_ROOT>\longread_validation'
FIG = os.path.join(ATTJ, '05_figures')
SR = os.path.join(ATTJ, '03_supporting_reads')
CT = os.path.join(ATTJ, '04_controls')

C = dict(smbu06='#2c7fb8', smbu08='#d95f02', pl2='#7570b3', acc='#e7298a', grey='#888888', light='#cccccc')
for sub in ['png', 'pdf', 'svg']:
    os.makedirs(os.path.join(FIG, sub), exist_ok=True)


def load_tsv(p):
    rows = []
    if not os.path.exists(p):
        return rows
    lines = open(p, encoding='utf-8').read().splitlines()
    hdr = lines[0].lstrip('#').split('\t')
    for l in lines[1:]:
        if not l or l.startswith('#'):
            continue
        rows.append(dict(zip(hdr, l.split('\t'))))
    return rows


V = {r['item']: r['value'] for r in load_tsv(os.path.join(ATTJ, '01_references', 'independent_verification.tsv'))}
S_pos = int(V['S_pos_0based']); B = 38719; A = 146; PL2 = 77437
j_att_mid = S_pos + A // 2

cand = load_tsv(os.path.join(SR, 'junction_candidates2.tsv'))
cls3 = load_tsv(os.path.join(SR, 'read_classification3.tsv'))
bs = load_tsv(os.path.join(CT, 'depth_ratio_bootstrap3.tsv'))
dep = load_tsv(os.path.join(CT, 'depth_region_stats3.tsv'))
sjc = load_tsv(os.path.join(CT, 'state_junction_counts.tsv'))

CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')
JP = {'attJ_chr08_state': ('attJ (chr08, excised state)', C['smbu08']),
      'dimer_forward': ('plasmid u1→u2 boundary', C['pl2']),
      'dimer_back_wrap': ('plasmid back boundary (origin)', C['acc']),
      'chr06_retained_state': ('retained (chr06)', C['smbu06'])}


def sjc_get(ds, jn, col):
    for r in sjc:
        if r['dataset'] == ds and r['junction'].startswith(jn):
            return r[col]
    return '?'


fig = plt.figure(figsize=(7.4, 12.0))
gs = gridspec.GridSpec(6, 1, height_ratios=[1.35, 1.5, 1.3, 1.0, 0.9, 0.85], hspace=0.68)

# ---------------- A 结构模型与混合状态 ----------------
ax = fig.add_subplot(gs[0]); ax.axis('off')
ax.set_xlim(0, 100); ax.set_ylim(0, 40)
ax.text(0.5, 39.2, 'A   Structural models and population states  (solid = read-supported; dashed = inference)',
        fontsize=7.2, fontweight='bold')
# retained chromosome (majority)
ax.plot([3, 47], [28.0, 28.0], color=C['smbu06'], lw=3.2, solid_capstyle='butt')
ax.add_patch(Rectangle((20, 27.3), 12.5, 1.4, facecolor=C['smbu06'], alpha=0.35, edgecolor='none'))
ax.text(26.2, 30.2, 'attL—X′ (38,573 bp)—attR', fontsize=5.4, ha='center', color=C['smbu06'])
ax.text(3, 33.6, 'retained state (majority of smbu08 reads, ~85%)', fontsize=5.8, color=C['smbu06'])
# excised chromosome (minority)
ax.plot([3, 27.5], [19.0, 19.0], color=C['smbu08'], lw=3.2, solid_capstyle='butt')
ax.plot([29, 47], [19.0, 19.0], color=C['smbu08'], lw=3.2, solid_capstyle='butt')
ax.plot(28.3, 19.0, marker='|', color='black', ms=9, mew=1.6)
ax.text(28.9, 20.9, 'attJ (146 bp, byte-identical to attR)', fontsize=5.4, color=C['smbu08'])
ax.text(3, 24.8, 'excised state (~15%; the delivered smbu08 chromosome assembly)', fontsize=5.8, color=C['smbu08'])
ax.add_patch(FancyArrowPatch((49, 28.0), (56, 28.0), arrowstyle='-|>', mutation_scale=8, lw=0.8, color=C['grey'], ls='--'))
ax.text(52.5, 30.2, 'excision (both states)', fontsize=4.8, ha='center', color=C['grey'])
# plasmid circle
cc = (80, 22.0); r = 6.0
ax.add_patch(Circle(cc, r, fill=False, edgecolor=C['pl2'], lw=2.6))
ax.plot(cc[0], cc[1] + r, marker='o', ms=3.4, color=C['pl2'])
ax.text(cc[0] + 1.0, cc[1] + r + 1.6, 'boundary A (read-supported)', fontsize=4.7, color=C['pl2'])
ax.plot(cc[0], cc[1] - r, marker='o', ms=3.4, color=C['acc'])
ax.text(cc[0] + 1.0, cc[1] - r - 2.4, 'boundary B / origin (read-supported)', fontsize=4.7, color=C['acc'])
ax.text(80, 32.8, 'Plasmid2 77,437 bp: two tandem units of attL+X', fontsize=5.4, ha='center', color=C['pl2'])
ax.text(80, 8.6, 'tandem-dimer representation; a monomer circle\nat 2× copy number is sequence-equivalent',
        fontsize=5.0, ha='center', color='#555555')
ax.text(2.5, 3.0, 'Read evidence: 41 reads strongly place the attJ state (chr08) over the retained model;\n'
                  '1,883 reads span the plasmid boundary (mapQ≥20); retained:excised chromosome counts ≈ 5.6:1.',
        fontsize=5.2, color='#333333')

# ---------------- B 跨 junction read 示意 ----------------
ax2 = fig.add_subplot(gs[1]); ax2.axis('off')
ax2.set_xlim(0, 100); ax2.set_ylim(0, 46)
ax2.text(0.5, 45.0, 'B   Junction-spanning Nanopore reads  (single primary alignment, mapQ ≥ 20, ≥ 1 kb anchors both sides)',
         fontsize=7.2, fontweight='bold')
rowsB = [
    ('attJ (chr08)', C['smbu08'], sjc_get('smbu08', 'excised_left', 'loose_reads'), sjc_get('smbu08', 'excised_left', 'strict_reads'),
     [c for c in cand if c['jtype'] == 'attJ_chr08_state' and c['pass_strict'] == 'True']),
    ('plasmid boundary (×2 on circle)', C['pl2'], '1,883', '992',
     [c for c in cand if c['jtype'] == 'dimer_forward' and c['pass_strict'] == 'True']),
    ('retained junction (chr06)', C['smbu06'], sjc_get('smbu08', 'retained_left', 'loose_reads'), sjc_get('smbu08', 'retained_left', 'strict_reads'),
     [c for c in cand if c['jtype'] == 'chr06_retained_state' and c['pass_strict'] == 'True']),
]
y = 41
for lab, col, nloose, nstrict, rows in rowsB:
    ax2.text(0.5, y, '%s  — loose reads: %s   strict reads: %s' % (lab, nloose, nstrict), fontsize=6.0, fontweight='bold', color=col)
    y -= 3.0
    rows = sorted(rows, key=lambda c: -int(c['AS'] or 0))
    for c in rows[:3]:
        al = int(c['a_left']); ar = int(c['a_right']); total = max(al + ar, 1)
        rl = int(c['read_len']) if c['read_len'] else 5000
        span = min(rl, 60000) / 60000 * 44
        xc = 50
        xl = xc - span * (al / total); xr = xc + span * (ar / total)
        ax2.plot([xl, xr], [y, y], color=col, lw=1.4, solid_capstyle='butt')
        ax2.plot([xl, xl], [y - 0.7, y + 0.7], color=col, lw=1.4)
        ax2.plot([xr, xr], [y - 0.7, y + 0.7], color=col, lw=1.4)
        ax2.plot(xc, y, marker='|', color='black', ms=8, mew=1.5)
        ax2.text(xr + 1.2, y, '%s  len=%s  anchors %s|%s  mapQ=%s  id=%s' %
                 (c['read'][:12], c['read_len'], al, ar, c['mapq'], c['ident']), fontsize=4.0, va='center')
        y -= 2.6
    y -= 3.6
ax2.plot([50, 50], [1.0, y + 5], color='black', lw=0.4, ls=':')
ax2.text(50, 0.3, 'junction position', fontsize=4.6, ha='center')

# ---------------- C IGV 风格 ----------------
ax3 = fig.add_subplot(gs[2])
ax3.set_title('C   Read alignments across key junctions (±2 kb window of the 10-kb junction references), IGV-style',
              fontsize=7.2, fontweight='bold', loc='left')
order = ['attJ_chr08_state', 'dimer_forward', 'chr06_retained_state']
yy = 0
yticks = []
for jt in order:
    lab, col = JP[jt]
    jabs = 5000
    rel = [c for c in cand if c['jtype'] == jt and c['primary'] == 'True' and c['covered'] == 'True' and float(c['ident'] or 0) >= 0.9]
    rel.sort(key=lambda c: -(int(c['a_left'] or 0) + int(c['a_right'] or 0)))
    y0 = yy
    for c in rel[:12]:
        ops = [(int(a), b) for a, b in CIG_RE.findall(c['cigar'])]
        cur = int(c['r0'])
        for n, op in ops:
            if op in 'M=X':
                a = cur - jabs; b = a + n
                if b > -2100 and a < 2100:
                    ax3.plot([max(a, -2100), min(b, 2100)], [yy, yy], lw=1.5, color=col, solid_capstyle='butt')
                cur += n
            elif op in 'DN':
                cur += n
        yy += 1
    yticks.append(((y0 + yy) / 2 - 0.5, lab))
    yy += 3
    ax3.axhline(yy - 1.8, color='#dddddd', lw=0.5)
ax3.axvline(0, color='black', lw=0.7, ls='--')
ax3.set_xlim(-2100, 2100); ax3.set_ylim(-1, yy + 1)
ax3.set_yticks([t[0] for t in yticks]); ax3.set_yticklabels([t[1] for t in yticks], fontsize=5.2)
ax3.set_xlabel('position relative to junction (bp)', fontsize=6)
ax3.tick_params(labelsize=5.5)

# ---------------- D 深度 ----------------
ax4 = fig.add_subplot(gs[3])
want = ['chr08_whole', 'pl1_whole', 'pl2_whole', 'attJ_±2kb(chr08)', 'attR_±2kb(chr06)',
        'retained_X_internal(chr06)', 'dimer_forward_±2kb(pl2)', 'dimer_back_wrap_±2kb(pl2)']
short = {'chr08_whole': 'chr08\nwhole', 'pl1_whole': 'Plasmid1\nwhole', 'pl2_whole': 'Plasmid2\nwhole',
         'attJ_±2kb(chr08)': 'attJ\n±2kb', 'attR_±2kb(chr06)': 'attR\n±2kb',
         'retained_X_internal(chr06)': 'element X\n(chr06)', 'dimer_forward_±2kb(pl2)': 'pl2 boundary A\n±2kb',
         'dimer_back_wrap_±2kb(pl2)': 'pl2 origin\n±2kb'}
depd = {d['region']: d for d in dep}
bsd = {b['region']: b for b in bs}
sel = [k for k in want if k in depd]
xs = np.arange(len(sel))
vals = [float(depd[k]['depth_primary_mean']) for k in sel]
cols = [C['grey'] if k == 'chr08_whole' else (C['smbu08'] if 'chr08' in k or 'chr06' in k else C['pl2']) for k in sel]
ax4.bar(xs, vals, color=cols, edgecolor='black', lw=0.4)
for i, k in enumerate(sel):
    b = bsd.get(k)
    lab = ('%.2f×\n(%.2f–%.2f)' % (float(b['ratio_primary_algA']), float(b['bootstrap_CI95_low']), float(b['bootstrap_CI95_high']))) if b else '1.00×'
    ax4.text(i, vals[i] * 1.03, lab, fontsize=4.2, ha='center')
ax4.set_ylim(0, max(vals) * 1.34)
ax4.set_xticks(xs); ax4.set_xticklabels([short[k] for k in sel], fontsize=5.0)
ax4.set_ylabel('Nanopore depth (primary, ×)', fontsize=6)
ax4.set_title('D   Region depth (primary alignments) and ratio versus chr08 (bootstrap CI95, read resampling)',
              fontsize=7.2, loc='left')
ax4.tick_params(labelsize=5.5)

# ---------------- E 模型支持 ----------------
ax5 = fig.add_subplot(gs[4])
classes = [('attJ_chr08_state', 'attJ state\n(chr08)', C['smbu08']),
           ('dimer_forward', 'dimer boundary\n(read-spanning)', C['pl2']),
           ('dimer_back_wrap', 'origin boundary\n(read-spanning)', C['acc']),
           ('chr06_retained_state', 'retained state\n(chr06)', C['smbu06']),
           ('negative_control', 'negative\ncontrols', '#bbbbbb')]
cnt = defaultdict(Counter)
for r in cls3:
    for jt in r['jtype'].split(','):
        cnt[jt][r['verdict']] += 1
xs = np.arange(len(classes))
strong = [cnt[k]['strong'] for k, _, _ in classes]
amb = [cnt[k]['ambiguous'] for k, _, _ in classes]
ax5.bar(xs, strong, color=[c for _, _, c in classes], edgecolor='black', lw=0.4, label='strong (model-distinguishing)')
ax5.bar(xs, amb, bottom=strong, color=[c for _, _, c in classes], alpha=0.4, edgecolor='black', lw=0.3,
        label='ambiguous (ΔAS<10 vs best alternative)')
for i in range(len(classes)):
    ax5.text(i, strong[i] + amb[i] + 8, str(strong[i] + amb[i]), ha='center', fontsize=5)
ax5.set_xticks(xs); ax5.set_xticklabels([lab for _, lab, _ in classes], fontsize=5.0)
ax5.set_ylabel('candidate reads', fontsize=6)
ax5.set_title('E   Per-read best-model assignment (pre-registered: strong = ΔAS≥10, identity≥90%, coverage≥30%)',
              fontsize=7.2, loc='left')
ax5.legend(fontsize=4.6, frameon=False, loc='upper right')
ax5.tick_params(labelsize=5.5)

# ---------------- F 染色体状态计数 ----------------
ax6 = fig.add_subplot(gs[5])
groups = [('retained junction\n(chr06 flank|attL)', 'retained_left'), ('excised junction\n(chr08 flank|attJ)', 'excised_left')]
x = np.arange(len(groups)); w = 0.36
v08 = [int(sjc_get('smbu08', g[1], 'loose_reads')) for g in groups]
v06 = [int(sjc_get('smbu06_control', g[1], 'loose_reads')) for g in groups]
ax6.bar(x - w / 2, v08, w, color=C['smbu08'], edgecolor='black', lw=0.4, label='smbu08')
ax6.bar(x + w / 2, v06, w, color=C['smbu06'], edgecolor='black', lw=0.4, label='smbu06 (pure-retained control)')
for i in range(len(groups)):
    ax6.text(i - w / 2, v08[i] + 4, str(v08[i]), ha='center', fontsize=5)
    ax6.text(i + w / 2, v06[i] + 4, str(v06[i]), ha='center', fontsize=5)
ax6.set_xticks(x); ax6.set_xticklabels([g[0] for g in groups], fontsize=5.2)
ax6.set_ylabel('junction-spanning reads\n(loose, mapQ≥20)', fontsize=6)
ax6.set_title('F   Chromosome state counts: smbu08 contains both states (%d:%d ≈ 5.6:1 retained:excised); '
              'smbu06 control is 24:1' % (v08[0], v08[1]), fontsize=7.2, loc='left')
ax6.legend(fontsize=5, frameon=False)
ax6.tick_params(labelsize=5.5)

for fmt, dpi in [('png', 600), ('pdf', None), ('svg', None)]:
    fig.savefig(os.path.join(FIG, fmt, 'FigS3_junction_reads.%s' % fmt), dpi=dpi, bbox_inches='tight')
print('figS3 done')
