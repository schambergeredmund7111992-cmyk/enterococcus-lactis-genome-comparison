# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""fig10_ani.py — Figure 10：ANI 热图（含参考论文使用的比较株 CX 2-6_2、IDCC 2105）
ANIb（1020 bp 片段，≥70% id & ≥70% 长度，双向平均）；计算一次后缓存 TSV。"""
import glob
import os
import sys
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK, IN06, IN08, BLASTN, MAKEDB, read_fasta, run, open_log
from fig_style import plt, savefig, C
import numpy as np

OUT = os.path.join(WORK, '09_figures')
FA = os.path.join(WORK, '12_functional_annotation')
os.makedirs(FA, exist_ok=True)
TMP = os.path.join(FA, 'tmp_ani')
os.makedirs(TMP, exist_ok=True)
LOG = open_log(os.path.join(FA, 'logs', 'fig10_ani.log'))
CACHE = os.path.join(FA, 'summary', 'ani_matrix.tsv')
FRAG = 1020

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
    ('KCTC21015', find_fna('GCF_015767715.1') or find_fna('GCA_015767715.1')),
    ('E_faecium_64-3', find_fna('GCF_001298485.1')),
    ('E_faecalis_LD33', find_fna('GCF_001598635.1')),
    ('L_lactis_MG1363', WORK + r'\02_taxonomy\refs\GCF_000009425.1\ncbi_dataset\data\GCF_000009425.1\GCF_000009425.1_ASM942v1_genomic.fna'),
]
genomes = {}
for name, path in SOURCES:
    if path and os.path.exists(path):
        genomes[name] = path
    else:
        LOG.write('MISSING source: %s (%s)\n' % (name, path))
names = list(genomes)
LOG.write('genomes: %s\n' % names)


def frag_write(path, out):
    n = 0
    with open(out, 'w') as f:
        for nm, s in read_fasta(path):
            for i in range(0, len(s) - FRAG + 1, FRAG):
                n += 1
                f.write('>%s_f%d\n%s\n' % (nm.split()[0], n, s[i:i + FRAG]))
    return n


def ani_one(qf, db):
    r = run([BLASTN, '-query', qf, '-db', db, '-outfmt', '6 qseqid pident length qlen',
             '-evalue', '1e-5', '-max_target_seqs', '5', '-num_threads', '8'], log=None)
    best = {}
    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        q, pid, alen, qlen = line.split('\t')
        if q not in best:
            best[q] = (float(pid), int(alen), int(qlen))
    tot = sum(1 for _ in open(qf)) // 2
    passed = [(p, a) for p, a, q in best.values() if p >= 70 and a >= 0.7 * q]
    if not passed:
        return None, 0, tot
    tl = sum(a for _, a in passed)
    return sum(p * a for p, a in passed) / tl, len(passed), tot


if os.path.exists(CACHE):
    LOG.write('reuse cached matrix\n')
else:
    frags = {}
    for n in names:
        qf = os.path.join(TMP, n + '_frags.fa')
        nf = frag_write(genomes[n], qf)
        frags[n] = qf
        db = os.path.join(TMP, n + '_db')
        run([MAKEDB, '-in', genomes[n], '-dbtype', 'nucl', '-out', db], log=None, check=True)
        LOG.write('%s frags=%d\n' % (n, nf))
    with open(CACHE, 'w', encoding='utf-8') as f:
        f.write('g1\tg2\tani\tnpass_ab\tnpass_ba\n')
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                a_ab, na, ta = ani_one(frags[a], os.path.join(TMP, b + '_db'))
                a_ba, nb, tb = ani_one(frags[b], os.path.join(TMP, a + '_db'))
                vals = [v for v in (a_ab, a_ba) if v]
                ani = sum(vals) / len(vals) if vals else 0.0
                f.write('%s\t%s\t%.2f\t%d\t%d\n' % (a, b, ani, na, nb))
                f.flush()
                LOG.write('%s vs %s: %.2f%%\n' % (a, b, ani))

# 绘图
M = np.eye(len(names)) * 100
with open(CACHE, encoding='utf-8') as f:
    next(f)
    for line in f:
        g1, g2, ani, _, _ = line.rstrip('\n').split('\t')
        i, j = names.index(g1), names.index(g2)
        M[i, j] = M[j, i] = float(ani)

fig, ax = plt.subplots(figsize=(6.4, 5.6))
im = ax.imshow(M, cmap='RdYlBu_r', vmin=80, vmax=100)
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=38, ha='right', fontsize=6)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=6)
for i in range(len(names)):
    for j in range(len(names)):
        v = M[i, j]
        ax.text(j, i, '%.2f' % v if i != j else '—', ha='center', va='center', fontsize=5.2,
                color='black' if v < 96 else 'white')
ax.set_title('Pairwise ANIb (%) — project strains, 194, and comparison genomes\n(including the strains used in the reference-style comparison: IDCC 2105, CX 2-6_2)', fontsize=7.5)
cb = fig.colorbar(im, ax=ax, shrink=0.8)
cb.set_label('ANIb (%)', fontsize=6)
cb.ax.tick_params(labelsize=5.5)
fig.text(0.5, 0.005, 'ANIb: 1,020-bp fragments, ≥70% identity and ≥70% fragment length, bidirectional average (JSpeciesWS-style). Values were cached in ani_matrix.tsv.',
         fontsize=5.6, ha='center')
savefig(fig, os.path.join(OUT, 'Fig10_ani', 'Fig10'))
print('fig10 done')
