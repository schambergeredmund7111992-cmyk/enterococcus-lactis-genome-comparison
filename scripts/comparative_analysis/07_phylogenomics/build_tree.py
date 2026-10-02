# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
build_tree.py — 模块 07：核心基因超级矩阵与最大似然树
1) 每标记 MAFFT --auto
2) 列过滤：去除 gap 比例 >20% 的列
3) 标记保留：覆盖 ≥80% 基因组；基因组保留：出现在 ≥80% 保留标记
4) 超级矩阵 + IQ-TREE（LG+G，UFBoot 1000）
输出：trees/marker_aln/, trees/supermatrix.{phy,partition}, trees/core_tree.*
"""
import os
import sys
import glob

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, MAFFT, IQTREE, read_fasta, write_fasta, run, open_log

OUT = os.path.join(WORK, '07_phylogenomics')
MKD = os.path.join(OUT, 'markers')
TRD = os.path.join(OUT, 'trees')
ALD = os.path.join(TRD, 'marker_aln')
os.makedirs(ALD, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'build_tree.log'))

marker_files = sorted(glob.glob(os.path.join(MKD, 'marker_*.faa')))
LOG.write('marker files: %d\n' % len(marker_files))


def align_one(mf):
    tag = os.path.basename(mf)[7:-4]
    out = os.path.join(ALD, tag + '.aln.faa')
    if os.path.exists(out) and os.path.getsize(out) > 100:
        return tag, 'cached'
    r = run(['cmd', '/c', MAFFT, '--quiet', '--auto', mf], log=None)
    if r.returncode != 0 or not r.stdout.startswith('>'):
        return tag, 'FAIL rc=%d' % r.returncode
    with open(out, 'w') as f:
        f.write(r.stdout)
    return tag, 'ok'


from concurrent.futures import ThreadPoolExecutor
n_ok = 0
with ThreadPoolExecutor(max_workers=6) as ex:
    for tag, status in ex.map(align_one, marker_files):
        if status in ('ok', 'cached'):
            n_ok += 1
        else:
            LOG.write('MAFFT %s: %s\n' % (tag, status))
LOG.write('aligned markers: %d/%d\n' % (n_ok, len(marker_files)))

aligned = {}
for mf in marker_files:
    tag = os.path.basename(mf)[7:-4]
    out = os.path.join(ALD, tag + '.aln.faa')
    if os.path.exists(out) and os.path.getsize(out) > 100:
        aligned[tag] = read_fasta(out)
LOG.write('loaded aligned: %d\n' % len(aligned))

# 列过滤（gap>20% 去除）
def trim(recs):
    L = len(recs[0][1])
    keep = []
    for i in range(L):
        gaps = sum(1 for _, s in recs if s[i] == '-')
        if gaps / len(recs) <= 0.2:
            keep.append(i)
    return [('>'.join(n.split('|')[:1]) if False else n, ''.join(s[i] for i in keep)) for n, s in recs]


trimmed = {tag: trim(recs) for tag, recs in aligned.items()}
# 标记保留（≥80% 基因组）
N_GEN = 209
sel_markers = [tag for tag, recs in trimmed.items() if len(recs) >= 0.8 * N_GEN]
LOG.write('selected markers (>=80%% coverage of %d genomes): %d\n' % (N_GEN, len(sel_markers)))

# 基因组保留（≥80% 标记）
from collections import Counter
cnt = Counter()
for tag in sel_markers:
    for n, _ in trimmed[tag]:
        cnt[n.split('|')[0]] += 1
keep_genomes = sorted(g for g, c in cnt.items() if c >= 0.8 * len(sel_markers))
LOG.write('kept genomes: %d\n' % len(keep_genomes))

# 超级矩阵
seqs = {g: [] for g in keep_genomes}
part_lines = []
pos = 0
for tag in sel_markers:
    recs = dict((n.split('|')[0], s) for n, s in trimmed[tag])
    L = len(next(iter(recs.values())))
    for g in keep_genomes:
        seqs[g].append(recs.get(g, '-' * L))
    pos += L
    part_lines.append('LG+G, %s = %d-%d' % (tag, pos - L + 1, pos))
nchar = pos
LOG.write('supermatrix: %d genomes x %d aa\n' % (len(keep_genomes), nchar))

phy = os.path.join(TRD, 'supermatrix.fa')
with open(phy, 'w') as f:
    for g in keep_genomes:
        f.write('>%s\n%s\n' % (g, ''.join(seqs[g])))
with open(os.path.join(TRD, 'supermatrix.partition'), 'w') as f:
    f.write('\n'.join(part_lines) + '\n')

# IQ-TREE（用 FASTA，避免 PHYLIP 解析歧义）
prefix = os.path.join(TRD, 'core_tree')
r = run([IQTREE, '-s', phy, '-m', 'LG+G', '-B', '1000', '-T', '8', '--prefix', prefix,
         '-redo'], log=LOG, check=True)
LOG.write('IQ-TREE done: %s\n' % (prefix + '.treefile'))
print('build_tree done')
