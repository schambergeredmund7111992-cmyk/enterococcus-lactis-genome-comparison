# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
copy_number_corrected.py — 元件相对丰度（修正串联重复双计）与 bootstrap CI
问题：Plasmid2 为两近同源单元串联，元件 read 会对 pl2 产生 2 条比对（u1/u2），
      逐比对计深度会把元件丰度放大 2×；chr06 参考上还有第 3 条比对（保留态模型位点）。
方法：以 read 为单位；某 read 对某参考的贡献 = 其在该参考上的比对长度 × 1/(该 read 在该参考上的比对数)。
      A = Σ(元素 read 对 pl2 的贡献)（每个 read 至多贡献一个 read 长度）
      B = Σ(染色体 read 对 chr08 的贡献)
      元件/染色体碱基丰度比 = (A/L_pl2) / (B/L_chr08)；另给出按 38,719 bp 单体型与 77,437 bp 二聚体
      两种解释下的"分子拷贝数/染色体当量"。
bootstrap：对 read 重采样 1000×（对 A、B 联合重采样以保留相关性）。
输出：04_controls/copy_number_corrected.tsv
"""
import gzip
import os
import random
import re
from collections import defaultdict

ATTJ = r'<PROJECT_ROOT>\longread_validation'
MAP = os.path.join(ATTJ, '02_mapping')
OUT = os.path.join(ATTJ, '04_controls')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'copy_number_corrected.log'), 'a', encoding='utf-8')

CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')
L_CHR = 2638187
L_DIMER = 77437
L_MONOMER = 38719

# read -> ref -> 比对长度列表
aln = defaultdict(lambda: defaultdict(list))
with gzip.open(os.path.join(MAP, 'full_db.sam.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[0] == '@':
            continue
        p = line.rstrip('\n').split('\t')
        ref = p[2]
        if ref not in ('chr08_full', 'pl2_asassembled', 'chr06_full'):
            continue
        ops = [(int(a), b) for a, b in CIG_RE.findall(p[5])]
        rl = sum(n for n, o in ops if o in 'M=X')
        if rl > 0:
            aln[p[0]][ref].append(rl)
LOG.write('reads with alignments to chr08/pl2/chr06: %d\n' % len(aln))

# 元素 read：对 pl2 有比对的 read；染色体 read：对 chr08 有比对的 read（两组可重叠——保留真实重叠）
elem = {}
chrom = {}
both = 0
for rd, refs in aln.items():
    a = sum(x / len(refs['pl2_asassembled']) for x in refs.get('pl2_asassembled', []))
    b = sum(x / len(refs['chr08_full']) for x in refs.get('chr08_full', []))
    if 'pl2_asassembled' in refs:
        elem[rd] = a
    if 'chr08_full' in refs:
        chrom[rd] = b
    if 'pl2_asassembled' in refs and 'chr08_full' in refs:
        both += 1
LOG.write('element reads=%d chrom reads=%d overlap=%d\n' % (len(elem), len(chrom), both))

A = sum(elem.values())
B = sum(chrom.values())
r_bases = (A / L_DIMER) / (B / L_CHR)
r_dimer = (A / L_DIMER) / (B / L_CHR)
r_monomer = (A / L_MONOMER) / (B / L_CHR)


def boot(n=1000, seed=17):
    rng = random.Random(seed)
    ek = list(elem.items()); ck = list(chrom.items())
    ne, nc = len(ek), len(ck)
    out = []
    for _ in range(n):
        a = sum(ek[rng.randrange(ne)][1] for _ in range(ne))
        b = sum(ck[rng.randrange(nc)][1] for _ in range(nc))
        if b > 0:
            out.append((a / L_DIMER) / (b / L_CHR))
    out.sort()
    return out[int(n * 0.025)], out[int(n * 0.975)]

lo, hi = boot()
rows = [
    ('element_bases_A', A, 'Σ read×pl2 比对长度 / 该 read 的 pl2 比对数'),
    ('chromosome_bases_B', B, 'Σ read×chr08 比对长度 / 该 read 的 chr08 比对数'),
    ('ratio_A_per_base_B_per_base', round(r_bases, 3), '（A/77,437）/(B/2,638,187)：按元件 77,437 bp 计'),
    ('bootstrap_CI95_low', round(lo, 3), 'read 重采样 1000×'),
    ('bootstrap_CI95_high', round(hi, 3), ''),
    ('molecule_copies_per_chromosome_dimer_interpretation', round(r_dimer, 2), '若为 77,437 bp 二聚体分子'),
    ('molecule_copies_per_chromosome_monomer_interpretation', round(r_monomer, 2), '若为 38,719 bp 单体分子'),
]
with open(os.path.join(OUT, 'copy_number_corrected.tsv'), 'w', encoding='utf-8') as f:
    f.write('# 元件相对丰度（修正串联重复双计）；高丰度不等于自主复制\n')
    f.write('item\tvalue\tdetail\n')
    for r_ in rows:
        f.write('\t'.join(str(x) for x in r_) + '\n')
print('copy_number_corrected done')
for r_ in rows:
    print(r_)
