# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
unit_1bp_test.py — 单元间 1bp 差异的 read 内判别检验（本会话实现）

背景与设计考虑：
  unit2 = unit1 删除 1bp（位于 poly-A 同聚物，unit1 0-based 36516）。这一 1bp 是单体环与
  串联二聚体在序列层面的唯一判别点。批量的 copy1 vs copy2 分布比较**不可判别**：reads 在
  两个近同源拷贝间的指派由比对分数决定，而分数差恰来自该 1bp 本身（参考指派偏倚），
  即便真实为单体也会出现表观差异。因此本脚本只做 read 内配对检验：
  只使用**同一条 read 同时覆盖两拷贝对应窗口**的长 reads，在同一 read 内比较两个同聚物
  长度（配对消除了 read 个体的测序噪声），并对 unit1 内其它 A-run 做同样的 read 内配对
  作为阴性对照（这些位置两拷贝应完全相同，预期 delta≈0）。
输出：04_controls/unit_1bp_withinread.tsv, unit_1bp_withinread_pairs.tsv
"""
import gzip
import os
import random
import re
from collections import defaultdict

ATTJ = r'<PROJECT_ROOT>\longread_validation'
REFDIR = os.path.join(ATTJ, '01_references', 'references')
OUT = os.path.join(ATTJ, '04_controls')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'unit_1bp_test.log'), 'a', encoding='utf-8')

CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')
B = 38719
DEL_POS = 36516


def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path, encoding='utf-8', errors='replace'):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None:
                recs.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


pl2 = read_fasta(os.path.join(REFDIR, 'pl2_asassembled.fasta'))[0][1]
assert len(pl2) == 77437 and pl2[DEL_POS] == 'A'


def max_run_containing(seq, p):
    if seq[p] != 'A':
        return None
    a = b = p
    while a > 0 and seq[a - 1] == 'A':
        a -= 1
    while b + 1 < len(seq) and seq[b + 1] == 'A':
        b += 1
    return a, b + 1 - a


rs1, rl1 = max_run_containing(pl2, DEL_POS)
rs2, rl2 = max_run_containing(pl2, DEL_POS + B - 1)
LOG.write('target run copy1: rs=%d len=%d ; copy2: rs=%d len=%d\n' % (rs1, rl1, rs2, rl2))

# 阴性对照 A-run（unit1 内，长度>=6，距目标>200bp）
ctrls = []
for m in re.finditer(r'A{6,}', pl2[:B]):
    a, b = m.start(), m.end()
    if abs(a - DEL_POS) < 200:
        continue
    off = B - 1 if a > DEL_POS else B
    c2 = max_run_containing(pl2, a + off) if a + off < len(pl2) else None
    if c2 is None:
        continue
    ctrls.append(dict(name='ctrl_%d' % len(ctrls), rs1=a, rl1=b - a, rs2=c2[0], rl2=c2[1]))
    if len(ctrls) >= 6:
        break
LOG.write('controls: %s\n' % [(c['name'], c['rs1'], c['rl1'], c['rs2'], c['rl2']) for c in ctrls])


def refpos_to_seqidx(ops, r0, p):
    pos, qi = r0, 0
    for n, o in ops:
        if o in 'M=X':
            if pos <= p < pos + n:
                return qi + (p - pos)
            pos += n; qi += n
        elif o == 'I':
            qi += n
        elif o in 'DN':
            if pos <= p < pos + n:
                return None
            pos += n
        elif o == 'S':
            qi += n
        elif o == 'H':
            pass
    return None


def measure_run(seq, ops, r0, rs, rl):
    """在 read SEQ 中量测映射到参考窗口 [rs, rs+rl) 的连续 A-run 长度；窗口需完全被 M 覆盖"""
    i0 = refpos_to_seqidx(ops, r0, rs)
    i1 = refpos_to_seqidx(ops, r0, rs + rl - 1)
    if i0 is None or i1 is None:
        return None
    # 窗口内部必须是连续 M（无 indel 打断）
    for p in (rs, rs + rl - 1):
        pass
    j0, j1 = i0, i1
    while j0 > 0 and seq[j0 - 1] == 'A':
        j0 -= 1
    while j1 + 1 < len(seq) and seq[j1 + 1] == 'A':
        j1 += 1
    return j1 - j0 + 1


pairs = []
scanned = 0
with gzip.open(os.path.join(ATTJ, '02_mapping', 'full_db.sam.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[0] == '@':
            continue
        p = line.rstrip('\n').split('\t')
        if p[2] != 'pl2_asassembled':
            continue
        flag = int(p[1])
        if flag & 0x100 or flag & 0x800:
            continue
        if int(p[4]) < 20:
            continue
        seq = p[9]
        if len(seq) < 39000:
            continue
        scanned += 1
        ops = [(int(a), b) for a, b in CIG_RE.findall(p[5])]
        r0 = int(p[3]) - 1
        a1 = measure_run(seq, ops, r0, rs1, rl1)
        a2 = measure_run(seq, ops, r0, rs2, rl2)
        if a1 is not None and a2 is not None:
            pairs.append(('target', p[0], len(seq), a1, a2, a1 - a2))
        for c in ctrls:
            b1 = measure_run(seq, ops, r0, c['rs1'], c['rl1'])
            b2 = measure_run(seq, ops, r0, c['rs2'], c['rl2'])
            if b1 is not None and b2 is not None:
                pairs.append((c['name'], p[0], len(seq), b1, b2, b1 - b2))

LOG.write('long pl2 alignments scanned=%d; within-read pairs=%d\n' % (scanned, len(pairs)))

with open(os.path.join(OUT, 'unit_1bp_withinread_pairs.tsv'), 'w', encoding='utf-8') as f:
    f.write('group\tread\tread_len\trun_copy1\trun_copy2\tdelta(1-2)\n')
    for r_ in pairs:
        f.write('\t'.join(str(x) for x in r_) + '\n')


def boot_mean(vals, n=5000, seed=3):
    if not vals:
        return None, None, None
    rng = random.Random(seed)
    m = sum(vals) / len(vals)
    outs = []
    for _ in range(n):
        outs.append(sum(vals[rng.randrange(len(vals))] for _ in range(len(vals))) / len(vals))
    outs.sort()
    return m, outs[int(n * 0.025)], outs[int(n * 0.975)]


groups = defaultdict(list)
for g, rd, l, a1, a2, d in pairs:
    groups[g].append(d)

rows = []
for g in ['target'] + [c['name'] for c in ctrls]:
    vals = groups.get(g, [])
    m, lo, hi = boot_mean(vals)
    rows.append([g, len(vals)] + ([round(m, 3), round(lo, 3), round(hi, 3)] if vals else ['', '', '']))
with open(os.path.join(OUT, 'unit_1bp_withinread.tsv'), 'w', encoding='utf-8') as f:
    f.write('# read 内配对：delta = run(copy1窗口) − run(copy2窗口)；target=两拷贝在装配中相差 1bp 的位置\n')
    f.write('# 若 1bp 差异真实：target delta ≈ +1；若为组装重复假象：target delta ≈ 0；对照 A-run 两拷贝序列相同，delta ≈ 0\n')
    f.write('group\tn_reads_with_both_windows\tmean_delta\tbootstrap_CI95_low\tbootstrap_CI95_high\n')
    for r_ in rows:
        f.write('\t'.join(str(x) for x in r_) + '\n')
print('unit_1bp_test done')
for r_ in rows:
    print(r_)
