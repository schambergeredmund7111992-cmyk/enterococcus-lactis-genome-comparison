# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
unit_1bp_test2.py — 1bp 差异 read 内配对检验（改进版：双参考帧 + 近邻对照）
相对于 unit_1bp_test.py：
  1) 同时使用 pl2_asassembled 与 dimer_back_internal 两个线性化帧，覆盖跨原点（回接）的长 reads
  2) 阴性对照 A-run 选择在目标 ±25 kb 内（同一批 read 可同时提供对照配对）
  3) 允许 secondary 比对（同一 read 只保留每帧最佳 AS 的一条），mapq≥20，read_len≥38,800
输入：02_mapping/panel_db.sam.gz（含 pl2_asassembled 与 dimer_back_internal）
输出：04_controls/unit_1bp_withinread_pairs2.tsv, unit_1bp_withinread2.tsv
"""
import gzip
import os
import random
import re
from collections import defaultdict

ATTJ = r'<PROJECT_ROOT>\longread_validation'
REFDIR = os.path.join(ATTJ, '01_references', 'references')
OUT = os.path.join(ATTJ, '04_controls')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'unit_1bp_test2.log'), 'a', encoding='utf-8')

CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')
B = 38719
DEL_POS = 36516
PL2 = 77437


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
assert len(pl2) == PL2


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
LOG.write('windows pl2-frame: W1=%d+%d W2=%d+%d\n' % (rs1, rl1, rs2, rl2))

# 对照：unit1 内、距目标 ≤25kb 的 A-run（≥6bp），排除目标 ±300
ctrls = []
for m in re.finditer(r'A{6,}', pl2[:B]):
    a, b = m.start(), m.end()
    if abs(a - DEL_POS) < 300 or abs(a - DEL_POS) > 25000:
        continue
    off = B - 1 if a > DEL_POS else B
    if a + off >= PL2:
        continue
    c2 = max_run_containing(pl2, a + off)
    if c2 is None:
        continue
    ctrls.append(dict(name='ctrl_%d' % len(ctrls), rs1=a, rl1=b - a, rs2=c2[0], rl2=c2[1]))
    if len(ctrls) >= 8:
        break
LOG.write('controls: %s\n' % [(c['name'], c['rs1'], c['rl1'], c['rs2'], c['rl2']) for c in ctrls])

FRAMES = [('pl2_asassembled', 0), ('dimer_back_internal', B)]   # (ref, rotation offset)
def fw(p, rot):
    return (p - rot) % PL2


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


def measure(seq, ops, r0, rs, rl):
    i0 = refpos_to_seqidx(ops, r0, rs)
    i1 = refpos_to_seqidx(ops, r0, rs + rl - 1)
    if i0 is None or i1 is None:
        return None
    j0, j1 = i0, i1
    while j0 > 0 and seq[j0 - 1] == 'A':
        j0 -= 1
    while j1 + 1 < len(seq) and seq[j1 + 1] == 'A':
        j1 += 1
    return j1 - j0 + 1


# read -> frame -> dict(group -> best (AS, a1, a2)), dedupe by read+frame+group keeping max AS
best = defaultdict(dict)
scanned = 0
with gzip.open(os.path.join(ATTJ, '02_mapping', 'panel_db.sam.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[0] == '@':
            continue
        p = line.rstrip('\n').split('\t')
        ref = p[2]
        frame = None
        for fname, rot in FRAMES:
            if ref == fname:
                frame = (fname, rot)
        if frame is None:
            continue
        if int(p[4]) < 20:
            continue
        seq = p[9]
        if seq == '*' or len(seq) < 38800:
            continue
        scanned += 1
        ops = [(int(a), b) for a, b in CIG_RE.findall(p[5])]
        r0 = int(p[3]) - 1
        asv = 0
        for t in p[11:]:
            if t.startswith('AS:i:'):
                asv = int(t[5:])
        fname, rot = frame
        a1 = measure(seq, ops, r0, fw(rs1, rot), rl1)
        a2 = measure(seq, ops, r0, fw(rs2, rot), rl2)
        if a1 is not None and a2 is not None:
            key = (p[0], 'target')
            cur = best[fname].get(key)
            if cur is None or asv > cur[0]:
                best[fname][key] = (asv, a1, a2)
        for c in ctrls:
            b1 = measure(seq, ops, r0, fw(c['rs1'], rot), c['rl1'])
            b2 = measure(seq, ops, r0, fw(c['rs2'], rot), c['rl2'])
            if b1 is not None and b2 is not None:
                key = (p[0], c['name'])
                cur = best[fname].get(key)
                if cur is None or asv > cur[0]:
                    best[fname][key] = (asv, b1, b2)

LOG.write('long panel alignments scanned=%d\n' % scanned)

# 合并两帧；同一 read 的 target 只取一次（优先 AS 高）
pairs = []
seen_target = {}
for fname, keyed in best.items():
    for (rd, grp), (asv, a1, a2) in keyed.items():
        if grp == 'target':
            if rd not in seen_target or asv > seen_target[rd][0]:
                seen_target[rd] = (asv, fname, a1, a2)
        else:
            pairs.append((grp, rd, fname, a1, a2, a1 - a2))
for rd, (asv, fname, a1, a2) in seen_target.items():
    pairs.append(('target', rd, fname, a1, a2, a1 - a2))

with open(os.path.join(OUT, 'unit_1bp_withinread_pairs2.tsv'), 'w', encoding='utf-8') as f:
    f.write('group\tread\tframe\trun_copy1\trun_copy2\tdelta(1-2)\n')
    for r_ in sorted(pairs):
        f.write('\t'.join(str(x) for x in r_) + '\n')


def boot(vals, n=5000, seed=5):
    if not vals:
        return None, None, None
    rng = random.Random(seed)
    m = sum(vals) / len(vals)
    o = []
    for _ in range(n):
        o.append(sum(vals[rng.randrange(len(vals))] for _ in range(len(vals))) / len(vals))
    o.sort()
    return m, o[int(n * 0.025)], o[int(n * 0.975)]


groups = defaultdict(list)
for g, rd, fr, a1, a2, d in pairs:
    groups[g].append(d)
rows = []
for g in ['target'] + [c['name'] for c in ctrls]:
    vals = groups.get(g, [])
    m, lo, hi = boot(vals)
    rows.append([g, len(vals)] + ([round(m, 3), round(lo, 3), round(hi, 3)] if vals else ['', '', '']))
with open(os.path.join(OUT, 'unit_1bp_withinread2.tsv'), 'w', encoding='utf-8') as f:
    f.write('# read 内配对（双帧合并）：delta = run(copy1窗口) − run(copy2窗口)\n')
    f.write('# 真实 1bp 差异 => target delta ≈ +1；组装重复假象 => ≈ 0；对照 A-run 应 ≈ 0\n')
    f.write('group\tn_reads\tmean_delta\tbootstrap_CI95_low\tbootstrap_CI95_high\n')
    for r_ in rows:
        f.write('\t'.join(str(x) for x in r_) + '\n')
print('unit_1bp_test2 done')
for r_ in rows:
    print(r_)
print('target pairs:')
for r_ in sorted(pairs):
    if r_[0] == 'target':
        print('  ', r_)
