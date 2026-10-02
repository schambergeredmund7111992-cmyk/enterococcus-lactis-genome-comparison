# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
unit_1bp_test3.py — 1bp 差异 read 内配对检验（泛化版：从任意参考独立定位两个窗口）
思路：read 内配对需要把两个同聚物窗口映射回 read SEQ；每个窗口只需**某一条**比对能定位它
（不必同一比对同时覆盖两窗口）。因此把 panel 全部含窗口的参考建表，从任意参考收集
W1/W2 的 read 内量测，按 read 配对。
输出：04_controls/unit_1bp_withinread_pairs3.tsv, unit_1bp_withinread3.tsv
"""
import gzip
import os
import random
import re
from collections import defaultdict

ATTJ = r'<PROJECT_ROOT>\longread_validation'
REFDIR = os.path.join(ATTJ, '01_references', 'references')
OUT = os.path.join(ATTJ, '04_controls')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'unit_1bp_test3.log'), 'a', encoding='utf-8')

CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')
B = 38719
PL2 = 77437
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
V = {}
for line in open(os.path.join(ATTJ, '01_references', 'independent_verification.tsv'), encoding='utf-8').read().splitlines()[1:]:
    p = line.split('\t')
    V[p[0]] = p[1]
S_pos = int(V['S_pos_0based'])


def max_run_containing(seq, p):
    if p >= len(seq) or seq[p] != 'A':
        return None
    a = b = p
    while a > 0 and seq[a - 1] == 'A':
        a -= 1
    while b + 1 < len(seq) and seq[b + 1] == 'A':
        b += 1
    return a, b + 1 - a


RS1, RL1 = max_run_containing(pl2, DEL_POS)          # e.g. 36509, 8
RS2, RL2 = max_run_containing(pl2, DEL_POS + B - 1)  # e.g. 75228, 7
LOG.write('ref runs: W1(pl2)=%d+%d W2(pl2)=%d+%d\n' % (RS1, RL1, RS2, RL2))

# ref -> (w1_pos or None, w2_pos or None)，配合各参考的从 pl2 坐标映射函数
def idfn(p):
    return p
def rot(p):
    return (p - B) % PL2
def u1m(p):
    return p if p < B else None
def u2m(p):
    return p - B if p >= B else None
def repB(p):
    return (p - 19359) % B if p < B else None
def jfwd(p):
    return p - 33719 if 33719 <= p < 43719 else None
def jback(p):
    if p >= PL2 - 5000:
        return p - (PL2 - 5000)
    if p < 5000:
        return p + 5000
    return None
def retained(p):
    # chr06 坐标：element 起点 S_pos 对应 j=5000
    return p - S_pos + 5000

REFS = {
    'pl2_asassembled': (idfn, RS1, RS2),
    'dimer_back_internal': (rot, RS1, RS2),
    'unit1': (u1m, RS1, None),
    'unit2': (u2m, None, RS2),
    'monomer_circle_repA': (u1m, RS1, None),
    'monomer_circle_repB': (repB, RS1, None),
    'unit1_to_unit2_junction': (jfwd, RS1, RS2),
    'unit2_to_unit1_circular_junction': (jback, RS1, RS2),
    'retained_window': (retained, RS1, None),
}


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


m1 = {}   # read -> (AS, run)
m2 = {}
with gzip.open(os.path.join(ATTJ, '02_mapping', 'panel_db.sam.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[0] == '@':
            continue
        p = line.rstrip('\n').split('\t')
        ref = p[2]
        if ref not in REFS:
            continue
        if int(p[4]) < 20:
            continue
        seq = p[9]
        if seq == '*':
            continue
        asv = 0
        for t in p[11:]:
            if t.startswith('AS:i:'):
                asv = int(t[5:])
        fn, w1, w2 = REFS[ref]
        ops = [(int(a), b) for a, b in CIG_RE.findall(p[5])]
        r0 = int(p[3]) - 1
        rd = p[0]
        if w1 is not None:
            q = fn(w1)
            if q is not None:
                v = measure(seq, ops, r0, q, RL1)
                if v is not None and (rd not in m1 or asv > m1[rd][0]):
                    m1[rd] = (asv, v, ref)
        if w2 is not None:
            q = fn(w2)
            if q is not None:
                v = measure(seq, ops, r0, q, RL2)
                if v is not None and (rd not in m2 or asv > m2[rd][0]):
                    m2[rd] = (asv, v, ref)

pairs = []
for rd in set(m1) & set(m2):
    pairs.append(('target', rd, m1[rd][2], m2[rd][2], m1[rd][1], m2[rd][1], m1[rd][1] - m2[rd][1]))
LOG.write('reads with W1 measurement=%d, W2=%d, both=%d\n' % (len(m1), len(m2), len(pairs)))

with open(os.path.join(OUT, 'unit_1bp_withinread_pairs3.tsv'), 'w', encoding='utf-8') as f:
    f.write('group\tread\tref_W1\tref_W2\trun_copy1\trun_copy2\tdelta(1-2)\n')
    for r_ in sorted(pairs):
        f.write('\t'.join(str(x) for x in r_) + '\n')


def boot(vals, n=5000, seed=9):
    if not vals:
        return None, None, None
    rng = random.Random(seed)
    m = sum(vals) / len(vals)
    o = []
    for _ in range(n):
        o.append(sum(vals[rng.randrange(len(vals))] for _ in range(len(vals))) / len(vals))
    o.sort()
    return m, o[int(n * 0.025)], o[int(n * 0.975)]


vals = [r_[6] for r_ in pairs]
m, lo, hi = boot(vals)
with open(os.path.join(OUT, 'unit_1bp_withinread3.tsv'), 'w', encoding='utf-8') as f:
    f.write('# read 内配对（泛化定位）；真实 1bp 差异 => delta ≈ +1；组装重复假象 => ≈ 0\n')
    f.write('group\tn_reads\tmean_delta\tbootstrap_CI95_low\tbootstrap_CI95_high\n')
    if vals:
        f.write('target\t%d\t%.3f\t%.3f\t%.3f\n' % (len(vals), m, lo, hi))
    else:
        f.write('target\t0\t\t\t\n')
print('unit_1bp_test3 done; target pairs=%d' % len(vals))
if vals:
    print('mean delta=%.3f CI=(%.3f, %.3f)' % (m, lo, hi))
for r_ in sorted(pairs):
    print('  ', r_)
