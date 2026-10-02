# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
verify_structure.py — 独立复核（第二实现，字节级，不依赖任何比对器的边界放置）
  1) 复制子身份：长度 + 记录数核验
  2) 假设无关地求解：chr06 中是否存在长度 38,719 bp 的片段 S，使 chr08 == chr06−S（字节级）
     —— 用 rc(unit1) 的 200bp 种子在 chr06 中定位 S，再全长度验证
  3) 由 S 的精确位置导出 attL/attJ/attR，做字节比较并报告最大相同窗口
  4) pl2 两单元与 1bp 缺失（纯字节）
  5) 旋转等价：unit1/unit2 vs rc(S)、rc(attL+X)、rc(X+attR)、rc(attL+X+attR+…) 候选
  6) minimap2（asm10）作为独立比对器交叉确认 chr06/chr08 共线性与缺口
  7) 已存参考 .fasta 与本次重算的字节一致性
"""
import os
import re
import subprocess
import hashlib

ATTJ = r'<PROJECT_ROOT>\longread_validation'
IN08 = r'<PROJECT_ROOT>\smbu08\smbu08'
IN06 = r'<PROJECT_ROOT>\smbu06\smbu06'
MM2 = r'<TOOLS_ROOT>\minimap2-win\minimap2-31-r1302-windows-x86_64-ucrt64\minimap2.exe'
MM2 = r'<TOOLS_ROOT>\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe'
OUT = os.path.join(ATTJ, '01_references', 'independent_verification.tsv')
REFDIR = os.path.join(ATTJ, '01_references', 'references')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'verify_structure.log'), 'a', encoding='utf-8')

rows = []
def rec(k, v, d=''):
    rows.append((k, str(v), d))
    LOG.write('%s\t%s\t%s\n' % (k, v, d))
    print('%s = %s  %s' % (k, v, d), flush=True)


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
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


def rc(s):
    return s.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def max_common(a, b):
    """从 0 起的最长相等前缀长度"""
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


# ---------- 1) 复制子身份 ----------
c6 = read_fasta(os.path.join(IN06, '2.Assembly', 'smbu06.Chromosome.fasta'))[0][1]
c8 = read_fasta(os.path.join(IN08, '2.Assembly', 'smbu08.Chromosome.fasta'))[0][1]
pl06 = read_fasta(os.path.join(IN06, '2.Assembly', 'smbu06.Plasmid.fasta'))
pl08 = read_fasta(os.path.join(IN08, '2.Assembly', 'smbu08.Plasmid.fasta'))
rec('pl06_records', len(pl06), str([(n, len(s)) for n, s in pl06]))
rec('pl08_records', len(pl08), str([(n, len(s)) for n, s in pl08]))
assert len(c6) == 2676906 and len(c8) == 2638187, 'chromosome length mismatch'
p2 = [s for n, s in pl08 if len(s) == 77437][0]
p1_08 = [s for n, s in pl08 if len(s) == 128837][0]
p1_06 = [s for n, s in pl06 if len(s) == 128837][0]
rec('replicon_identity', 'ok', 'chr06=2676906 chr08=2638187 pl2=77437 pl1=128837（按长度唯一选取）')
rec('chr06_minus_chr08_len', len(c6) - len(c8), '期望 38719')

# ---------- 2) 定位 S：chr06 中被删除的 38,719 bp（由 chr08 右翼唯一种子定位） ----------
B = 38719
u1, u2 = p2[:B], p2[B:]
seed_r = c8[-500:]
hits = [m.start() for m in re.finditer('(?=%s)' % seed_r, c6)]
rec('chr08_right500_hits_in_chr06', len(hits), 'positions=%s' % hits)
assert len(hits) == 1, 'right-flank seed not unique'
assert hits[0] == len(c8) - 500 + B, 'tail offset != length difference'   # 尾部 c8[x] == c6[x+B]
# 二分查找 S_pos：性质 P(x) := c8[x:x+1000] == c6[x+B:x+B+1000]，x>=S_pos 时为真
K = 1000
def P(x):
    return c8[x:x + K] == c6[x + B:x + B + K]
lo, hi = 0, len(c8) - K
assert not P(lo) and P(hi), 'P(x) not monotone as expected'
while hi - lo > 1:
    mid = (lo + hi) // 2
    if P(mid):
        hi = mid
    else:
        lo = mid
S_pos = hi
rec('S_pos_0based', S_pos, 'S 起点（0-based），二分定位：P(%d)=False, P(%d)=True' % (lo, hi))
# chr08 == chr06 - S ?（左侧应 0 差异；右侧允许记录所有 SNP 位点）
l_diff = [i for i in range(S_pos) if c8[i] != c6[i]]
rlen = len(c8) - S_pos
r_diff = [i for i in range(rlen) if c8[S_pos + i] != c6[S_pos + B + i]]
rec('chr08_vs_chr06_minus_S_left_diffs', len(l_diff), 'chr08[0:%d] vs chr06[0:%d]（字节级差异数）' % (S_pos, S_pos))
rec('chr08_vs_chr06_minus_S_right_diffs', len(r_diff),
    'chr08[%d:] vs chr06[%d:]（%d bp）差异位点(0-based, 相对右翼起点)=%s' % (S_pos, S_pos + B, rlen, r_diff[:20]))
for j in r_diff[:20]:
    rec('right_flank_snp_%d' % j,
        'chr08:%d(%s) vs chr06:%d(%s)' % (S_pos + j + 1, c8[S_pos + j], S_pos + B + j + 1, c6[S_pos + B + j]),
        'chr08 == chr06 删除 S 后 + 该 SNP')
assert len(l_diff) == 0, 'left flank not identical'
S = c6[S_pos:S_pos + B]
rec('S_pos_0based', S_pos, 'S = chr06[%d:%d]，1-based %d-%d，长度 %d' % (S_pos, S_pos + B, S_pos + 1, S_pos + B, B))

# ---------- 3) att 结构 ----------
attlen = 146
attL = S[:attlen]                       # chr06 side, front copy
attR = c6[S_pos + B:S_pos + B + attlen]  # chr06 side, back copy
attJ = c8[S_pos:S_pos + attlen]          # chr08 junction copy
rec('attL', 'chr06 %d-%d' % (S_pos + 1, S_pos + attlen), attL)
rec('attR', 'chr06 %d-%d' % (S_pos + B + 1, S_pos + B + attlen), attR)
rec('attJ', 'chr08 %d-%d' % (S_pos + 1, S_pos + attlen), attJ)
rec('attJ_eq_attR_byte', attJ == attR, '146bp 字节比较')
rec('attJ_eq_attL_byte', attJ == attL, '146bp 字节比较')
dl = [i for i in range(attlen) if attL[i] != attJ[i]]
dr = [i for i in range(attlen) if attL[i] != attR[i]]
rec('attL_vs_attJ_diffs', len(dl), 'positions(1-based)=%s' % ';'.join(str(i + 1) for i in dl))
rec('attL_vs_attR_diffs', len(dr), 'positions(1-based)=%s' % ';'.join(str(i + 1) for i in dr))
# 最大相同窗口（边界是否为恰好 146bp）
ext_r = attlen + max_common(c8[S_pos + attlen:S_pos + 400], c6[S_pos + B + attlen:S_pos + B + 400])
ext_l = attlen
while ext_l < 400 and S_pos - (ext_l - attlen) > 0 and c8[S_pos - (ext_l - attlen) - 1] == c6[S_pos + B - (ext_l - attlen) - 1]:
    ext_l += 1
rec('attJ_attR_identical_window', ext_r, '从 attJ 起点向右延伸的最大完全相同窗口（bp）')
rec('attJ_attR_identical_window_upstream', ext_l, '向两侧合计最大相同窗口（bp，含 146 基线）')

X = S[attlen:]
rec('X_internal_len', len(X), 'chr06 %d-%d' % (S_pos + attlen + 1, S_pos + B))
rec('element_total_len', len(S) + attlen, 'chr06 %d-%d = attL+X+attR' % (S_pos + 1, S_pos + B + attlen))

# ---------- 4) pl2 两单元与 1bp 差异 ----------
rec('unit1_len', len(u1), 'pl2[0:%d]' % B)
rec('unit2_len', len(u2), 'pl2[%d:]' % B)
k = next((i for i in range(len(u2)) if u1[i] != u2[i]), None)
assert k is not None and u2 == u1[:k] + u1[k + 1:], 'unit2 != unit1 with single-base deletion'
rec('unit2_is_unit1_minus_1bp', True, 'unit1 位置 %d (1-based) 缺失碱基 %s' % (k + 1, u1[k]))
rec('unit1_around_del', u1[k - 15:k + 16], '缺失位点上下文（unit1 坐标 %d-%d）' % (k - 14, k + 16))

# ---------- 5) 旋转等价 ----------
def rot_offsets(target, seq):
    if len(target) != len(seq):
        return None
    d = target + target
    out, st = [], 0
    while True:
        i = d.find(seq, st)
        if i < 0 or i >= len(target):
            break
        out.append(i)
        st = i + 1
    return out

for nm, seq in [('unit1', u1), ('unit2', u2)]:
    for tnm, tgt in [('rc(attL+X)', rc(attL + X)), ('rc(X+attR)', rc(X + attR)),
                     ('rc(element)', rc(S + attR)), ('rc(S)', rc(S))]:
        off = rot_offsets(tgt, seq)
        rec('%s_rot_eq_%s' % (nm, tnm), bool(off), 'offsets=%s' % off)
rec('unit2_rot_eq_unit1', bool(rot_offsets(u1, u2)), 'offsets=%s' % rot_offsets(u1, u2))

# ---------- 6) minimap2 交叉确认（独立比对器） ----------
tmp = os.path.join(REFDIR, '_v2_chr06.fa')
with open(tmp, 'w') as f:
    f.write('>chr06\n%s\n' % c6)
r = subprocess.run([MM2, '-x', 'asm10', '--eqx', '-c', '-t', '8', tmp,
                    os.path.join(IN08, '2.Assembly', 'smbu08.Chromosome.fasta')],
                   capture_output=True, text=True, encoding='utf-8', errors='replace')
LOG.write('minimap2 rc=%d\n' % r.returncode)
paf = [l.split('\t') for l in r.stdout.strip().split('\n') if l]
paf = [p for p in paf if int(p[1]) > 1000000 and p[4] == '+']
paf.sort(key=lambda x: -int(x[1]))
rec('mm2_long_alignments', len(paf),
    str([(x[5], 'q%s-%s' % (x[2], x[3]), 's%s-%s' % (x[7], x[8])) for x in paf]))
if paf:
    p_ = paf[0]
    rec('mm2_best_span', int(p_[3]) - int(p_[2]), 'query %s-%s, target %s-%s' % (p_[2], p_[3], p_[7], p_[8]))

# ---------- 7) 已存参考一致性 ----------
j_att_mid = S_pos + attlen // 2
expect = {
    'chr08_full': c8, 'chr06_full': c6, 'pl2_asassembled': p2, 'pl1_smbu08': p1_08,
    'attJ_junction': c8[j_att_mid - 5000:j_att_mid + 5000],
    'retained_window': c6[S_pos - 5000:S_pos + B + attlen + 5000],
    'unit1': u1, 'unit2': u2,
    'unit1_to_unit2_junction': p2[B - 5000:B + 5000],
    'unit2_to_unit1_circular_junction': p2[-5000:] + p2[:5000],
    'monomer_circle_repA': u1,
    'monomer_circle_repB': u1[len(u1) // 2:] + u1[:len(u1) // 2],
    'dimer_back_internal': p2[B:] + p2[:B],
}
for name, seq in expect.items():
    stored = read_fasta(os.path.join(REFDIR, name + '.fasta'))[0][1]
    eq = stored == seq
    rec('stored_ref_%s_matches' % name, eq, 'len stored=%d recompute=%d' % (len(stored), len(seq)))
    if not eq:
        dd = [i for i in range(min(len(stored), len(seq))) if stored[i] != seq[i]]
        rec('stored_ref_%s_diffpos' % name, len(dd), 'first %s' % dd[:20])

# 阴性对照来源核验
negmeta = {}
for line in open(os.path.join(ATTJ, '01_references', 'reference_construction.tsv'), encoding='utf-8').read().splitlines()[1:]:
    p_ = line.split('\t')
    negmeta[p_[0]] = p_[3]
for nm in ['neg_attJ_shuffled', 'neg_ctrl_chr08_1', 'neg_ctrl_chr08_2', 'neg_ctrl_chr08_3']:
    s = read_fasta(os.path.join(REFDIR, nm + '.fasta'))[0][1]
    if nm.startswith('neg_ctrl_chr08'):
        pos = int(re.search(r'pos (\d+)', negmeta[nm]).group(1))
        rec('neg_ref_%s_exact_window' % nm, s == c8[pos - 5000:pos + 5000], 'pos=%d' % pos)
    else:
        rec('neg_ref_%s_same_composition' % nm,
            sorted(s) == sorted(c8[j_att_mid - 5000:j_att_mid + 5000]), '与 attJ 窗口碱基组成一致')

# 共享质粒对照（旋转/链等价）
d_f = p1_08 + p1_08
off_f = d_f.find(p1_06)
rec('pl1_08_eq_pl1_06_rotation', off_f >= 0 or (rc(p1_08) + rc(p1_08)).find(p1_06) >= 0,
    '正向旋转偏移=%d' % off_f)

with open(OUT, 'w', encoding='utf-8') as f:
    f.write('item\tvalue\tdetail\n')
    for r_ in rows:
        f.write('\t'.join(r_) + '\n')
LOG.write('DONE\n')
print('verify_structure done ->', OUT)
